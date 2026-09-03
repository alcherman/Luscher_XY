"""Constrained XY-model Markov updates."""

import math
import numbers

import numpy as np

from supervillain.action.xy import XY
from supervillain.generator import Generator
from supervillain.generator._xy_kernels import sweep_parallel, sweep_serial
from supervillain.generator._xy_cluster import cluster_serial
from supervillain.generator.reference_implementation.xy import reference_sweep
from supervillain.generator.reference_implementation.xy_cluster import reference_cluster
from supervillain.h5 import ReadWriteable
from supervillain.lattice import Form


_TWO_PI = 2.0 * math.pi


def _flat_neighbours(lattice):
    """Return forward and backward periodic nearest-neighbour tables."""
    sites = np.arange(lattice.sites, dtype=np.int64).reshape(lattice.dims)
    plus = np.empty((lattice.D, lattice.sites), dtype=np.int64)
    minus = np.empty_like(plus)
    for mu in range(lattice.D):
        plus[mu] = np.roll(sites, -1, axis=mu).reshape(-1)
        minus[mu] = np.roll(sites, +1, axis=mu).reshape(-1)
    return plus, minus


def _flat_colours(lattice):
    """Pack the lattice's valid checkerboard groups into flat arrays."""
    groups = tuple(
        np.ravel_multi_index(color, lattice.dims).astype(np.int64, copy=False)
        for color in lattice.checkerboarding
    )
    lengths = np.fromiter((len(group) for group in groups), dtype=np.int64)
    offsets = np.empty(len(groups) + 1, dtype=np.int64)
    offsets[0] = 0
    np.cumsum(lengths, out=offsets[1:])
    color_sites = np.concatenate(groups) if groups else np.empty(0, dtype=np.int64)

    if len(color_sites) != lattice.sites or len(np.unique(color_sites)) != lattice.sites:
        raise RuntimeError("Lattice.checkerboarding did not partition its sites")
    return color_sites, offsets


def _same_lattice(left, right):
    """Compare lattice geometry across independently H5-restored objects."""
    return left is right or (
        getattr(left, "D", None) == getattr(right, "D", None)
        and getattr(left, "N", None) == getattr(right, "N", None)
    )


class SiteUpdate(ReadWriteable, Generator):
    r"""A colour-parallel constrained XY single-site Metropolis update.

    A sweep proposes one change at every site.  Changes and Metropolis uniforms
    are generated before entering a kernel, so a fixed RNG seed replays across
    serial and parallel Numba executions.
    """

    def __init__(self, action, interval_phi=None, seed=None, backend="numba",
                 parallel=True):
        if not isinstance(action, XY):
            raise ValueError("Need an XY action")
        if action.Lattice.N < 2:
            raise ValueError("SiteUpdate requires Lattice.N >= 2")
        if backend not in {"numba", "python"}:
            raise ValueError("backend must be 'numba' or 'python'")

        epsilon = float(action.epsilon)
        if interval_phi is None:
            interval_phi = min(epsilon, math.pi)
        interval_phi = float(interval_phi)
        if (not np.isfinite(interval_phi) or interval_phi <= 0.0
                or interval_phi > math.pi):
            raise ValueError("interval_phi must be finite and in (0, pi]")

        self.Action = action
        self.Lattice = action.Lattice
        self.kappa = float(action.kappa)
        self.epsilon = epsilon
        self.interval_phi = interval_phi
        self.backend = backend
        self.parallel = bool(parallel)
        self.rng = np.random.default_rng(seed)

        self.plus, self.minus = _flat_neighbours(self.Lattice)
        self.color_sites, self.offsets = _flat_colours(self.Lattice)

        self.accepted = 0
        self.proposed = 0
        self.sweeps = 0
        # Sum of per-sweep empirical acceptance fractions.  This mirrors the
        # existing generator report convention while needing no extra kernel
        # reduction beyond the required accepted-count return value.
        self.acceptance = 0.0

    def __str__(self):
        return "XY SiteUpdate"

    def _validate_phi(self, phi):
        if not isinstance(phi, Form):
            raise TypeError("phi must be a lattice Form")
        if phi.degree != 0 or not _same_lattice(phi.lattice, self.Lattice):
            raise ValueError("phi must be a 0-form on this action's lattice")
        if phi.shape != (1,) + self.Lattice.dims:
            raise ValueError("phi has an incompatible shape")
        raw = np.asarray(phi)
        if np.any(raw < 0.0) or np.any(raw >= _TWO_PI):
            raise ValueError("phi angles must lie in [0, 2*pi)")
        if not self.Action.valid({"phi": phi}):
            raise ValueError("phi is not a valid constrained XY configuration")

    def _validate_randomness(self, proposals, uniforms):
        volume = self.Lattice.sites
        proposals = np.asarray(proposals, dtype=np.float64)
        uniforms = np.asarray(uniforms, dtype=np.float64)
        if proposals.shape != (volume,) or uniforms.shape != (volume,):
            raise ValueError("proposals and uniforms must each have shape (Lattice.sites,)")
        if not np.all(np.isfinite(proposals)):
            raise ValueError("proposals must be finite")
        if not np.all(np.isfinite(uniforms)) or np.any(uniforms < 0.0) or np.any(uniforms >= 1.0):
            raise ValueError("uniforms must lie in [0, 1)")
        return np.ascontiguousarray(proposals), np.ascontiguousarray(uniforms)

    def sweep(self, phi, proposals, uniforms):
        """Replay one supplied-randomness sweep and return a fresh 0-form."""
        self._validate_phi(phi)
        proposals, uniforms = self._validate_randomness(proposals, uniforms)

        # The action accepts numeric Forms, including integer and strided
        # arrays.  Kernels need a mutable, flat float64 C-order view; make
        # that copy here so neither case can mutate or silently truncate the
        # caller's field.
        updated = Form(
            np.array(np.asarray(phi), dtype=np.float64, copy=True, order="C"),
            degree=0,
            lattice=phi.lattice,
        )
        raw_phi = np.asarray(updated).reshape(-1)
        if self.backend == "python":
            accepted = reference_sweep(
                raw_phi, self.plus, self.minus, self.color_sites, self.offsets,
                proposals, uniforms, self.kappa, self.epsilon,
            )
        else:
            kernel = sweep_parallel if self.parallel else sweep_serial
            accepted = kernel(
                raw_phi, self.plus, self.minus, self.color_sites, self.offsets,
                proposals, uniforms, self.kappa, self.epsilon,
            )

        accepted = int(accepted)
        self.accepted += accepted
        self.proposed += self.Lattice.sites
        self.sweeps += 1
        self.acceptance += accepted / self.Lattice.sites
        return updated

    def step(self, cfg):
        """Generate one reproducible volume sweep and replace only ``phi``."""
        volume = self.Lattice.sites
        proposals = self.rng.uniform(-self.interval_phi, self.interval_phi, volume)
        uniforms = self.rng.uniform(0.0, 1.0, volume)
        return cfg | {"phi": self.sweep(cfg["phi"], proposals, uniforms)}

    def inline_observables(self, steps):
        return {}

    def report(self):
        if not self.proposed:
            return "There were no XY single-site proposals."
        return (
            f"There were {self.accepted} XY single-site proposals accepted of "
            f"{self.proposed} proposed updates.\n"
            f"    {self.accepted / self.proposed:.6f} acceptance rate\n"
            f"    {self.acceptance / self.sweeps:.6f} average sweep acceptance rate."
        )


class ClusterUpdate(ReadWriteable, Generator):
    r"""One admissibility-aware Wolff cluster for the ferromagnetic XY model.

    A random spin-space reflection is independent of the spatial dimension.
    Bonds whose single-end reflection violates ``|d phi| < epsilon`` are
    activated with probability one; other bonds have activation probability
    ``1 - exp(-max(0, 2*kappa*(r.s_x)*(r.s_y)))``.  A completed cluster is
    flipped without a Metropolis accept/reject step.

    Every :meth:`step` performs exactly one cluster, irrespective of its
    size.  ``equivalent_sweeps`` reports accumulated cluster membership over
    volume; it is a work diagnostic, not an instruction to stop on a
    state-dependent cluster size.  Both backends grow clusters serially.

    Only XY actions with ``kappa >= 0`` are supported.  Use :class:`SiteUpdate`
    for negative couplings.  Near strict-cutoff floating-point boundaries,
    a numerically inadmissible reflected field raises ``FloatingPointError``
    before mutation or statistics updates.  No tolerance relaxes the cutoff.
    The numerical guards are not a proof of floating-point detailed balance.
    """

    def __init__(self, action, seed=None, backend="numba"):
        if not isinstance(action, XY):
            raise ValueError("Need an XY action")
        if action.Lattice.N < 2:
            raise ValueError("ClusterUpdate requires Lattice.N >= 2")
        if not np.isfinite(action.kappa) or action.kappa < 0.0:
            raise ValueError("ClusterUpdate requires kappa >= 0; use SiteUpdate for negative kappa")
        if backend not in {"numba", "python"}:
            raise ValueError("backend must be 'numba' or 'python'")

        self.Action = action
        self.Lattice = action.Lattice
        self.kappa = float(action.kappa)
        self.epsilon = float(action.epsilon)
        self.backend = backend
        self.rng = np.random.default_rng(seed)
        self.plus, self.minus = _flat_neighbours(self.Lattice)
        self.last_cluster_size = 0
        self.last_bond_tests = 0
        self.clusters = 0
        self.flipped_sites = 0
        self.bond_tests = 0

    def __str__(self):
        return "XY ClusterUpdate"

    @property
    def equivalent_sweeps(self):
        """Accumulated cluster membership divided by lattice volume."""
        return self.flipped_sites / self.Lattice.sites

    def _validate_phi(self, phi):
        if not isinstance(phi, Form):
            raise TypeError("phi must be a lattice Form")
        if phi.degree != 0 or not _same_lattice(phi.lattice, self.Lattice):
            raise ValueError("phi must be a 0-form on this action's lattice")
        if phi.shape != (1,) + self.Lattice.dims:
            raise ValueError("phi has an incompatible shape")
        if not self.Action.valid({"phi": phi}):
            raise ValueError("phi must be a canonical, admissible constrained XY configuration")

    def _validate_randomness(self, seed_site, reflection_angle, uniforms):
        if (not isinstance(seed_site, numbers.Integral)
                or isinstance(seed_site, (bool, np.bool_))
                or not 0 <= seed_site < self.Lattice.sites):
            raise ValueError("seed_site must be an integer in [0, Lattice.sites)")
        if (not isinstance(reflection_angle, numbers.Real)
                or isinstance(reflection_angle, (bool, np.bool_))
                or not np.isfinite(reflection_angle)):
            raise ValueError("reflection_angle must be a finite real number")
        draws = np.asarray(uniforms)
        if draws.shape != (self.Lattice.D, self.Lattice.sites):
            raise ValueError("uniforms must have shape (Lattice.D, Lattice.sites)")
        if (not np.issubdtype(draws.dtype, np.number) or np.iscomplexobj(draws)
                or not np.all(np.isfinite(draws))
                or np.any(draws < 0.0) or np.any(draws >= 1.0)):
            raise ValueError("uniforms must be finite real numbers in [0, 1)")
        alpha = float(reflection_angle) % _TWO_PI
        if alpha >= _TWO_PI:
            alpha = 0.0
        return int(seed_site), alpha, np.ascontiguousarray(draws, dtype=np.float64)

    def cluster(self, phi, seed_site, reflection_angle, uniforms):
        """Replay one supplied-randomness cluster and return a fresh 0-form.

        ``reflection_angle`` specifies the normal ``r=(cos(alpha),sin(alpha))``.
        The independent uniform for the positive bond ``x -> x+mu`` is
        ``uniforms[mu, x]``, with flat site index ``x``.  Neither the field nor
        the supplied random tape is mutated.
        """
        self._validate_phi(phi)
        seed_site, reflection_angle, uniforms = self._validate_randomness(
            seed_site, reflection_angle, uniforms,
        )
        updated = Form(
            np.array(np.asarray(phi), dtype=np.float64, copy=True, order="C"),
            degree=0,
            lattice=phi.lattice,
        )
        kernel = reference_cluster if self.backend == "python" else cluster_serial
        size, tests = kernel(
            np.asarray(updated).reshape(-1), self.plus, self.minus,
            seed_site, reflection_angle, uniforms, self.kappa, self.epsilon,
        )
        self.last_cluster_size = int(size)
        self.last_bond_tests = int(tests)
        self.clusters += 1
        self.flipped_sites += int(size)
        self.bond_tests += int(tests)
        return updated

    def step(self, cfg):
        """Perform exactly one reproducible cluster, replacing only ``phi``."""
        seed_site = int(self.rng.integers(self.Lattice.sites))
        reflection_angle = self.rng.uniform(0.0, _TWO_PI)
        uniforms = self.rng.uniform(0.0, 1.0, (self.Lattice.D, self.Lattice.sites))
        return cfg | {"phi": self.cluster(cfg["phi"], seed_site, reflection_angle, uniforms)}

    def inline_observables(self, steps):
        return {}

    def report(self):
        if not self.clusters:
            return "There were no XY cluster updates."
        return (
            f"There were {self.clusters} XY cluster updates.\n"
            f"    {self.flipped_sites / self.clusters:.6f} mean cluster size\n"
            f"    {self.equivalent_sweeps:.6f} equivalent sweeps (cluster sites / volume)\n"
            f"    {self.bond_tests} boundary bond tests."
        )


class WindingUpdate(ReadWriteable, Generator):
    r"""Offer one global unit-twist proposal for the constrained XY model.

    This supplements local updates, which can freeze a winding sector.  It is
    not by itself an ergodic sampler.
    """

    def __init__(self, action, seed=None):
        if not isinstance(action, XY):
            raise ValueError("Need an XY action")
        if action.Lattice.N < 2:
            raise ValueError("WindingUpdate requires Lattice.N >= 2")
        self.Action = action
        self.Lattice = action.Lattice
        self.rng = np.random.default_rng(seed)
        self.accepted = 0
        self.proposed = 0
        self.sweeps = 0
        self.acceptance = 0.0

    def __str__(self):
        return "XY WindingUpdate"

    def step(self, cfg):
        phi = cfg["phi"]
        if (not isinstance(phi, Form) or phi.degree != 0
                or not _same_lattice(phi.lattice, self.Lattice)):
            raise ValueError("phi must be a 0-form on this action's lattice")
        if np.any(np.asarray(phi) < 0.0) or np.any(np.asarray(phi) >= _TWO_PI):
            raise ValueError("phi angles must lie in [0, 2*pi)")
        if not self.Action.valid({"phi": phi}):
            raise ValueError("phi is not a valid constrained XY configuration")

        axis = int(self.rng.integers(self.Lattice.D))
        sign = 1 if self.rng.integers(2) else -1
        coordinates = np.arange(self.Lattice.N, dtype=float)
        shape = [1] * self.Lattice.D
        shape[axis] = self.Lattice.N
        twist = sign * _TWO_PI * coordinates.reshape(shape) / self.Lattice.N
        candidate = self.Action.canonical(phi + twist[np.newaxis, ...])

        acceptance = 0.0
        if self.Action.valid({"phi": candidate}):
            delta_s = float(self.Action(candidate) - self.Action(phi))
            acceptance = 1.0 if delta_s <= 0.0 else math.exp(-delta_s)

        accepted = self.rng.uniform(0.0, 1.0) < acceptance
        self.proposed += 1
        self.sweeps += 1
        self.acceptance += acceptance
        if accepted:
            self.accepted += 1
            return cfg | {"phi": candidate}
        return cfg | {"phi": phi.copy()}

    def inline_observables(self, steps):
        return {}

    def report(self):
        if not self.proposed:
            return "There were no XY winding proposals."
        return (
            f"There were {self.accepted} XY winding proposals accepted of "
            f"{self.proposed} proposed updates.\n"
            f"    {self.accepted / self.proposed:.6f} acceptance rate\n"
            f"    {self.acceptance / self.sweeps:.6f} average Metropolis acceptance probability."
        )
