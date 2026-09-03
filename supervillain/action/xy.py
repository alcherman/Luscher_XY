"""The compact XY action with a Luescher admissibility constraint."""

import numbers

import numpy as np

from supervillain.batch import Batch
from supervillain.configurations import Configurations
from supervillain.h5 import ReadWriteable
from supervillain.lattice import Form, Lattice


_TWO_PI = 2.0 * np.pi


class XY(ReadWriteable):
    r"""Compact XY action on a periodic hypercubic lattice.

    The field is an angle at each site and the action is

    .. math::

       S[\varphi] = -\kappa\sum_{x,\mu}
       \cos\!\left(\operatorname{wrap}
       [\varphi(x+\hat\mu)-\varphi(x)]\right).

    Each positive-direction bond is counted once, so the quadratic term is
    ``kappa/2 * sum(theta**2)`` after subtracting the constant cold action.
    The action accepts every finite real representative of a compact angle:
    adding independent multiples of ``2*pi`` at sites leaves it unchanged.
    :meth:`canonical` maps such representatives to the storage convention
    ``[0, 2*pi)``.  :meth:`valid` additionally requires that convention, so a
    configuration accepted for storage is both canonical and admissible.
    Admissibility requires every principal link difference to lie strictly
    inside ``(-epsilon, epsilon)``.  Invalid or non-finite fields have
    infinite action.
    """

    def __init__(self, lattice, kappa, epsilon=np.pi):
        if not isinstance(lattice, Lattice):
            raise TypeError(
                "XY requires a supervillain.lattice.Lattice, got "
                f"{type(lattice).__name__}"
            )
        if (not isinstance(lattice.D, numbers.Integral)
                or isinstance(lattice.D, (bool, np.bool_))
                or lattice.D < 1):
            raise ValueError("XY requires a lattice with integer D >= 1.")
        if (not isinstance(lattice.N, numbers.Integral)
                or isinstance(lattice.N, (bool, np.bool_))
                or lattice.N < 2):
            raise ValueError("XY requires a lattice with integer N >= 2.")
        if (not isinstance(kappa, numbers.Real)
                or isinstance(kappa, (bool, np.bool_))
                or not np.isfinite(kappa)):
            raise ValueError("kappa must be a finite real number.")
        if (not isinstance(epsilon, numbers.Real)
                or isinstance(epsilon, (bool, np.bool_))
                or not np.isfinite(epsilon)
                or epsilon <= 0.0 or epsilon > np.pi):
            raise ValueError("epsilon must be finite and satisfy 0 < epsilon <= pi.")

        self.Lattice = lattice
        self.kappa = float(kappa)
        self.epsilon = float(epsilon)

    def __str__(self):
        return f"XY({self.Lattice}, kappa={self.kappa}, epsilon={self.epsilon})"

    def _phi_form(self, phi):
        """Validate a scalar field and return it as a degree-zero Form."""
        try:
            array = np.asarray(phi)
        except Exception as exc:
            raise TypeError("phi must be convertible to a real-valued ndarray.") from exc

        expected_spatial = self.Lattice.dims
        expected_form = (1,) + expected_spatial
        if array.shape == expected_spatial:
            array = array[np.newaxis, ...]
        elif array.shape != expected_form:
            raise ValueError(
                f"phi must have shape {expected_spatial} or {expected_form}; "
                f"got {array.shape}."
            )

        if not np.issubdtype(array.dtype, np.number):
            raise TypeError("phi must have a numeric dtype.")
        if np.iscomplexobj(array):
            raise TypeError("phi must be real-valued.")
        if not np.isfinite(array).all():
            raise ValueError("phi must contain only finite values.")
        return Form(array, degree=0, lattice=self.Lattice)

    @staticmethod
    def _wrap(array):
        """Principal branch in [-pi, pi), without losing small differences."""
        raw = np.asarray(array)
        # ``(raw + pi) % (2*pi) - pi`` rounds raw+pi back to pi for
        # sufficiently small raw.  Preserving values already on the principal
        # branch keeps a strict cutoff meaningful even for tiny epsilon.
        wrapped = np.remainder(raw, _TWO_PI)
        wrapped = np.where(wrapped >= np.pi, wrapped - _TWO_PI, wrapped)
        return np.where((raw >= -np.pi) & (raw < np.pi), raw, wrapped)

    def canonical(self, phi):
        r"""Return a degree-zero copy whose angles are in ``[0, 2*pi)``."""
        field = self._phi_form(phi)
        # Work in float64.  For a tiny negative input, modulo can round to
        # exactly 2*pi; 0 is then the nearest representable point on the
        # circle and, unlike nextafter(2*pi, 0), creates no spurious ULP-sized
        # link across the branch cut.
        wrapped = np.mod(np.asarray(field, dtype=float), _TWO_PI)
        wrapped[wrapped >= _TWO_PI] = 0.0
        return Form(wrapped,
                    degree=0, lattice=self.Lattice)

    def links(self, phi):
        r"""Principal forward differences as a degree-one :class:`Form`."""
        field = self._phi_form(phi)
        # Convert before subtraction: unsigned angles otherwise underflow
        # on descending links, and narrow integer arithmetic can overflow.
        spatial = np.asarray(field, dtype=np.float64)[0]
        links = np.empty((self.Lattice.D,) + self.Lattice.dims, dtype=float)
        for direction in range(self.Lattice.D):
            links[direction] = self._wrap(
                np.roll(spatial, shift=-1, axis=direction) - spatial
            )
        return Form(links, degree=1, lattice=self.Lattice)

    def local(self, phi):
        r"""Per-link action contributions, with ``inf`` on inadmissible links."""
        links = self.links(phi)
        result = -self.kappa * np.cos(np.asarray(links))
        result[np.abs(np.asarray(links)) >= self.epsilon] = np.inf
        return Form(result, degree=1, lattice=self.Lattice)

    def __call__(self, phi, **kwargs):
        r"""Return the total action, or ``inf`` for an invalid configuration."""
        try:
            local = self.local(phi)
        except (TypeError, ValueError):
            return float("inf")
        if not np.isfinite(local).all():
            return float("inf")
        return float(np.sum(local))

    def valid(self, configuration):
        r"""Whether ``configuration['phi']`` is canonical and admissible.

        Unlike :meth:`__call__`, this method enforces the stored-angle
        convention ``0 <= phi < 2*pi`` as well as the strict link cutoff.
        Use :meth:`canonical` before placing an arbitrary finite angle field
        into a configuration.
        """
        try:
            if not isinstance(configuration, dict):
                return False
            field = self._phi_form(configuration["phi"])
            if (np.any(np.asarray(field) < 0.0)
                    or np.any(np.asarray(field) >= _TWO_PI)):
                return False
            links = self.links(field)
            return bool(np.isfinite(links).all()
                        and np.all(np.abs(np.asarray(links)) < self.epsilon))
        except (KeyError, TypeError, ValueError):
            return False

    def configurations(self, count):
        r"""Allocate ``count`` zeroed compact-angle configurations."""
        if (not isinstance(count, numbers.Integral)
                or isinstance(count, (bool, np.bool_)) or count < 0):
            raise ValueError("count must be a nonnegative integer.")
        return Configurations({
            "phi": Batch(int(count), cls=Form, degree=0, lattice=self.Lattice),
        })
