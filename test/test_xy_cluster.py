"""Independent graph, detailed-balance, and API checks for XY Wolff updates."""

import itertools
import math

import numpy as np
import pytest

from supervillain import Lattice, XY
from supervillain.action.villain import Villain
from supervillain.generator.xy import ClusterUpdate
from supervillain.generator.reference_implementation.xy import reference_action
from supervillain.generator.reference_implementation.xy_cluster import reference_cluster
from supervillain.lattice import Form


TWO_PI = 2.0 * math.pi


def _principal(delta):
    # Preserve tiny principal differences rather than adding pi first.
    delta = np.asarray(delta)
    reduced = np.remainder(delta, TWO_PI)
    return np.where((-math.pi <= delta) & (delta < math.pi), delta,
                    np.where(reduced >= math.pi, reduced - TWO_PI, reduced))


def _edge_probabilities(field, plus, alpha, kappa, epsilon):
    """Compute probabilities independently from the original pair weights."""
    reflected = (2 * alpha + math.pi - field) % TWO_PI
    probabilities = np.empty_like(plus, dtype=float)
    for mu in range(plus.shape[0]):
        for site, neighbor in enumerate(plus[mu]):
            new_delta = float(_principal(reflected[site] - field[neighbor]))
            if abs(new_delta) >= epsilon:
                probabilities[mu, site] = 1.0
            else:
                old_dot = math.cos(field[site] - field[neighbor])
                new_dot = math.cos(new_delta)
                probabilities[mu, site] = -math.expm1(-max(0, kappa * (old_dot-new_dot)))
    return probabilities


@pytest.mark.parametrize("D,N", [(1, 2), (1, 5), (2, 2), (2, 3), (3, 3), (4, 2)])
@pytest.mark.parametrize("epsilon", [0.7, 2.3, math.pi])
def test_independent_full_graph_oracle_matches_numba(D, N, epsilon):
    action = XY(Lattice(D, N), 0.81, epsilon)
    update = ClusterUpdate(action)
    rng = np.random.default_rng(871 + D + N)
    phi = action.canonical(rng.uniform(0.0, 0.25, action.Lattice.dims))
    for _ in range(8):
        seed_site = int(rng.integers(action.Lattice.sites))
        alpha = rng.uniform(0, TWO_PI)
        uniforms = rng.uniform(size=(D, action.Lattice.sites))
        oracle = np.asarray(phi).copy().reshape(-1)
        count, tests = reference_cluster(oracle, update.plus, update.minus,
                                        seed_site, alpha, uniforms, action.kappa, epsilon)
        result = update.cluster(phi, seed_site, alpha, uniforms)
        assert np.array_equal(np.asarray(result).reshape(-1), oracle)
        assert update.last_cluster_size == count
        assert update.last_bond_tests == tests
        assert action.valid({"phi": result})
        phi = result


@pytest.mark.parametrize("backend", ["python", "numba"])
def test_parallel_length_two_bonds_have_independent_activation(backend):
    action = XY(Lattice(1, 2), 0.5)
    updater = ClusterUpdate(action, backend=backend)
    field = action.canonical([0.0, 0.0])
    # Each parallel bond activates with p=1-exp(-2*kappa)=0.632...
    # Reject the first, accept the second: deduplicating neighbors is wrong.
    result = updater.cluster(field, 0, 0.3, np.array([[0.999, 0.0]]))
    assert updater.last_cluster_size == 2
    assert updater.last_bond_tests == 2
    assert np.allclose(result, math.pi + 0.6)


@pytest.mark.parametrize("epsilon", [math.pi/2, math.pi])
@pytest.mark.parametrize("backend", ["python", "numba"])
def test_kappa_zero_still_activates_strict_hard_bonds(epsilon, backend):
    action = XY(Lattice(2, 3), 0.0, epsilon)
    update = ClusterUpdate(action, backend=backend)
    phi = action.canonical(np.zeros(action.Lattice.dims))
    # Reflection of a zero spin is exactly epsilon.  Equality is forbidden,
    # including the measure-zero antipodal exception at epsilon=pi.
    alpha = (epsilon - math.pi)/2
    result = update.cluster(phi, 0, alpha, np.full((2, 9), 0.999))
    assert update.last_cluster_size == 9
    assert action.valid({"phi": result})
    assert np.allclose(result, epsilon)


@pytest.mark.parametrize("backend", ["python", "numba"])
def test_zero_coupling_without_forbidden_boundary_is_single_site(backend):
    action = XY(Lattice(3, 3), 0.0)
    update = ClusterUpdate(action, backend=backend)
    phi = action.canonical(np.full(action.Lattice.dims, 0.3))
    result = update.cluster(phi, 4, 0.7, np.zeros((3, 27)))
    assert update.last_cluster_size == 1
    assert np.count_nonzero(np.asarray(result) != np.asarray(phi)) == 1


@pytest.mark.parametrize("D,N,kappa,epsilon", [
    (1, 2, 0.7, math.pi),
    (1, 3, 0.8, 1.6),
    (2, 2, 0.6, math.pi),
    (2, 2, 0.9, 1.8),
    (2, 2, 0.0, 1.8),
])
def test_exact_embedded_state_transition_matrix_obeys_detailed_balance(D, N, kappa, epsilon):
    """Enumerate every bond outcome and seed, without stochastic tolerances.

    A fixed reflection leaves a finite orbit of 2**V Ising-like states.
    Enumerating its transition matrix tests the complete transition law,
    including hard constraints, seed-size factors, and duplicate N=2 links.
    Target weights come from the independent full cosine action.
    """
    action = XY(Lattice(D, N), kappa, epsilon)
    update = ClusterUpdate(action)
    volume = action.Lattice.sites
    alpha = 0.0
    base = np.linspace(0.3, 1.4, volume)
    reflected_base = math.pi - base
    all_bits = np.array(list(itertools.product([False, True], repeat=volume)))
    all_fields = np.where(all_bits, reflected_base, base)
    energies = np.array([reference_action(field.reshape(action.Lattice.dims), kappa, epsilon)
                         for field in all_fields])
    allowed = np.flatnonzero(np.isfinite(energies))
    index = {int(state): position for position, state in enumerate(allowed)}
    transition = np.zeros((len(allowed), len(allowed)))
    state_ids = 2 ** np.arange(volume-1, -1, -1)

    for state in allowed:
        field = all_fields[state]
        probabilities = _edge_probabilities(field, update.plus, alpha, kappa, epsilon)
        flat = probabilities.ravel()
        variable = np.flatnonzero((flat > 0) & (flat < 1))
        for choices in itertools.product([False, True], repeat=len(variable)):
            active = flat == 1
            active[variable] = choices
            probability = np.prod(np.where(active, flat, 1-flat))
            # Midpoints of the correct activation/rejection intervals avoid
            # artificial ulp-sensitive comparisons of p between formulas.
            uniforms = np.where(active, flat/2, (1+flat)/2).reshape(probabilities.shape)
            uniforms = np.minimum(uniforms, np.nextafter(1.0, 0.0))
            for seed_site in range(volume):
                result = update.cluster(action.canonical(field.reshape(action.Lattice.dims)),
                                        seed_site, alpha, uniforms)
                values = np.asarray(result).reshape(-1)
                bits = np.abs(values-reflected_base) < np.abs(values-base)
                target = int(bits @ state_ids)
                assert target in index, "Cluster escaped admissible reflection orbit"
                assert np.allclose(values, all_fields[target], atol=2e-15, rtol=0)
                transition[index[state], index[target]] += probability/volume

    weights = np.exp(-energies[allowed])
    weights /= weights.sum()
    flux = weights[:, None] * transition
    assert np.allclose(transition.sum(axis=1), 1, atol=2e-13, rtol=0)
    assert np.allclose(flux, flux.T, atol=2e-13, rtol=0)
    assert np.allclose(weights @ transition, weights, atol=2e-13, rtol=0)


def test_cluster_can_remove_global_winding_without_vortices():
    action = XY(Lattice(2, 8), 0.0, math.pi/2)
    update = ClusterUpdate(action)
    spatial = np.broadcast_to(np.arange(8)[:, None] * TWO_PI/8, (8, 8))
    phi = action.canonical(spatial)
    assert action.valid({"phi": phi})
    old_winding = _principal(np.roll(spatial, -1, axis=0)-spatial).sum(axis=0)/TWO_PI
    result = update.cluster(phi, 0, 0.0, np.full((2, 64), 0.5))
    after = np.asarray(result)[0]
    new_winding = _principal(np.roll(after, -1, axis=0)-after).sum(axis=0)/TWO_PI
    assert np.allclose(old_winding, 1)
    assert np.allclose(new_winding, 0)
    assert action.valid({"phi": result})
    assert 1 < update.last_cluster_size < action.Lattice.sites


def test_seeded_backends_replay_and_preserve_input_configuration():
    action = XY(Lattice(2, 3), 0.92, 2.1)
    python = ClusterUpdate(action, seed=913, backend="python")
    numba = ClusterUpdate(action, seed=913, backend="numba")
    phi = action.canonical(np.full((3, 3), 0.8))
    tag = object()
    a = b = {"phi": phi, "metadata": tag}
    for _ in range(25):
        snapshot = np.asarray(a["phi"]).copy()
        previous = a
        a, b = python.step(a), numba.step(b)
        assert a is not previous
        assert a["phi"] is not previous["phi"]
        assert np.array_equal(previous["phi"], snapshot)
        assert a["metadata"] is tag
        assert np.array_equal(a["phi"], b["phi"])
        assert python.last_cluster_size == numba.last_cluster_size
        assert python.last_bond_tests == numba.last_bond_tests
    assert python.clusters == numba.clusters == 25
    assert python.flipped_sites == numba.flipped_sites
    assert python.bond_tests == numba.bond_tests
    assert python.equivalent_sweeps == pytest.approx(python.flipped_sites/9)
    assert python.inline_observables(3) == {}
    assert isinstance(python.report(), str)


def test_only_ferromagnetic_xy_is_accepted():
    with pytest.raises((ValueError, TypeError), match="XY"):
        ClusterUpdate(Villain(Lattice(2, 3), 1.0))
    with pytest.raises((ValueError, TypeError)):
        ClusterUpdate(XY(Lattice(2, 3), -0.1))
    with pytest.raises((ValueError, TypeError)):
        ClusterUpdate(XY(Lattice(2, 3), 1.0), backend="gpu")


@pytest.mark.parametrize("bad_seed", [-1, 9, 0.5, True, np.nan])
def test_invalid_seed_site_rejected_without_counters_or_mutation(bad_seed):
    action = XY(Lattice(2, 3), 1.0)
    update = ClusterUpdate(action)
    phi = action.canonical(np.zeros((3, 3)))
    with pytest.raises((ValueError, TypeError)):
        update.cluster(phi, bad_seed, 0.0, np.full((2, 9), 0.5))
    assert np.all(np.asarray(phi) == 0)
    assert update.clusters == 0


@pytest.mark.parametrize("alpha", [np.nan, np.inf, -np.inf, True, 1j])
def test_invalid_reflection_angle_rejected(alpha):
    action = XY(Lattice(1, 3), 1.0)
    update = ClusterUpdate(action)
    with pytest.raises((ValueError, TypeError)):
        update.cluster(action.canonical([0, 0, 0]), 0, alpha, np.full((1, 3), 0.5))


@pytest.mark.parametrize("uniforms", [np.ones((1, 3)), np.full((1, 3), -0.1),
                                    np.full((1, 3), np.nan), np.zeros(3),
                                    np.zeros((1, 3), dtype=complex)])
def test_invalid_uniform_tape_rejected(uniforms):
    action = XY(Lattice(1, 3), 1.0)
    with pytest.raises((ValueError, TypeError)):
        ClusterUpdate(action).cluster(action.canonical([0, 0, 0]), 0, 0.0, uniforms)


def test_integer_and_strided_forms_are_copied_to_float64():
    action = XY(Lattice(2, 3), 0.0)
    update = ClusterUpdate(action)
    for raw in [np.zeros((1, 3, 3), dtype=int), np.zeros((1, 6, 3))[:, ::2, :]]:
        phi = Form(raw, degree=0, lattice=action.Lattice)
        result = update.cluster(phi, 0, 0.2, np.full((2, 9), 0.5))
        assert result.dtype == np.float64
        assert np.all(raw == 0)
        assert not np.shares_memory(result, raw)


@pytest.mark.parametrize("backend", ["python", "numba"])
@pytest.mark.parametrize("epsilon", [np.nextafter(0.0, 1.0), 1e-16, 1e-12])
def test_tiny_admissibility_stays_strict_for_cold_cluster(epsilon, backend):
    action = XY(Lattice(2, 3), 0.0, epsilon)
    update = ClusterUpdate(action, backend=backend)
    phi = action.canonical(np.zeros((3, 3)))
    result = update.cluster(phi, 0, 0.37, np.full((2, 9), 0.5))
    assert update.last_cluster_size == 9
    assert action.valid({"phi": result})


def test_invalid_configuration_is_rejected_before_cluster_growth():
    action = XY(Lattice(1, 3), 0.4, 0.5)
    update = ClusterUpdate(action)
    invalid = action.canonical([0.0, 0.5, 0.0])
    with pytest.raises((ValueError, TypeError)):
        update.cluster(invalid, 0, 0.0, np.full((1, 3), 0.5))
    with pytest.raises((ValueError, TypeError)):
        update.cluster(np.zeros((1, 3)), 0, 0.0, np.full((1, 3), 0.5))
    with pytest.raises((ValueError, TypeError)):
        update.cluster(Lattice(1, 4).zeros(0), 0, 0.0, np.full((1, 3), 0.5))


@pytest.mark.parametrize("backend", ["python", "numba"])
def test_reflection_roundoff_failure_is_explicit_and_transactional(backend):
    action = XY(Lattice(1, 3), 0.0, 8e-16)
    update = ClusterUpdate(action, backend=backend)
    phi = action.canonical([0.0, 7e-16, 0.0])
    before = np.asarray(phi).copy()
    assert action.valid({"phi": phi})
    # Reflection moves these angles from near zero to near four, where
    # representable spacing rounds the 7e-16 link up to 8.88e-16.
    with pytest.raises(FloatingPointError):
        update.cluster(phi, 0, 0.37, np.full((1, 3), 0.5))
    assert np.array_equal(phi, before)
    assert update.clusters == update.flipped_sites == update.bond_tests == 0


@pytest.mark.parametrize("backend", ["python", "numba"])
def test_extreme_finite_kappa_activates_ferromagnetic_bonds(backend):
    action = XY(Lattice(1, 3), np.finfo(float).max)
    update = ClusterUpdate(action, backend=backend)
    phi = action.canonical([0.0, 0.0, 0.0])
    with np.errstate(over="ignore"):
        result = update.cluster(phi, 0, 0.2, np.full((1, 3), 0.999))
    assert update.last_cluster_size == 3
    assert action.valid({"phi": result})


def test_random_tape_is_never_modified():
    action = XY(Lattice(2, 3), 1.1, 2.0)
    phi = action.canonical(np.zeros((3, 3)))
    tape = np.random.default_rng(39).uniform(size=(2, 9))
    before = tape.copy()
    for backend in ["python", "numba"]:
        ClusterUpdate(action, backend=backend).cluster(phi, 0, 0.1, tape)
        assert np.array_equal(tape, before)


@pytest.mark.parametrize("backend", ["python", "numba"])
def test_large_kappa_resolves_small_positive_bond_cost(backend):
    # Direct cos(0)-cos(2e-9) rounds to zero; multiplied by kappa this
    # incorrectly erases an order-one exponent. Its stable value is ~2.
    action = XY(Lattice(1, 3), 1e18, math.pi)
    update = ClusterUpdate(action, backend=backend)
    phi = action.canonical([0.0, 0.0, 0.0])
    update.cluster(phi, 0, math.pi/2 + 1e-9, np.full((1, 3), 0.5))
    assert update.last_cluster_size == 3


@pytest.mark.parametrize("backend", ["python", "numba"])
@pytest.mark.parametrize("dtype", [np.uint8, np.uint64, np.int8, np.int64])
def test_nonuniform_integer_angles_use_signed_float_link_differences(backend, dtype):
    action = XY(Lattice(1, 3), 1.0, 1.2)
    raw = np.array([[0, 1, 0]], dtype=dtype)
    phi = Form(raw, degree=0, lattice=action.Lattice)
    assert np.array_equal(action.links(phi), [[1.0, -1.0, 0.0]])
    assert action(phi) == pytest.approx(-(1 + 2 * math.cos(1)))
    assert action.valid({"phi": phi})
    tape = np.full((1, 3), 0.5)
    integer = ClusterUpdate(action, backend=backend)
    floating = ClusterUpdate(action, backend=backend)
    result = integer.cluster(phi, 0, 0.3, tape)
    expected = floating.cluster(action.canonical(phi), 0, 0.3, tape)
    assert np.array_equal(result, expected)
    assert np.array_equal(phi, raw)
    assert integer.last_cluster_size == floating.last_cluster_size
