#!/usr/bin/env python
r"""Independent acceptance, geometry, and serialization tests for compact XY.

The action oracle below deliberately uses ``np.roll`` on the spatial array.
It therefore does not share the action's Form/link-table implementation.
"""

import h5py as h5
import numba
import numpy as np
import pytest

import supervillain
from supervillain.action.xy import XY
from supervillain.generator.xy import SiteUpdate
from supervillain.generator._xy_kernels import sweep_parallel, sweep_serial
from supervillain.generator.reference_implementation.xy import (
    reference_action,
    reference_sweep,
)


TWO_PI = 2.0 * np.pi


def _spatial(phi):
    """Return a scalar field's spatial data, accepting a Form or raw array."""
    value = np.asarray(phi)
    return value[0] if value.ndim and value.shape[0] == 1 else value


def _principal(delta):
    """Principal representative in [-pi, pi); endpoint sign is immaterial to cos."""
    return (delta + np.pi) % TWO_PI - np.pi


def full_action(phi, kappa, epsilon=np.pi):
    """Slow, implementation-independent periodic positive-link oracle."""
    data = _spatial(phi)
    links = np.stack([_principal(np.roll(data, -1, axis=mu) - data)
                      for mu in range(data.ndim)])
    return -kappa * np.cos(links).sum()


def make_phi(lattice, seed=0):
    rng = np.random.default_rng(seed)
    return lattice.form(0) + rng.uniform(0.0, TWO_PI, size=lattice.dims)


@pytest.mark.parametrize('D,N', [(1, 2), (1, 5), (2, 3), (3, 3), (4, 2)])
def test_action_matches_independent_positive_link_oracle(D, N):
    lattice = supervillain.lattice.Lattice(D=D, N=N)
    action = XY(lattice, kappa=0.73)
    phi = make_phi(lattice, seed=100 * D + N)

    expected = full_action(phi, action.kappa)
    assert action(phi) == pytest.approx(expected, abs=2e-13)
    assert reference_action(_spatial(phi), action.kappa, action.epsilon) == pytest.approx(expected)
    assert np.asarray(action.links(phi)).shape == (D,) + lattice.dims
    assert action.local(phi).sum() == pytest.approx(expected, abs=2e-13)


@pytest.mark.parametrize('D,N', [(1, 2), (2, 3), (4, 2)])
def test_cold_normalization_and_canonical_form_shape(D, N):
    lattice = supervillain.lattice.Lattice(D=D, N=N)
    kappa = 1.25
    action = XY(lattice, kappa=kappa)
    cold_spatial = np.zeros(lattice.dims)

    canonical = action.canonical(cold_spatial)
    assert isinstance(canonical, supervillain.lattice.Form)
    assert canonical.degree == 0
    assert canonical.lattice is lattice
    assert canonical.shape == (1,) + lattice.dims
    assert action(canonical) == pytest.approx(-kappa * D * N**D)
    assert action.valid({'phi': canonical})


def test_n2_counts_both_positive_oriented_links():
    lattice = supervillain.lattice.Lattice(D=1, N=2)
    kappa, theta = 0.87, 0.41
    action = XY(lattice, kappa=kappa)
    phi = action.canonical(np.array([0.0, theta]))

    # x=0 -> 1 and x=1 -> 0 are both positive-oriented periodic links.
    assert action(phi) == pytest.approx(-2 * kappa * np.cos(theta))


def test_small_gradient_quadratic_coefficient_is_kappa_over_two():
    lattice = supervillain.lattice.Lattice(D=2, N=3)
    action = XY(lattice, kappa=1.3, epsilon=0.5)
    cold = lattice.zeros(0)
    perturbed = cold.copy()
    perturbed[0, 0, 0] = 1e-3
    theta_squared = np.asarray(action.links(perturbed)) ** 2
    coefficient = (action(perturbed) - action(cold)) / theta_squared.sum()
    assert coefficient == pytest.approx(action.kappa / 2, rel=2e-7)


def test_action_is_invariant_under_independent_site_2pi_lifts():
    lattice = supervillain.lattice.Lattice(D=3, N=3)
    action = XY(lattice, kappa=0.61, epsilon=2.7)
    phi = make_phi(lattice, seed=4)
    lifts = np.random.default_rng(5).integers(-3, 4, size=lattice.dims)
    lifted = _spatial(phi) + TWO_PI * lifts

    assert action.canonical(lifted).shape == phi.shape
    assert action(lifted) == pytest.approx(action(phi), abs=2e-13)
    assert np.allclose(action.canonical(lifted), action.canonical(phi), atol=2e-13)
    assert action.valid({'phi': lifted}) == action.valid({'phi': phi})


def test_validity_uses_strict_principal_link_cutoff_and_spin_domain():
    lattice = supervillain.lattice.Lattice(D=1, N=3)

    # The antipodal bond is invalid at epsilon=pi because the cutoff is strict.
    pi_action = XY(lattice, kappa=1.0, epsilon=np.pi)
    antipodal = pi_action.canonical(np.array([0.0, np.pi, 0.0]))
    assert not pi_action.valid({'phi': antipodal})

    epsilon = 0.7
    action = XY(lattice, kappa=1.0, epsilon=epsilon)
    exactly_at_cut = action.canonical(np.array([0.0, epsilon, 0.0]))
    safely_inside = action.canonical(np.array([0.0, epsilon - 1e-10, 0.0]))
    assert not action.valid({'phi': exactly_at_cut})
    assert action.valid({'phi': safely_inside})
    assert not action.valid({'phi': np.array([[0.0, TWO_PI, 0.0]])})
    assert not action.valid({'phi': np.array([[0.0, -1e-12, 0.0]])})


def test_tiny_epsilon_rejects_a_link_above_its_strict_cutoff():
    lattice = supervillain.lattice.Lattice(D=1, N=3)
    action = XY(lattice, kappa=1.0, epsilon=1e-16)
    phi = action.canonical(np.array([0.0, 2e-16, 0.0]))
    assert not action.valid({'phi': phi})


def _plaquette_circulation(phi):
    """Integer vorticity charge from principal link differences, independently."""
    data = _spatial(phi)
    d0 = _principal(np.roll(data, -1, axis=0) - data)
    d1 = _principal(np.roll(data, -1, axis=1) - data)
    return (d0 + np.roll(d1, -1, axis=0)
            - np.roll(d0, -1, axis=1) - d1) / TWO_PI


def test_vortex_requires_strictly_more_than_pi_over_two_link_cutoff():
    lattice = supervillain.lattice.Lattice(D=2, N=2)
    # Around one plaquette this increases by pi/2 on every directed edge.
    # Periodicity supplies the compensating antivortex elsewhere on the torus.
    phi = np.array([[0.0, np.pi / 2], [3 * np.pi / 2, np.pi]])
    charge = _plaquette_circulation(phi)
    assert set(np.rint(charge).astype(int).flat) == {-1, 1}

    above = XY(lattice, kappa=1.0, epsilon=np.pi / 2 + 1e-12)
    at = XY(lattice, kappa=1.0, epsilon=np.pi / 2)
    below = XY(lattice, kappa=1.0, epsilon=np.pi / 2 - 1e-12)
    assert above.valid({'phi': above.canonical(phi)})
    assert not at.valid({'phi': at.canonical(phi)})
    assert not below.valid({'phi': below.canonical(phi)})


@pytest.mark.parametrize('kappa', [np.nan, np.inf, -np.inf, True])
def test_invalid_kappa_is_rejected(kappa):
    lattice = supervillain.lattice.Lattice(D=2, N=3)
    with pytest.raises((TypeError, ValueError)):
        XY(lattice, kappa=kappa)


@pytest.mark.parametrize('epsilon', [0.0, -0.1, np.pi + 1e-12, np.nan, np.inf])
def test_invalid_epsilon_is_rejected(epsilon):
    lattice = supervillain.lattice.Lattice(D=2, N=3)
    with pytest.raises((TypeError, ValueError)):
        XY(lattice, kappa=1.0, epsilon=epsilon)


def test_invalid_lattice_and_configuration_shape_are_rejected():
    with pytest.raises(TypeError):
        XY(object(), kappa=1.0)

    action = XY(supervillain.lattice.Lattice(D=2, N=3), kappa=1.0)
    assert not action.valid({})
    assert not action.valid({'phi': np.zeros((2, 3, 3))})
    assert action.valid({'phi': np.zeros((3, 3))})


def _assert_color_partition_is_independent(update):
    lattice = update.Action.Lattice
    flat_colors = np.asarray(update.color_sites)
    offsets = np.asarray(update.offsets)
    volume = lattice.N**lattice.D
    assert offsets[0] == 0
    assert offsets[-1] == volume
    assert np.array_equal(np.sort(flat_colors), np.arange(volume))

    owner = np.empty(volume, dtype=int)
    for color, (lo, hi) in enumerate(zip(offsets[:-1], offsets[1:])):
        owner[flat_colors[lo:hi]] = color
    for site in range(volume):
        coords = np.unravel_index(site, lattice.dims)
        for axis in range(lattice.D):
            neighbor = list(coords)
            neighbor[axis] = (neighbor[axis] + 1) % lattice.N
            assert owner[site] != owner[np.ravel_multi_index(tuple(neighbor), lattice.dims)]


@pytest.mark.parametrize('D,N', [(1, 3), (2, 3), (3, 3), (2, 4)])
def test_site_coloring_covers_sites_without_periodic_neighbors(D, N):
    update = SiteUpdate(XY(supervillain.lattice.Lattice(D=D, N=N), 0.6), seed=7)
    _assert_color_partition_is_independent(update)


@pytest.mark.parametrize('D,N', [(1, 3), (2, 3), (3, 2)])
def test_fixed_random_tape_matches_reference_and_both_numba_kernels(D, N):
    lattice = supervillain.lattice.Lattice(D=D, N=N)
    action = XY(lattice, kappa=0.43, epsilon=2.2)
    update = SiteUpdate(action, interval_phi=0.9, seed=8, backend='numba', parallel=True)
    rng = np.random.default_rng(9)
    # The kernel contract starts from an admissible field.  A near-cold field
    # keeps this precondition true while the proposals still exercise rejects.
    phi = action.canonical(rng.uniform(0.0, 0.1, lattice.dims))
    proposals = rng.uniform(-update.interval_phi, update.interval_phi, lattice.N**lattice.D)
    uniforms = rng.uniform(0.0, 1.0, lattice.N**lattice.D)

    ref_phi = phi.copy()
    ref_accepted = reference_sweep(
        ref_phi.reshape(-1), update.plus, update.minus, update.color_sites, update.offsets,
        proposals, uniforms, action.kappa, action.epsilon)
    serial = phi.copy()
    serial_accepted = sweep_serial(
        serial.reshape(-1), update.plus, update.minus, update.color_sites, update.offsets,
        proposals, uniforms, action.kappa, action.epsilon)
    parallel = phi.copy()
    parallel_accepted = sweep_parallel(
        parallel.reshape(-1), update.plus, update.minus, update.color_sites, update.offsets,
        proposals, uniforms, action.kappa, action.epsilon)

    assert np.array_equal(serial, ref_phi)
    assert np.array_equal(parallel, ref_phi)
    assert serial_accepted == parallel_accepted == ref_accepted
    assert np.all((0.0 <= serial) & (serial < TWO_PI))
    assert action.valid({'phi': serial})


def test_step_copies_input_and_uses_generator_random_tape():
    lattice = supervillain.lattice.Lattice(D=2, N=3)
    action = XY(lattice, kappa=0.5)
    update = SiteUpdate(action, interval_phi=1.0, seed=12)
    phi = make_phi(lattice, seed=13)
    before = phi.copy()
    cfg = {'phi': phi, 'tag': 'kept'}

    stepped = update.step(cfg)
    assert stepped is not cfg
    assert stepped['phi'] is not cfg['phi']
    assert np.array_equal(phi, before)
    assert stepped['tag'] == 'kept'
    assert action.valid(stepped)


def test_sweep_accepts_integer_and_strided_forms_without_mutating_them():
    lattice = supervillain.lattice.Lattice(D=1, N=3)
    action = XY(lattice, kappa=0.5, epsilon=1.0)
    update = SiteUpdate(action, seed=14, parallel=False)
    zero_randomness = np.zeros(lattice.sites)

    integer = lattice.form(0, dtype=int)
    integer_before = integer.copy()
    integer_updated = update.sweep(integer, zero_randomness, zero_randomness)
    assert integer_updated.dtype == np.float64
    assert np.array_equal(integer, integer_before)
    assert action.valid({'phi': integer_updated})

    strided = action.canonical(np.zeros(lattice.dims))[:, ::-1]
    assert not strided.flags.c_contiguous
    strided_before = strided.copy()
    strided_updated = update.sweep(strided, zero_randomness, zero_randomness)
    assert strided_updated.flags.c_contiguous
    assert np.array_equal(strided, strided_before)
    assert action.valid({'phi': strided_updated})


def test_two_bond_metropolis_threshold_uses_the_full_action_delta():
    lattice = supervillain.lattice.Lattice(D=1, N=2)
    action = XY(lattice, kappa=1.3)
    update = SiteUpdate(action, interval_phi=1.0, seed=0, parallel=False)
    initial = action.canonical(np.array([0.0, 0.1]))
    proposals = np.array([0.6, 0.0])
    candidate = action.canonical(np.array([0.6, 0.1]))
    delta_s = action(candidate) - action(initial)
    assert delta_s == pytest.approx(2 * action.kappa * (np.cos(0.1) - np.cos(0.5)))
    assert delta_s > 0.0
    probability = np.exp(-delta_s)

    accepted = initial.copy()
    accepted_count = sweep_serial(
        accepted.reshape(-1), update.plus, update.minus, update.color_sites, update.offsets,
        proposals, np.array([0.5 * probability, 0.0]), action.kappa, action.epsilon)
    rejected = initial.copy()
    rejected_count = sweep_serial(
        rejected.reshape(-1), update.plus, update.minus, update.color_sites, update.offsets,
        proposals, np.array([(1.0 + probability) / 2.0, 0.0]), action.kappa, action.epsilon)
    # The zero proposal at the other site is always accepted; these counts
    # isolate the one uphill proposal at site zero.
    assert accepted_count == 2
    assert rejected_count == 1
    assert np.array_equal(accepted, candidate)
    assert np.array_equal(rejected, initial)


def _exact_relative_cosine(kappa, epsilon):
    theta = np.linspace(-epsilon, epsilon, 200_001)
    weight = np.exp(2 * kappa * np.cos(theta))
    return np.trapezoid(np.cos(theta) * weight, theta) / np.trapezoid(weight, theta)


@pytest.mark.parametrize('kappa,epsilon,expected', [
    (0.0, 0.8, np.sin(0.8) / 0.8),
    (1.7, 1.8, None),
])
def test_n2_relative_angle_matches_the_constrained_exact_integral(kappa, epsilon, expected):
    """A bounded sampler check: 2 sites have weight exp(2*kappa cos(theta))."""
    lattice = supervillain.lattice.Lattice(D=1, N=2)
    action = XY(lattice, kappa=kappa, epsilon=epsilon)
    update = SiteUpdate(action, interval_phi=epsilon, seed=310, parallel=False)
    phi = lattice.zeros(0)
    values = []
    warmup, retained, block = 2_000, 10_000, 200
    for sweep in range(warmup + retained):
        phi = update.step({'phi': phi})['phi']
        if sweep >= warmup:
            theta = _principal(float(phi[0, 1] - phi[0, 0]))
            values.append(np.cos(theta))
    values = np.asarray(values)
    blocks = values.reshape(-1, block).mean(axis=1)
    estimate = blocks.mean()
    standard_error = blocks.std(ddof=1) / np.sqrt(len(blocks))
    exact = _exact_relative_cosine(kappa, epsilon) if expected is None else expected
    assert abs(estimate - exact) < 6.0 * standard_error


def test_parallel_kernel_is_reproducible_across_numba_thread_counts():
    maximum = numba.config.NUMBA_NUM_THREADS
    if maximum < 2:
        pytest.skip("the active Numba runtime permits only one thread")
    lattice = supervillain.lattice.Lattice(D=2, N=3)
    action = XY(lattice, kappa=0.4, epsilon=1.4)
    update = SiteUpdate(action, interval_phi=0.8, seed=0, parallel=True)
    rng = np.random.default_rng(31)
    phi = action.canonical(rng.uniform(0.0, 0.1, lattice.dims))
    proposals = rng.uniform(-0.8, 0.8, lattice.sites)
    uniforms = rng.uniform(0.0, 1.0, lattice.sites)
    previous = numba.get_num_threads()
    try:
        numba.set_num_threads(1)
        single = phi.copy()
        accepted_single = sweep_parallel(
            single.reshape(-1), update.plus, update.minus, update.color_sites, update.offsets,
            proposals, uniforms, action.kappa, action.epsilon)
        numba.set_num_threads(min(4, maximum))
        many = phi.copy()
        accepted_many = sweep_parallel(
            many.reshape(-1), update.plus, update.minus, update.color_sites, update.offsets,
            proposals, uniforms, action.kappa, action.epsilon)
    finally:
        numba.set_num_threads(previous)
    assert accepted_single == accepted_many
    assert np.array_equal(single, many)


def test_xy_ensemble_h5_roundtrip_and_restart(tmp_path):
    lattice = supervillain.lattice.Lattice(D=1, N=3)
    action = XY(lattice, kappa=0.45, epsilon=2.5)
    generator = SiteUpdate(action, interval_phi=0.7, seed=22, backend='numba', parallel=False)
    ensemble = supervillain.Ensemble(action).generate(6, generator, start='cold')
    path = tmp_path / 'xy.h5'

    with h5.File(path, 'w') as handle:
        ensemble.to_h5(handle.create_group('ensemble'))
    with h5.File(path, 'r') as handle:
        restored = supervillain.Ensemble.from_h5(handle['ensemble'])

    assert isinstance(restored.Action, XY)
    assert restored.Action.kappa == action.kappa
    assert restored.Action.epsilon == action.epsilon
    assert np.array_equal(restored.configuration.phi, ensemble.configuration.phi)
    continued = supervillain.Ensemble.continue_from(restored, 2)
    assert len(continued) == 2
    assert all(restored.Action.valid(cfg) for cfg in continued.configuration)
