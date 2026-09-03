#!/usr/bin/env python
r"""Independent tests for XY measurements and durable scan checkpoints."""

import h5py
import json
import numpy as np
import pytest

import supervillain
from supervillain.action.xy import XY
from supervillain.xy import main, measure, run_chain


def _run(path, sweeps, *, resume=False, **overrides):
    parameters = dict(
        D=2, N=3, kappa=0.4, epsilon=2.4, sweeps=sweeps, warmup=1,
        seed=128, backend="numba", parallel=True, interval_phi=0.8,
        winding=True, checkpoint_every=2, resume=resume,
    )
    parameters.update(overrides)
    return run_chain(path, **parameters)


def _payload(path):
    with h5py.File(path, "r") as handle:
        return {name: handle[name][...] for name in handle}, dict(handle.attrs)


def _assert_payload_equal(first, second):
    first_data, first_attrs = _payload(first)
    second_data, second_attrs = _payload(second)
    assert first_data.keys() == second_data.keys()
    for name in first_data:
        assert np.array_equal(first_data[name], second_data[name]), name
    for name in ("schema", "parameters", "committed_rows", "analysis_start_row"):
        assert first_attrs[name] == second_attrs[name]


def test_uninterrupted_chain_equals_partial_then_resumed_chain(tmp_path):
    uninterrupted = _run(tmp_path / "uninterrupted.h5", sweeps=7)
    resumed = _run(tmp_path / "resumed.h5", sweeps=3)
    _run(resumed, sweeps=7, resume=True)

    _assert_payload_equal(uninterrupted, resumed)
    data, attrs = _payload(resumed)
    assert attrs["committed_rows"] == 8
    assert np.array_equal(data["sweep"], np.arange(8))


def test_resume_discards_uncommitted_rows_and_replays_them_exactly(tmp_path):
    expected = _run(tmp_path / "expected.h5", sweeps=6)
    interrupted = _run(tmp_path / "interrupted.h5", sweeps=3)

    # Simulate a process killed after every dataset was extended but before the
    # commit marker advanced.  The bytes are deliberately wrong: correct code
    # must truncate them and replay from the last committed RNG state.
    with h5py.File(interrupted, "r+") as handle:
        committed = int(handle.attrs["committed_rows"])
        for dataset in handle.values():
            old = dataset.shape[0]
            dataset.resize(old + 2, axis=0)
            dataset[old:] = dataset[committed - 1]
        assert handle["sweep"].shape[0] > committed

    _run(interrupted, sweeps=6, resume=True)
    _assert_payload_equal(expected, interrupted)


def test_resume_rejects_parameter_mismatch_and_truncated_committed_checkpoint(tmp_path):
    path = _run(tmp_path / "chain.h5", sweeps=3)
    with pytest.raises(ValueError, match="parameters"):
        _run(path, sweeps=4, kappa=0.41, resume=True)

    with h5py.File(path, "r+") as handle:
        handle.attrs.modify("committed_rows", int(handle.attrs["committed_rows"]) + 1)
    with pytest.raises(ValueError, match="Incomplete or corrupt checkpoint"):
        _run(path, sweeps=4, resume=True)


def test_seed_and_warmup_rows_are_retained_and_labeled(tmp_path):
    path = _run(tmp_path / "warmup.h5", sweeps=5, warmup=2, winding=False)
    data, attrs = _payload(path)

    assert len(data["sweep"]) == 6             # row zero is the exact seed
    assert np.array_equal(data["sweep"], np.arange(6))
    assert attrs["analysis_start_row"] == 3
    assert attrs["committed_rows"] == 6
    assert np.all(data["phi"] >= 0.0)
    assert np.all(data["phi"] < 2.0 * np.pi)


def test_measured_action_density_matches_the_corrected_action():
    action = XY(supervillain.lattice.Lattice(D=2, N=3), kappa=1.3, epsilon=1.4)
    phi = action.Lattice.zeros(0)
    phi[0, 0, 0] = 0.4
    measured = measure(action, phi)
    expected = -action.kappa * (14 + 4 * np.cos(0.4)) / 9
    assert measured["action_density"] == pytest.approx(expected)
    assert measured["action_density"] == pytest.approx(action(phi) / 9)


def test_requested_initial_winding_is_stored_in_seed_row(tmp_path):
    path = run_chain(
        tmp_path / "twisted.h5", D=2, N=5, kappa=0.4, epsilon=1.4,
        sweeps=1, warmup=0, seed=9, winding=False, initial_winding=[1, 0],
        checkpoint_every=1,
    )
    data, _ = _payload(path)
    assert np.array_equal(data["winding"][0, 0], np.ones(5, dtype=int))
    assert np.array_equal(data["winding"][0, 1], np.zeros(5, dtype=int))


def test_initial_winding_cannot_alias_to_a_different_sector(tmp_path):
    path = tmp_path / "aliased.h5"
    with pytest.raises(ValueError, match="initial winding"):
        run_chain(path, D=1, N=5, kappa=0.4, epsilon=1.4,
                  sweeps=1, warmup=0, initial_winding=[5])
    assert not path.exists()


def test_measure_detects_vortex_threshold_and_a_smooth_winding():
    vortex_lattice = supervillain.lattice.Lattice(D=2, N=2)
    vortex_phi = np.array([[0.0, np.pi / 2], [3 * np.pi / 2, np.pi]])
    permissive = XY(vortex_lattice, kappa=0.5, epsilon=np.pi / 2 + 1e-12)
    forbidden = XY(vortex_lattice, kappa=0.5, epsilon=np.pi / 2)
    measured = measure(permissive, vortex_phi)
    assert measured["vortex_density"] == pytest.approx(1.0)
    with pytest.raises(ValueError, match="inadmissible"):
        measure(forbidden, vortex_phi)

    winding_lattice = supervillain.lattice.Lattice(D=2, N=5)
    winding_action = XY(winding_lattice, kappa=0.5, epsilon=1.4)
    coordinates = np.indices(winding_lattice.dims)
    twist = winding_action.canonical(2 * np.pi * coordinates[0] / winding_lattice.N)
    smooth = measure(winding_action, twist)
    assert smooth["vortex_density"] == 0.0
    assert np.array_equal(smooth["winding"][0], np.ones(winding_lattice.N, dtype=int))
    assert np.array_equal(smooth["winding"][1], np.zeros(winding_lattice.N, dtype=int))


def test_xy_scan_cli_smoke_creates_a_chain(tmp_path):
    output = tmp_path / "scan"
    main([
        "--D", "1", "--N", "3", "--kappa", "0.5", "--epsilon", "2.5",
        "--sweeps", "2", "--warmup", "0", "--seed", "7", "--serial",
        "--no-winding", "--checkpoint-every", "1", "--output", str(output),
    ])
    path = output / "chain-00000.h5"
    assert path.exists()
    data, attrs = _payload(path)
    assert attrs["committed_rows"] == 3
    assert np.array_equal(data["sweep"], np.arange(3))


@pytest.mark.parametrize("backend", ["python", "numba"])
@pytest.mark.parametrize("winding", [False, True])
def test_cluster_checkpoint_replay_and_fixed_cycle_work(tmp_path, backend, winding):
    options = dict(update="cluster", clusters_per_step=3, interval_phi=None,
                   backend=backend, winding=winding, epsilon=1.4)
    expected = _run(tmp_path / "expected.h5", 7, **options)
    resumed = _run(tmp_path / "resumed.h5", 3, **options)
    with h5py.File(resumed, "r+") as handle:
        for dataset in handle.values():
            old = len(dataset)
            dataset.resize(old + 2, axis=0)
            dataset[old:] = dataset[old - 1]
    _run(resumed, 7, resume=True, **options)
    _assert_payload_equal(expected, resumed)
    data, attrs = _payload(expected)
    parameters = json.loads(attrs["parameters"])
    assert parameters["time_unit"] == "cluster_cycle"
    assert parameters["parallel"] is False
    assert "local_accepted" not in data
    assert data["cluster_sizes"].shape == (8, 3)
    assert np.all(data["cluster_sizes"][0] == 0)
    assert np.all((data["cluster_sizes"][1:] >= 1) & (data["cluster_sizes"][1:] <= 9))
    assert np.array_equal(data["clusters_total"], 3 * np.arange(8))
    assert np.array_equal(data["cluster_sites_total"], np.cumsum(data["cluster_sizes"].sum(axis=1)))
    assert np.array_equal(data["cluster_bond_tests_total"], np.cumsum(data["cluster_bond_tests"].sum(axis=1)))
    assert np.array_equal(data["cluster_equivalent_sweeps"], data["cluster_sites_total"] / 9)
    action = XY(supervillain.Lattice(D=2, N=3), 0.4, 1.4)
    assert all(action.valid({"phi": action.canonical(phi)}) for phi in data["phi"])


def test_cluster_scan_backends_share_the_same_random_stream(tmp_path):
    options = dict(update="cluster", clusters_per_step=2, interval_phi=None, winding=False)
    first = _run(tmp_path / "python.h5", 6, backend="python", **options)
    second = _run(tmp_path / "numba.h5", 6, backend="numba", **options)
    a, _ = _payload(first)
    b, _ = _payload(second)
    for name in a:
        assert np.array_equal(a[name], b[name]), name


@pytest.mark.parametrize("options", [
    dict(update="cluster", interval_phi=0.4),
    dict(update="cluster", interval_phi=None, kappa=-0.1),
    dict(update="cluster", interval_phi=None, clusters_per_step=0),
    dict(update="cluster", interval_phi=None, clusters_per_step=True),
    dict(update="metropolis", clusters_per_step=2),
    dict(update="unknown"),
])
def test_invalid_cluster_scan_parameters_fail_before_creating_file(tmp_path, options):
    path = tmp_path / "invalid.h5"
    with pytest.raises(ValueError):
        _run(path, 3, **options)
    assert not path.exists()


def test_cluster_scan_resume_rejects_a_changed_schedule(tmp_path):
    path = _run(tmp_path / "cluster.h5", 3, update="cluster", interval_phi=None)
    with pytest.raises(ValueError, match="parameters"):
        _run(path, 4, resume=True, update="cluster", interval_phi=None, clusters_per_step=2)
    with pytest.raises(ValueError, match="parameters"):
        _run(path, 4, resume=True)


def test_cluster_cli_smoke(tmp_path):
    main(["--D", "3", "--N", "2", "--kappa", "0", "--epsilon", "1.4",
          "--update", "cluster", "--clusters-per-step", "2", "--sweeps", "2",
          "--warmup", "0", "--no-winding", "--output", str(tmp_path)])
    data, attrs = _payload(tmp_path / "chain-00000.h5")
    assert attrs["committed_rows"] == 3
    assert data["cluster_sizes"].shape == (3, 2)


@pytest.mark.parametrize("backend", ["python", "numba"])
def test_cluster_ensemble_hdf5_restores_exact_continuation(tmp_path, backend):
    from supervillain.generator.xy import ClusterUpdate

    action = XY(supervillain.Lattice(2, 3), 1.1, 1.6)
    generator = ClusterUpdate(action, seed=319, backend=backend)
    ensemble = supervillain.Ensemble(action).generate(4, generator, start="cold")
    path = tmp_path / "ensemble.h5"
    with h5py.File(path, "w") as handle:
        ensemble.to_h5(handle.create_group("ensemble"))
    with h5py.File(path, "r") as handle:
        restored = supervillain.Ensemble.from_h5(handle["ensemble"])
    expected = supervillain.Ensemble.continue_from(ensemble, 5)
    actual = supervillain.Ensemble.continue_from(restored, 5)
    assert np.array_equal(expected.configuration.phi, actual.configuration.phi)
    for name in ("clusters", "flipped_sites", "bond_tests", "last_cluster_size",
                 "last_bond_tests", "equivalent_sweeps"):
        assert getattr(expected.generator, name) == getattr(actual.generator, name)


def test_cluster_resume_after_actual_mid_row_write_exception(tmp_path, monkeypatch):
    options = dict(update="cluster", clusters_per_step=2, interval_phi=None)
    expected = _run(tmp_path / "expected.h5", 7, **options)
    interrupted = tmp_path / "interrupted.h5"
    original = h5py.Dataset.__setitem__

    def write_with_fault(dataset, key, value):
        if dataset.name == "/cluster_sites_total" and key == 4:
            raise RuntimeError("injected mid-row write failure")
        return original(dataset, key, value)

    with monkeypatch.context() as fault:
        fault.setattr(h5py.Dataset, "__setitem__", write_with_fault)
        with pytest.raises(RuntimeError, match="injected mid-row"):
            _run(interrupted, 7, **options)
    with h5py.File(interrupted, "r") as handle:
        assert handle.attrs["committed_rows"] == 3
        assert len({len(dataset) for dataset in handle.values()}) == 2
    _run(interrupted, 7, resume=True, **options)
    _assert_payload_equal(expected, interrupted)


@pytest.mark.parametrize("backend", ["python", "numba"])
def test_hybrid_matches_explicit_fixed_composition_and_resumes(tmp_path, backend):
    from supervillain.generator.xy import ClusterUpdate, SiteUpdate, WindingUpdate

    options = dict(update="hybrid", clusters_per_step=2, backend=backend, epsilon=1.4)
    expected = _run(tmp_path / "expected.h5", 6, **options)
    partial = _run(tmp_path / "partial.h5", 3, **options)
    _run(partial, 6, resume=True, **options)
    _assert_payload_equal(expected, partial)
    data, attrs = _payload(expected)
    params = json.loads(attrs["parameters"])
    assert params["time_unit"] == "hybrid_cycle"
    action = XY(supervillain.Lattice(2, 3), 0.4, 1.4)
    streams = np.random.SeedSequence(128).spawn(3)
    local = SiteUpdate(action, interval_phi=0.8, seed=streams[0], backend=backend, parallel=True)
    cluster = ClusterUpdate(action, seed=streams[2], backend=backend)
    winding = WindingUpdate(action, seed=streams[1])
    cfg = {"phi": action.Lattice.zeros(0)}
    for cycle in range(1, 7):
        before = local.accepted
        cfg = local.step(cfg)
        assert data["local_accepted"][cycle] == local.accepted-before
        for index in range(2):
            cfg = cluster.step(cfg)
            assert data["cluster_sizes"][cycle,index] == cluster.last_cluster_size
        cfg = winding.step(cfg)
        assert np.array_equal(data["phi"][cycle], cfg["phi"])
        assert data["cluster_sites_total"][cycle] == cluster.flipped_sites
        state = json.loads(data["rng_state"][cycle])
        assert state["cluster"] == cluster.rng.bit_generator.state
        assert state["local"] == local.rng.bit_generator.state
        assert state["winding"] == winding.rng.bit_generator.state


@pytest.mark.parametrize("update", ["metropolis", "cluster", "hybrid"])
def test_recorded_timings_do_not_change_sampling_and_replay(tmp_path, update):
    options = dict(update=update, interval_phi=None)
    plain = _run(tmp_path / "plain.h5", 6, **options)
    timed = _run(tmp_path / "timed.h5", 3, record_timing=True, **options)
    _run(timed, 6, resume=True, record_timing=True, **options)
    a,_ = _payload(plain)
    b,_ = _payload(timed)
    for name in a:
        assert np.array_equal(a[name], b[name]), name
    for name in set(b)-set(a):
        assert name.endswith("_seconds")
        assert np.all(b[name]>=0) and np.isfinite(b[name]).all()
    assert np.all(b["cycle_seconds"][1:] >= b["update_seconds"][1:] + b["measurement_seconds"][1:])


def test_hybrid_cli_and_negative_coupling_guard(tmp_path):
    main(["--D","2","--N","3","--kappa","0.8","--update","hybrid",
          "--record-timing","--sweeps","2","--warmup","0","--serial",
          "--output",str(tmp_path)])
    data,_ = _payload(tmp_path / "chain-00000.h5")
    assert "local_accepted" in data and "cluster_sizes" in data and "update_seconds" in data
    with pytest.raises(ValueError,match="kappa"):
        _run(tmp_path / "invalid.h5",3,update="hybrid",kappa=-0.1)
    assert not (tmp_path / "invalid.h5").exists()


@pytest.mark.parametrize("update", ["metropolis", "cluster", "hybrid"])
@pytest.mark.parametrize("buffer_rows", [3, 20])
@pytest.mark.parametrize("kappa,epsilon", [(0.0, np.pi/2), (0.8, 1.8)])
def test_buffered_records_match_direct_writer_and_resume(tmp_path, update, buffer_rows, kappa, epsilon):
    options = dict(N=4, update=update, interval_phi=None, kappa=kappa, epsilon=epsilon,
                   checkpoint_every=5, record_timing=True)
    direct = _run(tmp_path / "direct.h5", 13, **options)
    buffered = _run(tmp_path / "buffered.h5", 7, buffer_rows=buffer_rows, **options)
    _run(buffered, 13, resume=True, buffer_rows=buffer_rows, **options)
    a, _ = _payload(direct)
    b, attrs = _payload(buffered)
    assert a.keys() == b.keys()
    assert attrs["committed_rows"] == 14
    for name in a:
        if not name.endswith("_seconds"):
            assert np.array_equal(a[name], b[name]), name
        else:
            assert np.isfinite(b[name]).all() and np.all(b[name] >= 0)
    assert np.all(b["cycle_seconds"][1:] >= b["update_seconds"][1:] + b["measurement_seconds"][1:])


@pytest.mark.parametrize("backend", ["python", "numba"])
def test_buffered_hybrid_replays_partial_batch_failure(tmp_path, monkeypatch, backend):
    options = dict(N=4, update="hybrid", backend=backend, buffer_rows=3, checkpoint_every=5)
    expected = _run(tmp_path / "expected.h5", 13, **options)
    interrupted = tmp_path / "interrupted.h5"
    original = h5py.Dataset.__setitem__

    def write_with_fault(dataset, key, value):
        if dataset.name == "/cluster_sites_total" and isinstance(key, slice) and key.start == 6:
            raise RuntimeError("injected partial batch failure")
        return original(dataset, key, value)

    with monkeypatch.context() as fault:
        fault.setattr(h5py.Dataset, "__setitem__", write_with_fault)
        with pytest.raises(RuntimeError, match="partial batch"):
            _run(interrupted, 13, **options)
    with h5py.File(interrupted) as handle:
        assert handle.attrs["committed_rows"] == 6
        assert len({len(dataset) for dataset in handle.values()}) == 2
    _run(interrupted, 13, resume=True, **options)
    _assert_payload_equal(expected, interrupted)
