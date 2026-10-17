"""Independent finite-product stationarity and joint-checkpoint checks."""
from pathlib import Path
import hashlib
import importlib.util
import itertools
import json

import h5py
import numpy as np
import pytest

P = Path(__file__).resolve().parents[1]/'experiments/2026-10-17-xy-3d-negative-kappa/replica-exchange'
_spec = importlib.util.spec_from_file_location('independent_pt_module', P/'pt.py')
pt = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(pt)


def _swap_matrix(states, kappas, C, edges):
    lookup = {tuple(s): i for i, s in enumerate(states)}
    matrix = np.zeros((len(states), len(states)))
    for index, state in enumerate(states):
        possibilities = [(list(state), 1.)]
        for a, b in edges:
            next_paths = []
            for x, probability in possibilities:
                logratio = pt.swap_log_ratio(kappas[a], kappas[b], C[x[a]], C[x[b]])
                y = list(x); y[a], y[b] = y[b], y[a]
                # Independent direct product-weight ratio, not reused formula.
                old_weight = np.exp(sum(k*C[s] for k, s in zip(kappas, x)))
                new_weight = np.exp(sum(k*C[s] for k, s in zip(kappas, y)))
                assert np.exp(logratio) == pytest.approx(new_weight/old_weight)
                accept = min(1., np.exp(logratio))
                next_paths += [(y, probability*accept), (x, probability*(1-accept))]
            possibilities = next_paths
        for x, probability in possibilities:
            matrix[index, lookup[tuple(x)]] += probability
    return matrix


def test_exact_individual_swap_product_detailed_balance():
    C = np.array([-1.4, .3, 1.9]); ks = [-.6, 0.]
    states = list(itertools.product(range(3), repeat=2))
    weights = np.array([np.exp(sum(k*C[s] for k, s in zip(ks, x))) for x in states])
    weights /= weights.sum()
    swap = _swap_matrix(states, ks, C, [(0, 1)])
    np.testing.assert_allclose(swap.sum(1), 1., atol=2e-15)
    np.testing.assert_allclose(weights[:, None]*swap, (weights[:, None]*swap).T, atol=2e-15)
    assert np.any((swap > 0) & (swap < 1))


def test_exact_alternating_exchange_and_local_kernel_stationarity():
    C = np.array([-.9, 1.2]); ks = [-.6, -.3, 0.]
    states = list(itertools.product(range(2), repeat=3))
    weights = np.array([np.exp(sum(k*C[s] for k, s in zip(ks, x))) for x in states]); weights /= weights.sum()
    kernels = []
    for k in ks:
        K = np.array([[0., .3], [.3, 0.]])*np.minimum(1., np.exp(k*(C[None, :]-C[:, None])))
        K[np.diag_indices(2)] = 1-K.sum(1); kernels.append(K)
    local = np.array([[np.prod([K[a, b] for K, a, b in zip(kernels, x, y)]) for y in states] for x in states])
    even = _swap_matrix(states, ks, C, [(0, 1)])
    odd = _swap_matrix(states, ks, C, [(1, 2)])
    for kernel in [local@even, local@odd, local@even@local@odd]:
        np.testing.assert_allclose(kernel.sum(1), 1., atol=2e-15)
        np.testing.assert_allclose(weights@kernel, weights, atol=2e-15)
    # Global composition is invariant but need not itself be reversible.
    whole = local@even@local@odd
    assert np.max(abs(weights[:, None]*whole-(weights[:, None]*whole).T)) > 1e-4


def test_odd_joint_epoch_and_unequal_uncommitted_dataset_tails_replay(tmp_path):
    kw = dict(N=4, epsilon=2.1, kappas=[-.5, -.25, 0.], warmup=2,
              seed=811, initial_winding=[1, 0, 0], checkpoint_every=3)
    whole, interrupted = tmp_path/'whole.h5', tmp_path/'interrupted.h5'
    pt.run(whole, steps=15, **kw)
    with pytest.raises(RuntimeError, match='Injected crash'):
        pt.run(interrupted, steps=15, inject_crash_after_flush=11, **kw)
    with h5py.File(interrupted, 'r+') as h:
        assert h.attrs['committed_rows'] == 10
        assert json.loads(h['rng_state'][9])['next_exchange_parity'] == 1
        h['magnetization_squared'].resize(17, axis=0)
        h['magnetization_squared'][12:] = -123.
        h['exchange_uniform'].resize(10, axis=0)
    pt.run(interrupted, steps=15, resume=True, **kw)
    with h5py.File(whole) as a, h5py.File(interrupted) as b:
        for key in a:
            if not key.endswith('_seconds'):
                np.testing.assert_array_equal(a[key][:], b[key][:], err_msg=key)


def test_source_mismatch_refusal_preserves_uncommitted_tail_bytes(tmp_path):
    file = tmp_path/'mismatch.h5'
    kw = dict(N=4, epsilon=2.1, kappas=[-.25, 0.], warmup=1,
              seed=819, initial_winding=[0, 0, 0], checkpoint_every=2)
    pt.run(file, steps=4, **kw)
    with h5py.File(file, 'r+') as h:
        sources = json.loads(h.attrs['source_hashes']); sources['replica_exchange_sha256'] = 'different'
        h.attrs['source_hashes'] = json.dumps(sources)
        h['magnetization_squared'].resize(7, axis=0)
    before = hashlib.sha256(file.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match='parameters/source mismatch'):
        pt.run(file, steps=6, resume=True, **kw)
    assert hashlib.sha256(file.read_bytes()).hexdigest() == before


def test_runner_uses_sibling_joint_audit_and_tiny_worker_closes(tmp_path, monkeypatch):
    import inspect
    import sys
    monkeypatch.syspath_prepend(str(P))
    spec = importlib.util.spec_from_file_location('independent_pt_runner', P/'runner.py')
    runner = importlib.util.module_from_spec(spec); spec.loader.exec_module(runner)
    assert Path(inspect.getsourcefile(runner.audit)).resolve() == (P/'audit.py').resolve()
    assert len(inspect.signature(runner.audit).parameters) == 3
    frozen_sources = runner.sources()
    point = runner.point(-.1, .65, 4, 2, steps=4)
    for directory in ['jobs', 'chains']:
        (tmp_path/directory).mkdir()
    manifest = tmp_path/'tiny-batch.json'
    manifest.write_text(json.dumps(dict(source_hashes=frozen_sources, points=[point])))
    monkeypatch.setattr(runner, 'P', tmp_path)
    monkeypatch.setattr(runner, 'sources', lambda: frozen_sources)
    monkeypatch.setattr(runner, 'route', lambda: None)
    runner.worker(manifest, point['id'])
    job = json.loads((tmp_path/'jobs'/f"{point['id']}.json").read_text())
    assert job['state'] == 'complete' and job['audit']['passed']
    assert job['audit']['field_configurations'] == 5*len(point['kappas'])
