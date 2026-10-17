"""Independent contour, transport and shared-time resampling checks."""
from pathlib import Path
import importlib.util
import sys

import h5py
import numpy as np
import pytest

P = Path(__file__).resolve().parents[1]/'experiments/2026-10-17-xy-3d-negative-kappa/replica-exchange'
sys.path.insert(0, str(P))
_spec = importlib.util.spec_from_file_location('independent_pt_analysis', P/'analyse.py')
a = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(a)


@pytest.mark.parametrize('method', ['linear', 'pchip'])
def test_root_knot_and_branch_ledger(method):
    ks = np.array([-1., -.5, 0.])
    assert a.root(ks, [.1, .25, .4], method=method) == -.5
    assert a.root(ks, [.25, .25, .4], method=method, all_roots=True) == []
    assert a.root(ks, [.1, np.nan, .4], method=method, all_roots=True) == []
    assert a.root(ks, [.1, .15, .2], method=method) is None
    grid = np.array([-1., -.6, -.3, 0.])
    roots = a.root(grid, [.1, .4, .1, .4], method=method, all_roots=True)
    assert len(roots) == 3
    assert a.root(grid, [.1, .4, .1, .4], method=method) == max(roots)
    assert all(grid[i] < r < grid[i+1] for i, r in enumerate(roots))


def test_counted_round_trips_are_full_target_bridge_target_cycles():
    walkers = np.array([[0, 1], [1, 0], [0, 1], [1, 0], [0, 1]])
    report = a.trips(walkers)
    assert report['round_trips'] == 3 and report['durations'] == [2, 2, 2]
    assert report['target_unique_walkers'] == report['bridge_unique_walkers'] == 2
    stuck = a.trips(np.tile([0, 1], (10, 1)))
    assert stuck['round_trips'] == 0


def test_diagnose_resamples_all_slot_moments_with_identical_time_blocks(tmp_path, monkeypatch):
    # Artificial exactly dependent slot/moment series isolate covariance;
    # these deliberately are not audited physical MC configurations.
    monkeypatch.setattr(a, 'P', tmp_path)
    monkeypatch.setattr(a.helper, 'tau_int', lambda values: .5)
    (tmp_path/'chains').mkdir()
    T, S, N = 121, 3, 4
    scales = np.array([1., 2., 3.])
    m2 = np.linspace(.02, .08, T)[:, None]*scales
    file = tmp_path/'chains/fake.h5'
    with h5py.File(file, 'w') as h:
        h['magnetization_squared'] = m2
        h['magnetization_fourth'] = m2*m2
        h['structure_zero'] = N**3*m2
        h['structure_min'] = np.repeat((N**3*m2/3)[:, :, None], 3, axis=2)
        h['winding'] = np.zeros((T, S, 3, N*N), dtype=int)
        h['mean_link_cosine'] = np.ones((T, S))*.7
        h['sweep'] = np.arange(T)
        h['exchange_accepted'] = np.zeros((T, S-1), dtype=int)
        h['exchange_attempted'] = np.zeros((T, S-1), dtype=int)
        h['walker_id'] = np.tile(np.arange(S), (T, 1))
        h['phi'] = np.zeros((T, S, 1, N, N, N))
    point = dict(id='fake', N=N, kappas=[-.5, -.25, 0.], warmup=0)
    mean, draws, report = a.diagnose(dict(point=point, generation_seconds=0., audit_seconds=0.), np.random.default_rng(315), B=200)
    assert draws is not None and report['blocks'] == 3
    for slot, scale in enumerate(scales):
        np.testing.assert_allclose(draws[:, slot, 0], draws[:, 0, 0]*scale, atol=2e-17)
        np.testing.assert_allclose(draws[:, slot, 1], draws[:, 0, 1]*scale**2, atol=2e-17)
        np.testing.assert_allclose(draws[:, slot, 2], draws[:, slot, 0]*N**3, atol=2e-15)
