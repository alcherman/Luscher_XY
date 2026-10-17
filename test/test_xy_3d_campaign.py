"""Focused tests for the independent 3D campaign audit."""
import importlib.util
from pathlib import Path

import h5py
import numpy as np

from supervillain.xy import run_chain


CAMPAIGN = Path(__file__).resolve().parents[1] / "experiments" / "2026-09-03-xy-3d-kappa-zero"
SPEC = importlib.util.spec_from_file_location("xy_3d_campaign_run", CAMPAIGN / "run.py")
run = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(run)


def test_3d_campaign_grid_has_two_distinct_starts():
    points = run.points()
    assert len(points) == 30
    assert sorted({point["epsilon"] for point in points}) == [2.4, 2.45, 2.5, 2.55, 2.6]
    assert sorted({point["L"] for point in points}) == [8, 12, 16]
    for epsilon in {point["epsilon"] for point in points}:
        for L in {point["L"] for point in points}:
            group = [point for point in points if point["epsilon"] == epsilon and point["L"] == L]
            assert [point["initial_winding"] for point in group] == [[0, 0, 0], [1, 0, 0]]


def test_3d_audit_recomputes_observables(tmp_path):
    point = run.point(250, 3, 0) | {"steps": 4}
    path = tmp_path / "tiny.h5"
    run_chain(path, D=3, N=3, kappa=0, epsilon=2.5, sweeps=4, warmup=1,
              seed=point["seed"], backend="numba", parallel=False,
              interval_phi=np.pi, winding=True, checkpoint_every=2,
              update="hybrid", clusters_per_step=4,
              initial_winding=[0, 0, 0], record_timing=True, buffer_rows=100)
    result = run.audit(path, point)
    assert result["passed"] and result["rows"] == 5
    with h5py.File(path, "r+") as f:
        f["vortex_density"][2] += 0.25
    try:
        run.audit(path, point)
    except AssertionError:
        pass
    else:
        raise AssertionError("audit accepted a corrupted 3D observable")
