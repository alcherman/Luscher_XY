"""Focused tests for the independent 4D campaign audit."""
import importlib.util
from pathlib import Path

import h5py
import numpy as np

from supervillain.xy import run_chain


CAMPAIGN = Path(__file__).resolve().parents[1] / "experiments" / "2026-09-03-xy-4d-kappa-zero"
SPEC = importlib.util.spec_from_file_location("xy_4d_campaign_run", CAMPAIGN / "run.py")
run = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(run)


def test_4d_campaign_grid_has_two_distinct_starts():
    points = run.points()
    assert len(points) == 24
    assert [round(value, 12) for value in sorted({point["epsilon_over_pi"] for point in points})] == [
        0.5, round(2 / 3, 12), round(5 / 6, 12), 1.0]
    assert sorted({point["L"] for point in points}) == [8, 12, 16]
    for ratio in {point["epsilon_over_pi"] for point in points}:
        for L in {point["L"] for point in points}:
            group = [point for point in points
                     if point["epsilon_over_pi"] == ratio and point["L"] == L]
            assert [point["initial_winding"] for point in group] == [
                [0, 0, 0, 0], [1, 0, 0, 0]]


def test_4d_refinement_divides_transition_bracket_evenly():
    points = run.refinement_points()
    assert len(points) == 18
    ratios = sorted({point["epsilon_over_pi"] for point in points})
    assert ratios == [7 / 8, 11 / 12, 23 / 24]
    complete_grid = [5 / 6] + ratios + [1]
    np.testing.assert_allclose(np.diff(complete_grid), np.full(4, 1 / 24), rtol=0, atol=2e-16)


def test_4d_near_critical_refinement():
    points = run.near_critical_points()
    assert len(points) == 12
    assert sorted({point["epsilon_over_pi"] for point in points}) == [0.84, 0.85]
    assert sorted({point["L"] for point in points}) == [8, 12, 16]


def test_4d_transition_refinement():
    points = run.transition_refinement_points()
    assert len(points) == 12
    assert sorted({point["epsilon_over_pi"] for point in points}) == [0.86, 0.87]
    assert sorted({point["L"] for point in points}) == [8, 12, 16]


def test_4d_large_volume_targets_bracket_endpoints():
    points = run.large_volume_points()
    assert len(points) == 4
    assert {point["L"] for point in points} == {20}
    assert sorted({point["epsilon_over_pi"] for point in points}) == [0.85, 7 / 8]
    assert [point["initial_winding"] for point in points] == [
        [0, 0, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0], [1, 0, 0, 0]]


def test_4d_audit_recomputes_observables(tmp_path):
    point = run.point(2 / 3, 3, 0) | {"steps": 4}
    path = tmp_path / "tiny.h5"
    run_chain(path, D=4, N=3, kappa=0, epsilon=2 * np.pi / 3,
              sweeps=4, warmup=1, seed=point["seed"], backend="numba",
              parallel=False, interval_phi=np.pi, winding=True,
              checkpoint_every=2, update="hybrid", clusters_per_step=4,
              initial_winding=[0, 0, 0, 0], record_timing=True, buffer_rows=100)
    result = run.audit(path, point)
    assert result["passed"] and result["rows"] == 5
    with h5py.File(path, "r+") as f:
        f["structure_min"][2, 3] += 0.25
    try:
        run.audit(path, point)
    except AssertionError:
        pass
    else:
        raise AssertionError("audit accepted a corrupted 4D observable")
