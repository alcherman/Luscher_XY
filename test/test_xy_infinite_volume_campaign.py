"""Independent audit and schedule/checkpoint invariants for the October campaign."""
from datetime import datetime
import importlib.util
from pathlib import Path
import signal
import subprocess
import sys
import time

import h5py
import numpy as np
import pytest

from supervillain.xy import run_chain

PATH = Path(__file__).resolve().parents[1] / "experiments/2026-10-13-xy-infinite-volume/run.py"
SPEC = importlib.util.spec_from_file_location("xy_infinite_volume_run", PATH)
run = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(run)


@pytest.mark.parametrize("when,expected", [
    (datetime(2026, 10, 13, 8, 59, tzinfo=run.CENTRAL), 18),
    (datetime(2026, 10, 13, 9, 0, tzinfo=run.CENTRAL), 16),
    (datetime(2026, 10, 13, 16, 59, tzinfo=run.CENTRAL), 16),
    (datetime(2026, 10, 13, 17, 0, tzinfo=run.CENTRAL), 18),
    (datetime(2026, 10, 17, 12, 0, tzinfo=run.CENTRAL), 18),
])
def test_office_hours_reserve(when, expected):
    assert run.worker_limit(when) == expected


@pytest.mark.parametrize("D", [2, 3, 4])
def test_audit_rejects_corrupt_observable_and_resume_is_exact(tmp_path, D):
    spec = run.point(D, 100, 302, 3, 1, 8)
    base = dict(D=D, N=3, kappa=spec["beta"], epsilon=spec["epsilon"], warmup=2,
                seed=spec["seed"], backend="numba", parallel=False, interval_phi=np.pi,
                winding=True, initial_winding=spec["initial_winding"], checkpoint_every=2,
                update="hybrid", clusters_per_step=4, record_timing=True, buffer_rows=100)
    direct, resumed = tmp_path / "direct.h5", tmp_path / "resumed.h5"
    run_chain(direct, sweeps=8, **base)
    run_chain(resumed, sweeps=4, **base)
    run_chain(resumed, sweeps=8, resume=True, **base)
    assert run.audit(direct, spec)["passed"]
    assert run.audit(resumed, spec)["passed"]
    with h5py.File(direct) as f, h5py.File(resumed) as g:
        for key in f:
            if not key.endswith("_seconds"):
                np.testing.assert_array_equal(f[key][:], g[key][:], err_msg=key)
    with h5py.File(resumed, "r+") as f:
        f["vortex_density"][3] += .25
    with pytest.raises(AssertionError):
        run.audit(resumed, spec)


def test_worker_termination_recovers_committed_checkpoint(tmp_path):
    """Exercise the actual worker journal through a separate terminated process."""
    spec = run.point(2, 100, 302, 3, 1, 10000)
    spec["warmup"] = 100
    code = """
import importlib.util, json, signal, sys
from pathlib import Path
spec = importlib.util.spec_from_file_location('worker_run', sys.argv[1])
run = importlib.util.module_from_spec(spec)
spec.loader.exec_module(run)
run.OUT = Path(sys.argv[2])
def terminate(signum, frame):
    raise KeyboardInterrupt('termination test')
signal.signal(signal.SIGTERM, terminate)
run.trajectory(json.loads(sys.argv[3]))
"""
    args = [sys.executable, "-u", "-c", code, str(PATH), str(tmp_path), __import__("json").dumps(spec)]
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    transcript = []
    try:
        for line in proc.stdout:
            transcript.append(line)
            if "checkpoint" in line:
                proc.send_signal(signal.SIGTERM)
                break
        rest, _ = proc.communicate(timeout=30)
        transcript.append(rest)
    finally:
        if proc.poll() is None:
            proc.kill()
    (tmp_path / "interruption.log").write_text("".join(transcript))
    assert proc.returncode != 0
    job = __import__("json").loads((tmp_path / "jobs" / f"{spec['id']}.json").read_text())
    assert job["state"] == "interrupted"
    resumed = subprocess.run(args, capture_output=True, text=True, timeout=60)
    assert resumed.returncode == 0, resumed.stdout + resumed.stderr
    chain = tmp_path / "chains" / f"{spec['id']}.h5"
    assert run.audit(chain, spec)["passed"]
    direct = tmp_path / "direct.h5"
    run_chain(direct, D=2, N=3, kappa=spec["beta"], epsilon=spec["epsilon"],
              sweeps=10000, warmup=100, seed=spec["seed"], backend="numba",
              parallel=False, interval_phi=np.pi, winding=True,
              initial_winding=spec["initial_winding"], checkpoint_every=100,
              update="hybrid", clusters_per_step=4, record_timing=True, buffer_rows=100)
    with h5py.File(direct) as f, h5py.File(chain) as g:
        for key in f:
            if not key.endswith("_seconds"):
                np.testing.assert_array_equal(f[key][:], g[key][:], err_msg=key)
