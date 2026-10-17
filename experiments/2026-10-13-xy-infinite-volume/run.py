"""Frozen-batch XY runner, detached from Codex with launchd.

All fields, RNG states, observables, transient rows and logs are retained.
The coordinator limits running single-thread workers by Central office hours.
"""
from __future__ import annotations

import os
for key in ("NUMBA_NUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"

import argparse
from datetime import datetime
import fcntl
import hashlib
import json
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import traceback
from zoneinfo import ZoneInfo

import h5py
import numpy as np

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
sys.path.insert(0, str(ROOT))
from supervillain.xy import run_chain

CENTRAL = ZoneInfo("America/Chicago")
DEADLINE = datetime(2026, 10, 16, 15, tzinfo=CENTRAL)
RESERVE = 25 * 2**30
LABEL = "org.luscherxy.phase-boundary-20261013"


def now():
    return datetime.now(CENTRAL).isoformat()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, prefix=path.name + ".", suffix=".tmp", delete=False) as stream:
        temporary = Path(stream.name)
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024**2), b""):
            digest.update(block)
    return digest.hexdigest()


def source_hashes():
    paths = [Path(__file__), ROOT / "supervillain/xy.py", ROOT / "supervillain/action/xy.py"]
    paths += [ROOT / "supervillain/generator" / name for name in
              ("xy.py", "_xy_kernels.py", "_xy_cluster.py", "reference_implementation/xy.py", "reference_implementation/xy_cluster.py")]
    paths += [ROOT / ".python-version", ROOT / "requirements-xy.lock.txt"]
    return {str(p.relative_to(ROOT)): sha256(p) for p in paths}


def worker_limit(when, cores=18):
    when = when.astimezone(CENTRAL)
    office = when.weekday() < 5 and 9 <= when.hour < 17
    return max(1, cores - 2) if office else cores


def point(D, ratio_centi, beta_milli, N, replica, steps):
    identifier = f"d{D}_e{ratio_centi:03d}_b{beta_milli:04d}_N{N:03d}_r{replica:02d}_s{steps:06d}"
    seed = int(np.random.SeedSequence([2026101301, D, ratio_centi, beta_milli, N, replica]).generate_state(1, dtype=np.uint64)[0])
    winding = [0] * D
    if replica % 2:
        winding[(replica // 2) % D] = 1
    return dict(id=identifier, D=D, epsilon_over_pi=ratio_centi / 100,
                epsilon=ratio_centi * np.pi / 100, beta=beta_milli / 1000,
                N=N, replica=replica, seed=seed, initial_winding=winding,
                steps=steps, warmup=min(1000, steps // 4), clusters_per_cycle=4)


def estimated_bytes(spec):
    # Dense float64 fields and every-loop winding, plus conservative HDF5,
    # scalar, RNG and timing overhead. No compression/thinning is assumed.
    D, N = spec["D"], spec["N"]
    per_row = 8 * (N**D + D * N**(D - 1)) + 4096
    return int(1.15 * (spec["steps"] + 1) * per_row + 16 * 2**20)


def freeze(name, points, purpose):
    if any(p["D"] not in (2, 3, 4, 5) or p["N"] < 3 or p["steps"] < 1 for p in points):
        raise ValueError("Invalid simulation specification")
    if any(p["D"] == 4 and p["N"] > 24 for p in points):
        raise ValueError("4D boxes above N=24 require a separately reviewed budget")
    if any(p["D"] == 5 and p["N"] > 12 for p in points):
        raise ValueError("Optional 5D boxes above N=12 require a separately reviewed budget")
    if len({p["id"] for p in points}) != len(points):
        raise ValueError("Duplicate trajectory IDs")
    body = dict(name=name, purpose=purpose, points=points, source_hashes=source_hashes(),
                estimated_bytes=sum(map(estimated_bytes, points)), created=now())
    path = OUT / "batches" / f"{name}.json"
    if path.exists():
        existing = json.loads(path.read_text())
        if any(existing[k] != body[k] for k in ("points", "source_hashes", "purpose")):
            raise ValueError("Existing frozen batch differs")
        return existing
    atomic_json(path, body)
    return body


def audit(path, spec):
    """Independent full-field recomputation, generalized from prior 3D/4D audits."""
    D, N, steps = spec["D"], spec["N"], spec["steps"]
    volume = N**D
    spatial = tuple(range(1, D + 1))
    with h5py.File(path, "r") as f:
        pars = json.loads(f.attrs["parameters"])
        expected = dict(D=D, N=N, kappa=spec["beta"], epsilon=spec["epsilon"],
                        seed=spec["seed"], update="hybrid", parallel=False,
                        clusters_per_step=4, initial_winding=spec["initial_winding"],
                        winding=True, buffer_rows=100)
        for key, value in expected.items():
            assert pars[key] == value, (key, pars[key], value)
        assert f.attrs["committed_rows"] == steps + 1
        assert all(len(ds) == steps + 1 for ds in f.values())
        np.testing.assert_array_equal(f["sweep"][:], np.arange(steps + 1))
        coordinates = np.indices((N,) * D)
        initial = np.mod(sum(2 * np.pi * w * coordinates[a] / N for a, w in enumerate(spec["initial_winding"])), 2 * np.pi)
        np.testing.assert_array_equal(f["phi"][0, 0], initial)
        sizes = f["cluster_sizes"][:]
        assert np.all(sizes[0] == 0)
        assert np.all((sizes[1:] >= 1) & (sizes[1:] <= volume))
        np.testing.assert_array_equal(f["clusters_total"][:], 4 * np.arange(steps + 1))
        np.testing.assert_array_equal(f["cluster_sites_total"][:], np.cumsum(sizes.sum(1)))
        np.testing.assert_array_equal(f["cluster_bond_tests_total"][:], np.cumsum(f["cluster_bond_tests"][:].sum(1)))
        np.testing.assert_allclose(f["cluster_equivalent_sweeps"][:], np.cumsum(sizes.sum(1)) / volume, rtol=0, atol=1e-12)
        assert np.all((f["local_accepted"][:] >= 0) & (f["local_accepted"][:] <= volume))
        for key in f:
            if key.endswith("_seconds"):
                v = f[key][:]
                assert np.isfinite(v).all() and np.all(v >= 0)
        for state in f["rng_state"][:]:
            assert set(json.loads(state)) == {"local", "cluster", "winding"}
        for first in range(0, steps + 1, 16):
            stop = min(first + 16, steps + 1)
            field = f["phi"][first:stop, 0]
            assert np.isfinite(field).all() and np.all((field >= 0) & (field < 2 * np.pi))
            links = np.stack([np.remainder(np.roll(field, -1, axis=a + 1) - field + np.pi, 2 * np.pi) - np.pi for a in range(D)], axis=1)
            assert np.all(np.abs(links) < spec["epsilon"])
            spin = np.exp(1j * field)
            magnetization = spin.mean(axis=spatial)
            m2 = np.abs(magnetization)**2
            cosine = np.cos(links).mean(axis=tuple(range(1, D + 2)))
            vortex_sum = np.zeros(stop - first)
            for a in range(D):
                for b in range(a + 1, D):
                    curl = links[:, a] + np.roll(links[:, b], -1, axis=a + 1) - np.roll(links[:, a], -1, axis=b + 1) - links[:, b]
                    vortex_sum += np.abs(np.rint(curl / (2 * np.pi))).mean(axis=spatial)
            winding = np.stack([np.rint(links[:, a].sum(axis=a + 1) / (2 * np.pi)).reshape(stop - first, -1) for a in range(D)], axis=1)
            minimum = []
            for a in range(D):
                shape = [1] * (D + 1)
                shape[a + 1] = N
                phase = np.exp(-2j * np.pi * np.arange(N) / N).reshape(shape)
                minimum.append(np.abs((spin * phase).sum(axis=spatial))**2 / volume)
            values = dict(magnetization_real=magnetization.real, magnetization_imag=magnetization.imag,
                          magnetization_abs=np.abs(magnetization), magnetization_squared=m2,
                          magnetization_fourth=m2 * m2, mean_link_cosine=cosine,
                          action_density=-D * spec["beta"] * cosine,
                          vortex_density=vortex_sum / (D * (D - 1) / 2), winding=winding,
                          structure_zero=volume * m2, structure_min=np.stack(minimum, axis=1),
                          max_link_angle=np.abs(links).max(axis=tuple(range(1, D + 2))))
            for key, value in values.items():
                np.testing.assert_allclose(f[key][first:stop], value, rtol=2e-12, atol=3e-13, err_msg=key)
    return dict(passed=True, rows=steps + 1, bytes=Path(path).stat().st_size, sha256=sha256(path))


def trajectory(spec):
    for folder in ("chains", "jobs", "logs"):
        (OUT / folder).mkdir(exist_ok=True)
    identifier = spec["id"]
    job = OUT / "jobs" / f"{identifier}.json"
    path = OUT / "chains" / f"{identifier}.h5"
    with (OUT / "jobs" / f"{identifier}.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        status = json.loads(job.read_text()) if job.exists() else {}
        if status.get("state") == "complete":
            return
        status.update(state="running", point=spec, pid=os.getpid(), updated=now(), paused=False)
        atomic_json(job, status)
        started = time.perf_counter()
        try:
            print(now(), "start/resume", identifier, flush=True)
            targets = [t for t in range(1000, spec["steps"], 1000) if t > spec["warmup"]] + [spec["steps"]]
            committed = 0
            if path.exists():
                with h5py.File(path) as f:
                    committed = int(f.attrs["committed_rows"]) - 1
            for target in targets:
                if target <= committed:
                    continue
                run_chain(path, D=spec["D"], N=spec["N"], kappa=spec["beta"], epsilon=spec["epsilon"],
                          sweeps=target, warmup=spec["warmup"], seed=spec["seed"],
                          backend="numba", parallel=False, interval_phi=np.pi,
                          winding=True, initial_winding=spec["initial_winding"],
                          checkpoint_every=100, resume=path.exists(), update="hybrid",
                          clusters_per_step=4, record_timing=True, buffer_rows=100)
                print(now(), "checkpoint", identifier, target, "cycles", flush=True)
            generated = time.perf_counter()
            print(now(), "generation complete; auditing", identifier, flush=True)
            result = audit(path, spec)
            status.update(state="complete", audit=result, updated=now(),
                          attempt_generation_seconds=generated - started,
                          audit_seconds=time.perf_counter() - generated)
            atomic_json(job, status)
            print(now(), "closed", identifier, result, flush=True)
        except BaseException as error:
            status.update(state="interrupted" if isinstance(error, KeyboardInterrupt) else "failed", updated=now(), error=traceback.format_exc())
            atomic_json(job, status)
            traceback.print_exc()
            raise


def alive(pid, identifier):
    if not isinstance(pid, int):
        return False
    result = subprocess.run(["/bin/ps", "-p", str(pid), "-o", "command="], capture_output=True, text=True)
    return result.returncode == 0 and str(OUT / "run.py") in result.stdout and "worker" in result.stdout and identifier in result.stdout


def coordinate():
    OUT.mkdir(exist_ok=True)
    with (OUT / "coordinator.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        config = json.loads((OUT / "control.json").read_text())
        cores = min(config["cores"], os.cpu_count() or 1)
        print(now(), "coordinator", os.getpid(), "cores", cores, flush=True)
        child_handles = {}
        while True:
            config = json.loads((OUT / "control.json").read_text())
            points = {}
            for manifest_path in sorted((OUT / "batches").glob("*.json")):
                manifest = json.loads(manifest_path.read_text())
                if manifest["source_hashes"] != source_hashes():
                    raise RuntimeError("Source changed after freezing " + str(manifest_path))
                for p in manifest["points"]:
                    if p["id"] in points and points[p["id"]] != p:
                        raise ValueError("Conflicting frozen trajectory")
                    points[p["id"]] = p
            states, live, pending = {}, [], []
            for identifier, p in points.items():
                job = OUT / "jobs" / f"{identifier}.json"
                s = json.loads(job.read_text()) if job.exists() else dict(state="pending", point=p)
                if s["state"] == "running" and alive(s.get("pid"), identifier):
                    live.append((identifier, s))
                elif s["state"] == "running":
                    s.update(state="interrupted", updated=now())
                    atomic_json(job, s)
                if s["state"] in ("pending", "interrupted"):
                    pending.append(p)
                states[identifier] = s
            for identifier, (proc, log) in list(child_handles.items()):
                if proc.poll() is not None:
                    log.close()
                    del child_handles[identifier]
            when = datetime.now(CENTRAL)
            limit = min(config.get("max_workers", cores), worker_limit(when, cores))
            # Pause surplus process groups at 9 a.m. without killing writes.
            # Resume the same in-memory state when slots become available.
            for index, (identifier, s) in enumerate(live):
                should_pause = index >= limit
                if should_pause != s.get("paused", False):
                    os.killpg(s["pid"], signal.SIGSTOP if should_pause else signal.SIGCONT)
                    s.update(paused=should_pause, updated=now())
                    atomic_json(OUT / "jobs" / f"{identifier}.json", s)
                    print(now(), "pause" if should_pause else "resume", identifier, flush=True)
            active = sum(not s.get("paused", False) for _, s in live)
            free = shutil.disk_usage(OUT).free
            remaining_live = sum(max(0, estimated_bytes(points[i]) - (OUT / "chains" / f"{i}.h5").stat().st_size) if (OUT / "chains" / f"{i}.h5").exists() else estimated_bytes(points[i]) for i, _ in live)
            reserve = max(RESERVE, config.get("disk_reserve_bytes", RESERVE))
            permitted = when < DEADLINE and config.get("launch_new", True)
            for p in pending:
                if active >= limit or not permitted:
                    break
                if p["D"] == 5 and not config.get("enable_d5", False):
                    continue
                requirement = estimated_bytes(p)
                if free - remaining_live - requirement < reserve:
                    continue
                identifier = p["id"]
                (OUT / "logs").mkdir(exist_ok=True)
                log = (OUT / "logs" / f"{identifier}.log").open("a")
                proc = subprocess.Popen([sys.executable, "-u", str(OUT / "run.py"), "worker", "--id", identifier],
                                        stdout=log, stderr=subprocess.STDOUT, cwd=ROOT,
                                        start_new_session=True, env=os.environ.copy())
                child_handles[identifier] = (proc, log)
                s = dict(state="running", point=p, pid=proc.pid, paused=False, updated=now())
                atomic_json(OUT / "jobs" / f"{identifier}.json", s)
                states[identifier] = s
                remaining_live += requirement
                active += 1
                print(now(), "launched", identifier, "pid", proc.pid, flush=True)
            summary = dict(updated=now(), pid=os.getpid(), worker_limit=limit,
                           active=active, paused=sum(s.get("paused", False) for _, s in live),
                           total=len(points), complete=sum(s["state"] == "complete" for s in states.values()),
                           failed=[i for i, s in states.items() if s["state"] == "failed"],
                           pending=sum(s["state"] in ("pending", "interrupted") for s in states.values()),
                           disk_free_bytes=free, reserve_bytes=reserve, deadline=DEADLINE.isoformat(),
                           launch_new=permitted, state="deadline" if when >= DEADLINE else "running")
            atomic_json(OUT / "status.json", summary)
            if when >= DEADLINE and not live and not child_handles:
                print(now(), "deadline reached; coordinator finished", flush=True)
                return
            time.sleep(10)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("freeze-calibration", "worker", "coordinate"))
    parser.add_argument("--id")
    args = parser.parse_args()
    if args.mode == "freeze-calibration":
        points = [point(D, 100, beta, N, replica, 4000)
                  for D, betas in [(3, [440, 455, 470]), (4, [290, 302, 314])]
                  for beta in betas for N in [6, 8, 12] for replica in range(2)]
        print(json.dumps(freeze("calibration", points, "Ordinary XY endpoint pilot; estimator calibration and actual resource measurement, not final phase-boundary inference."), indent=2))
    elif args.mode == "worker":
        def terminate(signum, frame):
            raise KeyboardInterrupt("Worker termination requested; resume from committed checkpoint")
        signal.signal(signal.SIGTERM, terminate)
        matches = [p for manifest in (OUT / "batches").glob("*.json")
                   for p in json.loads(manifest.read_text())["points"] if p["id"] == args.id]
        if not matches or any(p != matches[0] for p in matches):
            raise ValueError("Unknown or conflicting worker specification")
        trajectory(matches[0])
    else:
        coordinate()


if __name__ == "__main__":
    main()
