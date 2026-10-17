"""Checkpointed three-dimensional kappa=0 admissible-XY scout."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import sys
from time import perf_counter
import traceback

os.environ["NUMBA_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MPLBACKEND"] = "Agg"
os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/luscher-mpl")

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
sys.path.insert(0, str(ROOT))

import h5py
import numpy as np
from supervillain.xy import _provenance, run_chain


D = 3
KAPPA = 0.0
EPSILON_CENTI = [240, 245, 250, 255, 260]
SIZES = [8, 12, 16]
STEPS = 10_000
REPLICAS = 2
ROOT_SEED = 2026090401
CLUSTERS_PER_CYCLE = 4


def now():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temporary.replace(path)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def point(epsilon_centi, L, replica):
    seed = int(np.random.SeedSequence(
        [ROOT_SEED, epsilon_centi, L, replica]
    ).generate_state(1, dtype=np.uint64)[0])
    initial_winding = [0, 0, 0] if replica == 0 else [1, 0, 0]
    return dict(
        id=f"e{epsilon_centi:03d}_L{L:03d}_r{replica:02d}",
        D=D, L=L, kappa=KAPPA,
        epsilon_centi=epsilon_centi, epsilon=epsilon_centi / 100,
        replica=replica, seed=seed, steps=STEPS,
        initial_winding=initial_winding,
        clusters_per_cycle=CLUSTERS_PER_CYCLE,
    )


def points():
    return [point(e, L, replica) for e in EPSILON_CENTI
            for L in SIZES for replica in range(REPLICAS)]


def principal(delta):
    reduced = np.remainder(delta, 2 * np.pi)
    return np.where((-np.pi <= delta) & (delta < np.pi), delta,
                    np.where(reduced >= np.pi, reduced - 2 * np.pi, reduced))


def audit(path, spec):
    """Recompute every physical observable from every saved three-dimensional field."""
    steps, L = spec["steps"], spec["L"]
    volume = L ** D
    spatial_axes = tuple(range(1, D + 1))
    with h5py.File(path, "r") as f:
        parameters = json.loads(f.attrs["parameters"])
        expected_parameters = dict(D=D, N=L, kappa=KAPPA, epsilon=spec["epsilon"],
                                   seed=spec["seed"], update="hybrid",
                                   clusters_per_step=spec["clusters_per_cycle"],
                                   initial_winding=spec["initial_winding"])
        for key, value in expected_parameters.items():
            assert parameters[key] == value, (key, parameters[key], value)
        assert parameters["winding"] and not parameters["parallel"]
        assert parameters["buffer_rows"] == 100
        assert f.attrs["committed_rows"] == steps + 1
        assert all(len(dataset) == steps + 1 for dataset in f.values())
        np.testing.assert_array_equal(f["sweep"][:], np.arange(steps + 1))

        coordinates = np.indices((L,) * D)
        initial = np.mod(sum(2 * np.pi * winding * coordinates[axis] / L
                             for axis, winding in enumerate(spec["initial_winding"])), 2 * np.pi)
        np.testing.assert_array_equal(f["phi"][0, 0], initial)

        cluster_sizes = f["cluster_sizes"][:]
        assert np.all(cluster_sizes[0] == 0)
        assert np.all((cluster_sizes[1:] >= 1) & (cluster_sizes[1:] <= volume))
        summed_sizes = cluster_sizes.sum(axis=1)
        np.testing.assert_array_equal(
            f["clusters_total"][:], spec["clusters_per_cycle"] * np.arange(steps + 1))
        np.testing.assert_array_equal(f["cluster_sites_total"][:], np.cumsum(summed_sizes))
        np.testing.assert_allclose(
            f["cluster_equivalent_sweeps"][:], np.cumsum(summed_sizes) / volume,
            rtol=0, atol=1e-12)
        np.testing.assert_array_equal(
            f["cluster_bond_tests_total"][:],
            np.cumsum(f["cluster_bond_tests"][:].sum(axis=1)))
        assert np.all((f["local_accepted"][:] >= 0) & (f["local_accepted"][:] <= volume))
        for key in f:
            if key.endswith("_seconds"):
                values = f[key][:]
                assert np.isfinite(values).all() and np.all(values >= 0)
        for state in f["rng_state"][:]:
            assert set(json.loads(state)) == {"local", "cluster", "winding"}

        for start in range(0, steps + 1, 64):
            stop = min(start + 64, steps + 1)
            field = f["phi"][start:stop, 0]
            assert np.isfinite(field).all()
            assert np.all((field >= 0) & (field < 2 * np.pi))
            links = np.stack([
                principal(np.roll(field, -1, axis=axis + 1) - field)
                for axis in range(D)
            ], axis=1)
            assert np.all(np.abs(links) < spec["epsilon"])

            spin = np.exp(1j * field)
            magnetization = spin.mean(axis=spatial_axes)
            m2 = np.abs(magnetization) ** 2
            mean_cosine = np.cos(links).mean(axis=tuple(range(1, D + 2)))

            vortex_sum = np.zeros(stop - start)
            plaquette_orientations = 0
            for mu in range(D):
                for nu in range(mu + 1, D):
                    curl = (links[:, mu]
                            + np.roll(links[:, nu], -1, axis=mu + 1)
                            - np.roll(links[:, mu], -1, axis=nu + 1)
                            - links[:, nu])
                    vortex_sum += np.abs(np.rint(curl / (2 * np.pi))).mean(axis=spatial_axes)
                    plaquette_orientations += 1
            vortex_density = vortex_sum / plaquette_orientations

            winding = np.stack([
                np.rint(links[:, axis].sum(axis=axis + 1) / (2 * np.pi)).reshape(stop - start, -1)
                for axis in range(D)
            ], axis=1)
            structure_min = []
            for axis in range(D):
                shape = [1] * (D + 1)
                shape[axis + 1] = L
                phase = np.exp(-2j * np.pi * np.arange(L) / L).reshape(shape)
                structure_min.append(np.abs((spin * phase).sum(axis=spatial_axes)) ** 2 / volume)
            structure_min = np.stack(structure_min, axis=1)

            expected = dict(
                magnetization_real=magnetization.real,
                magnetization_imag=magnetization.imag,
                magnetization_abs=np.abs(magnetization),
                magnetization_squared=m2,
                magnetization_fourth=m2 * m2,
                mean_link_cosine=mean_cosine,
                action_density=-D * KAPPA * mean_cosine,
                vortex_density=vortex_density,
                winding=winding,
                structure_zero=volume * m2,
                structure_min=structure_min,
                max_link_angle=np.abs(links).max(axis=tuple(range(1, D + 2))),
            )
            for key, value in expected.items():
                np.testing.assert_allclose(f[key][start:stop], value,
                                           rtol=2e-12, atol=3e-13, err_msg=key)
    return dict(passed=True, rows=steps + 1, bytes=path.stat().st_size,
                sha256=sha256(path))


def run_point(spec):
    path = OUT / "chains" / f"{spec['id']}.h5"
    status_path = OUT / "jobs" / f"{spec['id']}.json"
    status = json.loads(status_path.read_text()) if status_path.exists() else dict(point=spec, segments=[])
    if status.get("state") == "complete":
        return status
    with (OUT / "logs" / f"{spec['id']}.log").open("a", buffering=1) as log, \
            redirect_stdout(log), redirect_stderr(log):
        try:
            status.setdefault("started_utc", now())
            status.update(state="running", pid=os.getpid(), updated_utc=now())
            atomic_json(status_path, status)
            done = 0
            if path.exists():
                with h5py.File(path, "r") as f:
                    done = int(f.attrs["committed_rows"]) - 1
            for target in range(1000, spec["steps"] + 1, 1000):
                if done >= target:
                    continue
                before = perf_counter()
                load_before = list(os.getloadavg())
                run_chain(
                    path, D=D, N=spec["L"], kappa=KAPPA, epsilon=spec["epsilon"],
                    sweeps=target, warmup=100, seed=spec["seed"], backend="numba",
                    parallel=False, interval_phi=np.pi, winding=True,
                    checkpoint_every=100, resume=path.exists(), update="hybrid",
                    clusters_per_step=spec["clusters_per_cycle"],
                    initial_winding=spec["initial_winding"], record_timing=True,
                    buffer_rows=100,
                )
                elapsed = perf_counter() - before
                with h5py.File(path, "r") as f:
                    assert f.attrs["committed_rows"] == target + 1
                    digest = hashlib.sha256(
                        f["phi"][-1].tobytes() + f["rng_state"][-1]
                    ).hexdigest()
                status["segments"].append(dict(
                    first=done + 1, last=target, seconds=elapsed,
                    load_before=load_before, load_after=list(os.getloadavg()),
                    field_rng_sha256=digest, finished_utc=now(),
                ))
                done = target
                status.update(cycles=done, updated_utc=now())
                atomic_json(status_path, status)
                print(f"{spec['id']}: {done}/{spec['steps']}, segment {elapsed:.3f}s", flush=True)
            status.update(state="auditing", updated_utc=now())
            atomic_json(status_path, status)
            before = perf_counter()
            result = audit(path, spec)
            status.update(
                state="complete", audit=result,
                audit_seconds=perf_counter() - before,
                generation_seconds=sum(segment["seconds"] for segment in status["segments"]),
                updated_utc=now(),
            )
            status.pop("error", None)
            atomic_json(status_path, status)
            print(f"{spec['id']}: complete; full-field 3D audit passed", flush=True)
            return status
        except BaseException:
            status.update(state="failed", error=traceback.format_exc(), updated_utc=now())
            atomic_json(status_path, status)
            traceback.print_exc()
            raise


def execute(batch, workers=4):
    for folder in ("chains", "jobs", "logs", "batches"):
        (OUT / folder).mkdir(exist_ok=True)
    with (OUT / "runner.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        manifest_path = OUT / "batches" / f"{batch}.json"
        manifest = json.loads(manifest_path.read_text())
        assert manifest["model_source_sha256"] == _provenance()["source_sha256"]
        assert manifest["runner_sha256"] == sha256(Path(__file__))
        all_points = manifest["points"]
        pending, completed = [], 0
        for spec in all_points:
            job = OUT / "jobs" / f"{spec['id']}.json"
            if job.exists() and json.loads(job.read_text()).get("state") == "complete":
                completed += 1
            else:
                pending.append(spec)
        projected = sum((spec["steps"] + 1) * (8 * spec["L"] ** D + 3500)
                        for spec in pending)
        if shutil.disk_usage(OUT).free < projected + 25 * 2**30:
            raise RuntimeError("Insufficient space for raw records plus 25 GiB reserve")
        runtime = dict(
            batch=batch, started_utc=now(), pid=os.getpid(), workers=workers,
            completed=completed, total=len(all_points), scheduled_this_launch=len(pending),
            failures=[], state="running", load_at_start=list(os.getloadavg()),
            disk_free_bytes=shutil.disk_usage(OUT).free, provenance=_provenance(),
        )
        runtime_path = OUT / f"{batch}-status.json"
        atomic_json(runtime_path, runtime)
        print(f"Starting {len(pending)} pending chains; {completed} already complete", flush=True)
        with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as pool:
            futures = {pool.submit(run_point, spec): spec for spec in pending}
            for future in as_completed(futures):
                spec = futures[future]
                try:
                    result = future.result()
                    runtime["completed"] += 1
                    print(f"{runtime['completed']}/{len(all_points)} {spec['id']}: "
                          f"generation {result['generation_seconds']:.1f}s; audit passed", flush=True)
                except Exception as error:
                    runtime["failures"].append(dict(id=spec["id"], error=repr(error)))
                    print(f"FAILED {spec['id']}: {error!r}", flush=True)
                runtime.update(updated_utc=now(), load=list(os.getloadavg()),
                               disk_free_bytes=shutil.disk_usage(OUT).free)
                atomic_json(runtime_path, runtime)
        runtime.update(finished_utc=now(),
                       state="complete" if runtime["completed"] == len(all_points) else "partial")
        atomic_json(runtime_path, runtime)
        if runtime["failures"]:
            raise RuntimeError(f"{len(runtime['failures'])} failed jobs")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batch")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if not 1 <= args.workers <= 8:
        raise ValueError("workers must be 1..8")
    execute(args.batch, args.workers)


if __name__ == "__main__":
    main()
