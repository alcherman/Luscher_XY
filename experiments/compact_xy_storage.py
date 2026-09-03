"""Compact completed XY trajectories without changing observable time series.

This is a destructive archival migration requested on 2026-09-03.  Each source
file is replaced only after the compact temporary file passes a value-by-value
comparison.  Progress is append-only and resumable.
"""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
from datetime import datetime, timezone

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parent
LOG = ROOT / "storage-compaction.jsonl"
PILOT_MANIFEST = ROOT / "2026-09-02-xy-bkt-refine" / "pilot-deletion-manifest.json"
SPECIAL_SWEEPS = {0, 100, 101, 500, 501, 1000, 1001, 2000, 2001, 4000, 4001}
DENSE_RELATIVE_PATHS = {
    "2026-09-02-xy-bkt-scout/chains/e0500_k0000_L016.h5",
    "2026-09-02-xy-bkt-refine/chains/e0565_k0000_L032_r000.h5",
    "2026-09-02-xy-bkt-refine/chains/e1000_k1120_L064_r000.h5",
    "2026-09-03-xy-bkt-low-kappa/chains/e0580_k0180_L064_r000.h5",
    "2026-09-03-xy-bkt-kappa-zero/chains/e0550_k0000_L064_r000.h5",
}


def now():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path: Path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    with temporary.open("rb") as stream:
        os.fsync(stream.fileno())
    temporary.replace(path)


def sha256(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def append_event(event):
    with LOG.open("a") as stream:
        stream.write(json.dumps(event, sort_keys=True) + "\n")
        stream.flush()
        os.fsync(stream.fileno())


def completed_paths():
    if not LOG.exists():
        return set()
    done = set()
    for line in LOG.read_text().splitlines():
        event = json.loads(line)
        if event.get("event") == "compacted":
            done.add(event["path"])
    return done


def acquire_runner_locks():
    handles = []
    for path in sorted(ROOT.glob("*/runner.lock")) + sorted(ROOT.glob("*/edges*-runner.lock")):
        handle = path.open("a+")
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(f"An experiment runner holds {path}") from error
        handles.append(handle)
    return handles


def job_for(path: Path):
    candidate = path.parent.parent / "jobs" / f"{path.stem}.json"
    return json.loads(candidate.read_text()) if candidate.exists() else None


def original_digest(path: Path):
    job = job_for(path)
    if job and job.get("state") == "complete" and job.get("audit", {}).get("passed"):
        audit = job["audit"]
        if audit["bytes"] != path.stat().st_size:
            raise RuntimeError(f"Size no longer matches closure audit: {path}")
        return audit["sha256"], "closure_audit"
    return sha256(path), "computed_during_compaction"


def equal(a, b):
    return np.array_equal(np.asarray(a), np.asarray(b))


def assert_attrs_equal(source, target):
    assert set(source.attrs) <= set(target.attrs)
    for key in source.attrs:
        if not equal(source.attrs[key], target.attrs[key]):
            raise AssertionError(f"Attribute mismatch: {source.name}:{key}")


def create_dataset(target, name, data, source_dataset, sparse_sweeps=None):
    data = np.asarray(data)
    kwargs = {}
    if data.ndim and data.dtype.kind != "O" and data.size:
        kwargs.update(compression="gzip", compression_opts=4, shuffle=True, fletcher32=True)
        rows = min(100, len(data))
        kwargs["chunks"] = (rows,) + data.shape[1:]
    dataset = target.create_dataset(name, data=data, dtype=source_dataset.dtype, **kwargs)
    for key, value in source_dataset.attrs.items():
        dataset.attrs[key] = value
    if sparse_sweeps is not None:
        dataset.attrs["storage_compaction"] = "sparse_exact_rows"
        dataset.attrs["retained_sweeps"] = np.asarray(sparse_sweeps, dtype=np.int64)
        dataset.attrs["original_rows"] = source_dataset.shape[0]
    return dataset


def copy_compact(source_path: Path, target_path: Path, dense: bool, old_digest: str):
    with h5py.File(source_path, "r") as source, h5py.File(target_path, "w") as target:
        for key, value in source.attrs.items():
            target.attrs[key] = value
        sweeps = np.asarray(source["sweep"])
        if dense:
            indices = np.arange(len(sweeps))
        else:
            keep = (sweeps % 100 == 0) | np.isin(sweeps, list(SPECIAL_SWEEPS))
            keep[0] = True
            keep[-1] = True
            indices = np.flatnonzero(keep)
        selected_sweeps = sweeps[indices]

        for name, dataset in source.items():
            if name in {"phi", "rng_state"}:
                values = dataset[:] if dense else dataset[indices]
                create_dataset(target, name, values, dataset,
                               None if dense else selected_sweeps)
            else:
                create_dataset(target, name, dataset[:], dataset)

        target.attrs["storage_compaction_version"] = 1
        target.attrs["storage_compaction_utc"] = now()
        target.attrs["original_file_sha256"] = old_digest
        target.attrs["original_file_bytes"] = source_path.stat().st_size
        target.attrs["field_retention"] = (
            "complete representative trajectory" if dense else
            "exact phi and RNG at cycles divisible by 100, cut-boundary candidates, and endpoints"
        )
        target.flush()

    # Re-open both files and compare all retained values and original metadata.
    with h5py.File(source_path, "r") as source, h5py.File(target_path, "r") as target:
        assert_attrs_equal(source, target)
        assert set(source) == set(target)
        sweeps = np.asarray(source["sweep"])
        for name, source_dataset in source.items():
            target_dataset = target[name]
            assert source_dataset.dtype == target_dataset.dtype
            assert_attrs_equal(source_dataset, target_dataset)
            if name in {"phi", "rng_state"} and not dense:
                retained = np.asarray(target_dataset.attrs["retained_sweeps"])
                indices = np.flatnonzero(np.isin(sweeps, retained))
                if not equal(source_dataset[indices], target_dataset[:]):
                    raise AssertionError(f"Sparse value mismatch: {source_path}:{name}")
            elif not equal(source_dataset[:], target_dataset[:]):
                raise AssertionError(f"Value mismatch: {source_path}:{name}")


def pilot_paths():
    campaign = ROOT / "2026-09-02-xy-bkt-refine"
    production = campaign / "batches" / "production.json"
    choices = json.loads(production.read_text())["selection"]["choices"]
    choice_by_id = {choice["pilot_id"]: choice for choice in choices}
    result = []
    for job_path in sorted((campaign / "jobs").glob("*.json")):
        job = json.loads(job_path.read_text())
        point = job.get("point", {})
        if point.get("replica", 0) < 100:
            continue
        chain = campaign / "chains" / f"{point['id']}.h5"
        if chain.exists():
            analysis = campaign / "analysis" / f"{point['id']}.json"
            bootstrap = campaign / "analysis" / f"{point['id']}_bootstrap.npz"
            choice = choice_by_id.get(point["id"])
            if not (choice and "pilot_U" in choice and "center" in choice
                    and job.get("audit", {}).get("passed")
                    and (campaign / "logs" / f"{point['id']}.log").exists()):
                raise RuntimeError(f"Pilot archival products incomplete: {point['id']}")
            actual = sha256(chain)
            if actual != job["audit"]["sha256"]:
                raise RuntimeError(f"Pilot hash mismatch: {chain}")
            result.append(dict(path=str(chain.relative_to(ROOT)), bytes=chain.stat().st_size,
                               sha256=actual, job=str(job_path.relative_to(ROOT)),
                               selection_record=str(production.relative_to(ROOT)),
                               pilot_U=choice["pilot_U"], selection_center=choice["center"],
                               analysis=str(analysis.relative_to(ROOT)) if analysis.exists() else None,
                               bootstrap=str(bootstrap.relative_to(ROOT)) if bootstrap.exists() else None))
    return result


def delete_pilots():
    if PILOT_MANIFEST.exists():
        record = json.loads(PILOT_MANIFEST.read_text())
        for item in record["files"]:
            path = ROOT / item["path"]
            if path.exists():
                if sha256(path) != item["sha256"]:
                    raise RuntimeError(f"Pilot changed after manifest: {path}")
                path.unlink()
        record["state"] = "complete"
        record.setdefault("completed_utc", now())
        atomic_json(PILOT_MANIFEST, record)
        return record

    files = pilot_paths()
    if len(files) != 22:
        raise RuntimeError(f"Expected 22 complete pilot chains, found {len(files)}")
    record = dict(state="prepared", created_utc=now(), files=files,
                  bytes=sum(item["bytes"] for item in files),
                  reason="Raw design-pilot trajectories removed at user request; manifests, jobs, logs, analysis summaries and bootstrap samples retained.")
    atomic_json(PILOT_MANIFEST, record)
    for item in files:
        (ROOT / item["path"]).unlink()
    record["state"] = "complete"
    record["completed_utc"] = now()
    atomic_json(PILOT_MANIFEST, record)
    return record


def compact_all():
    dense_missing = [path for path in DENSE_RELATIVE_PATHS if not (ROOT / path).exists()]
    if dense_missing:
        raise FileNotFoundError(dense_missing)
    done = completed_paths()
    paths = sorted(ROOT.glob("*/chains/*.h5"))
    paths = [path for path in paths if path.name not in {""}]
    print(f"Found {len(paths)} surviving trajectories; {len(done)} already compacted", flush=True)
    for number, path in enumerate(paths, 1):
        relative = str(path.relative_to(ROOT))
        if relative in done:
            continue
        # A crash can occur after the verified atomic replacement but before
        # its append-only completion event.  The replacement carries enough
        # provenance to recognize and record that state without touching it.
        with h5py.File(path, "r") as possible_compact:
            compact_version = possible_compact.attrs.get("storage_compaction_version")
            if compact_version == 1:
                event = dict(event="compacted", utc=now(), path=relative,
                             dense_field_history=relative in DENSE_RELATIVE_PATHS,
                             original_bytes=int(possible_compact.attrs["original_file_bytes"]),
                             compact_bytes=path.stat().st_size,
                             original_sha256=str(possible_compact.attrs["original_file_sha256"]),
                             original_digest_source="recovered_from_verified_replacement",
                             compact_sha256=sha256(path), recovered_after_interruption=True)
                append_event(event)
                done.add(relative)
                print(f"{number}/{len(paths)} recovered completed replacement {relative}", flush=True)
                continue
        temporary = path.with_name(f".{path.name}.compact.tmp")
        if temporary.exists():
            temporary.unlink()
        old_bytes = path.stat().st_size
        old_digest, digest_source = original_digest(path)
        dense = relative in DENSE_RELATIVE_PATHS
        copy_compact(path, temporary, dense, old_digest)
        new_digest = sha256(temporary)
        new_bytes = temporary.stat().st_size
        with temporary.open("rb") as stream:
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        directory_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
        append_event(dict(event="compacted", utc=now(), path=relative,
                          dense_field_history=dense, original_bytes=old_bytes,
                          compact_bytes=new_bytes, original_sha256=old_digest,
                          original_digest_source=digest_source,
                          compact_sha256=new_digest))
        done.add(relative)
        print(f"{number}/{len(paths)} {relative}: {old_bytes/2**20:.1f} -> {new_bytes/2**20:.1f} MiB", flush=True)


def self_test():
    import tempfile
    source_path = next((ROOT / "2026-09-02-xy-binder/chains").glob("*.h5"))
    with tempfile.TemporaryDirectory() as folder:
        target = Path(folder) / "compact.h5"
        digest = sha256(source_path)
        copy_compact(source_path, target, False, digest)
        with h5py.File(source_path, "r") as source, h5py.File(target, "r") as compact:
            assert len(compact["phi"]) < len(source["phi"])
            assert len(compact["magnetization_squared"]) == len(source["magnetization_squared"])
            assert compact.attrs["original_file_sha256"] == digest
    print("Self-test passed", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    lock = (ROOT / "storage-compaction.lock").open("a+")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    runner_locks = acquire_runner_locks()
    if not LOG.exists():
        append_event(dict(event="started", utc=now(), version=1,
                          field_rule="cycles divisible by 100, cut-boundary candidates, endpoints",
                          dense_representatives=sorted(DENSE_RELATIVE_PATHS)))
    pilot_record = delete_pilots()
    print(f"Pilot removal recorded: {len(pilot_record['files'])} files, {pilot_record['bytes']/2**30:.3f} GiB", flush=True)
    compact_all()
    events=[json.loads(line) for line in LOG.read_text().splitlines()]
    rows=[event for event in events if event.get('event')=='compacted']
    atomic_json(ROOT/'storage-compaction-summary.json',dict(state='complete',completed_utc=now(),files=len(rows),
        dense_representatives=sum(row['dense_field_history'] for row in rows),
        original_bytes=sum(row['original_bytes'] for row in rows),compact_bytes=sum(row['compact_bytes'] for row in rows),
        pilot_bytes_removed=pilot_record['bytes'],log_sha256=sha256(LOG)))
    print("Compaction complete", flush=True)
    _ = runner_locks  # keep lock handles alive through completion


if __name__ == "__main__":
    main()
