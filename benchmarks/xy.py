#!/usr/bin/env python
"""Reproducible, bounded timing check for the constrained XY site update.

This is a microbenchmark only.  It measures the update implementation on a
fixed random tape; it does not produce physics measurements or an MC stream.
Each finished repeat is atomically checkpointed in ``docs/xy-benchmark.json``.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import time

import numba
import numpy as np

# Permit ``.venv/bin/python benchmarks/xy.py`` from any current directory.
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from supervillain.action.xy import XY
from supervillain.generator.xy import SiteUpdate
from supervillain.lattice import Lattice


DEFAULT_OUTPUT = ROOT / "docs" / "xy-benchmark.json"
SOURCE_FILES = (
    "benchmarks/xy.py",
    "supervillain/action/xy.py",
    "supervillain/generator/xy.py",
    "supervillain/generator/_xy_kernels.py",
    "supervillain/generator/reference_implementation/xy.py",
)


def _digest(field):
    raw = np.ascontiguousarray(np.asarray(field)).view(np.uint8)
    return hashlib.sha256(raw).hexdigest()


def _source_digests():
    return {
        relative: hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        for relative in SOURCE_FILES
    }


def _atomic_write(path, record):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent,
                                     delete=False) as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = Path(handle.name)
    os.replace(temporary, path)


def _tapes(volume, sweeps, seed):
    rng = np.random.default_rng(seed)
    return [
        (rng.uniform(-0.2, 0.2, volume), rng.uniform(0.0, 1.0, volume))
        for _ in range(sweeps)
    ]


def _run(update, initial, tapes):
    field = initial
    started = time.perf_counter()
    for proposals, uniforms in tapes:
        field = update.sweep(field, proposals, uniforms)
    return {
        "seconds": time.perf_counter() - started,
        "accepted": int(update.accepted),
        "digest": _digest(field),
        "load_average": list(os.getloadavg()),
    }


def _label(implementation):
    backend, parallel = implementation
    return f"{backend}_{'parallel' if parallel else 'serial'}"


def _case(record, *, name, D, N, implementations, repeats, sweeps):
    lattice = Lattice(D, N)
    action = XY(lattice, kappa=0.43, epsilon=np.pi)
    initial = action.canonical(np.zeros(lattice.dims))
    tape_seed = 20260902 + 100 * D + N
    tapes = _tapes(lattice.sites, sweeps, seed=tape_seed)
    case = {
        "name": name,
        "D": D,
        "N": N,
        "volume": lattice.sites,
        "sweeps_per_repeat": sweeps,
        "action": {"kappa": action.kappa, "epsilon": action.epsilon},
        "random_tape": {
            "seed": tape_seed,
            "proposal_distribution": "uniform(-0.2, 0.2)",
            "uniform_distribution": "uniform(0.0, 1.0)",
        },
        "start_digest": _digest(initial),
        "implementations": {_label(implementation): [] for implementation in implementations},
    }
    record["cases"].append(case)
    _atomic_write(Path(record["output_path"]), record)

    # Compile outside all timed regions, including the parallel reduction.
    for implementation in implementations:
        SiteUpdate(action, backend=implementation[0], parallel=implementation[1]).sweep(
            initial, *tapes[0]
        )

    reference_digest = None
    for implementation in implementations:
        backend, parallel = implementation
        label = _label(implementation)
        for repeat in range(repeats):
            row = _run(
                SiteUpdate(action, backend=backend, parallel=parallel), initial, tapes
            )
            row["repeat"] = repeat + 1
            case["implementations"][label].append(row)
            if reference_digest is None:
                reference_digest = row["digest"]
            elif row["digest"] != reference_digest:
                raise RuntimeError(
                    f"fixed-tape parity mismatch in {name}: {implementation} repeat {repeat + 1}"
                )
            _atomic_write(Path(record["output_path"]), record)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--small-sweeps", type=int, default=1)
    parser.add_argument("--large-sweeps", type=int, default=20)
    args = parser.parse_args()
    if args.repeats < 1 or args.small_sweeps < 1 or args.large_sweeps < 1:
        raise ValueError("repeat and sweep counts must be positive")

    output = args.output.resolve()
    revision = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True,
        capture_output=True, check=False,
    ).stdout.strip()
    record = {
        "purpose": "warmed fixed-tape microbenchmark; not a production Monte Carlo run",
        "output_path": str(output),
        "script": str(Path(__file__).resolve()),
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "command": sys.argv,
        "git_revision": revision or None,
        "source_sha256": _source_digests(),
        "environment": {
            "python": sys.version,
            "numpy": np.__version__,
            "numba": numba.__version__,
            "numba_threads": numba.get_num_threads(),
            "machine": platform.machine(),
            "load_average": list(os.getloadavg()),
        },
        "cases": [],
    }
    _case(
        record, name="oracle-parity-d2", D=2, N=48,
        implementations=(("python", False), ("numba", False), ("numba", True)),
        repeats=args.repeats, sweeps=args.small_sweeps,
    )
    _case(
        record, name="oracle-parity-d3", D=3, N=12,
        implementations=(("python", False), ("numba", False), ("numba", True)),
        repeats=args.repeats, sweeps=args.small_sweeps,
    )
    _case(
        record, name="large-d3-serial-parallel", D=3, N=32,
        implementations=(("numba", False), ("numba", True)),
        repeats=args.repeats, sweeps=args.large_sweeps,
    )


if __name__ == "__main__":
    main()
