"""Measurements and a resumable parameter scan for the admissible XY model.

Run ``python -m supervillain.xy --help``. Every sweep, including warmup and
the initial configuration, is retained. No automatic equilibrium claim is made.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
from time import perf_counter

import h5py
import numba
import numpy as np

from supervillain.action.xy import XY
from supervillain.generator.xy import ClusterUpdate, SiteUpdate, WindingUpdate
from supervillain.lattice import Lattice


def measure(action, phi):
    """Raw moments, Fourier modes, plaquette vortices and all axial windings.

    ``winding[a]`` contains the circulation / 2pi on every straight periodic
    loop in direction a. In the presence of vortices these can differ between
    parallel loops. ``structure_min[a]`` is |sum exp(i phi-i p_a x)|²/V.
    Susceptibilities and Binder ratios must be formed from ensemble averages.
    """
    phi = action.canonical(phi)
    links = np.asarray(action.links(phi))
    if not np.all(np.abs(links) < action.epsilon):
        raise ValueError("Cannot measure an inadmissible configuration")
    L = action.Lattice
    spin = np.exp(1j * np.asarray(phi)[0])
    magnetization = spin.mean()
    m2 = float(abs(magnetization) ** 2)
    minimum = []
    winding = []
    for a in range(L.D):
        phase_shape = [1] * L.D
        phase_shape[a] = L.N
        phase = np.exp(-2j * np.pi * np.arange(L.N) / L.N).reshape(phase_shape)
        minimum.append(float(abs(np.sum(spin * phase)) ** 2 / L.sites))
        winding.append(np.rint(links[a].sum(axis=a) / (2 * np.pi)).astype(np.int64).ravel())
    vortices = 0.0
    plaquettes = 0
    for a in range(L.D):
        for b in range(a + 1, L.D):
            curl = (links[a] + np.roll(links[b], -1, axis=a)
                    - np.roll(links[a], -1, axis=b) - links[b])
            charge = np.rint(curl / (2 * np.pi))
            vortices += np.abs(charge).sum()
            plaquettes += charge.size
    mean_cos = float(np.cos(links).mean())
    return dict(
        action_density=-action.kappa * L.D * mean_cos,
        mean_link_cosine=mean_cos,
        magnetization_real=float(magnetization.real),
        magnetization_imag=float(magnetization.imag),
        magnetization_abs=float(abs(magnetization)),
        magnetization_squared=m2,
        magnetization_fourth=m2 * m2,
        structure_zero=L.sites * m2,
        structure_min=np.asarray(minimum),
        max_link_angle=float(np.abs(links).max()),
        vortex_density=float(vortices / plaquettes) if plaquettes else 0.0,
        winding=np.asarray(winding),
    )


def _integer(value, name, minimum):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def _provenance():
    root = Path(__file__).resolve().parents[1]
    try:
        revision = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True, stderr=subprocess.DEVNULL).strip()
        dirty = bool(subprocess.check_output(
            ["git", "status", "--porcelain"], cwd=root, text=True, stderr=subprocess.DEVNULL).strip())
    except (OSError, subprocess.CalledProcessError):
        revision, dirty = "unknown", True
    sources = [Path(__file__), root / "supervillain/action/xy.py",
               root / "supervillain/generator/xy.py", root / "supervillain/generator/_xy_kernels.py",
               root / "supervillain/generator/reference_implementation/xy.py",
               root / "supervillain/generator/_xy_cluster.py",
               root / "supervillain/generator/reference_implementation/xy_cluster.py"]
    digest = hashlib.sha256()
    for source in sources:
        digest.update(str(source.relative_to(root)).encode())
        digest.update(source.read_bytes())
    return dict(python=platform.python_version(), numpy=np.__version__, numba=numba.__version__,
                h5py=h5py.__version__, git_revision=revision, git_dirty=dirty,
                source_sha256=digest.hexdigest(), threads=numba.get_num_threads(),
                created_utc=datetime.now(timezone.utc).isoformat())


def run_chain(path, *, D, N, kappa, epsilon=np.pi, sweeps=1000, warmup=100,
              seed=0, backend="numba", parallel=True, interval_phi=None,
              winding=True, initial_winding=None, checkpoint_every=10, resume=False,
              update="metropolis", clusters_per_step=1, record_timing=False,
              buffer_rows=1):
    """Write/extend one independently seeded chain, returning its output path.

    ``sweeps`` is the total target number, including warmup; row 0 is the seed.
    ``warmup`` labels an initial analysis cut and does not delete data or certify
    thermalization. Resume restores the RNGs at the last committed checkpoint,
    discards only uncommitted rows, and requires matching model/run parameters
    and source/runtime versions. A crash can lose at most checkpoint_every
    newly computed sweeps; their random stream is replayed on resume.

    With ``update="cluster"``, each saved cycle contains the fixed number
    ``clusters_per_step`` of Wolff flips, followed by an optional winding move.
    ``sweeps`` and ``warmup`` then count these cycles, not volume sweeps.
    Every cluster size and bond-test count is retained. Cluster growth is
    serial; ``parallel`` applies only to Metropolis and is recorded as false
    for clusters. Size-dependent measurement stopping is deliberately avoided.
    ``update="hybrid"`` performs one Metropolis sweep before the fixed cluster
    batch. ``buffer_rows`` batches disk writes without thinning any records;
    every checkpoint flushes the pending buffer before advancing its marker.
    Optional per-cycle wall timings are diagnostics, not replayable RNG state;
    ``cycle_seconds`` excludes checkpoint flushes and buffered disk writes.
    """
    sweeps = _integer(sweeps, "sweeps", 1)
    warmup = _integer(warmup, "warmup", 0)
    checkpoint_every = _integer(checkpoint_every, "checkpoint_every", 1)
    buffer_rows = _integer(buffer_rows, "buffer_rows", 1)
    seed = _integer(seed, "seed", 0)
    if update not in {"metropolis", "cluster", "hybrid"}:
        raise ValueError("update must be 'metropolis', 'cluster', or 'hybrid'")
    has_metropolis = update in {"metropolis", "hybrid"}
    has_cluster = update in {"cluster", "hybrid"}
    clusters_per_step = _integer(clusters_per_step, "clusters_per_step", 1)
    if update == "metropolis" and clusters_per_step != 1:
        raise ValueError("clusters_per_step is only applicable to cluster updates")
    if update == "cluster" and interval_phi is not None:
        raise ValueError("interval_phi is only applicable to Metropolis updates")
    if warmup >= sweeps:
        raise ValueError("warmup must be smaller than sweeps")
    action = XY(Lattice(D=D, N=N), kappa, epsilon)
    initial_winding = [0] * D if initial_winding is None else list(initial_winding)
    if len(initial_winding) != D or any(isinstance(w, (bool, np.bool_)) or not isinstance(w, (int, np.integer)) for w in initial_winding):
        raise ValueError("initial_winding must contain D integers")
    if any(abs(2 * np.pi * int(w) / N) >= action.epsilon for w in initial_winding):
        raise ValueError("The requested initial winding is not admissible")
    # Prepare and validate the seed before creating a file. The explicit slope
    # bound also excludes windings that would alias to a different integer on
    # this finite lattice (for example w=N producing an apparent cold field).
    cfg = {"phi": action.Lattice.zeros(0)}
    coords = np.indices(action.Lattice.dims)
    cfg["phi"][0] = np.mod(sum(2 * np.pi * w * coords[a] / N for a, w in enumerate(initial_winding)), 2 * np.pi)
    if not action.valid(cfg):
        raise ValueError("The requested initial winding is not admissible")
    streams = np.random.SeedSequence(seed).spawn(3 if update == "hybrid" else 2)
    if update == "cluster":
        local = ClusterUpdate(action, seed=streams[0], backend=backend)
        parallel = False
    else:
        local = SiteUpdate(action, interval_phi=interval_phi, seed=streams[0], backend=backend, parallel=parallel)
    cluster = (ClusterUpdate(action, seed=streams[2], backend=backend) if update == "hybrid"
               else local if update == "cluster" else None)
    twist = WindingUpdate(action, seed=streams[1]) if winding else None
    parameters = dict(D=int(D), N=int(N), kappa=float(kappa), epsilon=float(epsilon),
                      warmup=warmup, seed=seed, backend=backend, parallel=bool(parallel),
                      interval_phi=float(local.interval_phi) if has_metropolis else None,
                      update=update, clusters_per_step=clusters_per_step,
                      time_unit={"cluster":"cluster_cycle", "hybrid":"hybrid_cycle",
                                 "metropolis":"metropolis_sweep"}[update],
                      record_timing=bool(record_timing),
                      buffer_rows=buffer_rows,
                      winding=bool(winding),
                      initial_winding=[int(w) for w in initial_winding])
    metadata = _provenance()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if resume and not path.exists():
        raise FileNotFoundError(path)
    with h5py.File(path, "r+" if resume else "x") as f:
        if resume:
            if f.attrs.get("schema") != "supervillain.xy.v2" or json.loads(f.attrs["parameters"]) != parameters:
                raise ValueError("Resume parameters do not match this XY chain")
            original = json.loads(f.attrs["provenance"])
            for key in ("source_sha256", "python", "numpy", "numba", "h5py"):
                if original.get(key) != metadata[key]:
                    raise ValueError(f"Resume provenance mismatch: {key}")
            committed = int(f.attrs["committed_rows"])
            if committed < 1 or any(ds.shape[0] < committed for ds in f.values()):
                raise ValueError("Incomplete or corrupt checkpoint")
            if committed > sweeps + 1:
                raise ValueError("Requested sweeps precede the saved checkpoint")
            for ds in f.values():
                ds.resize(committed, axis=0)
            cfg = {"phi": action.canonical(f["phi"][committed - 1])}
            state = json.loads(f["rng_state"][committed - 1])
            local.rng.bit_generator.state = state["local"]
            if update == "hybrid":
                cluster.rng.bit_generator.state = state["cluster"]
            if has_cluster:
                cluster.clusters = int(f["clusters_total"][committed - 1])
                cluster.flipped_sites = int(f["cluster_sites_total"][committed - 1])
                cluster.bond_tests = int(f["cluster_bond_tests_total"][committed - 1])
            if twist is not None:
                twist.rng.bit_generator.state = state["winding"]
            first_sweep = committed
        else:
            f.attrs.update(schema="supervillain.xy.v2", parameters=json.dumps(parameters),
                           provenance=json.dumps(metadata), committed_rows=0,
                           analysis_start_row=warmup + 1)
            first_sweep = 1

        pending = []

        def flush_pending():
            if not pending:
                return
            first = pending[0]["sweep"]
            stop = pending[-1]["sweep"] + 1
            for name in pending[0]:
                values = [row[name] for row in pending]
                arr = np.asarray(values, dtype=object if name == "rng_state" else None)
                if name not in f:
                    dtype = h5py.string_dtype("utf-8") if name == "rng_state" else arr.dtype
                    f.create_dataset(name, shape=(0,) + arr.shape[1:],
                                     maxshape=(None,) + arr.shape[1:], dtype=dtype, chunks=True)
                ds = f[name]
                ds.resize(stop, axis=0)
                ds[first:stop] = arr
            pending.clear()

        def append(sweep, accepted, twist_accepted, cluster_sizes=None, cluster_tests=None,
                   timings=None, cycle_start=None):
            measure_start = perf_counter() if record_timing else None
            record = measure(action, cfg["phi"])
            if record_timing:
                record.update(timings or dict(metropolis_seconds=0.0, cluster_seconds=0.0,
                                              winding_seconds=0.0, update_seconds=0.0))
                record["measurement_seconds"] = perf_counter() - measure_start
                record["cycle_seconds"] = 0.0
            state = {"local": local.rng.bit_generator.state}
            if update == "hybrid":
                state["cluster"] = cluster.rng.bit_generator.state
            if twist is not None:
                state["winding"] = twist.rng.bit_generator.state
            record.update(phi=np.asarray(cfg["phi"]), sweep=sweep,
                          winding_accepted=twist_accepted,
                          rng_state=json.dumps(state))
            if has_metropolis:
                record["local_accepted"] = accepted
            if has_cluster:
                record.update(
                    cluster_sizes=np.zeros(clusters_per_step, dtype=np.int64) if cluster_sizes is None else cluster_sizes,
                    cluster_bond_tests=np.zeros(clusters_per_step, dtype=np.int64) if cluster_tests is None else cluster_tests,
                    clusters_total=cluster.clusters, cluster_sites_total=cluster.flipped_sites,
                    cluster_bond_tests_total=cluster.bond_tests,
                    cluster_equivalent_sweeps=cluster.equivalent_sweeps,
                )
            if buffer_rows > 1:
                # Detach array data from generator-owned storage before the next update.
                record = {name: value.copy() if isinstance(value, np.ndarray) else value
                          for name, value in record.items()}
                if record_timing and cycle_start is not None:
                    record["cycle_seconds"] = perf_counter() - cycle_start
                pending.append(record)
                if len(pending) >= buffer_rows:
                    flush_pending()
                return
            for name, value in record.items():
                if name not in f:
                    arr = np.asarray(value)
                    dtype = h5py.string_dtype("utf-8") if name == "rng_state" else arr.dtype
                    f.create_dataset(name, shape=(0,) + arr.shape, maxshape=(None,) + arr.shape,
                                     dtype=dtype, chunks=True)
                ds = f[name]
                ds.resize(sweep + 1, axis=0)
                ds[sweep] = value
            if record_timing and cycle_start is not None:
                f["cycle_seconds"][sweep] = perf_counter() - cycle_start

        def commit():
            flush_pending()
            f.flush()
            f.attrs.modify("committed_rows", f["sweep"].shape[0])
            f.flush()

        if not resume:
            append(0, 0, 0)
            commit()
        for sweep in range(first_sweep, sweeps + 1):
            cycle_start = perf_counter() if record_timing else None
            mark = cycle_start
            timings = {} if record_timing else None
            sizes = tests = None
            accepted = 0
            if has_metropolis:
                before = local.accepted
                cfg = local.step(cfg)
                accepted = local.accepted - before
            if record_timing:
                now = perf_counter()
                timings["metropolis_seconds"] = now - mark if has_metropolis else 0.0
                mark = now
            if has_cluster:
                sizes = np.empty(clusters_per_step, dtype=np.int64)
                tests = np.empty(clusters_per_step, dtype=np.int64)
                for index in range(clusters_per_step):
                    cfg = cluster.step(cfg)
                    sizes[index] = cluster.last_cluster_size
                    tests[index] = cluster.last_bond_tests
            if record_timing:
                now = perf_counter()
                timings["cluster_seconds"] = now - mark if has_cluster else 0.0
                mark = now
            twist_accepted = 0
            if twist is not None:
                before_twist = twist.accepted
                cfg = twist.step(cfg)
                twist_accepted = twist.accepted - before_twist
            if record_timing:
                now = perf_counter()
                timings["winding_seconds"] = now - mark if twist is not None else 0.0
                timings["update_seconds"] = now - cycle_start
            append(sweep, accepted, twist_accepted, sizes, tests, timings, cycle_start)
            if sweep % checkpoint_every == 0:
                commit()
        commit()
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--D", type=int, required=True)
    parser.add_argument("--N", type=int, nargs="+", required=True)
    parser.add_argument("--kappa", type=float, nargs="+", required=True)
    parser.add_argument("--epsilon", type=float, nargs="+", default=[np.pi], help="radians, 0 < epsilon <= pi")
    parser.add_argument("--sweeps", type=int, default=1000, help="saved updates including warmup; fixed cluster cycles with --update cluster")
    parser.add_argument("--warmup", type=int, default=100)
    parser.add_argument("--replicas", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--backend", choices=("python", "numba"), default="numba")
    parser.add_argument("--update", choices=("metropolis", "cluster", "hybrid"), default="metropolis")
    parser.add_argument("--record-timing", action="store_true", help="record per-cycle wall timings")
    parser.add_argument("--buffer-rows", type=int, default=1, help="batch writes while retaining every cycle; checkpoints flush the buffer")
    parser.add_argument("--clusters-per-step", type=int, default=1, help="fixed cluster flips per saved cycle (cluster updates only)")
    parser.add_argument("--serial", action="store_true")
    parser.add_argument("--threads", type=int)
    parser.add_argument("--interval-phi", type=float)
    parser.add_argument("--no-winding", action="store_true")
    parser.add_argument("--initial-winding", type=int, nargs="+")
    parser.add_argument("--checkpoint-every", type=int, default=10)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    _integer(args.replicas, "replicas", 1)
    _integer(args.seed, "seed", 0)
    if args.threads is not None:
        numba.set_num_threads(args.threads)
    count = len(args.N) * len(args.kappa) * len(args.epsilon) * args.replicas
    seeds = np.random.SeedSequence(args.seed).spawn(count)
    index = 0
    for N in args.N:
        for kappa in args.kappa:
            for epsilon in args.epsilon:
                for replica in range(args.replicas):
                    path = args.output / f"chain-{index:05d}.h5"
                    print(f"{path}: D={args.D}, N={N}, kappa={kappa}, epsilon={epsilon}, replica={replica}", flush=True)
                    run_chain(path, D=args.D, N=N, kappa=kappa, epsilon=epsilon,
                              sweeps=args.sweeps, warmup=args.warmup,
                              seed=int(seeds[index].generate_state(1, dtype=np.uint64)[0]),
                              backend=args.backend, parallel=not args.serial,
                              update=args.update, clusters_per_step=args.clusters_per_step,
                              record_timing=args.record_timing,
                              buffer_rows=args.buffer_rows,
                              interval_phi=args.interval_phi, winding=not args.no_winding,
                              initial_winding=args.initial_winding, checkpoint_every=args.checkpoint_every,
                              resume=args.resume and path.exists())
                    index += 1


if __name__ == "__main__":
    main()
