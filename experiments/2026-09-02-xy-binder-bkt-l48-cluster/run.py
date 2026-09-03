"""D=2, L=48 cluster extension of the existing seven-coupling Binder grid."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
PREVIOUS = OUT.parent / '2026-09-02-xy-binder-bkt'
sys.path.insert(0, str(ROOT))
os.environ.setdefault('MPLBACKEND', 'Agg')
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/luscher-mpl')
os.environ.setdefault('NUMBA_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')

import h5py
import numpy as np
from supervillain.xy import run_chain, _provenance

KAPPAS = (.80, .95, 1.05, 1.12, 1.20, 1.35, 1.50)
STEPS = 10_000


def atomic_json(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def points():
    return [dict(L=48, kappa=k, seed=2026090501+i,
                 file=f'chains/L48_k{k:.2f}_cluster.h5') for i, k in enumerate(KAPPAS)]


def run_point(point):
    start = perf_counter()
    path = OUT / point['file']
    run_chain(path, D=2, N=48, kappa=point['kappa'], epsilon=np.pi,
              sweeps=STEPS, warmup=100, seed=point['seed'], backend='numba',
              update='cluster', clusters_per_step=1, parallel=False,
              winding=True, checkpoint_every=100, resume=path.exists())
    with h5py.File(path, 'r') as handle:
        assert handle.attrs['committed_rows'] == STEPS+1
        assert np.array_equal(handle['sweep'][:], np.arange(STEPS+1))
        assert np.array_equal(handle['clusters_total'][:], np.arange(STEPS+1))
        work = float(handle['cluster_equivalent_sweeps'][-1])
    return point | dict(seconds=perf_counter()-start, cluster_equivalent_sweeps=work,
                        completed_utc=datetime.now(timezone.utc).isoformat(),
                        bytes=path.stat().st_size)


def main():
    previous = json.loads((PREVIOUS / 'manifest-l32.json').read_text())
    assert tuple(p['kappa'] for p in previous['points']) == KAPPAS
    manifest = dict(D=2, L=48, epsilon=float(np.pi),
                    action='S=-kappa*sum_positive_links(cos(theta))',
                    cycles=STEPS, clusters_per_cycle=1,
                    global_update='one symmetric winding proposal per cycle',
                    backend='numba', initialization='cold phi=0',
                    provisional_warmup=100, analysis_cut='chosen after inspecting histories',
                    retention='cold seed and every cycle: full fields, all observables, RNG states, cluster work',
                    workers=3, points=points(),
                    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    prior_manifest_sha256=hashlib.sha256((PREVIOUS / 'manifest-l32.json').read_bytes()).hexdigest(),
                    prior_results_sha256=hashlib.sha256((PREVIOUS / 'results.json').read_bytes()).hexdigest())
    path = OUT / 'manifest.json'
    if path.exists() and json.loads(path.read_text()) != manifest:
        raise ValueError('Manifest differs; refusing to replace run provenance')
    atomic_json(path, manifest)
    status_path = OUT / 'status.json'
    status = (json.loads(status_path.read_text()) if status_path.exists() else
              dict(started_utc=datetime.now(timezone.utc).isoformat(), completed=[],
                   runtime_provenance=_provenance()))
    atomic_json(status_path, status)
    try:
        with ProcessPoolExecutor(max_workers=3) as pool:
            jobs = [pool.submit(run_point, p) for p in points()]
            for future in as_completed(jobs):
                result = future.result()
                status['completed'] = [p for p in status['completed'] if p['kappa'] != result['kappa']]
                status['completed'].append(result)
                atomic_json(status_path, status)
                print(f"{len(status['completed'])}/7 complete: L=48 kappa={result['kappa']:.2f}, "
                      f"{result['seconds']:.1f}s, cluster work={result['cluster_equivalent_sweeps']:.1f} volumes", flush=True)
    except BaseException as error:
        status['error'] = repr(error)
        atomic_json(status_path, status)
        raise
    status['finished_utc'] = datetime.now(timezone.utc).isoformat()
    atomic_json(status_path, status)


if __name__ == '__main__':
    main()
