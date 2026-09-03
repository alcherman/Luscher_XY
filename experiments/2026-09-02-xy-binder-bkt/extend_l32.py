"""Add seven L=32 chains using the original sampler and sweep budget."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import experiment as original

OUT = Path(__file__).resolve().parent


def points():
    return [dict(L=32, kappa=k, seed=2026090422 + i,
                 file=f'chains/L32_k{k:.2f}.h5')
            for i, k in enumerate(original.KAPPAS)]


def run():
    parent = json.loads((OUT / 'manifest.json').read_text())
    assert hashlib.sha256((OUT / 'experiment.py').read_bytes()).hexdigest() == parent['experiment_sha256']
    manifest = parent | dict(points=points(),
        extension='L=32, same seven couplings and 10,000-sweep budget',
        extension_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        parent_manifest_sha256=hashlib.sha256((OUT / 'manifest.json').read_bytes()).hexdigest())
    path = OUT / 'manifest-l32.json'
    if path.exists() and json.loads(path.read_text()) != manifest:
        raise ValueError('The existing L32 manifest differs; refusing to change provenance')
    original.atomic_json(path, manifest)
    status_path = OUT / 'status-l32.json'
    record = (json.loads(status_path.read_text()) if status_path.exists() else
              dict(started_utc=datetime.now(timezone.utc).isoformat(), completed=[]))
    original.atomic_json(status_path, record)
    with ProcessPoolExecutor(max_workers=3) as pool:
        jobs = [pool.submit(original.run_point, p) for p in points()]
        for future in as_completed(jobs):
            result = future.result()
            record['completed'] = [p for p in record['completed'] if p['kappa'] != result['kappa']]
            record['completed'].append(result)
            original.atomic_json(status_path, record)
            print(f"{len(record['completed'])}/7 complete: L=32 κ={result['kappa']:.2f} "
                  f"{result['seconds']:.1f}s", flush=True)
    record['finished_utc'] = datetime.now(timezone.utc).isoformat()
    original.atomic_json(status_path, record)


if __name__ == '__main__':
    run()
