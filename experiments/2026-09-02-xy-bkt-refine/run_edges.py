"""Run a frozen edge batch concurrently only when its chain IDs are disjoint.

The original runner keeps its global lock. This runner uses a per-batch lock
and explicitly proves file ownership does not overlap any other manifest.
"""
from run import (OUT,run_point,atomic_json,now,sha256,_provenance,fcntl,os,shutil,
                 ProcessPoolExecutor,as_completed,mp)
import argparse
import json
from pathlib import Path


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('batch');parser.add_argument('--workers',type=int,default=4);args=parser.parse_args()
    if not 1<=args.workers<=8:raise ValueError('workers must be 1..8')
    path=OUT/'batches'/f'{args.batch}.json';batch=json.loads(path.read_text());points=batch['points'];ids={p['id'] for p in points}
    assert len(ids)==len(points)
    assert batch['selection']['concurrent_runner_sha256']==sha256(Path(__file__))
    assert batch['runner_sha256']==sha256(OUT/'run.py') and batch['model_source_sha256']==_provenance()['source_sha256']
    all_points=list(points)
    for other in (OUT/'batches').glob('*.json'):
        if other==path:continue
        pp=json.loads(other.read_text())['points'];assert not(ids&{p['id'] for p in pp}),f'Overlapping output files with {other}'
        all_points+=pp
    for name in ['jobs','chains','logs']:(OUT/name).mkdir(exist_ok=True)
    with (OUT/f'{args.batch}-runner.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        pending=[];done=0;remaining_bytes=0
        for p in all_points:
            status_path=OUT/'jobs'/f"{p['id']}.json";s=json.loads(status_path.read_text()) if status_path.exists() else {}
            remaining_bytes+=max(0,p['steps']-s.get('cycles',0))*(8*p['L']**2+2600)
            if p['id'] in ids:
                if s.get('state')=='complete':done+=1
                else:pending.append(p)
        if shutil.disk_usage(OUT).free<remaining_bytes+25*1024**3:raise RuntimeError('Insufficient disk for ALL pending batches plus reserve')
        status_path=OUT/f'{args.batch}-status.json'
        runtime=dict(batch=args.batch,state='running',started_utc=now(),pid=os.getpid(),workers=args.workers,total=len(points),completed=done,
            failures=[],provenance=_provenance(),runner_sha256=sha256(Path(__file__)),
            concurrency='Disjoint chain IDs verified against every frozen manifest',disk_free_bytes=shutil.disk_usage(OUT).free)
        atomic_json(status_path,runtime)
        print(f"Starting {len(pending)} disjoint edge chains with {args.workers} workers",flush=True)
        with ProcessPoolExecutor(max_workers=args.workers,mp_context=mp.get_context('spawn')) as pool:
            futures={pool.submit(run_point,p):p for p in pending}
            for future in as_completed(futures):
                p=futures[future]
                try:
                    r=future.result();runtime['completed']+=1
                    print(f"{runtime['completed']}/{len(points)}: {p['id']} complete; generation {r['generation_seconds']:.1f}s; full-field audit passed",flush=True)
                except Exception as exc:
                    runtime['failures'].append(dict(id=p['id'],error=repr(exc)));print(f"FAILED {p['id']}: {exc!r}",flush=True)
                runtime.update(updated_utc=now(),load=list(os.getloadavg()),disk_free_bytes=shutil.disk_usage(OUT).free)
                atomic_json(status_path,runtime)
        runtime.update(state='complete' if runtime['completed']==len(points) else 'partial',finished_utc=now());atomic_json(status_path,runtime)
        if runtime['failures']:raise RuntimeError('Edge batch failures; see logs')

if __name__=='__main__':main()
