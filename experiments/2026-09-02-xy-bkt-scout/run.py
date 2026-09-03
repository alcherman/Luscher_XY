"""Checkpointed scouting scan. Every generated field and measurement is retained."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import redirect_stdout, redirect_stderr
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

os.environ['NUMBA_NUM_THREADS']='1'
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
os.environ['MPLBACKEND']='Agg'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
sys.path.insert(0,str(ROOT))
import h5py
import numpy as np
from supervillain.xy import run_chain, _provenance

EPSILON_MILLI=[500,510,520,530,540,565,580]+list(range(600,1001,20))
KAPPA_MILLI=list(range(0,1201,200))
SIZES=[16,32,64]
STEPS=10_000
ROOT_SEED=2026092201


def now(): return datetime.now(timezone.utc).isoformat()


def atomic_json(path,value):
    path=Path(path)
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    temporary.replace(path)


def sha256(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(4*1024**2),b''): digest.update(block)
    return digest.hexdigest()


def points():
    # Fix seed identity in parameter order; scheduling order may change independently.
    result=[]
    for e in EPSILON_MILLI:
        for k in KAPPA_MILLI:
            for L in SIZES:
                seed=int(np.random.SeedSequence([ROOT_SEED,e,k,L]).generate_state(1,dtype=np.uint64)[0])
                result.append(dict(id=f'e{e:04d}_k{k:04d}_L{L:03d}',L=L,
                    epsilon_milli=e,epsilon_over_pi=e/1000,epsilon=float(np.pi*(e/1000)),
                    kappa_milli=k,kappa=k/1000,seed=seed))
    return result


def principal(a):
    r=a%(2*np.pi)
    r=np.where(r>=np.pi,r-2*np.pi,r)
    return np.where((a>=-np.pi)&(a<np.pi),a,r)


def audit(path,point):
    """Independently recompute all physical observables from every saved field."""
    with h5py.File(path,'r') as f:
        p=json.loads(f.attrs['parameters']);L=point['L'];V=L*L
        assert p['D']==2 and p['N']==L and p['kappa']==point['kappa'] and p['epsilon']==point['epsilon']
        assert p['seed']==point['seed'] and p['update']=='hybrid'
        assert p['buffer_rows']==100 and p['clusters_per_step']==1 and p['winding'] and not p['parallel']
        assert f.attrs['committed_rows']==STEPS+1
        assert all(len(ds)==STEPS+1 for ds in f.values())
        assert np.array_equal(f['sweep'][:],np.arange(STEPS+1)) and np.all(f['phi'][0]==0)
        sizes=f['cluster_sizes'][:,0]
        assert np.all((sizes[1:]>=1)&(sizes[1:]<=V)) and sizes[0]==0
        assert np.array_equal(f['clusters_total'][:],np.arange(STEPS+1))
        assert np.array_equal(f['cluster_sites_total'][:],np.cumsum(sizes))
        assert np.array_equal(f['cluster_equivalent_sweeps'][:],np.cumsum(sizes)/V)
        assert np.array_equal(f['cluster_bond_tests_total'][:],np.cumsum(f['cluster_bond_tests'][:,0]))
        assert np.all((f['local_accepted'][:]>=0)&(f['local_accepted'][:]<=V))
        for key in f:
            if key.endswith('_seconds'): assert np.isfinite(f[key][:]).all() and np.all(f[key][:]>=0)
        for state in f['rng_state'][:]: assert set(json.loads(state))=={'local','cluster','winding'}
        for start in range(0,STEPS+1,128):
            stop=min(start+128,STEPS+1);a=f['phi'][start:stop,0]
            assert np.isfinite(a).all() and np.all((a>=0)&(a<2*np.pi))
            dx,dy=principal(np.roll(a,-1,axis=1)-a),principal(np.roll(a,-1,axis=2)-a)
            assert np.all(np.abs(dx)<p['epsilon']) and np.all(np.abs(dy)<p['epsilon'])
            spin=np.exp(1j*a);m=spin.mean((1,2));m2=np.abs(m)**2
            cosine=(np.cos(dx).mean((1,2))+np.cos(dy).mean((1,2)))/2
            vortex=np.abs(np.rint((dx+np.roll(dy,-1,axis=1)-np.roll(dx,-1,axis=2)-dy)/(2*np.pi))).mean((1,2))
            if point['epsilon_milli']==500: assert np.all(vortex==0)
            winding=np.stack((np.rint(dx.sum(1)/(2*np.pi)),np.rint(dy.sum(2)/(2*np.pi))),axis=1)
            phase=np.exp(-2j*np.pi*np.arange(L)/L)
            minimum=np.stack((np.abs((spin*phase[None,:,None]).sum((1,2)))**2/V,
                              np.abs((spin*phase[None,None,:]).sum((1,2)))**2/V),axis=1)
            expected=dict(magnetization_real=m.real,magnetization_imag=m.imag,
                magnetization_abs=np.abs(m),magnetization_squared=m2,magnetization_fourth=m2*m2,
                mean_link_cosine=cosine,action_density=-2*p['kappa']*cosine,vortex_density=vortex,
                winding=winding,structure_zero=V*m2,structure_min=minimum,
                max_link_angle=np.maximum(np.abs(dx).max((1,2)),np.abs(dy).max((1,2))))
            for key,value in expected.items():
                np.testing.assert_allclose(f[key][start:stop],value,rtol=2e-12,atol=3e-13,err_msg=key)
    return dict(passed=True,rows=STEPS+1,bytes=path.stat().st_size,sha256=sha256(path))


def run_point(point):
    status_path=OUT/'jobs'/f"{point['id']}.json"
    path=OUT/'chains'/f"{point['id']}.h5"
    status=json.loads(status_path.read_text()) if status_path.exists() else dict(point=point,segments=[])
    if status.get('state')=='complete': return status
    with (OUT/'logs'/f"{point['id']}.log").open('a',buffering=1) as log, redirect_stdout(log), redirect_stderr(log):
        try:
            status.update(state='running',pid=os.getpid(),updated_utc=now())
            status.setdefault('started_utc',now())
            atomic_json(status_path,status)
            done=0
            if path.exists():
                with h5py.File(path,'r') as f: done=int(f.attrs['committed_rows'])-1
            for target in range(1000,STEPS+1,1000):
                if done>=target: continue
                before=perf_counter();load=list(os.getloadavg())
                run_chain(path,D=2,N=point['L'],kappa=point['kappa'],epsilon=point['epsilon'],
                    sweeps=target,warmup=100,seed=point['seed'],backend='numba',parallel=False,
                    interval_phi=np.pi,winding=True,checkpoint_every=100,resume=path.exists(),
                    update='hybrid',clusters_per_step=1,record_timing=True,buffer_rows=100)
                elapsed=perf_counter()-before
                with h5py.File(path,'r') as f:
                    assert f.attrs['committed_rows']==target+1
                    digest=hashlib.sha256(f['phi'][-1].tobytes()+f['rng_state'][-1]).hexdigest()
                status['segments'].append(dict(first=done+1,last=target,seconds=elapsed,
                    load_before=load,load_after=list(os.getloadavg()),field_rng_sha256=digest,finished_utc=now()))
                done=target
                status.update(cycles=done,updated_utc=now())
                atomic_json(status_path,status)
                print(f"{point['id']}: {done}/{STEPS}, segment {elapsed:.3f}s",flush=True)
            status.update(state='auditing',updated_utc=now())
            atomic_json(status_path,status)
            before=perf_counter();result=audit(path,point)
            status.update(state='complete',audit=result,audit_seconds=perf_counter()-before,
                generation_seconds=sum(s['seconds'] for s in status['segments']),updated_utc=now())
            status.pop('error',None)
            atomic_json(status_path,status)
            print(f"{point['id']}: complete, full-field audit passed",flush=True)
            return status
        except BaseException:
            status.update(state='failed',error=traceback.format_exc(),updated_utc=now())
            atomic_json(status_path,status)
            traceback.print_exc()
            raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--limit',type=int,help='run only the first pending points; preserve full manifest')
    args=parser.parse_args()
    if not 1<=args.workers<=8: raise ValueError('workers must be 1..8')
    for name in ('chains','jobs','logs'): (OUT/name).mkdir(exist_ok=True)
    with (OUT/'runner.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        all_points=points();assert len(all_points)==588
        manifest=dict(schema='xy-scout.v1',D=2,sizes=SIZES,epsilon_milli=EPSILON_MILLI,
            kappa_milli=KAPPA_MILLI,steps=STEPS,points=all_points,root_seed=ROOT_SEED,
            initialization='cold, phi=0',update='1 Metropolis sweep + 1 Wolff cluster + 1 winding proposal',
            backend='numba serial, one CPU thread per independent job',interval_phi=float(np.pi),
            checkpoint_every=100,buffer_rows=100,record_timing=True,
            retention='cold seed and every cycle: full fields, all measurements, RNG, work and timing',
            analysis_status='scouting; cuts selected separately after reviewing histories',
            source_sha256=_provenance()['source_sha256'],runner_sha256=sha256(Path(__file__)))
        manifest_path=OUT/'manifest.json'
        if manifest_path.exists() and json.loads(manifest_path.read_text())!=manifest:
            raise ValueError('Campaign manifest differs; refusing to overwrite')
        atomic_json(manifest_path,manifest)
        done=[]
        for path in (OUT/'jobs').glob('*.json'):
            status=json.loads(path.read_text())
            if status.get('state')=='complete': done.append(status['point']['id'])
        pending=[p for p in all_points if p['id'] not in done]
        # Prioritize the full kappa=0 axis, then larger couplings in parameter order.
        pending.sort(key=lambda p:(p['kappa_milli']!=0,p['kappa_milli'],p['epsilon_milli'],p['L']))
        if args.limit: pending=pending[:args.limit]
        runtime=dict(started_utc=now(),pid=os.getpid(),workers=args.workers,completed=len(done),
            total=len(all_points),scheduled_this_launch=len(pending),failures=[],
            load_at_start=list(os.getloadavg()),disk_free_bytes=shutil.disk_usage(OUT).free,
            provenance=_provenance())
        atomic_json(OUT/'status.json',runtime)
        print(f"Starting {len(pending)} pending chains, {args.workers} workers, {len(done)} already complete",flush=True)
        with ProcessPoolExecutor(max_workers=args.workers,mp_context=mp.get_context('spawn')) as pool:
            futures={pool.submit(run_point,p):p for p in pending}
            for future in as_completed(futures):
                point=futures[future]
                try:
                    result=future.result();runtime['completed']+=1
                    print(f"{runtime['completed']}/588: {point['id']} complete; generation {result['generation_seconds']:.1f}s; audit passed",flush=True)
                except Exception as exc:
                    runtime['failures'].append(dict(id=point['id'],error=repr(exc)))
                    print(f"FAILED {point['id']}: {exc!r}",flush=True)
                runtime.update(updated_utc=now(),load=list(os.getloadavg()),disk_free_bytes=shutil.disk_usage(OUT).free)
                atomic_json(OUT/'status.json',runtime)
        runtime.update(finished_utc=now(),state='complete' if runtime['completed']==588 else 'partial')
        atomic_json(OUT/'status.json',runtime)
        if runtime['failures']: raise RuntimeError(f"{len(runtime['failures'])} failed jobs; see individual logs")


if __name__=='__main__': main()
