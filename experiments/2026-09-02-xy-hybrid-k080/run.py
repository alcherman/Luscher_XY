"""Same-lattice Metropolis, hybrid and fixed-batch Wolff comparison at kappa=.8.

Runs sequentially with one CPU thread; checkpoints every 100 cycles and records
an atomic progress/timing entry after each 1000-cycle segment. Full raw fields,
observables, random states and per-cluster work are retained in every chain.
"""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
os.environ.setdefault('MPLBACKEND','Agg')
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
os.environ.setdefault('NUMBA_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')

import h5py
import numpy as np
from supervillain.xy import run_chain, _provenance

STEPS=10_000
ARMS=[dict(name='metropolis',update='metropolis',clusters_per_step=1,seed=2026091201),
      dict(name='hybrid',update='hybrid',clusters_per_step=1,seed=2026091202),
      dict(name='cluster100',update='cluster',clusters_per_step=100,seed=2026091203)]
LABELS={'metropolis':'Metropolis: 1 sweep', 'hybrid':'Hybrid: 1 sweep + 1 cluster',
        'cluster100':'Wolff: 100 clusters'}
COLORS={'metropolis':'#0072B2','hybrid':'#009E73','cluster100':'#CC79A7'}


def atomic_json(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    temporary.replace(path)


def main():
    manifest=dict(D=2,L=48,kappa=.8,epsilon=float(np.pi),steps=STEPS,arms=ARMS,
        schedule='fixed counts per cycle; one winding proposal follows every arm',
        initialization='cold; phi=0',metropolis_proposal_half_width=float(np.pi),
        backend='numba, serial, one thread',concurrency='arms run sequentially',
        cluster_batch_reason='fixed 100 chosen from prior mean cluster fraction .0104; no production adaptation',
        provisional_warmup=100,analysis_cut='selected after cold-start history inspection',
        save='cold seed and every cycle: all fields/observables/RNG states/work/timings',
        run_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    path=OUT/'manifest.json'
    if path.exists() and json.loads(path.read_text())!=manifest:
        raise ValueError('Existing experiment manifest differs')
    atomic_json(path,manifest)
    path=OUT/'status.json'
    status=(json.loads(path.read_text()) if path.exists() else dict(
        started_utc=datetime.now(timezone.utc).isoformat(),runtime=_provenance(),
        load_at_start=list(os.getloadavg()),segments=[],completed=[]))
    atomic_json(path,status)
    try:
        for arm in ARMS:
            chain=OUT/'chains'/f"{arm['name']}.h5"
            done=0
            if chain.exists():
                with h5py.File(chain,'r') as f:
                    done=int(f.attrs['committed_rows'])-1
            for target in range(1000,STEPS+1,1000):
                if done>=target: continue
                start=perf_counter()
                load=list(os.getloadavg())
                run_chain(chain,D=2,N=48,kappa=.8,epsilon=np.pi,sweeps=target,warmup=100,
                    seed=arm['seed'],backend='numba',parallel=False,winding=True,
                    update=arm['update'],clusters_per_step=arm['clusters_per_step'],
                    checkpoint_every=100,record_timing=True,resume=chain.exists())
                elapsed=perf_counter()-start
                with h5py.File(chain,'r') as f:
                    assert f.attrs['committed_rows']==target+1
                    digest=hashlib.sha256(f['phi'][-1].tobytes()+f['rng_state'][-1]).hexdigest()
                status['segments'].append(dict(arm=arm['name'],first_cycle=done+1,last_cycle=target,
                    seconds=elapsed,load_before=load,load_after=list(os.getloadavg()),
                    last_field_rng_sha256=digest,finished_utc=datetime.now(timezone.utc).isoformat()))
                done=target
                atomic_json(path,status)
                print(f"{arm['name']}: {target}/{STEPS} cycles, segment {elapsed:.2f}s",flush=True)
            result=arm|dict(file=f"chains/{arm['name']}.h5",bytes=chain.stat().st_size,
                           seconds=sum(s['seconds'] for s in status['segments'] if s['arm']==arm['name']))
            status['completed']=[a for a in status['completed'] if a['name']!=arm['name']]+[result]
            atomic_json(path,status)
        status['finished_utc']=datetime.now(timezone.utc).isoformat()
        atomic_json(path,status)
    except BaseException as exc:
        status['error']=repr(exc)
        atomic_json(path,status)
        raise


if __name__=='__main__':
    main()
