"""Bounded negative-kappa reconnaissance; frozen native Metropolis kernels."""
import os
for key in ['NUMBA_NUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import sys,json,time,fcntl,subprocess,traceback,concurrent.futures,shutil,hashlib
P=Path(__file__).resolve().parent;ROOT=P.parents[1];OLD=ROOT/'experiments/2026-10-13-xy-infinite-volume'
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(OLD))
import run as legacy
import storage_admission as storage
from supervillain.xy import run_chain
import h5py,numpy as np
from scipy.optimize import brentq
from tree_level_prediction import v1,endpoint

def atomic(p,d):legacy.atomic_json(p,d)
def read(p):return json.loads(Path(p).read_text())
def hashes():return dict(native=legacy.source_hashes(),scout=legacy.sha256(__file__),audit=legacy.sha256(P/'audit.py'),analysis=legacy.sha256(P/'analyse.py'))
def make_points():
    ps=[]
    for ki,k in enumerate([0.,-.1,-.25,-.5]):
        tc=brentq(lambda q:5*v1(k,np.pi*q)-1,.01,1);center=tc-endpoint(3)+.80076
        for q in [round(center+x,3) for x in [-.025,0,.025]]:
            for n in [8,12]:
                for r in [0,1]:
                    seed=int(np.random.SeedSequence([2026101701,ki,round(q*1000),n,r]).generate_state(1,dtype=np.uint64)[0])
                    ps.append(dict(id=f'negative3_k{ki}_q{round(q*1000):04d}_N{n:03d}_r{r}',D=3,N=n,beta=k,kappa=k,epsilon_over_pi=q,epsilon=np.pi*q,seed=seed,initial_winding=[r,0,0],steps=4000,warmup=1000,update='metropolis'))
    return ps

def work(p,manifest_name="scout.json"):
    from audit import audit
    identifier=p['id'];job=P/'jobs'/f'{identifier}.json';path=P/'chains'/f'{identifier}.h5'
    with (P/'jobs'/f'{identifier}.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if job.exists() and read(job)['state']=='complete':return read(job)
        if read(P/'batches'/manifest_name)['source_hashes']!=hashes():raise RuntimeError('Frozen source mismatch')
        atomic(job,dict(state='running',point=p,pid=os.getpid(),started=legacy.now()))
        start=time.perf_counter()
        try:
            run_chain(path,D=3,N=p['N'],kappa=p['kappa'],epsilon=p['epsilon'],sweeps=p['steps'],warmup=p['warmup'],seed=p['seed'],backend='numba',parallel=False,interval_phi=np.pi,winding=True,initial_winding=p['initial_winding'],checkpoint_every=100,resume=path.exists(),update='metropolis',record_timing=True,buffer_rows=100)
            generation=time.perf_counter()-start;start=time.perf_counter();a=audit(path,p)
            j=dict(state='complete',point=p,generation_seconds=generation,audit_seconds=time.perf_counter()-start,audit=a,updated=legacy.now());atomic(job,j);return j
        except BaseException:
            atomic(job,dict(state='failed',point=p,error=traceback.format_exc(),updated=legacy.now()));raise

def prepare():
    points=make_points();target=Path('/Volumes/Office_SSD/Luscher_XY/2026-10-17-xy-3d-negative-kappa/chains')
    if not Path('/Volumes/Office_SSD').is_mount():raise RuntimeError('SSD absent')
    target.mkdir(parents=True,exist_ok=True)
    if not (P/'chains').exists():(P/'chains').symlink_to(target)
    if (P/'chains').resolve()!=target:raise RuntimeError('Unexpected chain route')
    ceiling=sum(storage.estimated_full_bytes(x) for x in points)
    admission=storage.acquire()
    try:
        active=0;copies=0
        for folder in [OLD,ROOT/'experiments/2026-10-16-xy-3d-spin-followup',P]:
            storage.OUT=folder;copies+=storage.reserved_remaining()
            for jf in (folder/'jobs').glob('*.json'):
                j=read(jf)
                if j.get('state') in ['starting','running'] and 'point' in j:
                    f=folder/'chains'/f"{j['point']['id']}.h5";active+=max(0,storage.estimated_full_bytes(j['point'])-(f.stat().st_size if f.exists() else 0))
        free=shutil.disk_usage(target).free
        if free-ceiling-active-copies<25*2**30:raise RuntimeError('Insufficient full storage budget')
    finally:admission.close()
    dest=P/'batches/scout.json'
    if dest.exists():raise FileExistsError(dest)
    body=dict(points=points,source_hashes=hashes(),full_ceiling_bytes=ceiling,physical_free_bytes=free,active_reserved_bytes=active,copy_reserved_bytes=copies,created=legacy.now(),scope='Finite-volume cost/coverage scout, no infinite-volume or mixing claim',workers=8)
    atomic(dest,body);return body

def main():
    with (P/'manager.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        review=read(P/'review.json')
        if not review.get('approved') or review['source_hashes']!=hashes():raise RuntimeError('Independent frozen design review required')
        manifest=read(P/'batches/scout.json');points=manifest['points'];pilot=next(p for p in points if p['N']==12 and p['kappa']==-.25 and p['initial_winding']==[0,0,0])
        awake=subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(os.getpid())])
        try:
            atomic(P/'status.json',dict(state='cost_pilot',updated=legacy.now()))
            with (P/'logs'/f"{pilot['id']}.log").open('a') as out:
                subprocess.run([sys.executable,'-u',__file__,'worker',pilot['id']],cwd=ROOT,stdout=out,stderr=subprocess.STDOUT,check=True)
            j=read(P/'jobs'/f"{pilot['id']}.json");cost=j['generation_seconds']+j['audit_seconds'];atomic(P/'cost-pilot.json',dict(job=j,ceiling_all_worker_seconds=cost*48,qualification='One trajectory cost, not a mixing or precision guarantee'))
            if cost>120:raise RuntimeError('Pilot exceeds bounded scout allowance; review cost before more sampling')
            atomic(P/'status.json',dict(state='running_scout',pilot_seconds=cost,updated=legacy.now()))
            # Each subprocess owns a full unfiltered trajectory log.
            def launch(point):
                with (P/'logs'/f"{point['id']}.log").open('a') as out:subprocess.run([sys.executable,'-u',__file__,'worker',point['id']],cwd=ROOT,stdout=out,stderr=subprocess.STDOUT,check=True)
                return point['id']
            with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
                for i,identifier in enumerate(pool.map(launch,[p for p in points if p['id']!=pilot['id']]),2):
                    atomic(P/'status.json',dict(state='running_scout',complete=i,total=len(points),last_id=identifier,updated=legacy.now()))
            atomic(P/'status.json',dict(state='all_audited_analysing',complete=48,updated=legacy.now()))
            with (P/'logs/analysis.log').open('a') as out:subprocess.run([sys.executable,'-u',str(P/'analyse.py')],cwd=ROOT,stdout=out,stderr=subprocess.STDOUT,check=True)
            atomic(P/'status.json',dict(state='complete_preliminary_finite_volume',complete=48,updated=legacy.now()))
        except BaseException:
            atomic(P/'status.json',dict(state='failed_needs_attention',error=traceback.format_exc(),updated=legacy.now()));raise
        finally:
            if awake.poll() is None:awake.terminate()
if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='prepare':print(json.dumps(prepare(),default=str)[:300])
    elif len(sys.argv)>1 and sys.argv[1]=='worker':
        manifest_name=sys.argv[3] if len(sys.argv)>3 else 'scout.json'
        work(next(p for p in read(P/'batches'/manifest_name)['points'] if p['id']==sys.argv[2]),manifest_name)
    else:main()
