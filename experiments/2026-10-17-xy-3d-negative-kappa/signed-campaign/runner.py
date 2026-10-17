"""Finite, source-pinned signed-cluster campaign; no old campaign controls modified."""
import os
for key in ['NUMBA_NUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import sys,json,fcntl,time,subprocess,concurrent.futures,traceback,shutil,plistlib
P=Path(__file__).resolve().parent;ROOT=P.parents[2];OLD=ROOT/'experiments/2026-10-13-xy-infinite-volume';sys.path[:0]=[str(P.parent),str(ROOT),str(OLD)]
import signed_cluster as signed
import run as legacy
import storage_admission as storage
import h5py,numpy as np
TARGET=Path('/Volumes/Office_SSD/Luscher_XY/2026-10-17-xy-3d-negative-kappa/signed-chains')
def read(p):return json.loads(Path(p).read_text())
def atomic(p,x):legacy.atomic_json(p,x)
def sources():return dict(native=legacy.source_hashes(),extension=signed.definition(),runner=legacy.sha256(__file__))
def route():
 if not Path('/Volumes/Office_SSD').is_mount() or (P/'chains').resolve()!=TARGET:raise RuntimeError('Dedicated SSD route unavailable')
 actual=plistlib.loads(subprocess.check_output(['/usr/sbin/diskutil','info','-plist','/Volumes/Office_SSD']))['VolumeUUID']
 if actual!='53650819-819D-42D6-BB15-24760A94AE0B':raise RuntimeError('Unexpected SSD volume')
def point(k,q,n,r,steps=4000):
 km=int(round(-k*1000));qm=int(round(q*1000000));seed=int(np.random.SeedSequence([2026101703,km,qm,n,r,steps]).generate_state(1,dtype=np.uint64)[0])
 return dict(id=f'signed3_km{km:04d}_q{qm:07d}_N{n:03d}_r{r:02d}_s{steps:06d}',D=3,N=n,beta=k,kappa=k,epsilon_over_pi=q,epsilon=q*np.pi,seed=seed,initial_winding=[r%2,0,0],replica=r,steps=steps,warmup=min(1000,steps//4),clusters_per_cycle=4)
def freeze(name,ps,purpose):
 route();path=P/'batches'/f'{name}.json'
 if path.exists():raise FileExistsError(path)
 lock=storage.acquire()
 try:
  copies=0
  for folder in [OLD,ROOT/'experiments/2026-10-16-xy-3d-spin-followup',P.parent,P]:
   storage.OUT=folder;copies+=storage.reserved_remaining()
   for jfile in (folder/'jobs').glob('*.json'):
    if read(jfile).get('state') in ['starting','running']:raise RuntimeError('Other known active campaign needs storage reconciliation')
  ceil=sum(storage.estimated_full_bytes(p) for p in ps);free=shutil.disk_usage(TARGET).free
  if free-copies-ceil<25*2**30:raise RuntimeError('Full batch reserve failed')
  atomic(path,dict(name=name,points=ps,purpose=purpose,source_hashes=sources(),created=legacy.now(),full_ceiling_bytes=ceil,physical_free_bytes=free,copy_reserved_bytes=copies,max_workers=8))
 finally:lock.close()
 folder=P/'source-snapshots'/name;folder.mkdir()
 for src in [Path(__file__),P.parent/'signed_cluster.py',P.parent/'_signed_cluster_kernel.py']:shutil.copyfile(src,folder/src.name)
 return path

def worker(manifest,identifier):
 batch=read(manifest)
 if batch['source_hashes']!=sources():raise RuntimeError('Frozen signed production source changed')
 p=next(v for v in batch['points'] if v['id']==identifier);route();job=P/'jobs'/f'{identifier}.json';file=P/'chains'/f'{identifier}.h5'
 with (P/'jobs'/f'{identifier}.lock').open('a+') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  if job.exists() and read(job)['state']=='complete':return
  atomic(job,dict(state='running',point=p,pid=os.getpid(),manifest=str(manifest),started=legacy.now()))
  try:
   start=time.perf_counter()
   signed.run_chain(file,D=3,N=p['N'],kappa=p['kappa'],epsilon=p['epsilon'],sweeps=p['steps'],warmup=p['warmup'],seed=p['seed'],backend='numba',parallel=False,interval_phi=np.pi,winding=True,initial_winding=p['initial_winding'],checkpoint_every=100,resume=file.exists(),update='hybrid',clusters_per_step=4,record_timing=True,buffer_rows=100)
   generation=time.perf_counter()-start;start=time.perf_counter();audit=legacy.audit(file,p)
   with h5py.File(file) as h:
    assert json.loads(h.attrs['provenance'])['signed_cluster_extension']==signed.definition()
   atomic(job,dict(state='complete',point=p,generation_seconds=generation,audit_seconds=time.perf_counter()-start,audit=audit,updated=legacy.now()))
  except BaseException:
   atomic(job,dict(state='failed',point=p,error=traceback.format_exc(),updated=legacy.now()));raise

def run(manifest):
 with (P/'manager.lock').open('a+') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);route();batch=read(manifest);review=read(Path(manifest).with_suffix('.review.json'))
  if not review.get('approved') or review['source_hashes']!=sources() or review['points']!=batch['points']:raise RuntimeError('Concrete signed-batch review required')
  awake=subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(os.getpid())]);ps=batch['points']
  try:
   def launch(p):
    with (P/'logs'/f"{p['id']}.log").open('a') as out:subprocess.run([sys.executable,'-u',__file__,'worker',str(manifest),p['id']],cwd=ROOT,stdout=out,stderr=subprocess.STDOUT,check=True)
    return p['id']
   atomic(P/'status.json',dict(state='running',manifest=str(manifest),total=len(ps),updated=legacy.now()))
   with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
    futures=[pool.submit(launch,p) for p in ps]
    for i,future in enumerate(concurrent.futures.as_completed(futures),1):atomic(P/'status.json',dict(state='running',complete=i,total=len(ps),last=future.result(),manifest=str(manifest),updated=legacy.now()))
   atomic(P/'status.json',dict(state='all_full_audits_passed_pending_analysis',complete=len(ps),manifest=str(manifest),updated=legacy.now()))
  except BaseException:
   atomic(P/'status.json',dict(state='failed_needs_attention',error=traceback.format_exc(),manifest=str(manifest),updated=legacy.now()));raise
  finally:
   if awake.poll() is None:awake.terminate()
if __name__=='__main__':
 if sys.argv[1]=='worker':worker(Path(sys.argv[2]),sys.argv[3])
 elif sys.argv[1]=='prepare-pilot':
  TARGET.mkdir(parents=True,exist_ok=True)
  if not (P/'chains').exists():(P/'chains').symlink_to(TARGET)
  ps=[point(k,q,n,r) for k,qs in [(-.25,[.65]),(-.5,[.57,.60])] for q in qs for n in [8,12] for r in [10,11]]
  print(freeze('signed-sector-mixing-pilot',ps,'Validate observed mixing of zero and unit-winding starts before refining negative continuation. No infinite-volume claim.'))
 elif sys.argv[1]=='run':run(Path(sys.argv[2]).resolve())
 else:raise ValueError('Unknown operation')
