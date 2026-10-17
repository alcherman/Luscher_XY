"""Bounded replica-exchange batches; per-ladder workers are single threaded."""
import importlib.util
from pathlib import Path
_ps=importlib.util.spec_from_file_location('asymptote_pt',Path(__file__).with_name('pt.py'));pt=importlib.util.module_from_spec(_ps);_ps.loader.exec_module(pt)
from pathlib import Path
import os,sys,json,fcntl,time,subprocess,concurrent.futures,traceback,shutil,plistlib
import numpy as np
import storage_admission as storage
import importlib.util
_audit_spec=importlib.util.spec_from_file_location('joint_replica_audit',Path(__file__).with_name('audit.py'));_audit_module=importlib.util.module_from_spec(_audit_spec);_audit_spec.loader.exec_module(_audit_module);audit=_audit_module.audit
P=Path(__file__).resolve().parent;ROOT=P.parents[2];TARGET=Path('/Volumes/Office_SSD/Luscher_XY/2026-10-17-xy-3d-negative-kappa/sign-sector-chains')
def atomic(p,x):pt.legacy.atomic_json(p,x)
def read(p):return json.loads(Path(p).read_text())
def sources():return dict(producer=pt.sources(),runner=pt.legacy.sha256(__file__),audit=pt.legacy.sha256(P/'audit.py'))
def route():
 if not Path('/Volumes/Office_SSD').is_mount() or (P/'chains').resolve()!=TARGET:raise RuntimeError('Dedicated PT SSD route unavailable')
 if plistlib.loads(subprocess.check_output(['/usr/sbin/diskutil','info','-plist','/Volumes/Office_SSD']))['VolumeUUID']!='53650819-819D-42D6-BB15-24760A94AE0B':raise RuntimeError('Unexpected SSD volume')
def point(k,q,n,r,steps=2000,spacing=.025):
 seed=int(np.random.SeedSequence([2026101709,round(-k*1000),round(q*1e6),n,r,steps]).generate_state(1,dtype=np.uint64)[0])
 return dict(id=f'scaled3_km{round(-k*1000):04d}_q{round(q*1e6):07d}_N{n:03d}_r{r:02d}_s{steps:06d}',N=n,D=3,target_kappa=k,kappas=pt.ladder(k,spacing).tolist(),epsilon=q*np.pi,epsilon_over_pi=q,steps=steps,warmup=min(1000,steps//4),seed=seed,initial_winding=[r%2,0,0],checkpoint_every=100)
def ceiling(p):
 slots=len(p['kappas']);base=storage.estimated_full_bytes(p)
 return int(slots*base+1.25*(p['steps']+1)*slots*4096+64*2**20)
def freeze(name,ps,purpose):
 route();path=P/'batches'/f'{name}.json'
 if path.exists():raise FileExistsError(path)
 if len({p['id'] for p in ps})!=len(ps) or len({p['seed'] for p in ps})!=len(ps):raise ValueError('Duplicate point/seed')
 if any(p['N'] not in [6,8] or p['steps']>16000 or p['target_kappa']>0 or p['target_kappa']< -16 for p in ps):raise ValueError('Unreviewed production scale')
 lock=storage.acquire()
 try:
  copies=0
  folders=[pt.OLD,ROOT/'experiments/2026-10-16-xy-3d-spin-followup',P.parent,P.parent/'signed-campaign',P.parent/'tempered-campaign',P.parent/'replica-exchange',P.parent/'asymptote',P.parent/'scaled-asymptote',P.parent/'strong-limit',P]
  for folder in folders:
   storage.OUT=folder;copies+=storage.reserved_remaining()
   for jf in (folder/'jobs').glob('*.json'):
    if read(jf).get('state') in ['starting','running']:raise RuntimeError('Reconcile known active campaign before fresh PT commitment')
  total=sum(ceiling(p) for p in ps);free=shutil.disk_usage(TARGET).free
  if free-total-copies<25*2**30:raise RuntimeError('Whole joint-history storage reserve failed')
  atomic(path,dict(name=name,purpose=purpose,points=ps,source_hashes=sources(),physical_free_bytes=free,full_ceiling_bytes=total,copy_reserved_bytes=copies,created=pt.legacy.now()))
 finally:lock.close()
 folder=P/'source-snapshots'/name;folder.mkdir()
 for src in [Path(__file__),P/'pt.py',P/'audit.py',P/'sector_moves.py',P/'base_audit.py',P.parent/'signed_cluster.py',P.parent/'_signed_cluster_kernel.py']:shutil.copyfile(src,folder/src.name)
 return path

def worker(manifest,identifier):
 batch=read(manifest)
 if batch['source_hashes']!=sources():raise RuntimeError('Frozen joint producer source mismatch')
 p=next(x for x in batch['points'] if x['id']==identifier);route();job=P/'jobs'/f'{identifier}.json';file=P/'chains'/f'{identifier}.h5'
 with (P/'jobs'/f'{identifier}.lock').open('a+') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  if job.exists() and read(job)['state']=='complete':return
  atomic(job,dict(state='running',point=p,pid=os.getpid(),started=pt.legacy.now(),manifest=str(manifest)))
  try:
   start=time.perf_counter();pt.run(file,N=p['N'],epsilon=p['epsilon'],kappas=p['kappas'],steps=p['steps'],warmup=p['warmup'],seed=p['seed'],initial_winding=p['initial_winding'],checkpoint_every=p['checkpoint_every'],initial_kind=p.get('initial_kind','winding'),sector_updates=p.get('sector_updates',True),resume=file.exists());generation=time.perf_counter()-start;start=time.perf_counter();a=audit(file,p,batch['source_hashes']['producer'])
   atomic(job,dict(state='complete',point=p,generation_seconds=generation,audit_seconds=time.perf_counter()-start,audit=a,updated=pt.legacy.now()))
  except BaseException:
   atomic(job,dict(state='failed',point=p,error=traceback.format_exc(),updated=pt.legacy.now()));raise

def run(manifest):
 with (P/'manager.lock').open('a+') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);route();batch=read(manifest);review=read(Path(manifest).with_suffix('.review.json'))
  if not review.get('approved') or review['source_hashes']!=sources() or review['points']!=batch['points']:raise RuntimeError('Independent concrete PT batch review required')
  ps=batch['points'];awake=subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(os.getpid())])
  try:
   def launch(p):
    with (P/'logs'/f"{p['id']}.log").open('a') as out:subprocess.run([sys.executable,'-u',__file__,'worker',str(manifest),p['id']],cwd=ROOT,stdout=out,stderr=subprocess.STDOUT,check=True)
    return p['id']
   atomic(P/'status.json',dict(state='running_joint_ladders',total=len(ps),manifest=str(manifest),updated=pt.legacy.now()))
   with concurrent.futures.ThreadPoolExecutor(max_workers=min(8,len(ps))) as pool:
    futures=[pool.submit(launch,p) for p in ps]
    for i,f in enumerate(concurrent.futures.as_completed(futures),1):atomic(P/'status.json',dict(state='running_joint_ladders',complete=i,total=len(ps),last=f.result(),manifest=str(manifest),updated=pt.legacy.now()))
   atomic(P/'status.json',dict(state='all_joint_field_swap_audits_passed_pending_analysis',complete=len(ps),manifest=str(manifest),updated=pt.legacy.now()))
  except BaseException:
   atomic(P/'status.json',dict(state='failed_needs_attention',error=traceback.format_exc(),updated=pt.legacy.now()));raise
  finally:
   if awake.poll() is None:awake.terminate()
if __name__=='__main__':
 if sys.argv[1]=='worker':worker(Path(sys.argv[2]),sys.argv[3])
 elif sys.argv[1]=='run':run(Path(sys.argv[2]).resolve())
 elif sys.argv[1]=='prepare-pilot':
  TARGET.mkdir(parents=True,exist_ok=True)
  if not (P/'chains').exists():(P/'chains').symlink_to(TARGET)
  ps=[point(k,q,12,rep) for k,q in [(-.25,.65),(-.5,.57)] for rep in [30,31]]
  print(freeze('joint-replica-exchange-mixing-pilot',ps,'Four-ladder fixed-epsilon bridge test; verify round trips and target-start convergence before any phase-boundary sampling.'))
 else:raise ValueError('Unknown action')
