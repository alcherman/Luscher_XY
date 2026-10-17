"""Finite SSD follow-up runner. Frozen model/update sources remain unchanged."""
import os
for k in ['NUMBA_NUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,time,fcntl,subprocess,signal,shutil,traceback,plistlib,uuid
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'experiments/2026-10-16-xy-3d-spin-followup';LEGACY=ROOT/'experiments/2026-10-13-xy-infinite-volume';TARGET=Path('/Volumes/Office_SSD/Luscher_XY/2026-10-16-xy-3d-spin-followup/chains')
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(LEGACY))
import run as frozen
import storage_admission as storage
import retain_fields as retain
sys.path.insert(0,str(OUT))
import online_spin as online
import retention_online
import h5py,numpy as np

RESERVE=25*2**30

def atomic(path,value):frozen.atomic_json(path,value)
def sources():return dict(kernels=frozen.source_hashes(),online=frozen.sha256(OUT/'online_spin.py'),runner=frozen.sha256(OUT/'campaign.py'),retention_policy=frozen.sha256(OUT/'retention_online.py'))
def read(path):return json.loads(Path(path).read_text())
def route():
    if not Path('/Volumes/Office_SSD').is_mount() or not TARGET.is_dir() or (OUT/'chains').resolve()!=TARGET:raise RuntimeError('Verified dedicated SSD route unavailable')
    expected=read(OUT/'control.json')['ssd_volume_uuid']
    actual=plistlib.loads(subprocess.check_output(['/usr/sbin/diskutil','info','-plist','/Volumes/Office_SSD']))['VolumeUUID']
    if actual!=expected:raise RuntimeError('SSD volume identity changed')
def live(pid,needle):
    try:
        os.kill(int(pid),0);command=subprocess.check_output(['/bin/ps','-p',str(pid),'-o','command='],text=True)
        return str(OUT/'campaign.py') in command and needle in command
    except (ProcessLookupError,subprocess.CalledProcessError):return False

def extra_audit(path,spec):
    maximum=0.
    with h5py.File(path,'r') as h:
        if json.loads(h.attrs['provenance']).get('spin_measurement_extension')!=online.definition():raise ValueError('Unpinned online measurement')
        for first in range(0,spec['steps']+1,16):
            stop=min(first+16,spec['steps']+1);phi=h['phi'][first:stop,0];real=np.empty((len(phi),2,3));imag=np.empty_like(real)
            for j,r in enumerate([spec['N']//4,spec['N']//2]):
                for axis in range(3):
                    delta=np.roll(phi,-r,axis=axis+1)-phi;real[:,j,axis]=np.cos(delta).mean((1,2,3));imag[:,j,axis]=np.sin(delta).mean((1,2,3))
            for key,expected in [('spin_correlation_real',real),('spin_correlation_imag',imag)]:
                got=h[key][first:stop];np.testing.assert_allclose(got,expected,rtol=2e-12,atol=3e-13,err_msg=key);maximum=max(maximum,float(np.max(abs(got-expected))))
            assert np.max(abs(h['spin_correlation_imag'][first:stop,1]))<1e-12
    return dict(passed=True,rows=spec['steps']+1,max_abs_error=maximum,method='Independent direct cos/sin field-difference audit on every row and direction')

def worker(manifest,identifier,token):
    batch=read(manifest)
    if batch['source_hashes']!=sources():raise RuntimeError('Frozen batch source mismatch')
    p=next(p for p in batch['points'] if p['id']==identifier);route();job=OUT/'jobs'/f'{identifier}.json';path=OUT/'chains'/f'{identifier}.h5'
    with (OUT/'jobs'/f'{identifier}.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX);old=read(job) if job.exists() else {}
        if old.get('state')=='complete':return
        if old.get('launch_token')!=token:
            print('Superseded launch claim; no simulation or archive writes',identifier,flush=True);return
        state=dict(state='running',launch_token=token,point=p,pid=os.getpid(),updated=frozen.now(),manifest=str(manifest),source_hashes=batch['source_hashes'],attempt_history=old.get('attempt_history',[])+([old] if old else []));atomic(job,state);start=time.perf_counter()
        try:
            committed=-1
            if path.exists():
                with h5py.File(path,'r') as h:committed=int(h.attrs['committed_rows'])-1
            kwargs=dict(D=3,N=p['N'],kappa=0.,epsilon=p['epsilon'],warmup=p['warmup'],seed=p['seed'],backend='numba',parallel=False,interval_phi=np.pi,winding=True,initial_winding=p['initial_winding'],checkpoint_every=100,update='hybrid',clusters_per_step=4,record_timing=True,buffer_rows=100)
            producer=online.run_chain if p['online'] else frozen.run_chain
            targets=list(range(1000,p['steps'],1000))+[p['steps']]
            for target in targets:
                if target<=committed:continue
                route();producer(path,sweeps=target,resume=path.exists(),**kwargs);print(frozen.now(),'checkpoint',identifier,target,flush=True)
            generation=time.perf_counter()-start;mark=time.perf_counter();audit=frozen.audit(path,p)
            if p['online']:audit['online_spin']=extra_audit(path,p)
            state.update(state='complete',audit=audit,generation_seconds=generation,audit_seconds=time.perf_counter()-mark,updated=frozen.now());atomic(job,state);print(frozen.now(),'closed',identifier,state['generation_seconds'],state['audit_seconds'],flush=True)
        except BaseException as error:
            state.update(state='interrupted' if isinstance(error,KeyboardInterrupt) else 'failed',error=traceback.format_exc(),updated=frozen.now());atomic(job,state);raise

def do_retention(identifier):
    storage.OUT=OUT;retain.OUT=OUT
    route();j=read(OUT/'jobs'/f'{identifier}.json')
    if j['state']!='complete' or not j['audit'].get('online_spin',{}).get('passed'):raise RuntimeError('Full original/online closure audits required')
    result=retention_online.compact(OUT/'chains'/f'{identifier}.h5',j['audit'])
    print(json.dumps(result),flush=True)

def available_for(points):
    """Admission caller holds shared lock; reserve active growth before fresh free."""
    growth=0
    for jpath in (OUT/'jobs').glob('*.json'):
        j=read(jpath)
        if j['state'] in ['starting','running']:
            path=OUT/'chains'/f"{j['point']['id']}.h5";size=path.stat().st_size if path.exists() else 0;growth+=max(0,storage.estimated_full_bytes(j['point'])-size)
    copies=storage.reserved_remaining();free=shutil.disk_usage(TARGET).free;requirements=sum(storage.estimated_full_bytes(p) for p in points)
    return dict(physical_free_bytes=free,active_remaining_bytes=growth,copy_remaining_bytes=copies,new_full_ceiling_bytes=requirements,reserve_bytes=RESERVE,admitted=free-growth-copies-requirements>=RESERVE)

def recover_starting(job_path):
    """Invalidate a launch nonce only while no worker owns the trajectory lock."""
    j=read(job_path);identifier=j['point']['id']
    with (OUT/'jobs'/f'{identifier}.lock').open('a+') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return j
        j=read(job_path)
        if j['state']=='starting':
            history=j.get('launch_recoveries',[])+[dict(token=j['launch_token'],updated=j.get('updated'))]
            j.update(state='interrupted',launch_token='invalidated-'+uuid.uuid4().hex,launch_recoveries=history,note='Recovered pre-PID launch claim under exclusive job lock; obsolete child must not run')
            atomic(job_path,j)
        elif j['state']=='running' and not live(j['pid'],'--worker '+identifier):
            j.update(state='interrupted',launch_token='invalidated-'+uuid.uuid4().hex,updated=frozen.now(),note='Owner disappeared; re-read under exclusive job lock before exact-checkpoint resume');atomic(job_path,j)
        return j

def admit_and_spawn(point,manifest):
    """Durable full-growth reservation and nonce precede Popen, under both locks."""
    identifier=point['id'];job=OUT/'jobs'/f'{identifier}.json';admission=storage.acquire()
    try:
        budget=available_for([point])
        if not budget['admitted']:return None,'space'
        with (OUT/'jobs'/f'{identifier}.lock').open('a+') as lock:
            try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            except BlockingIOError:return None,'owned'
            old=read(job) if job.exists() else {}
            if old.get('state') not in [None,'pending','interrupted']:return None,'already_claimed'
            token=uuid.uuid4().hex
            starting=dict(state='starting',point=point,launch_token=token,coordinator_pid=os.getpid(),updated=frozen.now(),attempt_history=old.get('attempt_history',[]),launch_recoveries=old.get('launch_recoveries',[]))
            atomic(job,starting)
            with (OUT/'logs'/f'{identifier}.log').open('a') as stream:
                try:proc=subprocess.Popen([sys.executable,'-u',str(OUT/'campaign.py'),'--worker',identifier,'--token',token,'--manifest',str(manifest)],cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
                except OSError:
                    starting.update(state='failed',error=traceback.format_exc());atomic(job,starting);raise
            atomic(job,dict(starting,state='running',pid=proc.pid,note='Full-growth reservation and nonce persisted before child creation'))
            return proc,'launched'
    finally:admission.close()

def coordinate(manifest,workers=14,compact=False):
    storage.OUT=OUT
    manifest=Path(manifest).resolve();batch=read(manifest);points=batch['points'];route()
    if batch['source_hashes']!=sources():raise RuntimeError('Batch source mismatch')
    lock=(OUT/'coordinator.lock').open('a+');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);children={};previous=None
    awake=subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(os.getpid())]);atomic(OUT/'coordinator-pid.json',dict(pid=os.getpid(),manifest=str(manifest),started=frozen.now()))
    try:
        while True:
            state={};running=[];pending=[];failed=[]
            for p in points:
                jpath=OUT/'jobs'/f"{p['id']}.json";j=read(jpath) if jpath.exists() else dict(state='pending',point=p)
                if j['state']=='starting':j=recover_starting(jpath)
                if j['state']=='running' and not live(j['pid'],'--worker '+p['id']):j=recover_starting(jpath)
                state[j['state']]=state.get(j['state'],0)+1
                if j['state'] in ['starting','running']:running.append(p)
                elif j['state'] in ['pending','interrupted']:pending.append(p)
                elif j['state']=='failed':failed.append(p['id'])
            for identifier,proc in list(children.items()):
                if proc.poll() is not None:
                    children.pop(identifier)
                    if proc.returncode>0 and (OUT/'jobs'/f'{identifier}.json').exists():
                        jp=OUT/'jobs'/f'{identifier}.json';jj=read(jp)
                        if jj['state'] not in ['complete','failed'] and not (jj['state']=='interrupted' and jj.get('error')):
                            jj.update(state='failed',error='Worker exited before publishing failure; inspect full log',exit_code=proc.returncode);atomic(jp,jj);raise RuntimeError('Worker startup failure: '+identifier)
            if failed:raise RuntimeError('Trajectory/audit failures; no new launches: '+str(failed))
            config=read(OUT/'control.json');allowed=config.get('launch_new',False)
            if state=={'complete':len(points)}:break
            blocked=[]
            for p in pending:
                if not allowed or len(running)>=min(workers,int(config.get('max_total_cpu_workers',14))):break
                proc,reason=admit_and_spawn(p,manifest)
                if reason=='space':blocked.append(p['id'])
                if proc is not None:children[p['id']]=proc;running.append(p)
            current=dict(updated=frozen.now(),pid=os.getpid(),manifest=str(manifest),job_states=state,running_ids=[p['id'] for p in running],storage_blocked_ids=blocked,launch_new=allowed,max_total_cpu_workers=workers,state='running')
            atomic(OUT/'status.json',current)
            if state!=previous:print(frozen.now(),state,flush=True);previous=state
            time.sleep(5)
        atomic(OUT/'status.json',dict(state='all_trajectories_audited',manifest=str(manifest),count=len(points),updated=frozen.now()))
        if compact:
            for p in points:
                if not p['online']:continue
                with (OUT/'logs'/f"retention-{p['id']}.log").open('a') as stream:subprocess.run([sys.executable,'-u',str(OUT/'campaign.py'),'--retain',p['id']],cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,check=True)
        atomic(OUT/'status.json',dict(state='complete',manifest=str(manifest),count=len(points),retention_requested=compact,updated=frozen.now()))
    finally:
        if awake.poll() is None:awake.terminate()
        lock.close()

if __name__=='__main__':
    import argparse
    a=argparse.ArgumentParser();a.add_argument('--manifest');a.add_argument('--worker');a.add_argument('--token');a.add_argument('--coordinate',action='store_true');a.add_argument('--workers',type=int,default=14);a.add_argument('--retain');a.add_argument('--compact',action='store_true');args=a.parse_args()
    if args.worker:
        def stop(signum,frame):raise KeyboardInterrupt('Safe checkpoint interruption requested')
        signal.signal(signal.SIGTERM,stop);worker(args.manifest,args.worker,args.token)
    elif args.retain:do_retention(args.retain)
    else:coordinate(args.manifest,args.workers,args.compact)
