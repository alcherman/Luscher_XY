"""Verified adaptive field retention. Observables/RNG histories remain complete.

phi retains its logical shape for native checkpoint replay. Unstored fields are
explicit NaN fill, never fabricated configurations. Exact retained indices and
original full-field closure audit are recorded. Replace only after verification.
"""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json,math,fcntl,time,shutil,hashlib,uuid,traceback
from pathlib import Path
import h5py,numpy as np
import run
import storage_admission
from analyze_calibration import helper
OUT=Path(__file__).resolve().parent

class StorageAdmissionDeferred(RuntimeError):pass

def append_event(record):
    with (OUT/'retention-events.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        with (OUT/'retention-events.jsonl').open('a') as f:
            f.write(json.dumps(record)+'\n');f.flush();os.fsync(f.fileno())
        sync_directory(OUT)

def sync_directory(path):
    fd=os.open(path,os.O_RDONLY);os.fsync(fd);os.close(fd)

def reconcile(path):
    ledger=OUT/'retention-records'/f'{path.stem}.json'
    record=json.loads(ledger.read_text())
    with h5py.File(path,'r') as h:
        assert len(h['sweep'])==record['original_fields'],'Archive was extended; require a new audited retention record. Resume from a working copy.'
    if not record.get('committed'):
        actual=run.sha256(path)
        if actual!=record['compact_sha256'] and actual==record.get('source_sha256') and record.get('previous_retention'):
            # A crash before the second atomic replacement left the previous
            # committed archive intact. Preserve the attempt and restore its ledger.
            attempts=OUT/'retention-attempts';attempts.mkdir(exist_ok=True)
            run.atomic_json(attempts/(path.stem+'-'+record['compact_sha256']+'.json'),record);sync_directory(attempts)
            record=record['previous_retention'];run.atomic_json(ledger,record);sync_directory(ledger.parent)
            return record
        assert actual==record['compact_sha256'],'Replaced file differs from verified prepared record'
        append_event(dict(id=path.stem,recovered=True,**record))
        record['committed']=True;run.atomic_json(ledger,record);sync_directory(ledger.parent)
    return record

def outside_evidence(path,audit):
    manifest=OUT/'outside-bracket-retention-design-20261014.json'
    if not manifest.exists():return None
    raw=manifest.read_bytes();design=json.loads(raw);entry=design['eligible'].get(path.stem)
    if entry is None:return None
    assert entry['audit_sha256']==audit['sha256'],'Outside-bracket manifest audit mismatch'
    return dict(manifest=str(manifest),manifest_sha256=hashlib.sha256(raw).hexdigest(),entry=entry)

def selection(h,evidence=None):
    if evidence is not None:
        n=len(h['sweep']);indices=sorted(set(range(0,n,100))|{0,n-1})
        return np.array(indices,dtype=np.int64),dict(mode='checkpoints_only_outside_brackets',checkpoint_interval=100,diagnostic_windows=[],outside_bracket_evidence=evidence)

    n=len(h['sweep']);p=json.loads(h.attrs['parameters']);cut=min(p['warmup']+1,n//2)
    m=h['magnetization_squared'][:];v=h['magnetization_fourth'][:];z=h['structure_zero'][:];f=h['structure_min'][:].mean(1)
    a=m[cut:];b=v[cut:];zz=z[cut:];ff=f[cut:]
    series={'m2':a,'m4':b,'binder_influence':b.mean()/a.mean()**3*(a-a.mean())-(b-b.mean())/(2*a.mean()**2),'xi_influence':(zz-zz.mean())/ff.mean()-zz.mean()/ff.mean()**2*(ff-ff.mean())}
    for name in ['mean_link_cosine','vortex_density','action_density']:series[name]=h[name][cut:]
    w=np.empty(n)
    for first in range(0,n,100):w[first:first+100]=(h['winding'][first:first+100].astype(float)**2).mean((1,2))
    series['winding2']=w[cut:];taus={}
    for name,s in series.items():
        for label,part in [('postwarmup',s),('latehalf',s[len(s)//2:])]:
            t=helper.tau_int(part) if len(part)>20 else None
            taus[name+'_'+label]=None if t is None or not np.isfinite(t) else float(t)
    tau=max([1.]+[v for v in taus.values() if v is not None]);stride=max(1,math.ceil(10*tau));width=max(64,math.ceil(50*tau))
    keep=set(range(0,n,100));keep.update([0,n-1]);keep.update(range(cut,n,stride));windows=[]
    for first in [0,cut,cut+(n-cut)//2-width//2,n-width]:
        first=max(0,min(first,n-1));last=min(n,first+width);keep.update(range(first,last));windows.append([first,last])
    return np.array(sorted(keep),dtype=np.int64),dict(tau_estimates=taus,max_measured_tau=tau,stride=stride,diagnostic_windows=windows,checkpoint_interval=100,independence='Spacing reduces measured correlation; does not certify independence for unmeasured/constant observables.')

def equal(a,b):
    a,b=np.asarray(a),np.asarray(b)
    return a.dtype==b.dtype and a.shape==b.shape and (np.array_equal(a,b) if a.dtype.kind in 'OUS' else a.tobytes()==b.tobytes())

def compact(path,audit,parallel=False):
    path=Path(path)
    copy_lock=(OUT/'retention-copy.lock').open('a+')
    fcntl.flock(copy_lock,fcntl.LOCK_SH if parallel else fcntl.LOCK_EX)
    (OUT/'jobs').mkdir(exist_ok=True)
    lock=(OUT/'jobs'/f'{path.stem}.lock').open('a+')
    fcntl.flock(lock,fcntl.LOCK_EX)
    try:
        control=OUT/'control.json'
        if control.exists() and not parallel:
            config=json.loads(control.read_text());config['retention_worker_active']=True
            config['max_workers']=run.worker_limit(run.datetime.now(run.CENTRAL),18)-1
            run.atomic_json(control,config)
        with h5py.File(path,'r') as h:
            previously_retained=h.attrs.get('adaptive_field_retention_version')==1
            evidence=outside_evidence(path,audit)
            prior=reconcile(path) if previously_retained else None
            if prior is not None and (evidence is None or prior['policy'].get('mode')=='checkpoints_only_outside_brackets'):return prior
            assert audit['passed'] and len(h['phi'])==len(h['sweep'])
        job=OUT/'jobs'/f'{path.stem}.json'
        if job.exists():assert json.loads(job.read_text())['state']=='complete','Cannot compact an unfinished trajectory'
        original=run.sha256(path);assert original==(prior['compact_sha256'] if prior is not None else audit['sha256']),'Full-field audit hash mismatch'
        old_bytes=path.stat().st_size;temp=path.with_name('.'+path.name+'.'+uuid.uuid4().hex+'.retention.tmp')
        with h5py.File(path,'r') as src:
            indices,policy=selection(src,evidence)
            if prior is not None:
                available=set(map(int,src['phi'].attrs['retained_indices']))
                assert set(map(int,indices))<=available,'Requested checkpoint was not retained previously'
                assert all(np.isfinite(src['phi'][int(i)]).all() for i in indices),'Missing checkpoint configuration'

            admission=storage_admission.acquire()
            try:
                # The trajectory lock excludes another owner for THIS archive.
                # Never clear a different active converter's reservation.
                storage_admission.release(path.stem)
                remaining=0
                for job_path in (OUT/'jobs').glob('*.json'):
                    j=json.loads(job_path.read_text())
                    if j['state']=='running':
                        live=OUT/'chains'/f"{j['point']['id']}.h5";remaining+=max(0,storage_admission.estimated_full_bytes(j['point'])-(live.stat().st_size if live.exists() else 0))
                field_bytes=src['phi'].size*src['phi'].dtype.itemsize
                estimate=int(1.15*(old_bytes-src['phi'].id.get_storage_size()+field_bytes*len(indices)/len(src['phi']))+128*2**20)
                # Snapshot all remaining growth before reading physical free.
                if not storage_admission.admit_retention_copy(path.stem,temp,estimate,remaining):
                    raise StorageAdmissionDeferred('Retention temporary copy would breach disk reserve')
            finally:admission.close()
            # A unique attempt path is created only AFTER its lease is admitted.
            with h5py.File(temp,'w') as dst:
                for k,v in src.attrs.items():dst.attrs[k]=v
                for name,d in src.items():
                    options=dict(shape=d.shape,maxshape=(None,)+d.shape[1:],dtype=d.dtype,track_order=True)
                    if name=='phi':options.update(chunks=(1,)+d.shape[1:],fillvalue=np.nan,compression='gzip',compression_opts=1,shuffle=True,fletcher32=True)
                    elif d.dtype.kind!='O':options.update(chunks=True,compression='gzip',compression_opts=1,shuffle=True,fletcher32=True)
                    else:options.update(chunks=True)
                    target=dst.create_dataset(name,**options)
                    for k,v in d.attrs.items():target.attrs[k]=v
                    if name=='phi':
                        for i in indices:target[int(i)]=d[int(i)]
                        target.attrs['retained_indices']=indices;target.attrs['original_rows']=len(d);target.attrs['missing_fields']='Explicit NaN fill; consult retained_indices. Resume exact final checkpoint from a working copy; production archive is immutable.'
                    else:
                        for first in range(0,len(d),100):target[first:first+100]=d[first:first+100]
                record=dict(original_sha256=prior['original_sha256'] if prior else original,original_bytes=prior['original_bytes'] if prior else old_bytes,original_full_field_audit=audit,retained_fields=len(indices),original_fields=len(src['phi']),policy=policy,created=run.now(),source_sha256=original,source_bytes=old_bytes)
                if prior is not None:record['previous_retention']=prior
                dst.attrs['adaptive_field_retention_version']=1;dst.attrs['retention_record']=json.dumps(record);dst.flush()
        with h5py.File(path,'r') as src,h5py.File(temp,'r') as dst:
            assert set(src)==set(dst)
            for key in src.attrs:
                if key!='retention_record':assert equal(src.attrs[key],dst.attrs[key]),key
            for name,d in src.items():
                other=dst[name];assert d.shape==other.shape and d.dtype==other.dtype
                for key in d.attrs:
                    if name!='phi' or key not in ['retained_indices','original_rows','missing_fields']:assert equal(d.attrs[key],other.attrs[key]),key
                for first in range(0,len(indices),16) if name=='phi' else range(0,len(d),100):
                    rows=indices[first:first+16] if name=='phi' else slice(first,first+100)
                    assert equal(d[rows],other[rows]),name
            assert all(i in indices for i in range(0,len(src['sweep']),100))
            np.testing.assert_array_equal(src['phi'][-1],dst['phi'][-1])
        with temp.open('rb') as f:os.fsync(f.fileno())
        record.update(compact_sha256=run.sha256(temp),compact_bytes=temp.stat().st_size,verified=True,committed=False)
        (OUT/'retention-records').mkdir(exist_ok=True);sync_directory(OUT)
        # A durable prepared record allows recovery after replace/event interruption.
        run.atomic_json(OUT/'retention-records'/f'{path.stem}.json',record)
        sync_directory(OUT/'retention-records')
        os.replace(temp,path)
        sync_directory(path.parent)
        append_event(dict(id=path.stem,**record))
        record['committed']=True;run.atomic_json(OUT/'retention-records'/f'{path.stem}.json',record);sync_directory(OUT/'retention-records')
        return record
    finally:
        admission=storage_admission.acquire()
        try:storage_admission.release(path.stem)
        finally:admission.close()
        lock.close();copy_lock.close()

def needs_conversion(job,eligible):
    identifier=job['point']['id'];ledger=OUT/'retention-records'/f'{identifier}.json'
    if not ledger.exists():return True
    record=json.loads(ledger.read_text())
    if not record.get('committed'):return True
    assert record['original_sha256']==job['audit']['sha256'],'Retention ledger/source audit mismatch'
    return identifier in eligible and record['policy'].get('mode')!='checkpoints_only_outside_brackets'

def parallel_convert(root,identifier,audit):
    global OUT
    OUT=Path(root);storage_admission.OUT=OUT
    try:
        r=compact(OUT/'chains'/f'{identifier}.h5',audit,parallel=True)
        return dict(id=identifier,status='committed',original_bytes=r['original_bytes'],compact_bytes=r['compact_bytes'],retained_fields=r['retained_fields'],pid=os.getpid())
    except StorageAdmissionDeferred as error:
        return dict(id=identifier,status='waiting_for_space',reason=str(error),pid=os.getpid())

def account_workers(workers):
    control=OUT/'control.json'
    if control.exists():
        config=json.loads(control.read_text());config['retention_worker_active']=workers>0
        config['retention_workers']=workers
        config['max_workers']=max(0,run.worker_limit(run.datetime.now(run.CENTRAL),18)-workers)
        run.atomic_json(control,config)

def daemon(workers=1):
    assert 1<=workers<=14,'Leave four cores of headroom on the 18-core host'
    from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
    import multiprocessing
    (OUT/'retention-records').mkdir(exist_ok=True)
    with (OUT/'retention-daemon.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        # EX waits for all surviving SH copies, so clearing stale leases is safe.
        with (OUT/'retention-copy.lock').open('a+') as copy_lock:
            fcntl.flock(copy_lock,fcntl.LOCK_EX)
            admission=storage_admission.acquire()
            try:
                storage_admission.clear_retention_reservations()
            finally:admission.close()
            stop=OUT/'retention-stop-request.json'
            if stop.exists():stop.unlink()
        # Spawn after releasing EX; do not inherit maintenance lock descriptors.
        futures={};inflight=set();deferred={};failed=[];converted=0;last_refresh=0;jobs=[];eligible={}
        with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as pool:
            while True:
                stopping=(OUT/'retention-stop-request.json').exists() or bool(failed)
                if not stopping and time.monotonic()-last_refresh>30:
                    import design_outside_bracket_retention as outside_design
                    outside_design.OUT=OUT;outside_design.refresh()
                    jobs=[json.loads(f.read_text()) for f in (OUT/'jobs').glob('*.json')]
                    eligible=json.loads((OUT/'outside-bracket-retention-design-20261014.json').read_text())['eligible']
                    jobs=sorted([j for j in jobs if j['state']=='complete' and j.get('audit',{}).get('passed')],key=lambda j:(j['point']['id'] not in eligible,-j['audit']['bytes']))
                    last_refresh=time.monotonic()
                if not stopping:
                    for j in jobs:
                        if len(futures)>=workers:break
                        identifier=j['point']['id']
                        if identifier in inflight or deferred.get(identifier,0)>time.monotonic():continue
                        if not needs_conversion(j,eligible):continue
                        future=pool.submit(parallel_convert,str(OUT),identifier,j['audit']);futures[future]=identifier;inflight.add(identifier)
                account_workers(len(futures))
                run.atomic_json(OUT/'retention-status.json',dict(updated=run.now(),pid=os.getpid(),state='draining' if stopping else ('converting' if futures else 'monitoring'),worker_limit=workers,active=len(futures),inflight=sorted(inflight),converted_this_run=converted,failures=failed,waiting_for_space=[i for i,t in deferred.items() if t>time.monotonic()]))
                if stopping and not futures:break
                if not futures:time.sleep(15);continue
                done,_=wait(futures,timeout=5,return_when=FIRST_COMPLETED)
                for future in done:
                    identifier=futures.pop(future);inflight.remove(identifier)
                    try:
                        result=future.result();print(run.now(),json.dumps(result),flush=True)
                        if result['status']=='waiting_for_space':deferred[identifier]=time.monotonic()+30
                        else:converted+=1;deferred.pop(identifier,None)
                    except BaseException as error:
                        failure=dict(id=identifier,error=repr(error),traceback=''.join(traceback.format_exception(error)),created=run.now());failed.append(failure);print('RETENTION FAILURE',json.dumps(failure),flush=True)
                        run.atomic_json(OUT/'retention-failures'/f'{identifier}-{os.getpid()}.json',failure)
            account_workers(0)
        if failed:raise RuntimeError(f'Parallel retention failed; preserved originals and journals: {failed}')

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--id');p.add_argument('--daemon',action='store_true');p.add_argument('--workers',type=int,default=1);a=p.parse_args();(OUT/'retention-records').mkdir(exist_ok=True)
    if a.daemon:daemon(a.workers)
    else:
        j=json.loads((OUT/'jobs'/f'{a.id}.json').read_text());print(compact(OUT/'chains'/f'{a.id}.h5',j['audit']),flush=True)
