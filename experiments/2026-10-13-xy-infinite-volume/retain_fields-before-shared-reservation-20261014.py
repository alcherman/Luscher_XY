"""Verified adaptive field retention. Observables/RNG histories remain complete.

phi retains its logical shape for native checkpoint replay. Unstored fields are
explicit NaN fill, never fabricated configurations. Exact retained indices and
original full-field closure audit are recorded. Replace only after verification.
"""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json,math,fcntl,time,shutil
from pathlib import Path
import h5py,numpy as np
import run
from analyze_calibration import helper
OUT=Path(__file__).resolve().parent

def sync_directory(path):
    fd=os.open(path,os.O_RDONLY);os.fsync(fd);os.close(fd)

def reconcile(path):
    ledger=OUT/'retention-records'/f'{path.stem}.json'
    record=json.loads(ledger.read_text())
    with h5py.File(path,'r') as h:
        assert len(h['sweep'])==record['original_fields'],'Archive was extended; require a new audited retention record. Resume from a working copy.'
    if not record.get('committed'):
        assert run.sha256(path)==record['compact_sha256'],'Replaced file differs from verified prepared record'
        with (OUT/'retention-events.jsonl').open('a') as f:
            f.write(json.dumps(dict(id=path.stem,recovered=True,**record))+'\n');f.flush();os.fsync(f.fileno())
        sync_directory(OUT)
        record['committed']=True;run.atomic_json(ledger,record);sync_directory(ledger.parent)
    return record

def selection(h):
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

def compact(path,audit):
    path=Path(path)
    (OUT/'jobs').mkdir(exist_ok=True)
    lock=(OUT/'jobs'/f'{path.stem}.lock').open('a+')
    fcntl.flock(lock,fcntl.LOCK_EX)
    try:
        control=OUT/'control.json'
        if control.exists():
            config=json.loads(control.read_text());config['retention_worker_active']=True
            config['max_workers']=run.worker_limit(run.datetime.now(run.CENTRAL),18)-1
            run.atomic_json(control,config)
        with h5py.File(path,'r') as h:
            if h.attrs.get('adaptive_field_retention_version')==1:return reconcile(path)
            assert audit['passed'] and len(h['phi'])==len(h['sweep'])
        job=OUT/'jobs'/f'{path.stem}.json'
        if job.exists():assert json.loads(job.read_text())['state']=='complete','Cannot compact an unfinished trajectory'
        original=run.sha256(path);assert original==audit['sha256'],'Full-field audit hash mismatch'
        old_bytes=path.stat().st_size;temp=path.with_name('.'+path.name+'.retention.tmp')
        with h5py.File(path,'r') as src,h5py.File(temp,'w') as dst:
            indices,policy=selection(src)
            remaining=0
            for job_path in (OUT/'jobs').glob('*.json'):
                j=json.loads(job_path.read_text())
                if j['state']=='running':
                    live=OUT/'chains'/f"{j['point']['id']}.h5";remaining+=max(0,run.estimated_bytes(j['point'])-(live.stat().st_size if live.exists() else 0))
            field_bytes=src['phi'].size*src['phi'].dtype.itemsize
            estimate=int(1.15*(old_bytes-src['phi'].id.get_storage_size()+field_bytes*len(indices)/len(src['phi']))+128*2**20)
            assert shutil.disk_usage(path.parent).free-remaining-estimate>=25*2**30,'Retention temporary copy would breach disk reserve'
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
            record=dict(original_sha256=original,original_bytes=old_bytes,original_full_field_audit=audit,retained_fields=len(indices),original_fields=len(src['phi']),policy=policy,created=run.now())
            dst.attrs['adaptive_field_retention_version']=1;dst.attrs['retention_record']=json.dumps(record);dst.flush()
        with h5py.File(path,'r') as src,h5py.File(temp,'r') as dst:
            assert set(src)==set(dst)
            for key in src.attrs:assert equal(src.attrs[key],dst.attrs[key]),key
            for name,d in src.items():
                other=dst[name];assert d.shape==other.shape and d.dtype==other.dtype
                for key in d.attrs:assert equal(d.attrs[key],other.attrs[key]),key
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
        with (OUT/'retention-events.jsonl').open('a') as f:f.write(json.dumps(dict(id=path.stem,**record))+'\n');f.flush();os.fsync(f.fileno())
        sync_directory(OUT)
        record['committed']=True;run.atomic_json(OUT/'retention-records'/f'{path.stem}.json',record);sync_directory(OUT/'retention-records')
        return record
    finally:lock.close()

def daemon():
    (OUT/'retention-records').mkdir(exist_ok=True)
    with (OUT/'retention-daemon.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        while True:
            jobs=[json.loads(f.read_text()) for f in (OUT/'jobs').glob('*.json')]
            jobs=sorted([j for j in jobs if j['state']=='complete' and j.get('audit',{}).get('passed')],key=lambda j:-j['audit']['bytes'])
            count=0
            for j in jobs:
                path=OUT/'chains'/f"{j['point']['id']}.h5"
                with h5py.File(path,'r') as h:
                    retained=h.attrs.get('adaptive_field_retention_version')==1
                if retained:
                    ledger=OUT/'retention-records'/f'{path.stem}.json';r=json.loads(ledger.read_text())
                    with h5py.File(path,'r') as h:assert len(h['sweep'])==r['original_fields'],'Archive extended without new audit/retention record'
                    if not r.get('committed'):compact(path,j['audit'])
                    continue
                r=compact(path,j['audit']);count+=1
                print(run.now(),path.stem,r['original_bytes'],r['compact_bytes'],r['retained_fields'],flush=True)
            control=OUT/'control.json'
            if control.exists():
                config=json.loads(control.read_text());config['max_workers']=18;config['retention_worker_active']=False;run.atomic_json(control,config)
            run.atomic_json(OUT/'retention-status.json',dict(updated=run.now(),state='monitoring',converted_this_pass=count))
            time.sleep(30)

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--id');p.add_argument('--daemon',action='store_true');a=p.parse_args();(OUT/'retention-records').mkdir(exist_ok=True)
    if a.daemon:daemon()
    else:
        j=json.loads((OUT/'jobs'/f'{a.id}.json').read_text());print(compact(OUT/'chains'/f'{a.id}.h5',j['audit']),flush=True)
