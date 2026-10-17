"""Finite detached finer-scan runner; immutable production kernels stay unchanged."""
from pathlib import Path
import os,sys,json,time,fcntl,subprocess,concurrent.futures
import measure as m
from design_fine_zero_beta import NAME

OUT=m.OUT

def states(points):
    result={};failed=[]
    for p in points:
        path=m.CAMPAIGN/'jobs'/f"{p['id']}.json";j=json.loads(path.read_text()) if path.exists() else {'state':'pending'}
        result[j['state']]=result.get(j['state'],0)+1
        if j['state']=='failed' or (j['state']=='complete' and not j.get('audit',{}).get('passed')):failed.append(p['id'])
    return result,failed

def coordinator_pid():
    f=m.CAMPAIGN/'status.json'
    pid=json.loads(f.read_text()).get('pid') if f.exists() else None
    if not isinstance(pid,int):return None
    r=subprocess.run(['/bin/ps','-p',str(pid),'-o','command='],capture_output=True,text=True)
    return pid if r.returncode==0 and str(m.CAMPAIGN/'coordinator_ssd.py') in r.stdout else None

def control(**changes):
    c=m.CAMPAIGN/'control.json';v=json.loads(c.read_text());v.update(changes);m.run.atomic_json(c,v)

def execute(tag,args,expected_ids=None):
    log=OUT/f'fine-scan-analysis-{tag}.log'
    with log.open('a') as f:r=subprocess.run([sys.executable,str(OUT/'zero_beta.py'),'--tag',tag]+args,cwd=m.run.ROOT,stdout=f,stderr=subprocess.STDOUT)
    if r.returncode:raise RuntimeError(f'Analysis failed: {tag}; see {log}')
    p=OUT/'zero-beta'/json.loads((OUT/'zero-beta'/f'latest-{tag}.json').read_text())['folder']/'analysis.json'
    if expected_ids is not None:
        actual={j['id'] for j in json.loads(p.read_text())['source_jobs']}
        if actual!=set(expected_ids):raise RuntimeError(f'Incomplete or changed source cohort: {tag}')
    return dict(tag=tag,path=str(p),sha256=m.run.sha256(p))

def main():
    with (OUT/'fine-scan-manager.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        ps=json.loads((m.CAMPAIGN/'batches'/f'{NAME}.json').read_text())['points']
        record=dict(pid=os.getpid(),created=m.run.now(),expected=len(ps),snapshots={})
        def publish(state,**extra):
            record.update(state=state,updated=m.run.now(),**extra);m.run.atomic_json(OUT/'fine-scan-status.json',record)
        awake=subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(os.getpid())])
        try:
            # Preserve all deferred jobs, and refuse to open a launch gate for unrelated work.
            allowed={p['id'] for p in ps}
            for f in (m.CAMPAIGN/'batches').glob('*.json'):
                for p in json.loads(f.read_text())['points']:
                    j=m.CAMPAIGN/'jobs'/f"{p['id']}.json";s=json.loads(j.read_text())['state'] if j.exists() else 'pending'
                    if s in ['pending','interrupted','running'] and p['id'] not in allowed:raise RuntimeError('Unrelated production commitment: '+p['id'])
            st,fail=states(ps)
            if fail:raise RuntimeError('Production/audit failures: '+str(fail))
            if st!={'complete':len(ps)}:
                control(launch_new=True,max_workers=14,current_scientific_priority='User-authorized October15 finer3D beta0 spin dimension scan; only82 reviewed stage1 jobs; no D2/D5/newN1024 work. Keep SSD originals.')
                if coordinator_pid() is None:
                    # The frozen coordinator itself owns its nonblocking global singleton lock.
                    log=(OUT/'fine-scan-coordinator.log').open('a')
                    child=subprocess.Popen([sys.executable,'-u',str(m.CAMPAIGN/'coordinator_ssd.py')],cwd=m.run.ROOT,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
                    time.sleep(2)
                    if child.poll() is not None:raise RuntimeError('Coordinator did not remain alive; inspect full log')
                previous=None
                while True:
                    st,fail=states(ps);publish('running_simulations',job_states=st,failed_ids=fail,coordinator_pid=coordinator_pid())
                    if st!=previous:print(m.run.now(),st,flush=True);previous=st
                    if fail:raise RuntimeError('Production/audit failures: '+str(fail))
                    if st=={'complete':len(ps)}:break
                    if m.run.datetime.now(m.run.CENTRAL)>=m.run.DEADLINE:raise RuntimeError('Prior Friday deadline reached; no new launches; safe checkpoints preserved')
                    if coordinator_pid() is None:raise RuntimeError('Coordinator stalled/exited')
                    time.sleep(15)
            control(launch_new=False,current_scientific_priority='Fine3D beta0 stage1 all audited; source-pinned analysis. No new launches.')
            pid=coordinator_pid()
            if pid is not None:os.kill(pid,15)
            publish('analysing')
            primary=execute('fine-primary',['--id-prefix','spinzero3_'],allowed);record['snapshots']['fine-primary']=primary;publish('analysing')
            variants={'fine-full':['--window','full'],'fine-cut2001':['--cut','2001'],'fine-cut4001':['--cut','4001'],'fine-block2':['--factor','2'],'fine-U20':['--target','.20'],'fine-U30':['--target','.30'],'fine-even':['--starts','even'],'fine-odd':['--starts','odd']}
            def one(item):
                tag,args=item;expected={p['id'] for p in ps if tag not in ['fine-even','fine-odd'] or p['replica']%2==int(tag=='fine-odd')}
                return execute(tag,['--source-analysis',primary['path']]+args,expected)
            with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
                for value in pool.map(one,variants.items()):record['snapshots'][value['tag']]=value;publish('analysing')
            value=execute('fine-pooled',[]);record['snapshots'][value['tag']]=value
            publish('snapshots_complete_pending_scientific_review')
            # Existing user-authorized retention only after all full closure audits.
            def retain(p):
                log=OUT/f"fine-retention-{p['id']}.log"
                with log.open('a') as stream:r=subprocess.run([sys.executable,str(m.CAMPAIGN/'retain_fields.py'),'--id',p['id']],cwd=m.run.ROOT,stdout=stream,stderr=subprocess.STDOUT)
                return p['id'],r.returncode
            failures=[]
            with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
                for identifier,code in pool.map(retain,ps):
                    if code:failures.append(identifier)
            publish('snapshots_complete_pending_scientific_review',retention_failed_ids=failures)
            if failures:raise RuntimeError('Retention errors, original/journal preserved: '+str(failures))
        except BaseException as error:
            control(launch_new=False);publish('failed_needs_attention',error=repr(error));raise
        finally:
            if awake.poll() is None:awake.terminate()

if __name__=='__main__':main()
