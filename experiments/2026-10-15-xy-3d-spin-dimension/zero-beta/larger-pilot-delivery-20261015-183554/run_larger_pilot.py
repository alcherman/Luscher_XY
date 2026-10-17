"""Finite restart-independent N48/N64 pilot, then qualified diagnostic snapshots."""
from pathlib import Path
import os,sys,json,time,fcntl,subprocess,concurrent.futures
import measure as m
from design_larger_zero_beta import NAME
from run_fine_scan import states,coordinator_pid,control,execute

OUT=m.OUT

def main():
    ps=json.loads((m.CAMPAIGN/'batches'/f'{NAME}.json').read_text())['points'];allowed={p['id'] for p in ps}
    record=dict(created=m.run.now(),pid=os.getpid(),expected_trajectories=len(ps),snapshots={})
    def publish(state,**extra):
        record.update(state=state,updated=m.run.now(),**extra);m.run.atomic_json(OUT/'larger-pilot-status.json',record)
    with (OUT/'larger-pilot-manager.lock').open('a+') as singleton:
        fcntl.flock(singleton,fcntl.LOCK_EX|fcntl.LOCK_NB)
        publish('waiting_for_previous_finite_manager')
        # Wait for the earlier manager to finish its final archival conversions.
        with (OUT/'fine-scan-manager.lock').open('a+') as study_lock:
            fcntl.flock(study_lock,fcntl.LOCK_EX)
            awake=subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(os.getpid())])
            try:
                for f in (m.CAMPAIGN/'batches').glob('*.json'):
                    for p in json.loads(f.read_text())['points']:
                        job=m.CAMPAIGN/'jobs'/f"{p['id']}.json";s=json.loads(job.read_text())['state'] if job.exists() else 'pending'
                        if s in ['pending','interrupted','running'] and p['id'] not in allowed:raise RuntimeError('Unrelated production commitment: '+p['id'])
                st,fail=states(ps)
                if fail:raise RuntimeError('Pilot contains simulation/audit failure: '+str(fail))
                if st!={'complete':len(ps)}:
                    control(launch_new=True,launch_allowed_ids=sorted(allowed),max_workers=14,current_scientific_priority='Authorized N48/N64 beta0 spin-dimension pilot only:28 reviewed20k jobs. No automatic extensions, longer runs or N96. Keep SSD originals.')
                    if coordinator_pid() is None:
                        stream=(OUT/'larger-pilot-coordinator.log').open('a');child=subprocess.Popen([sys.executable,'-u',str(m.CAMPAIGN/'coordinator_ssd.py')],cwd=m.run.ROOT,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
                        time.sleep(2)
                        if child.poll() is not None:raise RuntimeError('Coordinator exited on startup; inspect full log')
                    previous=None
                    while True:
                        st,fail=states(ps);publish('running_simulations',job_states=st,failed_ids=fail,coordinator_pid=coordinator_pid())
                        if st!=previous:print(m.run.now(),st,flush=True);previous=st
                        if fail:raise RuntimeError('Pilot simulation/audit failure: '+str(fail))
                        if st=={'complete':len(ps)}:break
                        if m.run.datetime.now(m.run.CENTRAL)>=m.run.DEADLINE:raise RuntimeError('Friday safety deadline reached; stopping new launches, preserving safe checkpoints')
                        if coordinator_pid() is None:raise RuntimeError('Coordinator stalled/exited')
                        running=[]
                        for p in ps:
                            f=m.CAMPAIGN/'jobs'/f"{p['id']}.json"
                            if f.exists():
                                j=json.loads(f.read_text())
                                if j['state']=='running':running.append(str(j['pid']))
                        if running:
                            with (OUT/'larger-pilot-resource-samples.log').open('a') as log:
                                log.write(m.run.now()+'\n');log.flush();subprocess.run(['/bin/ps','-p',','.join(running),'-o','pid=,rss=,pcpu=,command='],stdout=log,stderr=subprocess.STDOUT)
                        time.sleep(20)
                control(launch_new=False,launch_allowed_ids=[],current_scientific_priority='N48/N64 pilot all audited; scientific diagnostics only. No longer runs or extensions authorized.')
                pid=coordinator_pid()
                if pid is not None:os.kill(pid,15)
                publish('analysing')
                old=m.OUT/'zero-beta'/json.loads((m.OUT/'zero-beta/latest-fine-primary.json').read_text())['folder']/'analysis.json';baseline=json.loads(old.read_text())
                cohort=baseline['source_jobs']+[dict(id=p['id'],point=p,audit=json.loads((m.CAMPAIGN/'jobs'/f"{p['id']}.json").read_text())['audit']) for p in ps]
                cohort_path=OUT/'larger-pilot-analysis-cohort.json';m.run.atomic_json(cohort_path,dict(created=m.run.now(),previous_source=str(old),previous_sha256=m.run.sha256(old),source_jobs=cohort))
                expected={j['id'] for j in cohort}
                value=execute('larger-primary',['--source-analysis',str(cohort_path)],expected);record['snapshots'][value['tag']]=value;publish('analysing')
                source=value['path']
                variants={'larger-full':['--window','full'],'larger-cut4001':['--cut','4001'],'larger-cut8001':['--cut','8001'],'larger-block2':['--factor','2'],'larger-U20full':['--target','.20','--window','full'],'larger-U30full':['--target','.30','--window','full'],'larger-evenfull':['--starts','even','--window','full'],'larger-oddfull':['--starts','odd','--window','full'],'larger-only':['--id-prefix','spinpilot3_','--window','full'],'larger-only-late':['--id-prefix','spinpilot3_','--window','full','--cut','8001']}
                def one(item):
                    tag,args=item;ids=allowed if tag.startswith('larger-only') else {j['id'] for j in cohort if tag not in ['larger-evenfull','larger-oddfull'] or json.loads((m.CAMPAIGN/'jobs'/f"{j['id']}.json").read_text())['point']['replica']%2==int(tag=='larger-oddfull')}
                    return execute(tag,['--source-analysis',source]+args,ids)
                with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
                    for info in pool.map(one,variants.items()):record['snapshots'][info['tag']]=info;publish('analysing')
                for name in ['diagnose_larger_pilot.py','report_larger_pilot.py']:
                    with (OUT/(name.removesuffix('.py')+'.log')).open('a') as stream:
                        r=subprocess.run([sys.executable,str(OUT/name)],cwd=m.run.ROOT,stdout=stream,stderr=subprocess.STDOUT)
                    if r.returncode:raise RuntimeError(name+' failed; preserve all snapshots and inspect full log')
                publish('snapshots_complete_pending_scientific_review')
                def retain(p):
                    with (OUT/f"larger-retention-{p['id']}.log").open('a') as stream:r=subprocess.run([sys.executable,str(m.CAMPAIGN/'retain_fields.py'),'--id',p['id']],cwd=m.run.ROOT,stdout=stream,stderr=subprocess.STDOUT)
                    return p['id'],r.returncode
                failures=[]
                with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
                    for identifier,code in pool.map(retain,ps):
                        if code:failures.append(identifier)
                publish('snapshots_complete_pending_scientific_review',retention_failed_ids=failures)
                if failures:raise RuntimeError('Retention errors, originals/journals preserved: '+str(failures))
            except BaseException as error:
                control(launch_new=False,launch_allowed_ids=[]);publish('failed_needs_attention',error=repr(error));raise
            finally:
                if awake.poll() is None:awake.terminate()

if __name__=='__main__':main()
