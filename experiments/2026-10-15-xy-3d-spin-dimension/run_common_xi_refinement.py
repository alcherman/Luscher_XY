"""Finite restart-independent, reviewed xi crossing refinement and analysis."""
import concurrent.futures
import fcntl
import json
import os
import subprocess
import sys
import time
from pathlib import Path
import measure as m
from design_common_xi_refinement import NAME
from run_fine_scan import states,coordinator_pid,control

OUT=m.OUT


def analysis_file(path):
    path=Path(path)
    return path/'analysis.json' if path.is_dir() else path


def verify_no_unrelated(points):
    allowed={p['id'] for p in points}
    for f in (m.CAMPAIGN/'batches').glob('*.json'):
        for p in json.loads(f.read_text())['points']:
            job=m.CAMPAIGN/'jobs'/f"{p['id']}.json";state=json.loads(job.read_text())['state'] if job.exists() else 'pending'
            if state in ['pending','interrupted','running'] and p['id'] not in allowed:raise RuntimeError('Unrelated production commitment: '+p['id'])


def analyze(tag,args):
    with (OUT/f'common-xi-refined-{tag}.log').open('a') as stream:
        r=subprocess.run([sys.executable,str(OUT/'common_xi_quotients.py'),'--tag','refined-'+tag]+args,cwd=m.run.ROOT,stdout=stream,stderr=subprocess.STDOUT)
    if r.returncode:raise RuntimeError('Quotient analysis failed: '+tag)
    folder=OUT/json.loads((OUT/f'latest-common-xi-refined-{tag}.json').read_text())['folder']
    return dict(folder=str(folder),path=str(folder/'analysis.json'),sha256=m.run.sha256(folder/'analysis.json'))


def main():
    manifest=json.loads((m.CAMPAIGN/'batches'/f'{NAME}.json').read_text());ps=manifest['points'];allowed={p['id'] for p in ps}
    design=json.loads((OUT/'common-xi-refinement-design.json').read_text())
    provenance=json.loads((Path(design['source_ledger']).parent/'provenance.json').read_text())
    prior_analysis=analysis_file(provenance['source_analyses']['baseline']['path'])
    if m.run.sha256(prior_analysis)!=provenance['source_analyses']['baseline']['sha256']:raise RuntimeError('Previous scalar cohort source hash changed')
    expected={s['id'] for s in json.loads(prior_analysis.read_text())['source_jobs']}|allowed
    state=dict(created=m.run.now(),pid=os.getpid(),expected=len(ps),snapshots={})
    def publish(value,**extra):
        state.update(state=value,updated=m.run.now(),**extra);m.run.atomic_json(OUT/'common-xi-refinement-status.json',state)
    with (OUT/'common-xi-refinement-manager.lock').open('a+') as singleton:
        fcntl.flock(singleton,fcntl.LOCK_EX|fcntl.LOCK_NB)
        publish('waiting_for_previous_finite_manager')
        with (OUT/'fine-scan-manager.lock').open('a+') as study_lock:
            fcntl.flock(study_lock,fcntl.LOCK_EX)
            awake=subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(os.getpid())])
            try:
                verify_no_unrelated(ps);st,fail=states(ps)
                if fail:raise RuntimeError('Existing job/audit failure: '+str(fail))
                if st!={'complete':len(ps)}:
                    control(launch_new=True,launch_allowed_ids=sorted(allowed),max_workers=14,retention_workers=0,retention_worker_active=False,current_scientific_priority='User-authorized common-coupling xi/L epsilon refinement and24:48 bracket resolution only. Frozen120 trajectories; no larger boxes or unrelated jobs. SSD originals held.')
                    if coordinator_pid() is None:
                        stream=(OUT/'common-xi-refinement-coordinator.log').open('a')
                        proc=subprocess.Popen([sys.executable,'-u',str(m.CAMPAIGN/'coordinator_ssd.py')],cwd=m.run.ROOT,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
                        time.sleep(2)
                        if proc.poll() is not None:raise RuntimeError('Coordinator exited on startup')
                    previous_states=None
                    while True:
                        st,fail=states(ps);publish('running_simulations',job_states=st,failed_ids=fail,coordinator_pid=coordinator_pid())
                        if st!=previous_states:print(m.run.now(),st,flush=True);previous_states=st
                        if fail:raise RuntimeError('Simulation/audit failure: '+str(fail))
                        if st=={'complete':len(ps)}:break
                        if m.run.datetime.now(m.run.CENTRAL)>=m.run.DEADLINE:raise RuntimeError('Friday safety deadline reached; no new launches, preserve safe checkpoints')
                        if coordinator_pid() is None:raise RuntimeError('Coordinator stalled/exited')
                        time.sleep(15)
                control(launch_new=False,launch_allowed_ids=[],current_scientific_priority='All common-coupling refinement trajectories audited; scalar analysis only, no further launches.')
                pid=coordinator_pid()
                if pid is not None:os.kill(pid,15)
                publish('analysing')
                cohort=[]
                for identifier in sorted(expected):
                    j=json.loads((m.CAMPAIGN/'jobs'/f'{identifier}.json').read_text())
                    if j['state']!='complete' or not j.get('audit',{}).get('passed'):raise RuntimeError('Expected cohort is incomplete: '+identifier)
                    cohort.append(dict(id=identifier,point=j['point'],audit=j['audit']))
                pinned=OUT/'common-xi-refinement-analysis-cohort.json'
                m.run.atomic_json(pinned,dict(created=m.run.now(),previous_analysis=str(prior_analysis),previous_sha256=m.run.sha256(prior_analysis),source_jobs=cohort))
                baseline=analyze('primary',['--source-analysis',str(pinned)]);state['snapshots']['primary']=baseline;publish('analysing')
                cohort=json.loads(Path(baseline['path']).read_text());cohort_ids={s['id'] for s in cohort['source_jobs']}
                if cohort_ids!=expected:raise RuntimeError('Primary analysis cohort differs from reviewed prior+new sources')
                variants={'block4':['--factor','4'],'latecut':['--cut','4001'],'even':['--starts','even'],'odd':['--starts','odd']}
                def one(item):
                    tag,args=item;return tag,analyze(tag,['--source-analysis',baseline['path']]+args)
                with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
                    for tag,info in pool.map(one,variants.items()):state['snapshots'][tag]=info;publish('analysing')
                coverage=[]
                for tag,info in state['snapshots'].items():
                    a=json.loads(Path(info['path']).read_text())
                    for r in a['pairs']:
                        if r['N']==24 and r['window']=='local':coverage.append(dict(variant=tag,method=r['method'],status=r['status'],valid_fraction=r.get('valid_fraction',0.),epsilon_over_pi=r.get('epsilon_over_pi')))
                state['crossing_24_48_checks']=coverage
                state['crossing_24_48_covered']=all(r['status']=='conditional_quotient' for r in coverage if r['variant'] in ['primary','block4','latecut']) and len([r for r in coverage if r['variant'] in ['primary','block4','latecut']])==6
                with (OUT/'common-xi-refined-delivery-build.log').open('a') as stream:
                    r=subprocess.run([sys.executable,str(OUT/'report_common_xi.py'),baseline['folder'],state['snapshots']['block4']['folder']],cwd=m.run.ROOT,stdout=stream,stderr=subprocess.STDOUT)
                if r.returncode:raise RuntimeError('Quotient delivery build failed')
                delivery=OUT/json.loads((OUT/'latest-common-xi-delivery.json').read_text())['folder']
                m.run.atomic_json(delivery/'additional-snapshots.json',state['snapshots'])
                m.run.atomic_json(delivery/'crossing-24-48-coverage-checks.json',dict(checks=coverage,covered=state['crossing_24_48_covered'],scope='Conditional local linear/PCHIP crossing coverage, not exponent accuracy or mixing certification. If false, further measured refinement is still required.'))
                publish('snapshots_complete_pending_scientific_review' if state['crossing_24_48_covered'] else 'crossing_still_needs_measured_refinement',delivery=str(delivery))
                def retain(p):
                    with (OUT/f"common-xi-retention-{p['id']}.log").open('a') as stream:
                        r=subprocess.run([sys.executable,str(m.CAMPAIGN/'retain_fields.py'),'--id',p['id']],cwd=m.run.ROOT,stdout=stream,stderr=subprocess.STDOUT)
                    return p['id'],r.returncode
                failures=[]
                with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
                    for identifier,code in pool.map(retain,ps):
                        if code:failures.append(identifier)
                publish('snapshots_complete_pending_scientific_review' if state['crossing_24_48_covered'] else 'crossing_still_needs_measured_refinement',retention_failed_ids=failures)
                if failures:raise RuntimeError('Retention failure; originals/journals preserved: '+str(failures))
            except BaseException as error:
                control(launch_new=False,launch_allowed_ids=[]);publish('failed_needs_attention',error=repr(error));raise
            finally:
                if awake.poll() is None:awake.terminate()


if __name__=='__main__':main()
