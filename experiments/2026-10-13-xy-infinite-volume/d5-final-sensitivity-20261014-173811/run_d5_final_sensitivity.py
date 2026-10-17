"""Restart-independent, source-pinned D5 sensitivity snapshot orchestration."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[key]='1'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json,sys,time,fcntl,subprocess,re,shutil,argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import run
OUT=Path(__file__).resolve().parent

def frozen_points():
    return {p['id']:p for f in sorted((OUT/'batches').glob('*.json')) for p in json.loads(f.read_text())['points'] if p['D']==5}

def coverage(points):
    states={};failed=[]
    for ident in points:
        path=OUT/'jobs'/(ident+'.json');j=json.loads(path.read_text()) if path.exists() else dict(state='pending')
        state=j['state'];states[state]=states.get(state,0)+1
        if state=='failed' or (state=='complete' and not j.get('audit',{}).get('passed')):failed.append(ident)
    return states,failed

def compare_fit(baseline,variant):
    same=baseline['included_root_sizes']==variant['included_root_sizes'];a,b=baseline['result'],variant['result']
    exists=all('location_infinite' in x for x in [a,b])
    return dict(same_sizes=same,baseline_sizes=baseline['included_root_sizes'],variant_sizes=variant['included_root_sizes'],baseline_status=a['status'],variant_status=b['status'],central_shift=b['location_infinite']-a['location_infinite'] if same and exists else None,variant_MC95=b.get('MC95'),variant_p=b.get('p'),variant_joint_valid_fraction=b.get('joint_valid_fraction'),variant_flags=b.get('input_flags',[]),both_conditional_fits=same and a['status']==b['status']=='conditional_D5_size_fit',limitation='Central shift is descriptive; no paired-draw or independent-replication significance is inferred. Different/missing sizes are incomplete comparisons.')

def fit_key(f):return (f['axis'],f['q'],f['method'],f['minimum_N'],f['correction'])

def main(wait=False,parallel=3):
    with (OUT/'d5-final-sensitivity-manager.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        previous=None
        while True:
            points=frozen_points()
            if not points:raise RuntimeError('No frozen D5 cohort')
            states,failed=coverage(points)
            record=dict(updated=run.now(),pid=os.getpid(),state='waiting_for_audits',expected_trajectories=len(points),job_states=states,failed_ids=failed)
            run.atomic_json(OUT/'d5-final-sensitivity-manager-status.json',record)
            if failed:raise RuntimeError('D5 simulation/audit failure; refusing final cohort')
            if states=={'complete':len(points)}:break
            if states!=previous:print(run.now(),'waiting for full D5 closure',states,flush=True);previous=states
            if not wait:raise RuntimeError('D5 production still active; use --wait')
            if run.datetime.now(run.CENTRAL)>=run.DEADLINE:raise RuntimeError('Deadline reached before final cohort closed; deliver attained/partial separately')
            time.sleep(30)
        stamp=run.now()[:19].replace('-','').replace(':','').replace('T','-');folder=OUT/('d5-final-sensitivity-'+stamp);folder.mkdir()
        shutil.copyfile(__file__,folder/Path(__file__).name)
        manifest=dict(created=run.now(),pid=os.getpid(),source_ids=sorted(points),source_count=len(points),snapshots={},scope='Pinned D5 source cohort. Conditional MC and descriptive stability only; zeroDOF/rejected fits and missing root coverage remain explicit. Shared starts/axes/models are not independent evidence.')
        def publish(state):
            run.atomic_json(folder/'manifest.json',manifest)
            run.atomic_json(OUT/'d5-final-sensitivity-manager-status.json',dict(updated=run.now(),pid=os.getpid(),state=state,folder=folder.name,source_count=len(points),completed=list(manifest['snapshots'])))
        def execute(tag,args):
            log=OUT/'logs'/f'd5-final-sensitivity-{stamp}-{tag}.log'
            with log.open('w') as stream:
                result=subprocess.run([sys.executable,str(OUT/'analyze_d5_calibration.py')]+args,cwd=run.ROOT,stdout=stream,stderr=subprocess.STDOUT)
            found=re.findall(r'^'+re.escape(str(OUT))+r'/d5-calibration-analysis-\S+',log.read_text(),re.MULTILINE)
            path=Path(found[-1]) if found else None
            info=dict(returncode=result.returncode,log=str(log),folder=path.name if path else None)
            if result.returncode==0 and path:
                analysis=json.loads((path/'analysis.json').read_text());starts=analysis['settings']['starts'];expected={i for i,p in points.items() if starts=='all' or p['replica']%2==int(starts=='odd')}
                info.update(settings=analysis['settings'],analysis_sha256=run.sha256(path/'analysis.json'),source_ids_match_expected={j['id'] for j in analysis['source_jobs']}==expected)
            print(run.now(),'snapshot',tag,info,flush=True)
            return tag,info
        publish('primary')
        tag,primary=execute('primary',[]);manifest['snapshots'][tag]=primary
        if primary['returncode'] or not primary.get('source_ids_match_expected'):raise RuntimeError('Primary failed or its source cohort changed')
        baseline_path=OUT/primary['folder']/'analysis.json';baseline=json.loads(baseline_path.read_text());manifest['baseline']=primary['folder'];publish('variants')
        variants=[('full',['--window','full']),('cut2001',['--cut','2001']),('cut3001',['--cut','3001']),('block2',['--factor','2']),('U20',['--target','.20']),('U30',['--target','.30']),('even',['--starts','even']),('odd',['--starts','odd'])]
        common=['--no-promote','--source-analysis',str(baseline_path)]
        with ThreadPoolExecutor(max_workers=parallel) as pool:
            futures=[pool.submit(execute,tag,args+common) for tag,args in variants]
            for future in as_completed(futures):
                tag,info=future.result();manifest['snapshots'][tag]=info;publish('variants')
        comparisons=[];root_coverage=[]
        base_roots={(r['axis'],r['q'],r['N'],r['method']):r for r in baseline['roots']}
        for tag,info in manifest['snapshots'].items():
            if info['returncode'] or not info.get('source_ids_match_expected') or not info['folder']:continue
            analysis=json.loads((OUT/info['folder']/'analysis.json').read_text());fits={fit_key(f):f for f in analysis['size_fits']};roots={(r['axis'],r['q'],r['N'],r['method']):r for r in analysis['roots']}
            for key,r in base_roots.items():
                other=roots.get(key);root_coverage.append(dict(variant=tag,axis=key[0],q=key[1],N=key[2],method=key[3],status=other['status'] if other else 'missing',valid_fraction=other.get('valid_fraction') if other else None,flags=other.get('flags') if other else [],baseline_status=r['status']))
            for f in baseline['size_fits']:
                key=fit_key(f)
                comparisons.append(dict(variant=tag,axis=key[0],q=key[1],method=key[2],minimum_N=key[3],correction=key[4],comparison=compare_fit(f,fits[key]) if key in fits else dict(status='missing_variant_fit',same_sizes=False)))
        run.atomic_json(folder/'comparisons.json',dict(created=run.now(),baseline=manifest['baseline'],root_coverage=root_coverage,fit_comparisons=comparisons,scope=manifest['scope']))
        complete=len(manifest['snapshots'])==9 and all(v['returncode']==0 and v.get('source_ids_match_expected') for v in manifest['snapshots'].values())
        manifest['status']='snapshots_complete_pending_scientific_review' if complete else 'incomplete_snapshots';publish(manifest['status'])
        run.atomic_json(OUT/'latest-d5-sensitivity.json',dict(folder=folder.name,status=manifest['status'],baseline=manifest['baseline']))
        print(run.now(),manifest['status'],folder,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--wait',action='store_true');parser.add_argument('--parallel',type=int,default=3);args=parser.parse_args()
    try:main(args.wait,args.parallel)
    except BlockingIOError:
        print(run.now(),'Another final-sensitivity manager already holds the lock; no duplicate started',flush=True)
    except BaseException as error:
        failure=OUT/('d5-final-sensitivity-failure-'+run.now()[:19].replace(':','').replace('-','')+'-'+str(os.getpid())+'.json')
        run.atomic_json(failure,dict(created=run.now(),pid=os.getpid(),error=repr(error)))
        run.atomic_json(OUT/'d5-final-sensitivity-manager-status.json',dict(updated=run.now(),pid=os.getpid(),state='failed',failure=str(failure)));raise
