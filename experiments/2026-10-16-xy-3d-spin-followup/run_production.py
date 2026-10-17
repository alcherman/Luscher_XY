"""Finite, restartable production stages with fresh physical admission each stage."""
import os,sys,fcntl,subprocess,traceback
import campaign as c
import design
P=c.OUT
STAGES=[('production-inner',[.80068,.80084]),('production-outer',[.80060,.80092]),('production-center',[.80076])]

def points(qs):
    result=[design.point(n,q,r,100000,w,role='production') for n in [64,48,32,24] for q in qs for r,w in enumerate([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])]
    for point in result:point['warmup']=10000
    return result

def validate_batch(batch,ps):
    if batch['points']!=ps or batch['source_hashes']!=c.sources():raise RuntimeError('Existing stage differs from frozen production specification')

def remaining_commitments(exclude):
    """Unfinished frozen jobs outside this stage may not spend our reserve."""
    other={}
    for path in (P/'batches').glob('*.json'):
        for p in c.read(path)['points']:
            if p['id'] in exclude:continue
            jp=P/'jobs'/f"{p['id']}.json";j=c.read(jp) if jp.exists() else {}
            if j.get('state') not in ['complete','failed']:
                other[p['id']]=p
    return list(other.values())

def freeze_stage(name,ps):
    manifest=P/'batches'/f'{name}.json'
    others=remaining_commitments({p['id'] for p in ps})
    if others:raise RuntimeError('Other frozen unfinished commitments require explicit reconciliation before new stage: '+str([p['id'] for p in others]))
    if manifest.exists():validate_batch(c.read(manifest),ps);return manifest
    return design.freeze(name,ps,'Reviewed production: full-cadence G, common xi/L, four independent starts; 100000 cycles, 10000 warmup label, no claim of equilibrium or final Delta precision.')

def completed_stage(name,ps):
    manifest=P/'batches'/f'{name}.json';report=P/f'{name}-sampling.json'
    if not manifest.exists() or not report.exists():return False
    validate_batch(c.read(manifest),ps);saved=c.read(report)
    if saved.get('stop_flags') or {r['id'] for r in saved['rows']}!={p['id'] for p in ps} or sampling_stops(saved['rows']):return False
    for p in ps:
        jobpath=P/'jobs'/f"{p['id']}.json";ledgerpath=P/'retention-records'/f"{p['id']}.json";data=P/'chains'/f"{p['id']}.h5"
        if not jobpath.exists() or not ledgerpath.exists() or not data.exists():return False
        j=c.read(jobpath);r=c.read(ledgerpath)
        if j['state']!='complete' or not j.get('audit',{}).get('passed') or not j['audit'].get('online_spin',{}).get('passed'):return False
        if not r.get('committed') or not r.get('verified') or r['original_sha256']!=j['audit']['sha256'] or data.stat().st_size!=r['compact_bytes']:return False
    return True

def sampling_stops(rows):
    """Conservative descriptive screens, not calibrated hypothesis tests."""
    import numpy as np
    flags=[]
    for row in rows:
        stats=row['statistics'];values=np.array(stats['G_mean']+[stats['m2_mean'],stats['R_mean']],float)
        if not np.all(np.isfinite(values)) or np.any(values[:3]<=0):flags.append(dict(id=row['id'],reason='invalid_or_nonpositive_mean'))
        if row['available_blocks']<50:flags.append(dict(id=row['id'],reason='fewer_than_50_conservative_blocks'))
        se=np.array(stats.get('block_SE_Gq_Gh_m2',[float('nan')]*3));halves=np.array(stats['first_last_half_Gq_Gh_m2'])
        if np.any(~np.isfinite(se)) or np.any(se<=0):flags.append(dict(id=row['id'],reason='invalid_block_SE'))
        elif np.any(abs(halves[0]-halves[1])>10*se):flags.append(dict(id=row['id'],reason='large_half_shift_descriptive_5SE_screen'))
    groups={}
    for row in rows:groups.setdefault((row['N'],row['point']['epsilon_over_pi']),[]).append(row)
    for group in groups.values():
        for i,a in enumerate(group):
            for b in group[i+1:]:
                x=a['statistics'];y=b['statistics'];se=np.sqrt(np.array(x.get('block_SE_Gq_Gh_m2',[float('nan')]*3))**2+np.array(y.get('block_SE_Gq_Gh_m2',[float('nan')]*3))**2)
                difference=abs(np.array(x['G_mean']+[x['m2_mean']])-np.array(y['G_mean']+[y['m2_mean']]))
                if np.any(difference>5*se):flags.append(dict(ids=[a['id'],b['id']],reason='large_independent_start_shift_descriptive_5SE_screen'))
    return flags

def main():
    with (P/'production-manager.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        awake=subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(os.getpid())])
        state=dict(pid=os.getpid(),started=c.frozen.now(),completed_stages=[])
        def publish(value,**extra):
            state.update(state=value,updated=c.frozen.now(),**extra);c.atomic(P/'production-status.json',state)
        try:
            if not (P/'production-review.json').exists() or not c.read(P/'production-review.json').get('approved'):raise RuntimeError('Independent concrete design review required')
            review=c.read(P/'production-review.json')
            for filename,key in [('run_production.py','manager_sha256'),('design.py','design_sha256'),('analyze_pilot.py','analysis_sha256')]:
                if review.get(key)!=c.frozen.sha256(P/filename):raise RuntimeError('Production support changed since independent review: '+filename)
            import shutil
            snapshot=P/'source-snapshots/reviewed-production-support';snapshot.mkdir(exist_ok=True)
            for filename in ['run_production.py','design.py','analyze_pilot.py']:
                dest=snapshot/filename
                if dest.exists() and c.frozen.sha256(dest)!=c.frozen.sha256(P/filename):raise RuntimeError('Existing production support snapshot differs')
                if not dest.exists():shutil.copyfile(P/filename,dest)
            if c.read(P/'pilot-status.json')['state']!='pilot_complete_pending_main_design_review':raise RuntimeError('Audited pilot and analysis prerequisite')
            for name,qs in STAGES:
                ps=points(qs)
                if completed_stage(name,ps):
                    state['completed_stages'].append(name);publish('recovered_verified_completed_stage',stage=name);continue
                manifest=freeze_stage(name,ps)
                config=c.read(P/'control.json');config['launch_new']=True;c.atomic(P/'control.json',config)
                publish('generating_and_auditing',stage=name,trajectories=len(ps),manifest=str(manifest))
                with (P/'logs'/f'{name}-coordinator.log').open('a') as stream:
                    subprocess.run([sys.executable,'-u',str(P/'campaign.py'),'--coordinate','--manifest',str(manifest),'--workers','9'],cwd=c.ROOT,stdout=stream,stderr=subprocess.STDOUT,check=True)
                config=c.read(P/'control.json');config['launch_new']=False;c.atomic(P/'control.json',config)
                jobs=[c.read(P/'jobs'/f"{p['id']}.json") for p in ps]
                if any(j['state']!='complete' or not j['audit']['passed'] or not j['audit'].get('online_spin',{}).get('passed') for j in jobs):raise RuntimeError('Full closure audit prerequisite')
                publish('retaining_audited_fields',stage=name)
                for p in ps:
                    with (P/'logs'/f"retention-{p['id']}.log").open('a') as stream:
                        subprocess.run([sys.executable,'-u',str(P/'campaign.py'),'--retain',p['id']],cwd=c.ROOT,stdout=stream,stderr=subprocess.STDOUT,check=True)
                    ledger=c.read(P/'retention-records'/f"{p['id']}.json")
                    if not ledger.get('committed'):raise RuntimeError('Retention journal incomplete')
                publish('analysing_stage_sampling',stage=name)
                from analyze_pilot import diagnose
                rows=[diagnose(j) for j in jobs]
                stops=sampling_stops(rows)
                c.atomic(P/f'{name}-sampling.json',dict(created=c.frozen.now(),rows=rows,stop_flags=stops,scope='Per-start diagnostics and conservative advance screens only. No formal multiple-testing significance, final exponent or calibrated systematic-error inference.'))
                if stops:raise RuntimeError('Sampling diagnostics require scientific review before next stage: '+str(stops))
                state['completed_stages'].append(name);publish('stage_complete',stage=name)
            publish('production_complete_pending_correlated_exponent_analysis')
        except BaseException:
            config=c.read(P/'control.json');config['launch_new']=False;c.atomic(P/'control.json',config)
            publish('failed_needs_attention',error=traceback.format_exc());raise
        finally:
            if awake.poll() is None:awake.terminate()
if __name__=='__main__':main()
