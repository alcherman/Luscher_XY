"""D5 pilot coverage and periodic zero-mode size checks; no precision certification."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json,shutil,argparse
from pathlib import Path
import numpy as np
import run
import analyze_calibration as binder
from threshold_calibration import monotone_crossing,monotone_values
from d5_scaling import fit
OUT=Path(__file__).resolve().parent

def finite(v):
    if isinstance(v,dict):return {k:finite(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [finite(x) for x in v]
    if isinstance(v,np.ndarray):return finite(v.tolist())
    if isinstance(v,(float,np.floating)):return float(v) if np.isfinite(v) else None
    if isinstance(v,np.integer):return int(v)
    return v

def half_diagnostics(chains,points,cut,summary):
    """Descriptive first/last retained-half Binder check at measured block width."""
    rows=[]
    for a,p,d in zip(chains,points,summary['diagnostics']):
        a=a[cut:];width=d['block_width'];values=[];errors=[];counts=[]
        for half in np.array_split(a,2):
            m2,m4=half[:,:2].T;u=1-m4.mean()/(2*m2.mean()**2)
            influence=m4.mean()/m2.mean()**3*(m2-m2.mean())-(m4-m4.mean())/(2*m2.mean()**2)
            blocks=len(half)//width;counts.append(blocks);values.append(float(u))
            errors.append(float(influence[:blocks*width].reshape(blocks,width).mean(1).std(ddof=1)/np.sqrt(blocks)) if blocks>=2 else None)
        se=float(np.hypot(*errors)) if all(e is not None for e in errors) else None
        z=float((values[1]-values[0])/se) if se is not None and se>0 else None
        rows.append(dict(id=p['id'],Binder_halves=values,linearized_block_SE=errors,blocks=counts,z=z,
            flag='Binder_half_shift' if z is not None and abs(z)>4 else ('insufficient_half_blocks' if z is None else None)))
    return rows

def local_indices(y,B,target,sign):
    """Data-selected bracket/3SD anchors; weak anchors retain extra neighbors."""
    se=np.nanstd(B,axis=0,ddof=1);weights=1/np.maximum(se,1e-12)**2
    pooled=monotone_values(sign*np.asarray(y),weights)
    hits=np.flatnonzero((pooled[:-1]<=sign*target)&(pooled[1:]>sign*target))
    if len(hits)!=1:return np.arange(len(y))
    i=int(hits[0]);positions=np.arange(len(y));fraction=np.isfinite(B).mean(0)
    left=np.flatnonzero((positions<=i)&(pooled<=sign*target-3*se)&(fraction>=.99))
    right=np.flatnonzero((positions>=i+1)&(pooled>=sign*target+3*se)&(fraction>=.99))
    first=int(left[-1]) if len(left) else max(0,i-1)
    last=int(right[0]) if len(right) else min(len(y)-1,i+2)
    return np.arange(first,last+1)

def crossing_record(x,y,B,target=.25,sign=1,method='linear',window='local'):
    indices=local_indices(y,B,target,sign) if window=='local' else np.arange(len(x))
    x=np.asarray(x);y=np.asarray(y);B=np.asarray(B)
    weights=1/np.maximum(np.nanstd(B[:,indices],axis=0,ddof=1),1e-12)**2
    central=monotone_crossing(x[indices],sign*y[indices],sign*target,method,weights)
    draws=np.array([monotone_crossing(x[indices],sign*b[indices],sign*target,method,weights) for b in B])
    valid=np.isfinite(draws)
    full_weights=1/np.maximum(np.nanstd(B,axis=0,ddof=1),1e-12)**2
    full=monotone_crossing(x,sign*y,sign*target,method,full_weights)
    row=dict(method=method,target=target,location=float(central) if np.isfinite(central) else None,
        valid_fraction=float(valid.mean()),status='conditional_crossing' if np.isfinite(central) and valid.mean()>=.99 else 'unresolved',
        MC95=np.quantile(draws[valid],[.025,.975]).tolist() if valid.sum()>1 else None,
        measured_window=x.tolist(),selected_window=x[indices].tolist(),selected_indices=indices.tolist(),
        full_window_location=float(full) if np.isfinite(full) else None,
        local_minus_full=float(central-full) if np.isfinite(central) and np.isfinite(full) else None)
    return row,draws

def main(cut=1001,factor=1,target=.25,window='local',starts='all',promote=True,source_analysis=None):
    folder=OUT/('d5-calibration-analysis-'+run.now()[:19].replace('-','').replace(':','').replace('T','-')+f'-cut{cut}-block{factor}-U{round(target*100)}-{window}-{starts}');folder.mkdir()
    for name in ['analyze_d5_calibration.py','d5_scaling.py','analyze_calibration.py','threshold_calibration.py']:shutil.copyfile(OUT/name,folder/name)
    jobs=[json.loads(f.read_text()) for f in (OUT/'jobs').glob('*.json')]
    jobs=[j for j in jobs if j['point']['D']==5 and j['state']=='complete' and j.get('audit',{}).get('passed')]
    if source_analysis:
        frozen=json.loads(Path(source_analysis).read_text());ids={j['id'] for j in frozen['source_jobs']}
        assert ids<={j['point']['id'] for j in jobs},'Frozen source jobs not all completed/audited'
        jobs=[j for j in jobs if j['point']['id'] in ids]
    sets={}
    if starts!='all':jobs=[j for j in jobs if j['point']['replica']%2==int(starts=='odd')]
    for j in jobs:
        p=j['point']
        if p.get('campaign_role')!='zero_beta_endpoint':sets.setdefault(('beta',p['epsilon_over_pi'],p['N'],p['beta']),[]).append(p)
        if p['beta']==0:sets.setdefault(('epsilon',0,p['N'],p['epsilon_over_pi']),[]).append(p)
    groups=[];skipped_groups=[];boots={};cache=OUT/'d5-binder-cache';cache.mkdir(exist_ok=True)
    for key,points in sorted(sets.items()):
        if len(points)<(2 if starts=='all' else 1):
            skipped_groups.append(dict(key=key,source_ids=[p['id'] for p in points],reason='too_few_starts'));continue
        points=sorted(points,key=lambda p:p['id'])
        axis,q,N,x=key;stem=cache/f'{axis}_q{q:.12f}_N{N}_x{x:.12f}_cut{cut}_block{factor}_{starts}'
        meta=stem.with_name(stem.name+'.json');npz=stem.with_name(stem.name+'.npz');ids=sorted(p['id'] for p in points)
        old=json.loads(meta.read_text()) if meta.exists() else None
        if old is None or old['ids']!=ids or old.get('helper_sha256')!=run.sha256(Path(binder.__file__)) or old.get('diagnostic_version')!=1:
            chains=[binder.read(p) for p in points]
            if any(len(a)<=cut+40 for a in chains):
                skipped_groups.append(dict(key=key,source_ids=ids,reason='insufficient_retained_rows'));continue
            summary,B=binder.summarize(chains,points,cut,draws=2000,factor=factor)
            if B is None:
                skipped_groups.append(dict(key=key,source_ids=ids,reason='insufficient_blocks',diagnostics=summary));continue
            summary['half_diagnostics']=half_diagnostics(chains,points,cut,summary)
            summary['flags']=sorted(set(summary['flags'])|{h['flag'] for h in summary['half_diagnostics'] if h['flag']})
            np.savez_compressed(npz,samples=B);run.atomic_json(meta,finite(dict(ids=ids,summary=summary,diagnostic_version=1,helper_sha256=run.sha256(Path(binder.__file__)))))
        else:summary=old['summary'];B=np.load(npz)['samples']
        groups.append(dict(axis=axis,q=q,N=N,x=x,summary=summary,source_ids=ids));boots[key]=B[:,0]
    roots=[];fits=[]
    for axis,q,N in sorted({(g['axis'],g['q'],g['N']) for g in groups}):
        selected=sorted([g for g in groups if (g['axis'],g['q'],g['N'])==(axis,q,N)],key=lambda g:g['x'])
        if len(selected)<2:continue
        x=np.array([g['x'] for g in selected]);y=np.array([g['summary']['values'][0] for g in selected]);B=np.column_stack([boots[axis,q,N,g['x']] for g in selected]);sign=1 if axis=='beta' else -1
        for method in ['linear','pchip']:
            row,draws=crossing_record(x,y,B,target,sign,method,window)
            row.update(axis=axis,q=q,N=N,
                flags=sorted({flag for i in row['selected_indices'] for flag in selected[i]['summary']['flags'] if flag!='invalid_ratio_bootstrap'}),
                all_measured_flags=sorted({flag for g in selected for flag in g['summary']['flags'] if flag!='invalid_ratio_bootstrap'}),
                source_ids=[p for i in row['selected_indices'] for p in selected[i]['source_ids']])
            roots.append(row);np.savez_compressed(folder/f'{axis}_q{q:.12f}_N{N}_{method}.npz',samples=draws)
    for axis,q in sorted({(r['axis'],r['q']) for r in roots}):
        for method in ['linear','pchip']:
            selected=sorted([r for r in roots if (r['axis'],r['q'],r['method'])==(axis,q,method) and r['status']=='conditional_crossing'],key=lambda r:r['N'])
            if not selected:continue
            B=np.column_stack([np.load(folder/f'{axis}_q{q:.12f}_N{r["N"]}_{method}.npz')['samples'] for r in selected])
            for minimum,correction in [(4,False),(6,False),(8,False),(4,True),(6,True)]:
                f=fit([r['N'] for r in selected],[r['location'] for r in selected],B,minimum,correction)
                if 'samples' in f:np.savez_compressed(folder/f'fit{len(fits):04d}.npz',samples=f.pop('samples'))
                f['input_flags']=sorted({flag for r in selected if r['N']>=minimum for flag in r['flags']})
                fits.append(dict(axis=axis,q=q,method=method,minimum_N=minimum,correction=correction,
                    included_root_sizes=[r['N'] for r in selected if r['N']>=minimum],
                    omitted_unresolved_sizes=[r['N'] for r in roots if (r['axis'],r['q'],r['method'])==(axis,q,method) and r['N']>=minimum and r['status']!='conditional_crossing'],result=f))
    run.atomic_json(folder/'analysis.json',finite(dict(created=run.now(),settings=dict(cut=cut,factor=factor,target=target,window=window,starts=starts,source_analysis=str(source_analysis) if source_analysis else None),groups=groups,skipped_groups=skipped_groups,roots=roots,size_fits=fits,source_jobs=[dict(id=j['point']['id'],audit=j['audit']) for j in jobs],scope='D5 calibration cost/coverage only. Locally bracketed interpolation; data-selected3SDwindows and full-windowcomparison are descriptive sensitivity choices, not coverage certification. PeriodicabovecriticalleadingN^-2.5, N^-3correction stresscheck; zeroDOFfits diagnostics only. Fullhistories and start/block flags preserved; conditionalbootstrap isnot totalconfidencecoverage.')))
    if promote:run.atomic_json(OUT/'latest-d5-calibration.json',dict(folder=folder.name,status='preliminary_not_final'))
    print(folder,flush=True);print('closed_jobs',len(jobs),'roots',[(r['axis'],r['q'],r['N'],r['status']) for r in roots if r['method']=='linear'],flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--cut',type=int,default=1001);parser.add_argument('--factor',type=int,default=1)
    parser.add_argument('--target',type=float,default=.25);parser.add_argument('--window',choices=['local','full'],default='local')
    parser.add_argument('--starts',choices=['all','even','odd'],default='all');parser.add_argument('--no-promote',action='store_true');parser.add_argument('--source-analysis')
    args=parser.parse_args();main(args.cut,args.factor,args.target,args.window,args.starts,not args.no_promote,args.source_analysis)

