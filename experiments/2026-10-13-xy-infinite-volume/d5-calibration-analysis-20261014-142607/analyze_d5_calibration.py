"""D5 pilot coverage and periodic zero-mode size checks; no precision certification."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json,shutil
from pathlib import Path
import numpy as np
import run
import analyze_calibration as binder
from threshold_calibration import monotone_crossing
from d5_scaling import fit
OUT=Path(__file__).resolve().parent

def finite(v):
    if isinstance(v,dict):return {k:finite(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [finite(x) for x in v]
    if isinstance(v,np.ndarray):return finite(v.tolist())
    if isinstance(v,(float,np.floating)):return float(v) if np.isfinite(v) else None
    if isinstance(v,np.integer):return int(v)
    return v

def main():
    folder=OUT/('d5-calibration-analysis-'+run.now()[:19].replace('-','').replace(':','').replace('T','-'));folder.mkdir()
    for name in ['analyze_d5_calibration.py','d5_scaling.py','analyze_calibration.py','threshold_calibration.py']:shutil.copyfile(OUT/name,folder/name)
    jobs=[json.loads(f.read_text()) for f in (OUT/'jobs').glob('*.json')]
    jobs=[j for j in jobs if j['point']['D']==5 and j['state']=='complete' and j.get('audit',{}).get('passed')]
    sets={}
    for j in jobs:
        p=j['point']
        if p.get('campaign_role')!='zero_beta_endpoint':sets.setdefault(('beta',p['epsilon_over_pi'],p['N'],p['beta']),[]).append(p)
        if p['beta']==0:sets.setdefault(('epsilon',0,p['N'],p['epsilon_over_pi']),[]).append(p)
    groups=[];boots={};cache=OUT/'d5-binder-cache';cache.mkdir(exist_ok=True)
    for key,points in sorted(sets.items()):
        if len(points)<2:continue
        axis,q,N,x=key;stem=cache/f'{axis}_q{q:.12f}_N{N}_x{x:.12f}'
        meta=stem.with_name(stem.name+'.json');npz=stem.with_name(stem.name+'.npz');ids=sorted(p['id'] for p in points)
        old=json.loads(meta.read_text()) if meta.exists() else None
        if old is None or old['ids']!=ids or old.get('helper_sha256')!=run.sha256(Path(binder.__file__)):
            summary,B=binder.summarize([binder.read(p) for p in points],points,1001,draws=2000)
            if B is None:continue
            np.savez_compressed(npz,samples=B);run.atomic_json(meta,finite(dict(ids=ids,summary=summary,helper_sha256=run.sha256(Path(binder.__file__)))))
        else:summary=old['summary'];B=np.load(npz)['samples']
        groups.append(dict(axis=axis,q=q,N=N,x=x,summary=summary,source_ids=ids));boots[key]=B[:,0]
    roots=[];fits=[]
    for axis,q,N in sorted({(g['axis'],g['q'],g['N']) for g in groups}):
        selected=sorted([g for g in groups if (g['axis'],g['q'],g['N'])==(axis,q,N)],key=lambda g:g['x'])
        if len(selected)<2:continue
        x=np.array([g['x'] for g in selected]);y=np.array([g['summary']['values'][0] for g in selected]);B=np.column_stack([boots[axis,q,N,g['x']] for g in selected]);weights=1/np.maximum(B.std(0,ddof=1),1e-12)**2;sign=1 if axis=='beta' else -1
        for method in ['linear','pchip']:
            central=monotone_crossing(x,sign*y,sign*.25,method,weights);draws=np.array([monotone_crossing(x,sign*b,sign*.25,method,weights) for b in B]);valid=np.isfinite(draws)
            row=dict(axis=axis,q=q,N=N,method=method,target=.25,location=float(central) if np.isfinite(central) else None,valid_fraction=float(valid.mean()),status='conditional_crossing' if np.isfinite(central) and valid.mean()>=.99 else 'unresolved',MC95=np.quantile(draws[valid],[.025,.975]).tolist() if valid.sum()>1 else None,measured_window=x.tolist(),flags=sorted({flag for g in selected for flag in g['summary']['flags'] if flag!='invalid_ratio_bootstrap'}))
            roots.append(row);np.savez_compressed(folder/f'{axis}_q{q:.12f}_N{N}_{method}.npz',samples=draws)
    for axis,q in sorted({(r['axis'],r['q']) for r in roots}):
        for method in ['linear','pchip']:
            selected=sorted([r for r in roots if (r['axis'],r['q'],r['method'])==(axis,q,method) and r['status']=='conditional_crossing'],key=lambda r:r['N'])
            if not selected:continue
            B=np.column_stack([np.load(folder/f'{axis}_q{q:.12f}_N{r["N"]}_{method}.npz')['samples'] for r in selected])
            for minimum,correction in [(4,False),(6,False),(4,True)]:
                f=fit([r['N'] for r in selected],[r['location'] for r in selected],B,minimum,correction)
                if 'samples' in f:np.savez_compressed(folder/f'fit{len(fits):04d}.npz',samples=f.pop('samples'))
                fits.append(dict(axis=axis,q=q,method=method,result=f))
    run.atomic_json(folder/'analysis.json',dict(created=run.now(),groups=groups,roots=roots,size_fits=fits,source_jobs=[dict(id=j['point']['id'],audit=j['audit']) for j in jobs],scope='D5 calibration cost/coverage only. Coarse interpolation and3sizes are notaccepted finalphaseboundaryprecision. PeriodicabovecriticalleadingN^-2.5, N^-3correction stresscheck; zeroDOFfits diagnostics only. Fullhistories and start/block flags preserved; conditionalbootstrap isnot totalconfidencecoverage.'))
    run.atomic_json(OUT/'latest-d5-calibration.json',dict(folder=folder.name,status='pilot_only'))
    print(folder,flush=True);print('closed_jobs',len(jobs),'roots',[(r['axis'],r['q'],r['N'],r['status']) for r in roots if r['method']=='linear'],flush=True)
if __name__=='__main__':main()
