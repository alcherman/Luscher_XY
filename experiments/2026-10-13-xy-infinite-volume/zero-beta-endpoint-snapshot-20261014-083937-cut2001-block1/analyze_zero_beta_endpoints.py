"""Exact beta-zero epsilon scans; conditional local crossings, no final acceptance."""
import os
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/luscher-mpl')
os.environ['MPLBACKEND'] = 'Agg'
import json, shutil
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import run
import analyze_calibration as binder
import d2_joint_statistics as joint
from threshold_calibration import monotone_crossing, size_fit

OUT = Path(__file__).resolve().parent

def epsilon_crossing(x, y, target, increasing, weights=None, method='linear'):
    """Reverse the observable, not epsilon, for decreasing epsilon scans."""
    sign = 1 if increasing else -1
    return monotone_crossing(x, sign*np.asarray(y), sign*target, method, weights)

def finite(value):
    if isinstance(value, dict): return {k:finite(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)): return [finite(v) for v in value]
    if isinstance(value, np.ndarray): return finite(value.tolist())
    if isinstance(value, (float,np.floating)): return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer): return int(value)
    return value

def main(cut=2001, factor=1):
    stamp=run.now()[:19].replace('-','').replace(':','').replace('T','-')
    folder=OUT/f'zero-beta-endpoint-snapshot-{stamp}-cut{cut}-block{factor}'
    folder.mkdir()
    for source in [Path(__file__),Path(joint.__file__),Path(binder.__file__),OUT/'threshold_calibration.py']:
        shutil.copyfile(source,folder/source.name)
    grouped={}; sources=[]
    for path in sorted((OUT/'jobs').glob('*.json')):
        job=json.loads(path.read_text()); p=job['point']
        if job['state']!='complete' or not job['audit']['passed'] or p['beta']!=0: continue
        # Preserve exact epsilon coordinates; never round to the positive-beta grid.
        if p['D']==2 and p.get('campaign_role')!='zero_beta_endpoint': continue
        grouped.setdefault((p['D'],p['N'],p['epsilon_over_pi']),[]).append(p)
        sources.append(dict(id=p['id'],audit=job['audit']))
    groups=[]; samples={}
    cache=OUT/'zero-beta-joint-cache';cache.mkdir(exist_ok=True)
    for (D,N,q), points in sorted(grouped.items()):
        if len(points)<2: continue
        key=dict(ids=sorted(p['id'] for p in points),cut=cut,factor=factor,
                 source_sha256=run.sha256(Path(joint.__file__) if D==2 else Path(binder.__file__)))
        stem=cache/f'd{D}_N{N}_e{q:.12f}_cut{cut}_block{factor}'
        old=json.loads(stem.with_suffix('.json').read_text()) if stem.with_suffix('.json').exists() else None
        if old and old['key']==key:
            summary=old['summary']; boot=np.load(stem.with_suffix('.npz'))['samples']
        else:
            if D==2: summary,boot=joint.summarize(points,cut=cut,factor=factor)
            else: summary,boot=binder.summarize([binder.read(p) for p in points],points,cut,draws=2000,factor=factor)
            if boot is None: continue
            np.savez_compressed(stem.with_suffix('.npz'),samples=boot)
            run.atomic_json(stem.with_suffix('.json'),finite(dict(key=key,summary=summary)))
        groups.append(dict(D=D,N=N,epsilon_over_pi=q,summary=summary,source_ids=key['ids']))
        samples[D,N,q]=boot
    observations=[]
    for g in groups:
        D,N,q=g['D'],g['N'],g['epsilon_over_pi']; b=samples[D,N,q]
        if D==2:
            observations.append(dict(D=D,N=N,q=q,observable='xi_over_L',value=g['summary']['values']['xi_over_L'],draws=b[:,1],flags=g['summary']['flags']))
            if (D,2*N,q) in samples:
                big=next(v for v in groups if (v['D'],v['N'],v['epsilon_over_pi'])==(D,2*N,q))
                observations.append(dict(D=D,N=N,q=q,observable='eta_eff',
                    value=-np.log2(big['summary']['values']['m2']/g['summary']['values']['m2']),
                    draws=-np.log2(samples[D,2*N,q][:,2]/b[:,2]),
                    flags=sorted(set(g['summary']['flags']+big['summary']['flags']))))
        else:
            observations.append(dict(D=D,N=N,q=q,observable='binder',value=g['summary']['values'][0],draws=b[:,0],flags=g['summary']['flags']))
    roots=[]; fits=[]
    for D,observable in sorted({(v['D'],v['observable']) for v in observations}):
        target={'binder':.25,'eta_eff':.25,'xi_over_L':.7506912}[observable]
        fig,ax=plt.subplots(figsize=(7,4.8)); selected=[]; root_draws=[]
        for N in sorted({v['N'] for v in observations if v['D']==D and v['observable']==observable}):
            rows=sorted([v for v in observations if (v['D'],v['observable'],v['N'])==(D,observable,N)],key=lambda v:v['q'])
            rows=[v for v in rows if v['value'] is not None and np.isfinite(v['value'])]
            if len(rows)<2: continue
            x=np.array([v['q'] for v in rows]);y=np.array([v['value'] for v in rows]); B=np.column_stack([v['draws'] for v in rows])
            se=np.nanstd(B,axis=0,ddof=1); weights=1/np.maximum(se,1e-12)**2
            ax.errorbar(x,y,yerr=se,marker='o',capsize=2,label=f'N={N}')
            for method in ['linear','pchip']:
                central=epsilon_crossing(x,y,target,observable=='eta_eff',weights,method)
                draws=np.array([epsilon_crossing(x,b,target,observable=='eta_eff',weights,method) if np.isfinite(b).all() else np.nan for b in B])
                valid=np.isfinite(draws); fraction=float(valid.mean())
                row=dict(D=D,N=N,observable=observable,target=target,interpolation=method,
                    epsilon_over_pi_location=central,valid_fraction=fraction,
                    status='conditional_crossing' if np.isfinite(central) and fraction>=.99 else 'unresolved',
                    MC_95=np.quantile(draws[valid],[.025,.975]).tolist() if valid.sum()>1 else None,
                    measured_epsilon_over_pi=x.tolist(),flags=sorted({f for r in rows for f in r['flags']}))
                roots.append(row)
                np.savez_compressed(folder/f'd{D}_{observable}_N{N}_{method}_epsilon_roots.npz',samples=draws)
                if method=='linear' and row['status']=='conditional_crossing': selected.append(row);root_draws.append(draws)
        ax.axhline(target,color='black',linestyle='--');ax.set(xlabel=r'$\epsilon/\pi$ at $\beta=0$',ylabel=observable,title=f'd={D}: finite-size endpoint scan (preliminary)');ax.legend();fig.tight_layout()
        for ext in ['png','pdf','svg']:fig.savefig(folder/f'd{D}-{observable}-zero-beta-scan.{ext}',dpi=160)
        plt.close(fig)
        if D in [3,4] and len(selected)>=3:
            B=np.column_stack(root_draws); valid=np.isfinite(B).all(1)
            if valid.mean()>=.99:
                for minimum in ([6,12,16] if D==3 else [10,12,16]):
                    for correction in [False,True]:
                        f=size_fit(D,[v['N'] for v in selected],[v['epsilon_over_pi_location'] for v in selected],B[valid],minimum,correction)
                        if 'intercept_bootstrap' not in f:continue
                        np.savez_compressed(folder/f'd{D}_minN{minimum}_corr{int(correction)}_epsilon_infinite.npz',samples=f.pop('intercept_bootstrap'))
                        f['epsilon_over_pi_infinite']=f.pop('beta_infinite'); f.update(D=D,dependent_variable='epsilon/pi',beta_fixed=0,acceptance='exploratory_not_final')
                        fits.append(f)
    result=finite(dict(created=run.now(),beta_fixed=0,endpoint_error_axis='horizontal epsilon/pi only',status='preliminary_requires_legacy_reuse_and_sensitivity_review',groups=groups,locations=roots,size_fits=fits,source_jobs=sources,cut=cut,block_factor=factor,
        limitations=['New/current chains only; September zero-beta archives still need exact-coordinate integration.','D2 limiting targets require logarithmic size corrections and shared eta/R covariance; no accepted infinite-volume endpoint yet.','Full-scan interpolation is preliminary; test local windows, cuts, blocks, starts and curvature before acceptance.','Conditional bootstrap intervals do not certify total finite-size/model coverage.']))
    run.atomic_json(folder/'analysis.json',result)
    run.atomic_json(OUT/'latest-zero-beta-endpoints.json',dict(folder=folder.name,status=result['status']))
    print(folder,flush=True)
    for row in roots:
        if row['interpolation']=='linear':print(row['D'],row['observable'],row['N'],row['epsilon_over_pi_location'],row['valid_fraction'],row['status'],flush=True)

if __name__=='__main__':main()
