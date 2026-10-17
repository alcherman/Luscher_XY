"""Snapshot finite-box crossings and conditional d3/d4 threshold extrapolations."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
os.environ['MPLBACKEND']='Agg'
import csv
import json
from pathlib import Path
import shutil
import numpy as np
import matplotlib.pyplot as plt
import run
import analyze_calibration as diagnostic
from threshold_calibration import crossing,size_fit,monotone_crossing,monotone_values
from scout import json_finite

OUT=Path(__file__).resolve().parent


def main():
    stamp=run.now()[:19].replace('-','').replace(':','').replace('T','-')
    folder=OUT/('threshold-curve-snapshot-'+stamp);folder.mkdir()
    helper_sha=run.sha256(Path(diagnostic.__file__))
    cache=OUT/'threshold-group-cache'/helper_sha[:12];cache.mkdir(parents=True,exist_ok=True)
    for source in [Path(__file__),Path(diagnostic.__file__),OUT/'threshold_calibration.py']:
        shutil.copyfile(source,folder/source.name)
    grouped={};jobs=[]
    for path in sorted((OUT/'jobs').glob('*.json')):
        j=json.loads(path.read_text());p=j['point']
        if j['state']=='complete' and j['audit']['passed'] and p['D'] in [3,4]:
            grouped.setdefault((p['D'],round(p['epsilon_over_pi']*100),p['N'],round(p['beta']*1000)),[]).append(p)
            jobs.append(dict(id=p['id'],audit=j['audit']))
    groups=[];samples={}
    for (D,q,N,beta),points in sorted(grouped.items()):
        if len(points)<2:continue
        points=sorted(points,key=lambda p:p['replica'])
        stem=cache/f'd{D}_e{q:03d}_N{N:03d}_b{beta:04d}_cut1001'
        key=dict(source_ids=[p['id'] for p in points],helper_sha256=helper_sha)
        old=json.loads(stem.with_suffix('.json').read_text()) if stem.with_suffix('.json').exists() else None
        if old and old['key']==key:
            summary=old['summary'];boot=np.load(stem.with_suffix('.npz'))['samples']
        else:
            summary,boot=diagnostic.summarize([diagnostic.read(p) for p in points],points,1001)
            if boot is None:continue
            # Undefined xi samples do not censor the separately used Binder samples.
            np.savez_compressed(stem.with_suffix('.npz'),samples=boot)
            run.atomic_json(stem.with_suffix('.json'),json_finite(dict(key=key,summary=summary)))
        groups.append(dict(D=D,q=q,N=N,beta=beta/1000,summary=summary,source_ids=key['source_ids']))
        samples[D,q,N,beta]=boot[:,0]
    locations=[];roots={}
    for D,q,N in sorted({(g['D'],g['q'],g['N']) for g in groups}):
        rows=sorted([g for g in groups if (g['D'],g['q'],g['N'])==(D,q,N)],key=lambda g:g['beta'])
        beta=np.array([g['beta'] for g in rows]);u=np.array([g['summary']['values'][0] for g in rows])
        B=np.column_stack([samples[D,q,N,round(b*1000)] for b in beta])
        sigma=B.std(0,ddof=1);weights=1/np.maximum(sigma,1e-12)**2
        pooled=monotone_values(u,weights)
        max_pool_shift=float(np.max(np.abs(pooled-u)/np.maximum(sigma,1e-12)))
        for threshold in [.20,.25,.30]:
            for method in ['linear','pchip']:
                central=monotone_crossing(beta,u,threshold,method,weights)
                row=dict(D=D,q=q,N=N,threshold=threshold,interpolation=method,status='unbracketed',beta_location=None)
                if np.isfinite(central):
                    boot=np.array([monotone_crossing(beta,b,threshold,method,weights) for b in B]);valid=np.isfinite(boot)
                    i=int(np.flatnonzero((pooled[:-1]<=threshold)&(pooled[1:]>threshold))[0])
                    row.update(beta_location=central,valid_fraction=float(valid.mean()),bracket=beta[i:i+2].tolist(),
                               flags=sorted({f for g in rows for f in g['summary']['flags'] if f!='invalid_ratio_bootstrap'}),
                               raw_linear_location=crossing(beta,u,threshold),monotone_pool_max_shift_SD=max_pool_shift,
                               status='conditional_crossing' if valid.mean()>=.99 else 'bootstrap_bracket_unresolved')
                    if valid.sum()>1:row['bootstrap_95']=np.quantile(boot[valid],[.025,.975]).tolist()
                    if max_pool_shift>3:row['flags'].append('nonmonotonic_data_exceeds_3sigma')
                    roots[D,q,N,threshold,method]=boot
                    np.savez_compressed(folder/f'd{D}_e{q:03d}_N{N:03d}_U{round(threshold*100):02d}_{method}_crossing.npz',samples=boot)
                locations.append(row)
    estimates=[];fits=[]
    for D,q in sorted({(g['D'],g['q']) for g in groups}):
        local=[]
        for threshold in [.20,.25,.30]:
            for method in ['linear','pchip']:
                selected=[r for r in locations if (r['D'],r['q'],r['threshold'],r['interpolation'])==(D,q,threshold,method) and r['status']=='conditional_crossing']
                if len(selected)<3:continue
                N=np.array([r['N'] for r in selected]);variants=[(int(N.min()),False)]
                if len(N)>=4:variants += [(int(N[1]),False),(int(N.min()),True)]
                for minimum,correction in variants:
                    B=np.column_stack([roots[D,q,n,threshold,method] for n in N]);valid=np.isfinite(B[:,N>=minimum]).all(1)
                    result=size_fit(D,N,[r['beta_location'] for r in selected],B[valid],minimum,correction)
                    if 'intercept_bootstrap' not in result:continue
                    np.savez_compressed(folder/f'd{D}_e{q:03d}_U{round(threshold*100):02d}_{method}_minN{minimum}_corr{int(correction)}_fit.npz',samples=result.pop('intercept_bootstrap'))
                    result.update(joint_valid_fraction=float(valid.mean()),input_flags=sorted({f for r in selected if r['N']>=minimum for f in r['flags']}))
                    fit=dict(D=D,q=q,threshold=threshold,interpolation=method,minimum_N=minimum,correction=correction,result=result)
                    local.append(fit);fits.append(fit)
        primary=[f for f in local if f['threshold']==.25 and f['interpolation']=='linear' and not f['correction']]
        row=dict(D=D,epsilon_over_pi=q/100,status='insufficient_bracketed_sizes',beta_infinite=None)
        if primary:
            base=primary[0]['result'];tested=[f['result'] for f in local]
            row.update(status='conditional_preview',beta_infinite=base['beta_infinite'],MC_95=base['bootstrap_95'],sizes=base['sizes'],
                       tested_interval_envelope=[min(f['bootstrap_95'][0] for f in tested),max(f['bootstrap_95'][1] for f in tested)],
                       tested_central_range=[min(f['beta_infinite'] for f in tested),max(f['beta_infinite'] for f in tested)],
                       chi2=base['chi2'],df=base['df'],p=base['p'],input_flags=base['input_flags'])
            if base['p'] is not None and base['p']<.01:row['status']='leading_shift_fit_rejected'
            if base['beta_infinite']<=0:row['status']='nonpositive_extrapolation'
        estimates.append(row)
    run.atomic_json(folder/'analysis.json',json_finite(dict(created=run.now(),groups=groups,locations=locations,fits=fits,estimates=estimates,source_jobs=jobs,
        analysis_sha256=run.sha256(Path(__file__)),scope='Conditional d3/d4 preview; primary fixed U=.25. Weighted monotone Binder response is a declared shape assumption, applied to central data and every draw. Pool shifts and unpooled central roots are retained. MC uses1000 stratified contiguous-block draws per coordinate,cut1001. Envelope unions tested interpolation,threshold,size-cut and correction intervals. No joint fit across epsilon; no inference from nominal rejected global roots. Low-block and bracket flags retained; not final validated accuracy. Cut/block checks require additional snapshots before final acceptance.')))
    with (folder/'transition-locations.csv').open('w') as stream:
        w=csv.writer(stream);w.writerow(['d','epsilon_over_pi','status','beta_infinite','MC95_low','MC95_high','tested_envelope_low','tested_envelope_high','sizes','chi2','dof'])
        for r in estimates:w.writerow([r['D'],r['epsilon_over_pi'],r['status'],r['beta_infinite'],*r.get('MC_95',['','']),*r.get('tested_interval_envelope',['','']),' '.join(map(str,r.get('sizes',[]))),r.get('chi2'),r.get('df')])
    def panel(ax,dimension=None):
        for D in ([dimension] if dimension else [3,4]):
            rows=[r for r in estimates if r['D']==D and r['beta_infinite'] is not None and r['beta_infinite']>0]
            if not rows:continue
            x=np.array([r['epsilon_over_pi'] for r in rows]);y=np.array([r['beta_infinite'] for r in rows]);color='tab:blue' if D==3 else 'tab:orange'
            env=np.array([r['tested_interval_envelope'] for r in rows]).T;mc=np.array([r['MC_95'] for r in rows]).T
            ax.errorbar(x,y,yerr=np.array([y-env[0],env[1]-y]),fmt='none',color=color,alpha=.35,capsize=3,lw=1.2)
            ax.errorbar(x,y,yerr=np.array([y-mc[0],mc[1]-y]),fmt='o',color=color,capsize=2,ms=4,label=f'd={D}')
            bad=np.array([r['status']!='conditional_preview' for r in rows])
            if bad.any():ax.scatter(x[bad],y[bad],marker='x',s=65,color='black',zorder=4,label='Leading shift fit rejected')
        ax.set_xlabel(r'$\epsilon/\pi$');ax.set_ylabel(r'Extrapolated $\beta_c$');ax.grid(alpha=.2);ax.legend(fontsize=9)
        ax.set_title('Conditional Binder-quarter preview'+('' if dimension is None else f', d={dimension}'))
    for dimension in [None,3,4]:
        fig,ax=plt.subplots(figsize=(8,5.5));panel(ax,dimension)
        fig.text(.5,.02,'Dark bars: conditional MC95. Pale bars: tested sensitivity envelope. Crosses: rejected size fits.\nUnbracketed slices omitted; these are preliminary extrapolations, not final accepted transition locations.',ha='center',fontsize=8)
        fig.tight_layout(rect=[0,.10,1,1]);stem='combined' if dimension is None else f'd{dimension}'
        for extension in ['png','pdf','svg']:fig.savefig(folder/f'{stem}-transition-preview.{extension}',dpi=170)
        plt.close(fig)
    run.atomic_json(OUT/'latest-threshold-preview.json',dict(created=run.now(),folder=folder.name))
    print('snapshot',folder,'groups',len(groups),flush=True)
    print(json.dumps(estimates),flush=True)


if __name__=='__main__':main()
