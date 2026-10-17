"""Frozen-source 2D eta_eff and xi/L crossings; no Binder criterion or new runs."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl');os.environ['MPLBACKEND']='Agg'
import json,shutil
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import run
from threshold_calibration import monotone_crossing,monotone_values
from scout import json_finite
OUT=Path(__file__).resolve().parent;EXP=OUT.parent
lookup={};provenance=[]
p=EXP/'2026-09-02-xy-bkt-scout';f=p/'results.json'
for r in json.loads(f.read_text()):
    v=r['point'];lookup[v['epsilon_milli'],v['kappa_milli'],v['L']]=dict(e=v['epsilon_milli'],b=v['kappa_milli'],N=v['L'],values=r['summary']['values'],retained=r['summary']['retained'],flags=r['flags'],source=str(f),samples=str(p/'analysis'/f"{v['id']}_bootstrap.npz"))
provenance.append(dict(path=str(f),sha256=run.sha256(f)))
for name in ['2026-09-02-xy-bkt-refine','2026-09-03-xy-bkt-low-kappa']:
    p=EXP/name;f=p/'pooled.json';audit=json.loads((p/'execution-audit.json').read_text());assert audit['passed']
    provenance.append(dict(path=str(f),sha256=run.sha256(f),audit=str(p/'execution-audit.json'),audit_sha256=run.sha256(p/'execution-audit.json')))
    for r in json.loads(f.read_text()):
        lookup[r['epsilon_milli'],r['kappa_milli'],r['L']]=dict(e=r['epsilon_milli'],b=r['kappa_milli'],N=r['L'],values=r['values'],retained=r['retained'],flags=r['flags'],individual_flags=r.get('individual_flags',{}),source=str(f),samples=str(p/'pooled'/f"{r['id']}_bootstrap.npz"))
# New eta/R batches use the same physical model and independent seed streams.
import d2_joint_statistics as joint
allowed={v['id'] for f in (OUT/'batches').glob('d2-eta-correlation-*.json') for v in json.loads(f.read_text())['points']}
newgroups={}
for f in (OUT/'jobs').glob('d2_*.json'):
    j=json.loads(f.read_text());v=j['point']
    if v['id'] in allowed and j['state']=='complete' and j['audit']['passed']:
        newgroups.setdefault((round(v['epsilon_over_pi']*1000),round(v['beta']*1000),v['N']),[]).append(v)
cache=OUT/'d2-joint-cache'/run.sha256(Path(joint.__file__))[:12];cache.mkdir(parents=True,exist_ok=True)
for key,points in sorted(newgroups.items()):
    if len(points)<2:continue
    e,b,N=key;stem=cache/f'e{e:04d}_b{b:04d}_N{N:03d}_cut2001';ids=sorted(v['id'] for v in points)
    meta=json.loads(stem.with_suffix('.json').read_text()) if stem.with_suffix('.json').exists() else None
    if not meta or meta['source_ids']!=ids:
        summary,draws=joint.summarize(points)
        if draws is None:continue
        np.savez_compressed(stem.with_suffix('.npz'),samples=draws)
        meta=dict(source_ids=ids,summary=summary,helper_sha256=run.sha256(Path(joint.__file__)),audits=[json.loads((OUT/'jobs'/f"{v['id']}.json").read_text())['audit'] for v in points])
        run.atomic_json(stem.with_suffix('.json'),json_finite(meta))
    v=dict(e=e,b=b,N=N,values=meta['summary']['values'],retained=meta['summary']['retained'],flags=meta['summary']['flags'],individual_flags={},source=str(stem.with_suffix('.json')),samples=str(stem.with_suffix('.npz')))
    if key in lookup:
        # Combine sufficient moments with retained-row weights, not derived ratios.
        old=lookup[key];weights=[old['retained'],v['retained']]
        def sufficient(q):
            q=np.asarray(q);m2,m4,R=q[...,2],q[...,3],q[...,1];fm=N*N*m2/(1+(2*N*np.sin(np.pi/N)*R)**2)
            return np.stack([m2,m4,N*N*m2,fm,q[...,4],q[...,5],q[...,6]],axis=-1)
        names=['binder','xi_over_L','m2','m4','link_cosine','vortex_density','winding_squared']
        means=np.average(np.array([sufficient([old['values'][k] for k in names]),sufficient([v['values'][k] for k in names])]),weights=weights,axis=0)
        samples=np.average(np.array([sufficient(np.load(old['samples'])['samples']),sufficient(np.load(v['samples'])['samples'])]),weights=weights,axis=0)
        merged=stem.with_name(stem.name+'_with_archived_pool');np.savez_compressed(merged.with_suffix('.npz'),samples=joint.quantities(samples,N))
        v.update(values=dict(zip(names,joint.quantities(means,N).tolist())),retained=sum(weights),flags=sorted(set(old['flags']+v['flags'])),individual_flags=old.get('individual_flags',{}),source=old['source']+';'+v['source'],samples=str(merged.with_suffix('.npz')))
    lookup[key]=v
    provenance.append(dict(path=str(stem.with_suffix('.json')),sha256=run.sha256(stem.with_suffix('.json')),sample_sha256=run.sha256(stem.with_suffix('.npz')),source_ids=ids))

lookup={k:v for k,v in lookup.items() if k[0]>=560 and k[0]%20==0}
folder=OUT/('d2-eta-correlation-existing-'+run.now()[:19].replace('-','').replace(':','').replace('T','-'));folder.mkdir()
shutil.copyfile(__file__,folder/Path(__file__).name);shutil.copyfile(OUT/'threshold_calibration.py',folder/'threshold_calibration.py');shutil.copyfile(Path(joint.__file__),folder/'d2_joint_statistics.py')
boots={k:np.load(v['samples'])['samples'] for k,v in lookup.items()};assert all(b.shape==(2000,7) for b in boots.values())
observations=[];locations=[];undefined_ratio_groups=[]
for e in sorted({k[0] for k in lookup}):
    curves=[]
    for N in sorted({k[2] for k in lookup}):
        allkeys=sorted(k for k in lookup if k[0]==e and k[2]==N)
        keys=[k for k in allkeys if lookup[k]['values']['xi_over_L'] is not None and np.isfinite(lookup[k]['values']['xi_over_L'])]
        undefined_ratio_groups.extend(lookup[k] for k in allkeys if k not in keys)
        rows=[lookup[k] for k in keys]
        if rows:
            curves.append(('xi_over_L',N,.7506912,np.array([k[1]/1000 for k in keys]),np.array([v['values']['xi_over_L'] for v in rows]),np.column_stack([boots[k][:,1] for k in keys]),rows))
    for small in sorted({k[2] for k in lookup if 2*k[2] in {v[2] for v in lookup}}):
        keys=sorted(k for k in lookup if k[0]==e and k[2]==small and (e,k[1],2*small) in lookup)
        if keys:
            eta=np.array([-np.log2(lookup[e,k[1],2*small]['values']['m2']/lookup[k]['values']['m2']) for k in keys])
            B=np.column_stack([-np.log2(boots[e,k[1],2*small][:,2]/boots[k][:,2]) for k in keys])
            rows=[dict(small=lookup[k],large=lookup[e,k[1],2*small]) for k in keys]
            curves.append(('eta_eff',small,.25,np.array([k[1]/1000 for k in keys]),eta,B,rows))
    for observable,N,target,beta,y,B,source_rows in curves:
        sigma=np.nanstd(B,axis=0,ddof=1);weights=1/np.maximum(sigma,1e-12)**2;sign=-1 if observable=='eta_eff' else 1
        pooled=monotone_values(sign*y,weights)
        hits=np.flatnonzero((pooled[:-1]<=sign*target)&(pooled[1:]>sign*target))
        # Keep the local bracket and one neighboring knot on either side.
        # Undefined xi ratios at remote massive-phase points must not censor
        # otherwise valid critical-neighborhood bootstrap draws.
        if len(hits)==1:
            i=int(hits[0]);finite_fraction=np.isfinite(B).mean(0);positions=np.arange(len(beta))
            low=np.flatnonzero((positions<=i)&(pooled<=sign*target-3*sigma)&(finite_fraction>=.99))
            high=np.flatnonzero((positions>=i+1)&(pooled>=sign*target+3*sigma)&(finite_fraction>=.99))
            first=int(low[-1]) if len(low) else max(0,i-1)
            last=int(high[0]) if len(high) else min(len(beta)-1,i+2)
            ix=np.arange(first,last+1)
            window_policy='Closest finite-bootstrap endpoints at least3measuredSD fromtarget, or adjacentknots if unavailable; data-selected sensitivity choice'
        else:
            ix=np.arange(len(beta));window_policy='No central bracket; all available knots retained as unresolved diagnostics'
        central=monotone_crossing(beta[ix],sign*y[ix],sign*target,weights=weights[ix])
        samples=np.array([monotone_crossing(beta[ix],sign*b[ix],sign*target,weights=weights[ix]) if np.isfinite(b[ix]).all() else np.nan for b in B])
        valid=np.isfinite(samples);row=dict(epsilon_over_pi=e/1000,observable=observable,N=N,pair=[N,2*N] if observable=='eta_eff' else None,target=target,window_policy=window_policy,bootstrap_window_beta=beta[ix].tolist(),beta_location=float(central) if np.isfinite(central) else None,valid_fraction=float(valid.mean()),status='unbracketed' if not np.isfinite(central) else 'conditional_crossing' if valid.mean()>=.99 else 'bootstrap_bracket_unresolved')
        stem=f'e{e:04d}_{observable}_N{N:03d}';np.savez_compressed(folder/f'{stem}.npz',samples=samples)
        if valid.sum()>1:row['MC95']=np.quantile(samples[valid],[.025,.975]).tolist()
        relevant=[source_rows[i] for i in ix]
        flat=[part for r in relevant for part in ([r['small'],r['large']] if observable=='eta_eff' else [r])]
        flags=sorted({flag for r in flat for flag in r.get('flags',[])})
        individual_flags={key:values for r in flat for key,values in r.get('individual_flags',{}).items()}
        full=monotone_crossing(beta,sign*y,sign*target,weights=weights)
        row.update(flags=flags,individual_flags=individual_flags,full_curve_central_location=float(full) if np.isfinite(full) else None,local_vs_full_shift=float(abs(central-full)) if np.isfinite(central) and np.isfinite(full) else None,pooling_max_shift_sd=float(np.max(np.abs(pooled-sign*y)/np.maximum(sigma,1e-12))))
        locations.append(row);observations.append(dict(epsilon_over_pi=e/1000,observable=observable,N=N,beta=beta.tolist(),values=y.tolist(),MC_sd=sigma.tolist(),sources=source_rows))
for obs in ['eta_eff','xi_over_L']:
    eps=sorted({r['epsilon_over_pi'] for r in observations if r['observable']==obs});fig,axes=plt.subplots((len(eps)+3)//4,4,figsize=(14,2.6*((len(eps)+3)//4)),squeeze=False)
    for ax,e in zip(axes.flat,eps):
        for r in observations:
            if r['observable']!=obs or r['epsilon_over_pi']!=e:continue
            label=f"{r['N']}:{2*r['N']}" if obs=='eta_eff' else str(r['N'])
            ax.errorbar(r['beta'],r['values'],yerr=r['MC_sd'],fmt='o-',ms=2,lw=.7,capsize=1,label=label)
        ax.axhline(.25 if obs=='eta_eff' else .7506912,color='black',lw=.8,ls=':');ax.set_title(f'epsilon/pi={e:.2f}',fontsize=9);ax.set_xlabel('beta',fontsize=8);ax.set_ylabel(obs,fontsize=8);ax.legend(fontsize=6);ax.grid(alpha=.2);ax.tick_params(labelsize=7)
    for ax in list(axes.flat)[len(eps):]:ax.set_visible(False)
    fig.suptitle(f'Existing 2D data: {obs}\nBars: cached conditional bootstrap SD. Dotted line: infinite-size critical target; finite-size logarithmic shifts remain.',fontsize=11);fig.tight_layout(rect=(0,0,1,.95))
    for ext in ['png','pdf','svg']:fig.savefig(folder/f'{obs}-existing-data.{ext}',dpi=150)
    plt.close(fig)
run.atomic_json(folder/'analysis.json',json_finite(dict(created=run.now(),locations=locations,observations=observations,undefined_ratio_groups=undefined_ratio_groups,source_provenance=provenance,new_eta_R_coordinate_groups=len(newgroups),selected_source_coordinate_counts={name:sum(name in v['source'] for v in lookup.values()) for name in ['2026-09-02-xy-bkt-scout','2026-09-02-xy-bkt-refine','2026-09-03-xy-bkt-low-kappa']},source_sha256=run.sha256(Path(__file__)),scope='Reuse audited September scout/refinement archives. Replacements prefer later longer pooled starts at identical coordinates. Cached2000draws preserve xi/m2 covariance across derived eta/R quantities. Targets eta=.25 and xi/L=.7506912 are infinite-size limits; finite-box crossings are not yet extrapolated BKT estimates. Weighted monotone interpolation uses the nearest finite-bootstrap anchors at least3MCSD away fromtarget, falling back to adjacentknots if unavailable. This expands too-tight existing neighborhoods to preserve root tails, without requiring remote undefined massive-phase ratios. Data-selected window margins and monotonicity are declared sensitivity choices, not calibrated coverage. No Binder used as phase criterion, no old files modified.',critical_reference='https://arxiv.org/pdf/1307.0485 equations1.4,3.5')))
run.atomic_json(OUT/'latest-d2-eta-correlation.json',dict(created=run.now(),folder=folder.name,status='finite_size_existing_data_reconnaissance'))
print('snapshot',folder,flush=True)
for e in sorted({r['epsilon_over_pi'] for r in locations}):print(e,[(r['observable'],r['N'],r['status'],r['beta_location'],r['valid_fraction']) for r in locations if r['epsilon_over_pi']==e],flush=True)
