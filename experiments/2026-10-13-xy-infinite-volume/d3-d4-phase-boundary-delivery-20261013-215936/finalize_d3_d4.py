"""Assemble conditional phase-boundary estimates and explicit sensitivity bars."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl');os.environ['MPLBACKEND']='Agg'
import argparse,csv,json,shutil
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import run
from threshold_calibration import size_fit,coordinate
from scout import json_finite
OUT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('snapshot',type=Path);parser.add_argument('--cut-snapshot',type=Path,action='append',default=[]);args=parser.parse_args()
paths=[args.snapshot.resolve()]+[p.resolve() for p in args.cut_snapshot]
data=[json.loads((p/'analysis.json').read_text()) for p in paths]
ids={j['id'] for j in data[0]['source_jobs']}
assert all({j['id'] for j in a['source_jobs']}==ids for a in data),'Sensitivity must use identical physical trajectories'
assert all(j['audit']['passed'] for j in data[0]['source_jobs'])
stamp=run.now()[:19].replace('-','').replace(':','').replace('T','-')
folder=OUT/('d3-d4-phase-boundary-delivery-'+stamp);folder.mkdir()
for source in [Path(__file__),OUT/'threshold_calibration.py',OUT/'analyze_calibration.py']:shutil.copyfile(source,folder/source.name)
allfits=[];roots_cache={}
for path,a in zip(paths,data):
    for D,q in sorted({(g['D'],g['q']) for g in a['groups']}):
        for threshold in a['thresholds']:
            for method in ['linear','pchip']:
                rows=sorted([r for r in a['locations'] if (r['D'],r['q'],r['threshold'],r['interpolation'],r['status'])==(D,q,threshold,method,'conditional_crossing')],key=lambda r:r['N'])
                if len(rows)<3:continue
                ns=np.array([r['N'] for r in rows]);key=(str(path),D,q,threshold,method)
                B=np.column_stack([np.load(path/f'd{D}_e{q:03d}_N{N:03d}_U{round(threshold*100):02d}_{method}_crossing.npz')['samples'] for N in ns])
                for minimum in ([6,8,12] if D==3 else [6,8,10,12]):
                    keep=ns>=minimum
                    if keep.sum()<3:continue
                    good=np.isfinite(B[:,keep]).all(1)
                    if good.mean()<.99:continue
                    for correction in [False,True]:
                        if correction and keep.sum()<4:continue
                        r=size_fit(D,ns,[r['beta_location'] for r in rows],B[good],minimum,correction)
                        samples=r.pop('intercept_bootstrap')
                        stem=f'v{len(allfits):04d}'
                        np.savez_compressed(folder/f'{stem}.npz',samples=samples)
                        r.update(joint_valid_fraction=float(good.mean()),input_flags=sorted({f for row in rows if row['N']>=minimum for f in row.get('flags',[])}))
                        allfits.append(dict(D=D,q=q,threshold=threshold,interpolation=method,minimum_N=minimum,correction=correction,cut=a.get('cut',1001),block_factor=a.get('block_factor',1),source=str(path),samples_file=stem+'.npz',result=r))
rows=[]
for D,q in sorted({(g['D'],g['q']) for g in data[0]['groups']}):
    primary=next((f for f in allfits if (f['D'],f['q'],f['threshold'],f['interpolation'],f['minimum_N'],f['correction'],f['source'])==(D,q,.25,'linear',6 if D==3 else 10,False,str(paths[0]))),None)
    row=dict(d=D,epsilon_over_pi=q/100,status='unresolved_branch_neighbor',beta=None,MC95_low=None,MC95_high=None,sensitivity_low=None,sensitivity_high=None)
    if primary:
        r=primary['result'];variants=[f for f in allfits if f['D']==D and f['q']==q];lo=min(f['result']['bootstrap_95'][0] for f in variants);hi=max(f['result']['bootstrap_95'][1] for f in variants)
        row.update(status='conditional_practical_estimate',beta=r['beta_infinite'],MC95_low=r['bootstrap_95'][0],MC95_high=r['bootstrap_95'][1],sensitivity_low=lo,sensitivity_high=hi,sizes=r['sizes'],chi2=r['chi2'],dof=r['df'],p=r['p'],flags=r['input_flags'],primary_samples=primary['samples_file'],variant_count=len(variants),max_sensitivity_halfwidth=max(r['beta_infinite']-lo,hi-r['beta_infinite']),target_002_within_tested_sensitivity=max(r['beta_infinite']-lo,hi-r['beta_infinite'])<=.02)
    rows.append(row)
metadata=dict(created=run.now(),rows=rows,variants=allfits,sources=[str(p) for p in paths],source_jobs=data[0]['source_jobs'],primary_policy='U=.25,local linear monotone crossing; d3 minimumN6,d4 minimumN10. Uniform d4 minimum chosen after pilot inspection, not preregistered. Fits remain conditional.',errorbar_policy='Dark bars conditional MC95; light bars union of conditional MC95 across tested thresholds .20/.25/.30,linear/PCHIP,size cuts,correction models,and provided cut/block snapshots. Include rejected diagnostic variants; union is descriptive, not a total confidence interval.',scope='D3/D4 only. D2 deferred; d5 disabled. Neighboring zero-beta slices unresolved rather than inferred absent from missing quarter crossing.',source_hashes={str(s):run.sha256(s) for s in [Path(__file__),OUT/'threshold_calibration.py',OUT/'analyze_calibration.py']})
run.atomic_json(folder/'phase-boundary-analysis.json',json_finite(metadata))
columns=['d','epsilon_over_pi','status','beta','MC95_low','MC95_high','sensitivity_low','sensitivity_high','sizes','chi2','dof','p','max_sensitivity_halfwidth','flags']
with (folder/'transition-locations.csv').open('w') as stream:
    writer=csv.DictWriter(stream,fieldnames=columns,extrasaction='ignore');writer.writeheader()
    for r in rows:writer.writerow({**r,'sizes':' '.join(map(str,r.get('sizes',[]))),'flags':';'.join(r.get('flags',[]))})
for dim in [3,4,None]:
    fig,ax=plt.subplots(figsize=(8,5))
    for D,color,marker in [(3,'tab:blue','o'),(4,'tab:orange','s')]:
        if dim is not None and dim!=D:continue
        selected=[r for r in rows if r['d']==D and r['beta'] is not None]
        x=np.array([r['epsilon_over_pi'] for r in selected]);y=np.array([r['beta'] for r in selected])
        outer=np.array([[r['beta']-r['sensitivity_low'],r['sensitivity_high']-r['beta']] for r in selected]).T
        mc=np.array([[r['beta']-r['MC95_low'],r['MC95_high']-r['beta']] for r in selected]).T
        ax.errorbar(x,y,yerr=np.maximum(0,outer),fmt='none',ecolor=color,alpha=.32,elinewidth=5,capsize=0,zorder=1)
        ax.errorbar(x,y,yerr=np.maximum(0,mc),fmt=marker+'-',color=color,ms=4,lw=1,capsize=3,label=f'd={D}',zorder=2)
    ax.set_xlabel(r'$\epsilon/\pi$');ax.set_ylabel(r'$\beta_c(\epsilon)$');ax.grid(alpha=.22);ax.legend()
    ax.set_title('XY phase boundary conditional estimates'+('' if dim is None else f' in d={dim}')+'\nFixed Binder U=1/4 with size extrapolation')
    fig.text(.5,.015,'Dark: conditional MC 95%. Light: tested sensitivity union, without total coverage guarantee.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.04,1,1));stem='combined-transition-locations' if dim is None else f'd{dim}-transition-locations'
    for ext in ['png','pdf','svg']:fig.savefig(folder/f'{stem}.{ext}',dpi=180)
    plt.close(fig)
for D in [3,4]:
    qs=sorted({r['epsilon_over_pi'] for r in rows if r['d']==D});fig,axes=plt.subplots((len(qs)+2)//3,3,figsize=(12,2.8*((len(qs)+2)//3)),squeeze=False)
    for ax,q in zip(axes.flat,qs):
        qcent=round(q*100);local=sorted([r for r in data[0]['locations'] if (r['D'],r['q'],r['threshold'],r['interpolation'])==(D,qcent,.25,'linear') and r['beta_location'] is not None],key=lambda r:r['N'])
        primary=next((r for r in rows if r['d']==D and r['epsilon_over_pi']==q),None)
        for r in local:
            x=float(coordinate(D,r['N']));y=r['beta_location'];lo,hi=r['bootstrap_95'];color='black' if r['N']>=(6 if D==3 else 10) else 'silver'
            ax.errorbar([x],[y],yerr=[[max(0,y-lo)],[max(0,hi-y)]],fmt='o',ms=3,color=color,capsize=2)
            offset=-13 if r['N']==20 and any(v['N']==24 for v in local) else 5
            ax.annotate(str(r['N']),(x,y),xytext=(3,offset),textcoords='offset points',fontsize=7,color=color)
        if primary and primary['beta'] is not None:
            f=next(f for f in allfits if f['samples_file']==primary['primary_samples']);r=f['result'];xx=np.linspace(0,max(float(coordinate(D,v['N'])) for v in local)*1.05,100)
            ax.plot(xx,r['coefficients'][0]+r['coefficients'][1]*xx,color='tab:blue',lw=1)
            ax.errorbar([0],[primary['beta']],yerr=[[max(0,primary['beta']-primary['MC95_low'])],[max(0,primary['MC95_high']-primary['beta'])]],fmt='s',color='tab:blue',ms=3,capsize=2)
            ax.text(.98,.96,f"chi2/dof={r['chi2']:.1f}/{r['df']}\nconditional estimate",transform=ax.transAxes,ha='right',va='top',fontsize=7,bbox=dict(facecolor='white',edgecolor='none',alpha=.9,pad=2))
        else:ax.text(.5,.5,'Unresolved branch neighbor\nNo bracketed quarter crossing',transform=ax.transAxes,ha='center',va='center',fontsize=8)
        ax.set_title(f'epsilon/pi={q:.2f}',fontsize=10);ax.set_xlabel('N^(-1/nu)' if D==3 else 'N^(-2) [ln N]^(-1/10)',fontsize=8);ax.set_ylabel('beta at U=1/4',fontsize=8);ax.grid(alpha=.2);ax.tick_params(labelsize=8)
    for ax in list(axes.flat)[len(qs):]:ax.set_visible(False)
    fig.suptitle(f'd={D} conditional finite-size extrapolations\nBars: conditional MC 95%. Labels: N. Gray boxes are excluded from the primary fit.',fontsize=11);fig.tight_layout(rect=(0,0,1,.95))
    for ext in ['png','pdf','svg']:fig.savefig(folder/f'd{D}-finite-size-extrapolations.{ext}',dpi=160)
    plt.close(fig)
counts={D:sum(j['id'].startswith(f'd{D}_') for j in data[0]['source_jobs']) for D in [3,4]};accepted=[r for r in rows if r['beta'] is not None]
rawbytes=sum((OUT/'chains'/f'{i}.h5').stat().st_size for i in ids)
lines=['# XY phase boundary in three and four dimensions','',f'The fixed Binder quarter estimator gives {sum(r["d"]==3 for r in accepted)} conditional 3D and {sum(r["d"]==4 for r in accepted)} conditional 4D transition locations on the requested epsilon grid. The analysis considers {counts[3]} audited 3D and {counts[4]} audited 4D trajectories, including endpoint and diagnostic runs. The largest tested sensitivity half-width is {max(r["max_sensitivity_halfwidth"] for r in accepted):.5f} in beta. These are attained estimates with explicit sensitivity bars; total error coverage at the requested 0.02 has not been established.','', '## Transition estimates','', '| d | epsilon/pi | beta | Conditional MC 95 percent interval | Tested sensitivity union |','|---|---|---|---|---|']
for r in rows:
    if r['beta'] is None:lines.append(f'| {r["d"]} | {r["epsilon_over_pi"]:.2f} | unresolved neighbor | | |')
    else:lines.append(f'| {r["d"]} | {r["epsilon_over_pi"]:.2f} | {r["beta"]:.5f} | [{r["MC95_low"]:.5f}, {r["MC95_high"]:.5f}] | [{r["sensitivity_low"]:.5f}, {r["sensitivity_high"]:.5f}] |')
lines += ['', '## Estimator and finite size analysis','', 'The convention is U = 1 - <m^4>/(2<m^2>^2). At each size, weighted monotone regression pools only nonmonotone measured Binder values, and local interpolation brackets U=0.25. The same procedure acts on each bootstrap draw. Monotonicity is a declared assumption, not an independent validation. The primary interpolation is linear; PCHIP is a sensitivity check. No rejected broad observable fit supplies a critical coupling.','', 'The leading shift coordinate is N^(-1/0.67169) in 3D and N^(-2)(ln N)^(-1/10) in 4D. The primary size cuts are N>=6 in 3D and N>=10 in 4D. The uniform 4D cut was chosen after inspecting the pilots and small-box residuals; its goodness-of-fit p values are exploratory and do not have a confirmatory selection interpretation. All-size fits, further size cuts, and a subleading correction remain in the sensitivity calculations. This uses modest-box extrapolation and does not prove the asymptotic law.','', 'Independent review of the earlier N24 holdout check found observed-minus-predicted quarter-root differences of -0.000034, -0.000224, and +0.000266 at 4D epsilon/pi=0.88,0.92,0.96, respectively, each with a conditional bootstrap interval including zero. The comparison supports small sensitivity within tested sizes, without certifying the thermodynamic limit.','', '## Monte Carlo uncertainty and tested sensitivity','', 'The default analysis discards cycles through 1000 and uses 1000 independent-coordinate, contiguous-block bootstrap draws, stratified by independent starts and preserving joint observables. Block widths respond to measured autocorrelations. Dark plot bars show the primary conditional MC 95 percent interval. Light bars are the union of intervals across the saved threshold, interpolation, minimum-size, correction, transient-cut and block-length variants. The saved variant grid is not full factorial: transient-cut and block-length checks use U=0.25. The light union includes rejected exploratory fits and is a descriptive stress test, not a calibrated confidence interval. Sensitivity and MC uncertainty are kept as separate numerical columns.','', 'Many groups retain fewer than 50 blocks per start. Initial start disagreement at 3D epsilon/pi=.92,N16,beta=.293 and 4D .98,N10,beta=.261 remains recorded. New independent 16000-cycle starts have descriptively agreeing late halves; window checks do not certify mixing. Removing those below-crossing coordinates leaves the central primary local-linear roots unchanged; a 4D PCHIP root shifts by 0.000190. Bootstrap influence from those deletions was not independently quantified.','', 'Missing positive-beta quarter crossings at 3D .80 and 4D .86 are reported as unresolved branch neighbors. Their zero-beta finite-size behavior is consistent with an ordered side, but does not alone certify absence of a positive thermodynamic transition. No horizontal uncertainty bars or artificial endpoint coordinates are introduced.','', '## Validation and retained evidence','', 'Every contributing trajectory passed independent field/observable recomputation and its closure audit. The frozen numerical sources remain unchanged. Exact checkpoint replay and focused XY tests were validated before production; nine focused threshold/calibration/storage tests passed during final analysis development. Correct generation is separate from demonstrated mixing or extrapolation accuracy.','', f'The contributing raw HDF5 files occupy {rawbytes/2**40:.3f} TiB. Raw histories, fields, RNG states, checkpoints, all earlier snapshots and full unfiltered logs remain retained on the verified Office_SSD route. D2 remains deferred at the user request and d5 remains disabled.','', 'Source snapshots:', *[f'- [{p.name}]({p}/analysis.json)' for p in paths], '', 'Reproducible outputs:', f'- [Numerical table]({folder}/transition-locations.csv)',f'- [Analysis and all tested variants]({folder}/phase-boundary-analysis.json)', f'- [Combined plot PDF]({folder}/combined-transition-locations.pdf)', '']
(folder/'xy-phase-boundary-d3-d4-scientific-report.md').write_text('\n'.join(lines))
run.atomic_json(OUT/'latest-d3-d4-delivery.json',dict(created=run.now(),folder=folder.name,status='conditional_estimates_pending_final_review'))
print('delivery',folder,flush=True);print('counts',counts,'largest_sensitivity_halfwidth',max(r['max_sensitivity_halfwidth'] for r in accepted),flush=True)
