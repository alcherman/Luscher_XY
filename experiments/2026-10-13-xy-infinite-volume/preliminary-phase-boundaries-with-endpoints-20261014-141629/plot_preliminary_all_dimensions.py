"""Preliminary infinite-volume curves, with horizontal beta-zero endpoint bars."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl');os.environ['MPLBACKEND']='Agg'
import json,csv,shutil
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import run
from extrapolate_d2_eta_correlation import fit_roots
OUT=Path(__file__).resolve().parent

def main():
    paths={k:OUT/json.loads((OUT/f).read_text())['folder'] for k,f in [('d2','latest-d2-eta-correlation.json'),('endpoint','latest-zero-beta-endpoints.json'),('d34','latest-d3-d4-delivery.json')]}
    a=json.loads((paths['d2']/'analysis.json').read_text());e=json.loads((paths['endpoint']/'analysis.json').read_text());d34=json.loads((paths['d34']/'phase-boundary-analysis.json').read_text())
    folder=OUT/('preliminary-phase-boundaries-with-endpoints-'+run.now()[:19].replace('-','').replace(':','').replace('T','-'));folder.mkdir();shutil.copyfile(__file__,folder/Path(__file__).name)
    rows=[];fits=[];unresolved=[]
    for q in sorted({r['epsilon_over_pi'] for r in a['locations']}):
        roots=[r for r in a['locations'] if r['epsilon_over_pi']==q and r['status']=='conditional_crossing']
        B=np.column_stack([np.load(paths['d2']/f"e{round(q*1000):04d}_{r['observable']}_N{r['N']:03d}.npz")['samples'] for r in roots])
        variants=[];base=None
        for minimum in [32,64]:
            for offset in [0,1,2,4]:
                for mode in ['joint','xi_over_L','eta_eff']:
                    f=fit_roots(roots,B,minimum,offset,offset,mode)
                    if 'samples' not in f:continue
                    np.savez_compressed(folder/f'd2_q{round(q*1000)}_min{minimum}_C{offset}_{mode}.npz',samples=f.pop('samples'))
                    f.update(epsilon_over_pi=q,minimum_N=minimum,log_offset=offset,mode=mode);fits.append(f);variants.append(f)
                    if minimum==32 and offset==0 and mode=='joint':base=f
        if base is None:
            unresolved.append(dict(d=2,epsilon_over_pi=q,reason='No identifiable shared-limit fit with >=99% common valid root draws'));continue
        rows.append(dict(d=2,epsilon_over_pi=q,beta=base['beta_infinite'],MC_low=base['MC95'][0],MC_high=base['MC95'][1],sensitivity_low=min(f['MC95'][0] for f in variants),sensitivity_high=max(f['MC95'][1] for f in variants),error_axis='vertical',status='exploratory_BKT_extrapolation',p=base['p']))
    for r in d34['rows']:
        if r['beta'] is None:unresolved.append(r);continue
        rows.append(dict(d=r['d'],epsilon_over_pi=r['epsilon_over_pi'],beta=r['beta'],MC_low=r['MC95_low'],MC_high=r['MC95_high'],sensitivity_low=r['sensitivity_low'],sensitivity_high=r['sensitivity_high'],error_axis='vertical',status=r['status'],p=r['p']))
    for D in [2,3,4]:
        candidates=[f for f in e['size_fits'] if f['D']==D and 'epsilon_over_pi_infinite' in f]
        if D==2:
            variants=[f for f in candidates if f['minimum_N'] in [32,64] and f['log_offset'] in [0,1,2,4]]
            base=next((f for f in variants if f['minimum_N']==32 and f['log_offset']==0 and f['mode']=='joint'),None);interval='MC95'
        else:
            variants=candidates;base=next((f for f in candidates if f['minimum_N']==(12 if D==3 else 10) and not f['correction']),None);interval='bootstrap_95'
        if base is None:unresolved.append(dict(d=D,beta=0,reason='No primary endpoint extrapolation'));continue
        rows.append(dict(d=D,epsilon_over_pi=base['epsilon_over_pi_infinite'],beta=0,MC_low=base[interval][0],MC_high=base[interval][1],sensitivity_low=min(f[interval][0] for f in variants),sensitivity_high=max(f[interval][1] for f in variants),error_axis='horizontal',status='preliminary_zero_beta_endpoint',p=base.get('p')))
    colors={2:'#2563a7',3:'#218463',4:'#b25629'}
    def panel(ax,D):
        vals=sorted([r for r in rows if r['d']==D and r['error_axis']=='vertical'],key=lambda r:r['epsilon_over_pi']);color=colors[D]
        if vals:
            x=np.array([r['epsilon_over_pi'] for r in vals]);y=np.array([r['beta'] for r in vals])
            light=np.array([[max(0,r['beta']-r['sensitivity_low']),max(0,r['sensitivity_high']-r['beta'])] for r in vals]).T
            dark=np.array([[max(0,r['beta']-r['MC_low']),max(0,r['MC_high']-r['beta'])] for r in vals]).T
            ax.errorbar(x,y,yerr=light,fmt='none',ecolor=color,alpha=.25,elinewidth=4,capsize=0,zorder=1)
            ax.errorbar(x,y,yerr=dark,fmt='o',color=color,markersize=4,capsize=2,elinewidth=1,zorder=2)
        endpoint=next((r for r in rows if r['d']==D and r['error_axis']=='horizontal'),None)
        if endpoint:
            x=endpoint['epsilon_over_pi']
            for lo,hi,lw,alpha in [(endpoint['sensitivity_low'],endpoint['sensitivity_high'],4,.25),(endpoint['MC_low'],endpoint['MC_high'],1,1)]:
                ax.errorbar([x],[0],xerr=[[max(0,x-lo)],[max(0,hi-x)]],fmt='none',ecolor=color,elinewidth=lw,alpha=alpha,capsize=3,zorder=3)
            ax.plot(x,0,marker='D',ms=5,color=color,zorder=4)
        ax.set_title(f'd = {D}'+('  (exploratory BKT extrapolation)' if D==2 else ''))
        ax.set_xlabel(r'$\epsilon/\pi$');ax.set_ylabel(r'$\beta_c$')
        ax.set_ylim(-.025 if D==2 else -.009,{2:1.24,3:.49,4:.33}[D]);ax.grid(alpha=.15)
        ax.set_xlim({2:.52,3:.775,4:.845}[D],1.015)
    footer='Dark: conditional MC 95% intervals. Light: descriptive model/size sensitivity; not a total confidence interval.\nZero-beta diamonds have horizontal bars only; fixed-epsilon points have vertical bars only.'
    fig,axes=plt.subplots(1,3,figsize=(15,4.9))
    for D,ax in zip([2,3,4],axes):panel(ax,D)
    fig.suptitle('XY phase boundaries — preliminary infinite-volume estimates',fontsize=14)
    fig.text(.5,.015,footer,ha='center',fontsize=9);fig.tight_layout(rect=[0,.105,1,.93])
    for ext in ['png','pdf','svg']:fig.savefig(folder/f'combined-phase-boundaries.{ext}',dpi=180)
    plt.close(fig)
    for D in [2,3,4]:
        fig,ax=plt.subplots(figsize=(7,5.2));panel(ax,D);fig.text(.5,.015,footer,ha='center',fontsize=7.5);fig.tight_layout(rect=[0,.105,1,1])
        for ext in ['png','pdf','svg']:fig.savefig(folder/f'd{D}-phase-boundary.{ext}',dpi=180)
        plt.close(fig)
    with (folder/'preliminary-transition-locations.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    run.atomic_json(folder/'analysis.json',dict(created=run.now(),rows=rows,d2_fits=fits,unresolved=unresolved,sources={k:str(v) for k,v in paths.items()},scope='User-requested preliminary plots. D2 baseline shared-limit GLS at minimumN32, fixed logoffset0; lightbars union minimumN32/64, offsets0/1/2/4 and joint/separate observable fits, including poor exploratory fits. D3/D4 positive branch reused unchanged from reviewed delivery. Endpoint primary D3 minimumN12,D4 minimumN10; size/correction unions remain preliminary. Endpoint local-window and cut/block comparisons are available, but largest D2 R512 and eta256:512 roots remain unresolved and excluded pending wider-anchor runs. Common-root comparisons do not validate omitted roots; size/model sensitivity remains preliminary. D2 epsilon1 baseline misses known critical coupling despite good fit, so its small MC bars must not be interpreted as validated precision. No calibrated total uncertainty coverage.',plot_error_axes='horizontal at beta0; vertical only elsewhere'))
    run.atomic_json(OUT/'latest-preliminary-all-dimensions.json',dict(folder=folder.name,status='preliminary_not_final'))
    print(folder,flush=True);print('points',len(rows),'unresolved',len(unresolved),flush=True)

if __name__=='__main__':main()
