"""Render existing-field spin correlations and qualified finite-pair quotients."""
import os
os.environ['MPLBACKEND']='Agg';os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-spin-correlation-mpl')
from pathlib import Path
import argparse,json,csv,shutil
import numpy as np
import matplotlib.pyplot as plt
import analyze_spin_correlations as a
import spin_correlation_quotients as sq

def pair_result(groups,draws,L,kind,window,folder):
    q=sorted(set(g['q'] for g in groups if g['L']==L and (window=='full' or g['q']<=.801)))
    s=[next(g for g in groups if g['L']==L and g['q']==x) for x in q];l=[next(g for g in groups if g['L']==2*L and g['q']==x) for x in q]
    central,n=sq.quotient(q,np.array([g['derived'] for g in s]),np.array([g['derived'] for g in l]),kind)
    B=np.full((a.DRAWS,5),np.nan);missing=multiple=invalid=0
    missing_groups=[(n,z) for n in [L,2*L] for z in q if (n,z) not in draws]
    input_flags=sorted(set(f for g in s+l for f in g['flags']))
    for k in range(a.DRAWS):
        if missing_groups:invalid+=1;continue
        x=np.array([draws[L,z][k] for z in q]);y=np.array([draws[2*L,z][k] for z in q])
        if not np.isfinite(x).all() or not np.isfinite(y).all():invalid+=1;continue
        val,n=sq.quotient(q,x,y,kind)
        if val is None:
            if n==0:missing+=1
            else:multiple+=1
        else:B[k]=val
    valid=np.isfinite(B).all(1);contrasts={}
    for label,i,j in [('quarter_minus_half',2,3),('quarter_minus_m2',2,4),('half_minus_m2',3,4)]:
        v=B[valid,i]-B[valid,j];contrasts[label]=dict(central=central[i]-central[j] if central is not None else None,SD=v.std(ddof=1) if len(v)>1 else None,MC95=np.quantile(v,[.025,.975]) if len(v)>1 else None)
    name=f'pair-{L}-{2*L}-{kind}-{window}';np.savez_compressed(folder/(name+'.npz'),samples=B)
    return dict(pair=[L,2*L],interpolation=kind,window=window,q_support=q,missing_draw_groups=missing_groups,input_flags=input_flags,central=central,status='conditional_finite_pair_estimate' if central is not None and valid.mean()>=.99 else 'unresolved_crossing_tails',valid_fraction=valid.mean(),invalid_input_draws=invalid,missing_crossing_draws=missing,multiple_crossing_draws=multiple,paired_contrasts=contrasts,columns=['q_cross','R_cross','Delta_Gquarter','Delta_Ghalf','Delta_m2'],SD=np.std(B[valid],axis=0,ddof=1) if valid.sum()>1 else None,MC95=np.quantile(B[valid],[.025,.975],axis=0).T if valid.sum()>1 else None,covariance=np.cov(B[valid].T) if valid.sum()>1 else None,samples=name+'.npz')

def main(args):
    source=Path(args.summary);data=json.loads((source/'analysis.json').read_text());groups=data['groups'];out=Path(args.output);out.mkdir(exist_ok=False)
    shutil.copyfile(__file__,out/Path(__file__).name);shutil.copyfile(sq.__file__,out/Path(sq.__file__).name)
    draws={(g['L'],g['q']):np.load(source/f"L{g['L']}-q{g['q']:.8f}-draws.npz")['samples'] for g in groups if g['status']!='insufficient_blocks'}
    pairs=[pair_result(groups,draws,L,kind,window,out) for L in [24,32] for kind in ['linear','pchip'] for window in ['narrow','full']]
    result=dict(source=str(source),source_sha256=a.ex.sha(source/'analysis.json'),settings=data['settings'],groups=groups,pairs=pairs,fixed_q_slope=data['fixed_q_slope'],scope='Existing retained regular checkpoints only. R and m2 use full scalar histories in shared time blocks when scalar_mode=dense_bins. R and log ensemble G are interpolated directly. Unique supported crossings required; root tails never clipped or branch-selected. Finite-pair Delta differs from asymptotic CFT dimension; no calibrated systematic error. Spatial sites/axes not independent samples; short chain block counts limit error assessment.')
    a.write(out/'analysis.json',result)
    with (out/'spin-correlation-measurements.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['L','epsilon_over_pi','checkpoint_fields_postcut','scalar_cycles_postcut','G_Lquarter','G_Lhalf','R_xi_over_L','logGquarter_stat_SD','logGhalf_stat_SD','flags'])
        for g in groups:w.writerow([g['L'],g['q'],g['samples_total'],g['scalar_cycles_total'],*g['raw_mean'][:2],g['derived'][0],*(g.get('SD') or [None]*4)[1:3],';'.join(g['flags'])])
    with (out/'finite-pair-dimensions.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['small_L','large_L','interpolation','window','status','valid_fraction','q_cross','Delta_quarter','Delta_half','Delta_m2','quarter_SD','half_SD','m2_SD'])
        for p in pairs:
            c=p['central'] if p['central'] is not None else [None]*5;sd=p['SD'] if p['SD'] is not None else [None]*5;w.writerow([*p['pair'],p['interpolation'],p['window'],p['status'],p['valid_fraction'],c[0],*c[2:],*sd[2:]])
    plt.rcParams.update({'font.size':11,'axes.grid':True,'grid.alpha':.2})
    fig,axes=plt.subplots(1,2,figsize=(11,4.4))
    for L in [24,32,48,64]:
        gg=sorted([g for g in groups if g['L']==L],key=lambda x:x['q']);q=np.array([g['q'] for g in gg]);raw=np.array([g['raw_mean'] for g in gg])
        for k,ax in enumerate(axes):
            errs=np.array([(g.get('raw_MC95') or [[np.nan,np.nan]]*6)[k] for g in gg]);line=ax.plot(q,raw[:,k],'o',label=f'L={L}')[0];ax.vlines(q,errs[:,0],errs[:,1],color=line.get_color())
            ax.set(xlabel=r'$\epsilon/\pi$',ylabel=rf'$G_L(L/{4 if k==0 else 2})$')
    axes[0].legend();fig.suptitle('3D, beta=0: long-distance spin correlations\nConditional block-bootstrap 95% bars; finite boxes, few blocks per start',fontsize=12);fig.tight_layout();save(fig,out/'spin-correlations-vs-epsilon')
    fig,ax=plt.subplots(figsize=(7,5));gg=sorted([g for g in groups if abs(g['q']-.80075)<1e-12],key=lambda x:x['L']);N=np.array([g['L'] for g in gg])
    for k,label in [(0,r'$r=L/4$'),(1,r'$r=L/2$')]:
        y=np.array([g['raw_mean'][k] for g in gg]);ci=np.array([(g.get('raw_MC95') or [[np.nan,np.nan]]*6)[k] for g in gg]);line=ax.plot(N,y,'o',label=label)[0];ax.vlines(N,ci[:,0],ci[:,1],color=line.get_color())
    slope=data['fixed_q_slope'];z=np.linspace(24,64,100)
    for k in [0,1]:
        if 'Delta' in slope:
            amp=np.exp(slope['log_amplitude_at_L32'][k]);ax.plot(z,amp*(z/32.)**(-2*slope['Delta']),ls='--',alpha=.6)
    title=rf"Fixed $\epsilon/\pi=0.80075$: $\Delta_{{eff}}={slope['Delta']:.3f}\pm{slope['SD']:.3f}$ (stat)" if 'Delta' in slope else 'Fixed epsilon/pi=0.80075: slope unresolved'
    ax.set(xscale='log',yscale='log',xlabel='L',ylabel=r'$G_L(r)$',title=title);ax.set_xticks(N,labels=[str(n) for n in N]);ax.legend();fig.text(.5,.01,'Common leading slope only; critical coupling and finite-size bias are not resolved by this fit.',ha='center',fontsize=9);fig.tight_layout(rect=(0,.04,1,1));save(fig,out/'fixed-epsilon-spin-scaling')
    fig,ax=plt.subplots(figsize=(8,4.5))
    for j,(kind,window) in enumerate([('linear','narrow'),('pchip','narrow'),('linear','full'),('pchip','full')]):
        for i,L in enumerate([24,32]):
            p=next(p for p in pairs if p['pair'][0]==L and p['interpolation']==kind and p['window']==window)
            if p['central'] is None or p['MC95'] is None:continue
            for k,marker in [(2,'o'),(3,'s')]:
                x=i*5+j+(k-2)*.18;y=p['central'][k];ci=p['MC95'][k];alpha=1 if p['valid_fraction']>=.99 else .35;ax.plot(x,y,marker,color=f'C{k-2}',alpha=alpha);ax.vlines(x,ci[0],ci[1],color=f'C{k-2}',alpha=alpha)
    ax.set_xticks([1.6,6.6],labels=['24:48','32:64']);ax.set(ylabel=r'Finite-pair $\Delta_{eff}$',xlabel='Box-size pair; within each: linear narrow, PCHIP narrow, linear full, PCHIP full',title='Long-distance quotients at common xi/L');ax.plot([],[],marker='o',ls='',color='C0',label='r=L/4');ax.plot([],[],marker='s',ls='',color='C1',label='r=L/2');ax.legend();fig.text(.5,.01,'Conditional MC 95% bars; pale = unresolved crossing tails. No asymptotic or total-error claim.',ha='center',fontsize=9);fig.tight_layout(rect=(0,.04,1,1));save(fig,out/'common-R-spin-dimensions')
    a.write(out/'status.json',dict(state='plots_and_quotients_complete_pending_independent_review',groups=len(groups),pairs=len(pairs)))
    for p in pairs:print(a.clean({k:p[k] for k in ['pair','interpolation','window','status','valid_fraction','central','SD']}),flush=True)

def save(fig,path):
    for ext in ['png','pdf','svg']:fig.savefig(str(path)+'.'+ext,dpi=180,bbox_inches='tight')
    plt.close(fig)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--summary',required=True);p.add_argument('--output',required=True);main(p.parse_args())
