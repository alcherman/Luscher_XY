"""Matched existing-data CV comparisons; no tail-selected raw variance claims."""
from pathlib import Path
import argparse,json,shutil,time
import numpy as np
import matplotlib.pyplot as plt
import analyze_spin_correlations as b
P=Path(__file__).resolve().parent

def fixed_samples(folder):
    d=json.loads((folder/'analysis.json').read_text());g=sorted([g for g in d['groups'] if abs(g['q']-.80075)<1e-12],key=lambda x:x['L'])
    if len(g)!=4:return None
    files=[folder/f"L{x['L']}-q{x['q']:.8f}-draws.npz" for x in g]
    if not all(f.exists() for f in files):return None
    B=np.stack([np.load(f)['samples'][:,1:3] for f in files],axis=1);valid=np.isfinite(B).all((1,2))
    if valid.sum()<2:return None
    C=np.stack([np.cov(B[valid,i].T) for i in range(4)]);W=np.linalg.inv(np.linalg.cholesky(C));X=np.zeros((4,2,3));X[:,0,0]=1;X[:,1,1]=1;X[:,:,2]=-2*np.log(np.array([x['L'] for x in g])[:,None]/32)
    A=np.einsum('nij,njk->nik',W,X).reshape(-1,3);T=np.linalg.pinv(A);out=np.full(len(B),np.nan);out[valid]=(np.einsum('nij,dnj->dni',W,B[valid]).reshape(valid.sum(),-1)@T.T)[:,-1]
    assert abs(out[valid].std(ddof=1)-d['fixed_q_slope']['SD'])<1e-12
    return out

def summary_change(x,y,centerx,centery):
    v=np.isfinite(x).all(1)&np.isfinite(y).all(1)
    if v.sum()<2:return dict(common_valid_fraction=float(v.mean()),status='insufficient_common_draws')
    dx=y[v]-x[v]
    return dict(common_valid_fraction=float(v.mean()),common_valid_draws=int(v.sum()),original_SD=x[v].std(0,ddof=1),CV_SD=y[v].std(0,ddof=1),SD_ratio=y[v].std(0,ddof=1)/x[v].std(0,ddof=1),central_change=np.array(centery)-np.array(centerx),change_SD=dx.std(0,ddof=1),change_MC95=np.quantile(dx,[.025,.975],axis=0).T,scope='Paired estimator changes on the same common-valid draws, not measured true bias. Root-conditioned Delta ratios are separate from unconditioned raw-G mean comparisons.')

def main(args):
    out=Path(args.output);out.mkdir(exist_ok=False);shutil.copyfile(__file__,out/Path(__file__).name);results=[];allgroups={};pairbanks={};fixed={}
    for tag in ['primary','cut4000','block4']:
        root=P/f'spin-control-variate-20261016-v2-{tag}';raw=json.loads((root/'comparison.json').read_text());allgroups[tag]=raw['groups'];A=json.loads((root/'original-quotients/analysis.json').read_text());C=json.loads((root/'CV-quotients/analysis.json').read_text());pairs=[]
        for a,c in zip(A['pairs'],C['pairs']):
            assert (a['pair'],a['interpolation'],a['window'])==(c['pair'],c['interpolation'],c['window'])
            x=np.load(root/'original-quotients'/a['samples'])['samples'];y=np.load(root/'CV-quotients'/c['samples'])['samples']
            v=np.isfinite(x).all(1)&np.isfinite(y).all(1)
            if v.any():np.testing.assert_array_equal(x[v,:2],y[v,:2])
            contrast=summary_change(x,y,a['central'],c['central'])
            r=dict(pair=a['pair'],interpolation=a['interpolation'],window=a['window'],original=a,CV=c,comparison=contrast);pairs.append(r)
            pairbanks[tag,a['pair'][0],a['interpolation'],a['window']]=[x,y]
        size_changes=[]
        for kind in ['linear','pchip']:
            for window in ['narrow','full']:
                x=pairbanks[tag,24,kind,window][1][:,2:4];y=pairbanks[tag,32,kind,window][1][:,2:4];v=np.isfinite(x).all(1)&np.isfinite(y).all(1)
                lo=next(r['CV']['central'][2:4] for r in pairs if r['pair'][0]==24 and r['interpolation']==kind and r['window']==window);hi=next(r['CV']['central'][2:4] for r in pairs if r['pair'][0]==32 and r['interpolation']==kind and r['window']==window)
                size_changes.append(dict(interpolation=kind,window=window,larger_minus_smaller_pair=summary_change(x,y,lo,hi),scope='Difference between disjoint size pairs; not a calibrated finite-size bias or systematic error.'))
        x=fixed_samples(root/'original');y=fixed_samples(root/'CV');fx=dict(status='insufficient_common_draws')
        if x is not None and y is not None:
            fx=summary_change(x[:,None],y[:,None],[A['fixed_q_slope']['Delta']],[C['fixed_q_slope']['Delta']]);np.savez_compressed(out/f'fixed-q-{tag}-paired.npz',original=x,CV=y)
        fixed[tag]=dict(original=A['fixed_q_slope'],CV=C['fixed_q_slope'],paired_change=fx)
        finite=[g for g in raw['groups'] if 'SD_ratio' in g];ratios=np.array([g['SD_ratio'] for g in finite])
        results.append(dict(variant=tag,pairs=pairs,CV_size_pair_changes=size_changes,fixed_q=fixed[tag],raw_gain_summary=dict(groups_attempted=len(raw['groups']),groups_with_finite_gain=len(finite),median_SD_ratio=np.median(ratios,axis=0),range_SD_ratio=[ratios.min(0),ratios.max(0)],count_SD_ratio_below_one=(ratios<1).sum(0),scope='Descriptive conditional bootstrap estimates, not independent repeated-run improvement.'),sources=[dict(path=str(root/f),sha256=b.ex.sha(root/f)) for f in ['comparison.json','original-quotients/analysis.json','CV-quotients/analysis.json']]))
    b.write(out/'analysis.json',dict(created=time.time(),variants=results,scope='Existing-data cross-fitted approximate CV test. Statistical gains are conditional bootstrap estimates; finite-size, stationarity and total-error coverage remain unresolved. No new MC.'))
    plt.rcParams.update({'font.size':11,'axes.grid':True,'grid.alpha':.2})
    fig,axes=plt.subplots(1,2,figsize=(10.8,4.2))
    for L in [24,32,48,64]:
        gg=sorted([g for g in allgroups['primary'] if g['L']==L],key=lambda x:x['q'])
        for j,ax in enumerate(axes):ax.plot([g['q'] for g in gg],[g['SD_ratio'][j] for g in gg],'o-',label=f'L={L}');ax.axhline(1,color='.5',lw=.5);ax.set(xlabel=r'$\epsilon/\pi$',ylabel='Estimated SD(CV) / SD(original)',title=f'r=L/{4 if j==0 else 2}')
    axes[0].legend();fig.suptitle('Raw correlation means: paired time-block bootstrap\nBefore positivity/crossing cuts; few blocks and two starts per coordinate',fontsize=12);fig.tight_layout();save(fig,out/'control-variate-mean-uncertainties')
    fig,ax=plt.subplots(figsize=(8,4.8));primary=results[0]
    for i,L in enumerate([24,32]):
        r=next(r for r in primary['pairs'] if r['pair'][0]==L and r['interpolation']=='linear' and r['window']=='narrow')
        x,y=pairbanks['primary',L,'linear','narrow'];v=np.isfinite(x).all(1)&np.isfinite(y).all(1)
        for name,B,marker,offset in [('original',x,'o',-.20),('CV',y,'s',.16)]:
            center=r[name]['central']
            for j,k in enumerate([2,3]):
                ci=np.quantile(B[v,k],[.025,.975]);pos=i+offset+j*.10;label=f'{name}, r=L/{4 if j==0 else 2}' if i==0 else None
                ax.plot(pos,center[k],marker,color=f'C{j}',label=label);ax.vlines(pos,*ci,color=f'C{j}')
    ax.set_xticks([0,1],labels=['24:48','32:64']);ax.set(xlim=(-.45,1.45),ylabel=r'Finite-pair $\Delta_{eff}$',title='Common xi/L: original and cross-fitted control variate');ax.legend(ncol=2,fontsize=9);fig.text(.5,.01,'Conditional MC95 bars, linear interpolation. Finite pairs; no asymptotic or total-error claim.',ha='center',fontsize=9);fig.tight_layout(rect=(0,.04,1,1));save(fig,out/'control-variate-common-R-dimensions')
    b.write(out/'status.json',dict(state='comparisons_complete_pending_independent_review',variants=3,new_MC_trajectories=0))

def save(fig,path):
    for ext in ['png','pdf','svg']:fig.savefig(str(path)+'.'+ext,dpi=180,bbox_inches='tight')
    plt.close(fig)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);main(p.parse_args())
