"""Cross-fitted approximate control variate on existing spin/scalar histories.

No simulations or archive writes. Coefficients use the OTHER independent start;
dense mean uncertainty and both cross-fitted coefficients are re-estimated in
shared physical-time-block bootstrap draws. Raw configuration regression is a
prespecified heuristic, not a time-mean variance optimum.
"""
from pathlib import Path
import argparse,json,hashlib,time,shutil
import numpy as np
import analyze_spin_correlations as base
P=Path(__file__).resolve().parent

def augmented_series(d,cut):
    times,a=base.aligned_series(d,cut,'dense_bins');cp=d['series'][d['sweep']>cut,2]
    return times,np.column_stack([a,cp,a[:,0]*cp,a[:,1]*cp,cp*cp])

def coefficients(mu):
    """Last4 columns: mCP, Gquarter*mCP, Ghalf*mCP, mCP^2."""
    mu=np.asarray(mu);var=mu[...,9]-mu[...,6]**2
    tol=64*np.finfo(float).eps*np.maximum(np.maximum(abs(mu[...,9]),mu[...,6]**2),np.finfo(float).tiny)
    good=np.isfinite(var)&(var>tol)
    with np.errstate(divide='ignore',invalid='ignore'):
        c=(mu[...,7:9]-mu[...,:2]*mu[...,6,None])/var[...,None]
    return np.where(good[...,None],c,np.nan)

def corrected_start(mu,c):
    out=np.asarray(mu)[...,:6].copy();out[...,:2]-=c*(mu[...,6]-mu[...,2])[...,None];return out

def crossfit(mu0,mu1,n0,n1):
    """Each mu is a central moment vector or a bank of resampled moments."""
    n0=np.asarray(n0);n1=np.asarray(n1);den=n0+n1
    c0=coefficients(mu1);c1=coefficients(mu0)
    x0=corrected_start(mu0,c0);x1=corrected_start(mu1,c1)
    cv=(x0*n0[...,None]+x1*n1[...,None])/den[...,None]
    original=(mu0[...,:6]*n0[...,None]+mu1[...,:6]*n1[...,None])/den[...,None]
    return original,cv,np.stack([c0,c1],axis=-2)

def raw_gain(original,cv):
    """Unlogged means; finite-paired mask, NEVER positivity/root conditioned."""
    valid=np.isfinite(original[:,:2]).all(1)&np.isfinite(cv[:,:2]).all(1)
    if valid.sum()<2:return dict(finite_pair_fraction=float(valid.mean()),status='insufficient_finite_draws')
    x=original[valid,:2];y=cv[valid,:2];sd0=x.std(0,ddof=1);sd1=y.std(0,ddof=1)
    return dict(finite_pair_fraction=float(valid.mean()),original_SD=sd0,CV_SD=sd1,SD_ratio=sd1/sd0,variance_ratio=(sd1/sd0)**2,paired_change_MC95=np.quantile(y-x,[.025,.975],axis=0).T,negative_mean_draw_counts=dict(original=(x<=0).sum(0),CV=(y<=0).sum(0)),scope='Same unlogged finite paired resamples; before positivity/crossing filtering. Descriptive conditional bootstrap uncertainty, not independent-run performance.')

def read_group(g,cut):
    rows=[];metas=[];files=[]
    for f in g['source_metadata']:
        f=Path(f);m=json.loads(f.read_text());data=f.with_suffix('.npz')
        if base.ex.sha(data)!=m['derived_data_sha256']:raise ValueError('Changed measured history')
        with np.load(data) as d:times,x=augmented_series(d,cut)
        rows.append(x);metas.append(m);files.append(dict(metadata=str(f),metadata_sha256=base.ex.sha(f),data=str(data),data_sha256=base.ex.sha(data)))
    order=np.argsort([m['id'] for m in metas]);return [rows[i] for i in order],[metas[i] for i in order],[files[i] for i in order]

def summarize(g,cut,factor):
    rows,metas,files=read_group(g,cut)
    if len(rows)!=2 or len(set(m['id'] for m in metas))!=2:raise ValueError('Exactly two distinct independent starts required')
    mu=[r.mean(0) for r in rows];counts=[len(r) for r in rows];orig,cv,c=crossfit(*mu,*counts)
    # Keep every original diagnostic/flag and add corrected-series blocking checks.
    diagnostics=[];widths=[]
    old={v['id']:v for v in g['diagnostics']}
    for i,r in enumerate(rows):
        res=r[:,:2]-c[i]*(r[:,6]-r[:,2])[:,None]
        tau={}
        for j,name in enumerate(['Gquarter_CV','Ghalf_CV']):
            tau[name+'_postcut']=base.tau(res[:,j]);tau[name+'_latehalf']=base.tau(res[len(res)//2:,j])
        w,cycles=base.block_width_cycles(tau,{},100,factor)
        w=max(w,old[metas[i]['id']]['blocking']['width_checkpoints']);widths.append(w)
        diagnostics.append(dict(id=metas[i]['id'],source_start_for_coefficient=metas[1-i]['id'],applied_coefficient=c[i],checkpoint_tau_CV=tau,previous_width_checkpoints=old[metas[i]['id']]['blocking']['width_checkpoints'],width_checkpoints=w,physical_width_cycles=100*w,samples=counts[i]))
    seed=int(hashlib.sha256(json.dumps([[m['id'] for m in metas],cut,factor,'spin-CV-crossfit-v1']).encode()).hexdigest()[:16],16);rng=np.random.default_rng(seed)
    banks=[];flags=list(g['flags'])
    for r,w,di in zip(rows,widths,diagnostics):
        b,info=base.bootstrap_blocks(r,w,rng);di['blocking']=info
        halves=[]
        for half in np.array_split(r,2):
            nb=len(half)//w
            if nb<2:halves.append(None);continue
            blockmeans=np.array([corrected_start(b.mean(0),di['applied_coefficient'])[:2] for b in np.array_split(half,nb)])
            halves.append(dict(mean=corrected_start(half.mean(0),di['applied_coefficient'])[:2],SE=blockmeans.std(0,ddof=1)/np.sqrt(nb),blocks=nb))
        di['CV_half_summaries_fixed_other_start_coefficient']=halves
        di['CV_half_scope']='Descriptive fixed-trained-coefficient half-block check; coefficient uncertainty not included here. Later-cut study and crossfit bootstrap remain separate.'
        if all(h is not None for h in halves):
            zz=(halves[1]['mean']-halves[0]['mean'])/np.maximum(np.hypot(halves[0]['SE'],halves[1]['SE']),1e-15);di['CV_half_z']=zz
            if np.max(abs(zz))>4:flags.append('CV_half_shift_fixed_coefficient')
        else:flags.append('insufficient_CV_half_blocks')
        if not np.isfinite(di['applied_coefficient']).all():flags.append('invalid_central_training_variance')
        if info['blocks']<50:flags.append('fewer_than_50_blocks_per_start')
        if b is None:flags.append('insufficient_CV_aware_blocks');banks.append(None)
        else:banks.append(b)
    common=dict(L=g['L'],q=g['q'],ids=[m['id'] for m in metas],flags=sorted(set(flags)),diagnostics=diagnostics,source_metadata=g['source_metadata'],samples_total=sum(counts),scalar_cycles_total=100*sum(counts),scalar_mode='dense_bins',original_sampling_diagnostics=g['diagnostics'],control=dict(formula='G_a-c_other*(m2CP_a-m2dense_a); exchange starts, row-weighted pool',coefficient_method='Other-start raw configuration covariance/variance, refitted in every bootstrap draw; heuristic not time-mean optimum',central_applied_coefficients=c,checkpoint_minus_dense_m2=[u[6]-u[2] for u in mu],coefficient_training_start_ids=[metas[1]['id'],metas[0]['id']]))
    if any(b is None for b in banks):
        variants=[dict(common,raw_mean=x,derived=base.derived(x,g['L']),status='insufficient_blocks',joint_valid_fraction=0.) for x in [orig,cv]]
        return variants,None,None,None,dict(status='insufficient_blocks',L=g['L'],q=g['q'],ids=[m['id'] for m in metas],diagnostics=diagnostics),files
    mus=[b[0]/b[1][:,None] for b in banks];O,C,coeffs=crossfit(*mus,banks[0][1],banks[1][1]);Z=[base.derived(v,g['L']) for v in [O,C]]
    variants=[]
    for name,x,B,z in zip(['original','crossfit_CV'],[orig,cv],[O,C],Z):
        valid=np.isfinite(z).all(1);central=base.derived(x,g['L'])
        variants.append(dict(common,estimator=name,raw_mean=x,derived=central,status='conditional_checkpoint_block_summary' if valid.mean()>=.99 and np.isfinite(central).all() else 'unresolved_joint_tails',joint_valid_fraction=float(valid.mean()),SD=z[valid].std(0,ddof=1) if valid.sum()>1 else None,MC95=np.quantile(z[valid],[.025,.975],axis=0).T if valid.sum()>1 else None,raw_MC95=np.quantile(B[np.isfinite(B).all(1)],[.025,.975],axis=0).T if np.isfinite(B).all(1).sum()>1 else None,invalid_draws=dict(R=int(np.sum(~np.isfinite(z[:,0]))),G_quarter=int(np.sum(~np.isfinite(z[:,1]))),G_half=int(np.sum(~np.isfinite(z[:,2]))),m2=int(np.sum(~np.isfinite(z[:,3]))))))
    start_arrays=[corrected_start(mus[i],coeffs[:,i])[:,:2] for i in range(2)]
    start_delta=start_arrays[1]-start_arrays[0];start_valid=np.isfinite(start_delta).all(1)
    center_delta=corrected_start(mu[1],c[1])[:2]-corrected_start(mu[0],c[0])[:2]
    start_diag=dict(central_difference=center_delta,finite_fraction=float(start_valid.mean()),SD=start_delta[start_valid].std(0,ddof=1) if start_valid.sum()>1 else None,MC95=np.quantile(start_delta[start_valid],[.025,.975],axis=0).T if start_valid.sum()>1 else None,scope='Raw CV mean start1-start0, with both other-start coefficients refitted; shared-training covariance retained. Exploratory stationarity diagnostic, not a mixing certificate.')
    for v in variants:v['CV_start_difference']=start_diag
    validcoeff=np.isfinite(coeffs).all((1,2));gain=raw_gain(O,C)
    gain.update(L=g['L'],q=g['q'],coefficient_invalid_draws=int((~validcoeff).sum()),coefficient_MC95=np.quantile(coeffs[validcoeff],[.025,.975],axis=0) if validcoeff.sum()>1 else None,joint_log_valid_fraction=float((np.isfinite(Z[0]).all(1)&np.isfinite(Z[1]).all(1)).mean()),central_original_G=orig[:2],central_CV_G=cv[:2])
    return variants,Z,O,C,gain,files

def main(args):
    out=Path(args.output);out.mkdir(exist_ok=False);shutil.copyfile(__file__,out/Path(__file__).name)
    tag=args.variant;settings=dict(primary=(1000,2),block4=(1000,4),cut4000=(4000,2));cut,factor=settings[tag]
    source=P/f'spin-correlation-grid-summary-20261016-{tag}';d=json.loads((source/'analysis.json').read_text())
    groups=[[],[]];banks=[{},{}];gains=[];provenance=[]
    for name in ['original','CV']:(out/name).mkdir()
    for g in d['groups']:
        summaries,Z,O,C,gain,files=summarize(g,cut,factor);gains.append(gain);provenance+=files
        for k,name in enumerate(['original','CV']):
            groups[k].append(summaries[k])
            if Z is not None:
                banks[k][g['L'],g['q']]=Z[k];np.savez_compressed(out/name/f"L{g['L']}-q{g['q']:.8f}-draws.npz",samples=Z[k])
        if Z is not None:np.savez_compressed(out/f"L{g['L']}-q{g['q']:.8f}-raw-paired.npz",original=O,CV=C)
        print(g['L'],g['q'],'SD_ratio',base.clean(gain.get('SD_ratio')),'blocks',[v['blocking']['blocks'] for v in summaries[0]['diagnostics']],flush=True)
    for k,name in enumerate(['original','CV']):
        selected=[g for g in groups[k] if abs(g['q']-.80075)<1e-12];B=[banks[k].get((g['L'],g['q'])) for g in selected]
        slope=base.fixed_q_slope(selected,np.stack(B,axis=1)) if len(selected)==4 and all(b is not None for b in B) else dict(status='incomplete_sizes')
        record=dict(created=time.time(),settings=dict(cut_cycles=cut,block_factor=factor,scalar_mode='dense_bins',estimator=name,draws=base.DRAWS,coefficient_refit=True,block_rule='Maximum frozen original and new corrected-series full/late-half width'),groups=groups[k],fixed_q_slope=slope,provenance=provenance,scope='Exploratory cross-fitted CV. Same time-block draws for original/CV; training uncertainty included. Conditional bootstrap, not certified systematic or performance gain.')
        base.write(out/name/'analysis.json',record)
    base.write(out/'comparison.json',dict(variant=tag,groups=gains,source_summary=str(source),source_summary_sha256=base.ex.sha(source/'analysis.json'),scope='Variance gains use unlogged paired means before positivity/root masks. Few blocks and only two MC starts limit assessment. No raw MC modifications.'))
    base.write(out/'status.json',dict(state='summaries_complete_pending_review',variant=tag,groups=len(groups[0])))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--variant',choices=['primary','block4','cut4000'],default='primary');p.add_argument('--output',required=True);main(p.parse_args())
