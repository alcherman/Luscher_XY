"""Disk-only pilot diagnostics; these fits do not certify a phase boundary."""
import os
for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS']:
    os.environ[name]='1'
os.environ['MPLBACKEND']='Agg'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import importlib.util
import json
from pathlib import Path
import sys

import h5py
import numpy as np
from scipy.stats import chi2
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(OUT))
import run
helper_spec=importlib.util.spec_from_file_location('old_xy_diagnostics',OUT.parent/'2026-09-02-xy-binder-bkt/experiment.py')
helper=importlib.util.module_from_spec(helper_spec);helper_spec.loader.exec_module(helper)


def ratios(means,N):
    m2,m4,f0,fm=means[...,0],means[...,1],means[...,2],means[...,3]
    with np.errstate(divide='ignore',invalid='ignore'):
        return np.stack([1-m4/(2*m2*m2),np.sqrt(f0/fm-1)/(2*N*np.sin(np.pi/N)),m2*m2/m4],axis=-1)


def read(point):
    with h5py.File(OUT/'chains'/f"{point['id']}.h5") as f:
        a=f['magnetization_squared'][:];b=f['magnetization_fourth'][:]
        f0=f['structure_zero'][:];fm=f['structure_min'][:].mean(1)
        w=np.empty(len(a))
        for first in range(0,len(a),100):
            stop=min(first+100,len(a))
            w[first:stop]=(f['winding'][first:stop].astype(float)**2).mean((1,2))
        return np.column_stack([a,b,f0,fm,f['mean_link_cosine'][:],f['vortex_density'][:],w])


def summarize(chains,points,cut,draws=1000,factor=1):
    arrays=[a[cut:] for a in chains];N=points[0]['N'];rng=np.random.default_rng([20261013,N,points[0]['D'],round(points[0]['beta']*1000),round(points[0].get('epsilon_over_pi',1)*1000),cut,factor])
    central=ratios(np.concatenate(arrays).mean(0),N)
    total=np.zeros((draws,7));length=np.zeros(draws);diagnostics=[];individual=[];individual_boot=[]
    for a,p in zip(arrays,points):
        m2,m4,f0,fm=a[:,:4].T
        influence=m4.mean()/m2.mean()**3*(m2-m2.mean())-(m4-m4.mean())/(2*m2.mean()**2)
        ri=(f0-f0.mean())/fm.mean()-f0.mean()/fm.mean()**2*(fm-fm.mean())
        monitor={'m2':m2,'m4':m4,'binder':influence,'xi_ratio':ri,'cosine':a[:,4],'vortex':a[:,5],'winding2':a[:,6]}
        taus={k:helper.tau_int(x) for k,x in monitor.items()}
        width=factor*max(20,int(np.ceil(10*max(t for t in taus.values() if t is not None))))
        nb=len(a)//width
        diagnostics.append(dict(id=p['id'],retained=len(a),block_width=width,blocks=nb,taus=taus))
        individual.append(ratios(a.mean(0),N))
        if nb<2:
            return dict(status='insufficient_blocks',diagnostics=diagnostics),None
        edges=np.linspace(0,len(a),nb+1,dtype=int);sums=np.add.reduceat(a,edges[:-1],axis=0);sizes=np.diff(edges)
        indices=rng.integers(nb,size=(draws,nb));sample_sum=sums[indices].sum(1);sample_n=sizes[indices].sum(1)
        total+=sample_sum;length+=sample_n;individual_boot.append(ratios(sample_sum/sample_n[:,None],N))
    samples=ratios(total/length[:,None],N);valid=np.isfinite(samples).all(1)
    flags=[]
    if not valid.all():flags.append('invalid_ratio_bootstrap')
    if any(d['blocks']<50 for d in diagnostics):flags.append('fewer_than_50_blocks_per_start')
    for key,index in [('binder',0),('xi_over_N',1),('Q',2)]:
        vals=np.asarray(individual)[:,index];errs=np.asarray(individual_boot)[:,:,index].std(1,ddof=1)
        eligible=errs>0
        if eligible.sum()>1:
            center=np.average(vals[eligible],weights=1/errs[eligible]**2)
            pv=float(chi2.sf(np.sum(((vals[eligible]-center)/errs[eligible])**2),eligible.sum()-1))
            if pv<.001:flags.append('start_disagreement_'+key)
    finite=samples[valid]
    return dict(status='pilot_diagnostics' if not flags else 'needs_diagnostic_review',values=central.tolist(),
                errors=finite.std(0,ddof=1).tolist(),covariance=np.cov(finite.T).tolist(),
                flags=flags,diagnostics=diagnostics,bootstrap_valid_fraction=float(valid.mean()),
                individual_values=np.asarray(individual).tolist(),cut=cut),samples


def fss_fit(rows,samples,D,observable=0,minimum_N=0,subleading=False):
    """Linearized finite-size fit with explicit identifiability and ratio tails."""
    eligible=[i for i,r in enumerate(rows) if r['D']==D and r['N']>=minimum_N and 'values' in r['summary']]
    if not eligible:return dict(status='no_data')
    Ns=np.array([rows[i]['N'] for i in eligible],dtype=float);betas=np.array([rows[i]['beta'] for i in eligible]);origin=betas.mean()
    if D==3:thermal=Ns**(1/.67169);corr=Ns**(-.789)
    else:thermal=Ns**2*np.log(Ns)**.1;corr=np.log(Ns)**(-.5)
    columns=[np.ones(len(Ns)),(betas-origin)*thermal,thermal,corr]
    if subleading:columns.append(np.log(np.log(Ns))/np.log(Ns) if D==4 else Ns**(-2*.789))
    X=np.column_stack(columns);scale=np.linalg.norm(X,axis=0)
    normalized=X/scale
    rank=int(np.linalg.matrix_rank(normalized));condition=float(np.linalg.cond(normalized))
    if rank<X.shape[1]:return dict(status='unidentifiable',rank=rank,parameters=X.shape[1],sizes=sorted(set(Ns.astype(int).tolist())))
    y=np.array([rows[i]['summary']['values'][observable] for i in eligible]);errors=np.array([rows[i]['summary']['errors'][observable] for i in eligible])
    if not np.isfinite(errors).all() or np.any(errors<=0):return dict(status='invalid_errors')
    transform=np.linalg.pinv(normalized/errors[:,None])/errors[None,:]
    coef=(transform@y)/scale
    beta_c=origin-coef[2]/coef[1]
    B=np.array([samples[i][:,observable] for i in eligible]).T
    coefficients=(B@transform.T)/scale
    with np.errstate(divide='ignore',invalid='ignore'):
        roots=origin-coefficients[:,2]/coefficients[:,1]
    valid=np.isfinite(roots)&(coefficients[:,1]>0)
    if valid.mean()<.99:return dict(status='unresolved_bootstrap',positive_slope_fraction=float(valid.mean()))
    statistic=float(np.sum(((y-X@coef)/errors)**2));df=len(y)-len(coef)
    pv=float(chi2.sf(statistic,df)) if df>0 else None
    reasons=[]
    if pv is not None and pv<.01:reasons.append('chi_square_rejects_linear_scaling_window')
    if not betas.min()<=beta_c<=betas.max():reasons.append('root_outside_sampled_beta_window')
    if not 0<coef[0]<(1 if D==4 else .5):reasons.append('unphysical_fitted_critical_ratio')
    status='rejected_pilot_fit' if reasons else 'conditional_pilot_fit'
    return dict(status=status,rejection_reasons=reasons,beta_c=float(beta_c),bootstrap_95=np.quantile(roots[valid],[.025,.975]).tolist(),
                chi2=statistic,df=df,p=pv,
                rank=rank,condition=condition,sizes=sorted(set(Ns.astype(int).tolist())),
                beta_interval=[float(betas.min()),float(betas.max())],parameters=coef.tolist(),
                input_flags=sorted({flag for i in eligible for flag in rows[i]['summary']['flags']}),
                limitation='Conditional ansatz on pilot boxes; not a validated infinite-volume critical coupling.')


def main():
    folder=OUT/'calibration-analysis';folder.mkdir(exist_ok=True)
    grouped={}
    for path in (OUT/'jobs').glob('*.json'):
        job=json.loads(path.read_text());p=job['point']
        if job['state']=='complete' and p['epsilon_over_pi']==1 and job['audit']['passed']:
            grouped.setdefault((p['D'],p['beta'],p['N']),[]).append(p)
    rows=[];boots=[]
    for (D,beta,N),points in sorted(grouped.items()):
        if len(points)<2:continue
        points=sorted(points,key=lambda p:p['replica']);chains=[read(p) for p in points]
        primary,bootstrap=summarize(chains,points,1001)
        row=dict(D=D,beta=beta,N=N,summary=primary,source_ids=[p['id'] for p in points],stability=[])
        if bootstrap is None:continue
        for cut,factor in [(501,1),(2001,1),(1001,2)]:
            other,_=summarize(chains,points,cut,factor=factor)
            row['stability'].append(dict(cut=cut,width_factor=factor,summary=other))
            if 'values' in other:
                diff=np.abs(np.array(other['values'])-primary['values'])
                se=np.hypot(other['errors'],primary['errors'])
                if np.any(diff>2*se):primary['flags'].append('cut_or_block_sensitive')
        fig,axes=plt.subplots(3,1,figsize=(8,6),sharex=True)
        for data,p in zip(chains,points):
            for ax,index in zip(axes,[0,4,6]):
                n=(len(data)-1)//50;x=data[1:1+50*n,index].reshape(n,50).mean(1)
                ax.plot(np.arange(n)*50+25,x,lw=.8,label=str(p['initial_winding']))
        for ax,label in zip(axes,['magnetization squared','mean link cosine','mean squared winding']):
            ax.axvline(1001,color='black',ls='--',lw=.8);ax.set_ylabel(label);ax.grid(alpha=.2)
        axes[0].legend(fontsize=7);axes[-1].set_xlabel('cycle');fig.suptitle(f'Pilot d={D}, beta={beta}, N={N}');fig.tight_layout()
        fig.savefig(folder/f'd{D}_b{int(beta*1000):04d}_N{N:03d}_histories.png',dpi=130);plt.close(fig)
        np.savez_compressed(folder/f'd{D}_b{int(beta*1000):04d}_N{N:03d}_bootstrap.npz',samples=bootstrap)
        rows.append(row);boots.append(bootstrap)
        print('analysed',D,beta,N,'flags',primary['flags'],flush=True)
    fits=[]
    for D in [3,4]:
        for minimum in [6,8]:
            for correction in [False,True]:
                observable=0 if D==3 else 2
                fits.append(dict(D=D,minimum_N=minimum,subleading=correction,observable='Binder U' if D==3 else 'Binder Q',result=fss_fit(rows,boots,D,observable,minimum,correction)))
    run.atomic_json(folder/'diagnostics.json',dict(created=run.now(),groups=rows,fits=fits,
        source_sha256=run.sha256(Path(__file__)),interpretation='Pilot diagnostics and conditional fits only; full fit validation and independent scientific review remain necessary.'))
    print('conditional_fits',json.dumps(fits),flush=True)


if __name__=='__main__':main()
