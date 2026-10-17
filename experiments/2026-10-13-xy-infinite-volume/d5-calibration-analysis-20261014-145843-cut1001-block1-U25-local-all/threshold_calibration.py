"""Directly bracketed Binder-quarter locations and simple size extrapolations."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
os.environ['MPLBACKEND']='Agg'
for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS']:
    os.environ[name]='1'
import csv
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.stats import chi2
import matplotlib.pyplot as plt
import run
from analyze_calibration import read,summarize

OUT=Path(__file__).resolve().parent


def crossing(beta,values,threshold=.25,method='linear'):
    """Interpolate inside a unique measured increasing bracket; never extrapolate."""
    beta=np.asarray(beta);values=np.asarray(values)
    hits=np.flatnonzero((values[:-1]<=threshold)&(values[1:]>threshold))
    if len(hits)!=1:return np.nan
    i=hits[0]
    if method=='linear':
        return float(beta[i]+(beta[i+1]-beta[i])*(threshold-values[i])/(values[i+1]-values[i]))
    if method!='pchip':raise ValueError(method)
    spline=PchipInterpolator(beta,values,extrapolate=False)
    return float(brentq(lambda b:float(spline(b)-threshold),beta[i],beta[i+1],xtol=1e-12))


def coordinate(D,N):
    if D not in [3,4]:raise ValueError('This size coordinate is only for d3/d4')
    N=np.asarray(N,dtype=float)
    return N**(-1/.67169) if D==3 else N**(-2)*np.log(N)**(-.1)


def monotone_values(values,weights=None):
    """Weighted pool-adjacent-violators regression; a declared shape assumption."""
    values=np.asarray(values,dtype=float)
    weights=np.ones(len(values)) if weights is None else np.asarray(weights,dtype=float)
    if np.any(weights<=0) or not np.isfinite(values).all() or not np.isfinite(weights).all():
        raise ValueError('Finite values and positive weights required')
    blocks=[]
    for i,(v,w) in enumerate(zip(values,weights)):
        blocks.append([i,i+1,w,w*v])
        while len(blocks)>1 and blocks[-2][3]/blocks[-2][2]>blocks[-1][3]/blocks[-1][2]:
            right=blocks.pop();left=blocks.pop()
            blocks.append([left[0],right[1],left[2]+right[2],left[3]+right[3]])
    result=np.empty(len(values))
    for first,last,weight,total in blocks:result[first:last]=total/weight
    return result


def monotone_crossing(beta,values,threshold=.25,method='linear',weights=None):
    values=monotone_values(values,weights)
    exact=np.flatnonzero(np.isclose(values,threshold,rtol=0,atol=1e-15))
    if len(exact)>1:
        if len(exact)==len(values):return np.nan
        return float((np.asarray(beta)[exact[0]]+np.asarray(beta)[exact[-1]])/2)
    return crossing(beta,values,threshold,method)


def size_fit(D,N,values,boot,minimum=0,correction=False):
    keep=np.asarray(N)>=minimum;N=np.asarray(N,dtype=float)[keep]
    values=np.asarray(values)[keep];boot=np.asarray(boot)[:,keep]
    x=coordinate(D,N);columns=[np.ones(len(N)),x]
    if correction:columns.append(x*(N**(-.789) if D==3 else np.log(N)**(-.5)))
    X=np.column_stack(columns)
    if len(N)<X.shape[1] or np.linalg.matrix_rank(X)<X.shape[1]:return {'status':'unidentifiable'}
    se=boot.std(0,ddof=1)
    scale=np.linalg.norm(X,axis=0);T=np.linalg.pinv(X/scale/se[:,None])/se[None,:]
    coef=(T@values)/scale;coefboot=(boot@T.T)/scale
    statistic=float(np.sum(((values-X@coef)/se)**2));df=len(N)-X.shape[1]
    return dict(status='zero_dof_provisional' if df==0 else 'conditional_size_fit',
                sizes=N.astype(int).tolist(),minimum_N=minimum,correction=correction,
                beta_infinite=float(coef[0]),bootstrap_95=np.quantile(coefboot[:,0],[.025,.975]).tolist(),
                coefficients=coef.tolist(),chi2=statistic,df=df,p=float(chi2.sf(statistic,df)) if df else None,
                intercept_bootstrap=coefboot[:,0])


def main(output_name='binder-quarter-calibration-2026-10-13'):
    folder=OUT/output_name;folder.mkdir(exist_ok=True)
    grouped={};jobs=[]
    for path in sorted((OUT/'jobs').glob('*.json')):
        j=json.loads(path.read_text());p=j['point']
        if j['state']=='complete' and p['epsilon_over_pi']==1 and j['audit']['passed']:
            grouped.setdefault((p['D'],p['N'],p['beta']),[]).append(p)
            jobs.append(j)
    groups=[];samples={}
    for (D,N,beta),points in sorted(grouped.items()):
        if len(points)!=2:continue
        points=sorted(points,key=lambda p:p['replica'])
        chains=[read(p) for p in points]
        summary,boot=summarize(chains,points,1001)
        if boot is None:continue
        # The Binder component remains usable when noisy xi samples are undefined.
        if not np.isfinite(boot[:,0]).all():raise RuntimeError('Nonfinite Binder samples')
        summary['binder_error_unfiltered']=float(boot[:,0].std(ddof=1))
        stability=[]
        for cut,factor in [(501,1),(2001,1),(1001,2)]:
            alternate,alternate_boot=summarize(chains,points,cut,factor=factor)
            stability.append(dict(cut=cut,block_factor=factor,summary=alternate))
            if alternate_boot is not None:
                np.savez_compressed(folder/f'd{D}_N{N:03d}_b{round(beta*1000):04d}_cut{cut}_factor{factor}_bootstrap.npz',samples=alternate_boot)
        group=dict(D=D,N=N,beta=beta,summary=summary,source_ids=[p['id'] for p in points],stability=stability)
        groups.append(group);samples[D,N,beta]=boot[:,0]
        np.savez_compressed(folder/f'd{D}_N{N:03d}_b{round(beta*1000):04d}_bootstrap.npz',samples=boot)
    rows=[];boot_locations={};fits=[]
    for D,N in sorted({(g['D'],g['N']) for g in groups}):
        g=sorted([a for a in groups if (a['D'],a['N'])==(D,N)],key=lambda a:a['beta'])
        beta=np.array([a['beta'] for a in g]);u=np.array([a['summary']['values'][0] for a in g])
        row=dict(D=D,N=N,sampled_beta=beta.tolist(),sampled_U=u.tolist(),
                 flags=sorted({f for a in g for f in a['summary']['flags']}))
        central=crossing(beta,u)
        if not np.isfinite(central):
            row.update(status='unbracketed',beta_location=None,
                       reason='Only one beta sampled' if len(beta)<2 else 'No unique interior Binder-quarter bracket')
            rows.append(row);continue
        hits=np.flatnonzero((u[:-1]<=.25)&(u[1:]>.25));i=int(hits[0])
        B=np.column_stack([samples[D,N,b] for b in beta])
        roots={method:np.array([crossing(beta,b,method=method) for b in B]) for method in ['linear','pchip']}
        valid=np.isfinite(roots['linear'])&np.isfinite(roots['pchip'])
        if valid.mean()<.99:
            row.update(status='unresolved_bootstrap_bracket',valid_fraction=float(valid.mean()));rows.append(row);continue
        row.update(status='directly_bracketed',beta_location=central,
                   bootstrap_95=np.quantile(roots['linear'][valid],[.025,.975]).tolist(),
                   bootstrap_sd=float(roots['linear'][valid].std(ddof=1)),bracket=beta[i:i+2].tolist(),
                   pchip_location=crossing(beta,u,method='pchip'),valid_fraction=float(valid.mean()))
        # Index alignment is retained across interpolation methods and sizes.
        boot_locations[D,N]=roots
        np.savez_compressed(folder/f'd{D}_N{N:03d}_crossings.npz',**roots)
        rows.append(row)
    for D in [3,4]:
        selected=[r for r in rows if r['D']==D and r['status']=='directly_bracketed']
        N=np.array([r['N'] for r in selected])
        configurations=[(int(N.min()),False)]
        if len(N)>2:configurations += [(int(N[1]),False),(int(N[2]),False),(int(N.min()),True),(int(N[1]),True)]
        for method in ['linear','pchip']:
            values=[r['beta_location'] if method=='linear' else r['pchip_location'] for r in selected]
            B=np.column_stack([boot_locations[D,n][method] for n in N]);valid=np.isfinite(B).all(1)
            for minimum,correction in configurations:
                result=size_fit(D,N,values,B[valid],minimum,correction)
                if 'intercept_bootstrap' in result:
                    np.savez_compressed(folder/f'd{D}_{method}_minN{minimum}_correction{int(correction)}_fit.npz',intercepts=result.pop('intercept_bootstrap'))
                result['input_flags']=sorted({f for r in selected if r['N']>=minimum for f in r['flags']})
                fits.append(dict(D=D,interpolation=method,minimum_N=minimum,correction=correction,result=result))
    cut_block_checks=[]
    for D in [3,4]:
        for cut,factor in [(501,1),(2001,1),(1001,2)]:
            selected=[];boots=[]
            for N in sorted({g['N'] for g in groups if g['D']==D}):
                gs=sorted([g for g in groups if g['D']==D and g['N']==N],key=lambda g:g['beta'])
                alternatives=[next(s['summary'] for s in g['stability'] if s['cut']==cut and s['block_factor']==factor) for g in gs]
                if any('values' not in a for a in alternatives):continue
                beta=[g['beta'] for g in gs];values=[a['values'][0] for a in alternatives]
                location=crossing(beta,values)
                if not np.isfinite(location):continue
                B=np.column_stack([np.load(folder/f'd{D}_N{N:03d}_b{round(b*1000):04d}_cut{cut}_factor{factor}_bootstrap.npz')['samples'][:,0] for b in beta])
                roots=np.array([crossing(beta,b) for b in B])
                selected.append(dict(N=N,beta_location=location,flags=sorted({f for a in alternatives for f in a['flags']})))
                boots.append(roots)
            if len(selected)>=2:
                B=np.column_stack(boots);valid=np.isfinite(B).all(1)
                result=size_fit(D,[r['N'] for r in selected],[r['beta_location'] for r in selected],B[valid])
                if 'intercept_bootstrap' in result:
                    np.savez_compressed(folder/f'd{D}_cut{cut}_factor{factor}_fit.npz',intercepts=result.pop('intercept_bootstrap'))
                cut_block_checks.append(dict(D=D,cut=cut,block_factor=factor,locations=selected,result=result,valid_fraction=float(valid.mean())))
    provenance=dict(created=run.now(),threshold=.25,epsilon_over_pi=1.,groups=groups,locations=rows,fits=fits,cut_block_checks=cut_block_checks,
        source_jobs=[dict(id=j['point']['id'],audit=j['audit']) for j in jobs],
        analysis_sha256=run.sha256(Path(__file__)),bootstrap_draws=1000,cut=1001,
        limitation='No histogram reweighting or extrapolation in beta. Linear and shape-preserving cubic interpolation are compared. Fit intervals are conditional Monte Carlo intervals, not total accuracy bounds. Any zero-degree-of-freedom variant has no goodness-of-fit test. Valid-draw fractions are recorded; cut/block results are conditional on retaining complete crossing vectors. No accepted full-campaign critical coupling is asserted.')
    from scout import json_finite
    run.atomic_json(folder/'analysis.json',json_finite(provenance))
    with (folder/'finite-size-locations.csv').open('w') as stream:
        writer=csv.writer(stream);writer.writerow(['d','N','U_threshold','status','beta_location','MC_95_low','MC_95_high','beta_bracket_low','beta_bracket_high','pchip_location'])
        for r in rows:writer.writerow([r['D'],r['N'],.25,r['status'],r.get('beta_location'),*r.get('bootstrap_95',['','']),*r.get('bracket',['','']),r.get('pchip_location')])
    plt.rcParams.update({'font.size':11})
    def panel(ax,D):
        selected=[r for r in rows if r['D']==D and r['status']=='directly_bracketed']
        base=next(f['result'] for f in fits if f['D']==D and f['interpolation']=='linear' and not f['correction'])
        x=coordinate(D,[r['N'] for r in selected]);y=np.array([r['beta_location'] for r in selected])
        lo=np.array([r['bootstrap_95'][0] for r in selected]);hi=np.array([r['bootstrap_95'][1] for r in selected])
        ax.errorbar(x,y,yerr=np.array([y-lo,hi-y]),fmt='o',capsize=4,color='tab:blue',label='Finite boxes: local linear crossing',zorder=3)
        flagged=np.array([bool(r['flags']) for r in selected])
        if flagged.any():ax.scatter(x[flagged],y[flagged],s=85,facecolors='none',edgecolors='black',lw=.8,zorder=4)
        for xi,yi,r in zip(x,y,selected):ax.annotate(f"N={r['N']}",(xi,yi),xytext=(5,7),textcoords='offset points',fontsize=9)
        grid=np.linspace(0,x.max()*1.06,201);c=base['coefficients']
        ax.plot(grid,c[0]+c[1]*grid,color='tab:blue',label='Leading size extrapolation')
        interval=base['bootstrap_95'];ax.errorbar([0],[c[0]],yerr=[[c[0]-interval[0]],[interval[1]-c[0]]],fmt='D',color='tab:blue',capsize=4)
        alternate=next(f['result'] for f in fits if f['D']==D and f['interpolation']=='pchip' and not f['correction'])
        c2=alternate['coefficients'];ax.plot(grid,c2[0]+c2[1]*grid,color='tab:orange',ls='--',label='Cubic crossing interpolation')
        ax.scatter(x,[r['pchip_location'] for r in selected],facecolors='none',edgecolors='tab:orange',marker='s',zorder=3)
        if D==4:
            cut=next(f['result'] for f in fits if f['D']==D and f['interpolation']=='linear' and f['minimum_N']==8 and not f['correction'])
            ax.plot(grid,cut['coefficients'][0]+cut['coefficients'][1]*grid,ls=':',color='tab:green',label='Leading fit, omit N=6')
        ax.axvline(0,color='gray',lw=.8);ax.grid(alpha=.2)
        ax.set_xlabel(r'$N^{-1/\nu}$, $\nu=0.67169$; infinite volume at 0' if D==3 else r'$N^{-2}(\log N)^{-1/10}$; infinite volume at 0')
        ax.set_ylabel(r'$\beta$ at $U=1/4$');ax.set_title(f'd={D}: '+('two-size provisional extrapolation' if len(selected)==2 else 'conditional finite-size extrapolation'))
        ax.legend(fontsize=8,loc='best')
        missing=[r['N'] for r in rows if r['D']==D and r['status']!='directly_bracketed']
        note='Two sizes cannot test fit quality.' if len(selected)==2 else f"Leading fit: χ²/dof={base['chi2']:.1f}/{base['df']}."
        if missing:note+=' Unbracketed N='+','.join(map(str,missing))+'.'
        ax.text(.5,-.24,note,ha='center',transform=ax.transAxes,fontsize=9)
    fig,axes=plt.subplots(1,2,figsize=(13,5.5))
    for ax,D in zip(axes,[3,4]):panel(ax,D)
    fig.suptitle(r'Binder-quarter calibration, $\epsilon/\pi=1$: directly bracketed crossings',fontsize=14)
    fig.text(.5,.015,'Vertical bars: conditional 95% MC intervals; size/interpolation sensitivity is separate. Black rings: diagnostic flags (few blocks in 4D).',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.10,1,.93])
    for extension in ['png','pdf','svg']:fig.savefig(folder/f'binder-quarter-size-extrapolations.{extension}',dpi=170)
    plt.close(fig)
    for D in [3,4]:
        fig,ax=plt.subplots(figsize=(7,5.5));panel(ax,D);fig.tight_layout(rect=[0,.08,1,.97])
        for extension in ['png','pdf','svg']:fig.savefig(folder/f'd{D}-binder-quarter-size-extrapolation.{extension}',dpi=170)
        plt.close(fig)
    print(json.dumps(dict(locations=rows,fits=fits),indent=2))


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='binder-quarter-calibration-2026-10-13')
    args=parser.parse_args();main(args.output)
