"""Predeclared local Taylor fits to noisy common-coupling xi/L scans."""
import argparse
import json
import shutil
from pathlib import Path
import numpy as np
from scipy.stats import chi2
import common_xi_quotients as x

WINDOWS={'inner':(.800625,.801000),'broad':(.800500,.801500)}


def response_fit(q,means,B,N,lo,hi,degree):
    """Marginal WLS for R and m2; shared draws carry their full covariance."""
    t=(q-(lo+hi)/2)/((hi-lo)/2);X=np.polynomial.polynomial.polyvander(t,degree)
    central=np.column_stack([x.ratio(means,N),means[:,0]])
    samples=np.stack([x.ratio(B,N),B[:,:,0]],axis=-1)
    if not np.isfinite(samples).all():return None,dict(status='invalid_input_ratio')
    coefs=[];draws=[];diagnostics=[]
    for j,name in enumerate(['xi_over_L','m2']):
        sd=samples[:,:,j].std(0,ddof=1)
        if not (sd>0).all():return None,dict(status='zero_input_variance')
        Z=X/sd[:,None]
        if np.linalg.matrix_rank(Z)<degree+1:return None,dict(status='rank_unresolved')
        T=np.linalg.pinv(Z)/sd[None,:];coef=T@central[:,j];bs=samples[:,:,j]@T.T
        residual=(central[:,j]-X@coef)/sd;stat=float(residual@residual);dof=len(q)-degree-1
        coefs.append(coef);draws.append(bs)
        diagnostics.append(dict(observable=name,chi2=stat,dof=dof,p=float(chi2.sf(stat,dof)) if dof else None,coefficients=coef.tolist(),SD_at_knots=sd.tolist()))
    return (np.array(coefs),np.stack(draws,axis=1)),dict(status='local_response_fit',N=N,knots=q.tolist(),degree=degree,response_diagnostics=diagnostics,scope='Marginal descriptive WLS residual tests; not a model-selection or total-coverage guarantee.')


def quotient(c1,c2,N1,N2,lo,hi):
    delta=np.polynomial.Polynomial(c1[0]-c2[0]);coef=delta.coef
    if np.max(abs(coef))<1e-14:return None,'coincident_responses'
    roots=[float(z.real) for z in delta.roots() if abs(z.imag)<1e-8 and -1-1e-9<=z.real<=1+1e-9]
    if len(roots)==0:return None,'no_root_in_window'
    if len(roots)>1:return None,'multiple_roots_in_window'
    t=np.clip(roots[0],-1,1);r1,m1=[np.polynomial.polynomial.polyval(t,v) for v in c1];r2,m2=[np.polynomial.polynomial.polyval(t,v) for v in c2]
    if not np.isfinite([r1,r2,m1,m2]).all() or min(r1,r2,m1,m2)<=0:return None,'nonphysical_root'
    q=(lo+hi)/2+t*(hi-lo)/2;D=-np.log(m2/m1)/(2*np.log(N2/N1))
    return np.array([q,r1,m1,m2,D]),'one_physical_root'


def pair_fit(q,mean1,B1,N1,mean2,B2,N2,lo,hi,degree):
    if len(q)<degree+2:return dict(status='insufficient_knots',knots=q.tolist()),None
    first,d1=response_fit(q,mean1,B1,N1,lo,hi,degree);second,d2=response_fit(q,mean2,B2,N2,lo,hi,degree)
    if first is None or second is None:return dict(status='response_fit_unresolved',responses=[d1,d2]),None
    center,reason=quotient(first[0],second[0],N1,N2,lo,hi)
    samples=np.full((len(B1),5),np.nan);counts={}
    for i,(a,b) in enumerate(zip(first[1],second[1])):
        value,why=quotient(a,b,N1,N2,lo,hi);counts[why]=counts.get(why,0)+1
        if value is not None:samples[i]=value
    valid=np.isfinite(samples).all(1);fraction=float(valid.mean())
    record=dict(status='conditional_local_quotient' if center is not None and fraction>=.99 else 'unresolved_local_quotient',central_root_status=reason,valid_fraction=fraction,bootstrap_root_events=counts,responses=[d1,d2],knots=q.tolist())
    if center is not None:record.update(epsilon_over_pi=float(center[0]),xi_over_L=float(center[1]),Delta_eff=float(center[4]))
    if valid.sum()>1:record.update(SD=float(samples[valid,4].std(ddof=1)),MC95=np.quantile(samples[valid],[.025,.975],axis=0).T.tolist())
    return record,samples


def run(status_path):
    status=json.loads(Path(status_path).read_text());folder=x.m.OUT/('local-xi-response-'+x.m.run.now()[:19].replace('-','').replace(':','').replace('T','-'));folder.mkdir()
    shutil.copyfile(__file__,folder/'local_xi_response.py');records=[];sources={}
    for variant,info in status['snapshots'].items():
        path=Path(info['path']);a=json.loads(path.read_text());sources[variant]=dict(path=str(path),sha256=x.m.run.sha256(path));point={s['id']:s['point'] for s in a['source_jobs']};curves={}
        for g in a['groups']:
            if 'mean' not in g['summary'] or not .80049<=g['q']<=.801501:continue
            s,B=x.summary([point[i] for i in g['summary']['ids']],a['settings']['cut'],a['settings']['factor'])
            if B is not None:curves.setdefault(g['N'],[]).append((g['q'],s,B))
        for N in sorted(curves):
            if 2*N not in curves:continue
            for window,(lo,hi) in WINDOWS.items():
                matched=[]
                for v in curves[N]:
                    if not lo-1e-14<=v[0]<=hi+1e-14:continue
                    other=[w for w in curves[2*N] if abs(w[0]-v[0])<=8*max(abs(np.spacing(w[0])),abs(np.spacing(v[0])))]
                    if len(other)==1:matched.append((v,other[0]))
                q=np.array([v[0][0] for v in matched])
                if not matched:continue
                m1=np.array([v[0][1]['mean'] for v in matched]);m2=np.array([v[1][1]['mean'] for v in matched]);B1=np.stack([v[0][2] for v in matched],axis=1);B2=np.stack([v[1][2] for v in matched],axis=1)
                flags=sorted({f for v in matched for w in v for f in w[1]['flags']})
                for degree in [1,2]:
                    r,B=pair_fit(q,m1,B1,N,m2,B2,2*N,lo,hi,degree);r.update(variant=variant,N=N,N2=2*N,window=window,range=[lo,hi],degree=degree,input_flags=flags)
                    records.append(r)
                    if B is not None:np.savez_compressed(folder/f'{variant}-{N}-{2*N}-{window}-degree{degree}.npz',samples=B)
    report=dict(created=x.m.run.now(),method='Local polynomial regression of ensemble xi/L and m2, a distinct estimator from raw-moment knot interpolation. Both observables share original coordinate/start bootstrap draws; same coordinate draws reused across pairs. Only pairwise matching physical knots are used.',windows=WINDOWS,records=records,sources=sources,scope='All predeclared degree/window/variant outcomes retained, including tail/multiplicity failures, fit residual diagnostics and start/cut differences. Good residual p or reduced MC SD is not accuracy or model-selection evidence. No imposed Delta/critical R, no numerical sigma_sys, no new MC.')
    x.m.run.atomic_json(folder/'analysis.json',report)
    x.m.run.atomic_json(x.m.OUT/'latest-local-xi-response.json',dict(folder=folder.name,status='pending_independent_review'))
    print(folder,flush=True)
    for r in records:print(json.dumps(r),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--status',default=str(x.m.OUT/'common-xi-refinement-status.json'));a=p.parse_args();run(a.status)
