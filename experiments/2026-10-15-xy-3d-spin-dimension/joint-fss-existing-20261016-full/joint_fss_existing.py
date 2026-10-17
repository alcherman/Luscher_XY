"""Source-pinned joint finite-size fits, existing 3D beta=0 MC only.

No simulation calls. Joint within-coordinate covariance, profiled critical
coupling, explicit WF-informed exponents, and descriptive held-out predictions.
"""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-joint-fss-mpl')
os.environ['MPLBACKEND']='Agg'
from pathlib import Path
import json,hashlib,time,argparse,shutil
import numpy as np
from scipy.optimize import minimize_scalar,least_squares
from scipy.stats import chi2
import common_xi_quotients as cq

OUT=Path(__file__).resolve().parent
QREF=.80075
QSCALE=1000.
LREF=32.
QBOUNDS=(-.5,1.) # epsilon/pi .80025 through .80175, not fixed to a crossing.
NU=.67169
OMEGA=.789
WINDOWS={'inner':(.800625,.801000),'broad':(.800500,.801500)}

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,obj):
    def clean(x):
        if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [clean(v) for v in x]
        if isinstance(x,np.ndarray):return clean(x.tolist())
        if isinstance(x,(np.floating,float)):return float(x) if np.isfinite(x) else None
        if isinstance(x,np.integer):return int(x)
        if isinstance(x,np.bool_):return bool(x)
        return x
    cq.m.run.atomic_json(path,clean(obj))

def blocks_whitener(B):
    """B has shape(draw, coordinate, observable). Never invent cross-coordinate covariance."""
    cov=np.stack([np.cov(B[:,i,:],rowvar=False,ddof=1) for i in range(B.shape[1])])
    ev=np.linalg.eigvalsh(cov)
    if not np.isfinite(cov).all() or np.any(ev<=0):raise ValueError('Singular/nonpositive coordinate covariance')
    W=np.linalg.inv(np.linalg.cholesky(cov))
    return W,cov

def whiten(W,a):return np.einsum('nij,njk->nik',W,a) if a.ndim==3 else np.einsum('nij,nj->ni',W,a)

def design(L,q,Q,model='leading',nu=NU,omega=OMEGA):
    ell=np.asarray(L,dtype=float)/LREF;q=np.asarray(q);x=(QSCALE*(q-QREF)-Q)*ell**(1/nu);u=ell**(-omega);dx=-ell**(1/nu)
    r=[np.ones(len(L)),x,x*x,u];rn=['Rstar','R_x','R_x2','R_correction']
    z=[np.ones(len(L)),-2*np.log(ell),x,x*x,u];zn=['log_amplitude','Delta','m_x','m_x2','m_correction']
    dr=[np.zeros(len(L)),dx,2*x*dx,np.zeros(len(L))];dz=[np.zeros(len(L)),np.zeros(len(L)),dx,2*x*dx,np.zeros(len(L))]
    if model=='mixed':
        r.append(x*u);rn.append('R_mixed');dr.append(dx*u)
        z.append(x*u);zn.append('m_mixed');dz.append(dx*u)
    elif model=='analytic':
        z.append(QSCALE*(q-QREF)-Q);zn.append('m_analytic_amplitude');dz.append(-np.ones(len(L)))
    elif model=='second':
        r.append(u*u);rn.append('R_second_correction');dr.append(np.zeros(len(L)))
        z.append(u*u);zn.append('m_second_correction');dz.append(np.zeros(len(L)))
    elif model=='cubic':
        r.append(x**3);rn.append('R_x3');dr.append(3*x*x*dx)
        z.append(x**3);zn.append('m_x3');dz.append(3*x*x*dx)
    elif model!='leading':raise ValueError(model)
    X=np.zeros((len(L),2,len(r)+len(z)));dX=np.zeros_like(X)
    X[:,0,:len(r)]=np.array(r).T;X[:,1,len(r):]=np.array(z).T
    dX[:,0,:len(r)]=np.array(dr).T;dX[:,1,len(r):]=np.array(dz).T
    return X,dX,rn+zn

def linear_solution(X,y):
    A=X.reshape(-1,X.shape[-1]);v=y.ravel();scale=np.linalg.norm(A,axis=0)
    if np.any(scale==0):return None,np.inf,0,np.inf
    z,_,rank,sv=np.linalg.lstsq(A/scale,v,rcond=1e-11)
    coef=z/scale;res=A@coef-v
    return coef,float(res@res),int(rank),float(sv[0]/sv[-1]) if sv[-1]>0 else np.inf

def solve_profile(L,q,y,W,model='leading',nu=NU,omega=OMEGA,grid=11,fixed_delta=None,full=False):
    cache={}
    def objective(Q):
        X,_,names=design(L,q,Q,model,nu,omega);A=whiten(W,X);target=whiten(W,y)
        if fixed_delta is not None:
            j=names.index('Delta');target=target-A[:,:,j]*fixed_delta;A=np.delete(A,j,axis=2)
        c,s,rank,condition=linear_solution(A,target)
        if fixed_delta is not None and c is not None:c=np.insert(c,j,fixed_delta)
        cache[float(Q)]=(c,s,rank,condition)
        return s
    qs=np.linspace(*QBOUNDS,grid);scores=np.array([objective(v) for v in qs]);candidates=[(scores[0],qs[0]),(scores[-1],qs[-1])]
    minima=[]
    for i in range(1,len(qs)-1):
        if scores[i]<=scores[i-1] and scores[i]<=scores[i+1]:
            r=minimize_scalar(objective,bounds=(qs[i-1],qs[i+1]),method='bounded',options={'xatol':1e-8})
            candidates.append((float(r.fun),float(r.x)));minima.append((float(r.fun),float(r.x)))
    s,Q=min(candidates);objective(Q);coef,score,rank,condition=cache[Q]
    X,dX,names=design(L,q,Q,model,nu,omega);pred=np.einsum('nip,p->ni',X,coef) if coef is not None else np.full_like(y,np.nan)
    res=dict(Q=Q,qc=QREF+Q/QSCALE,coefficients=coef,names=names,Delta=coef[names.index('Delta')] if coef is not None else np.nan,chi2=score,rank=rank,linear_parameters=len(names)-(fixed_delta is not None),condition=condition,prediction=pred,boundary=abs(Q-QBOUNDS[0])<1e-6 or abs(Q-QBOUNDS[1])<1e-6,multiple_near_minima=sum(v[0]-s<4 for v in minima)>1)
    if full:
        J=np.concatenate([np.einsum('nip,p->ni',dX,coef)[:,:,None],X],axis=2)
        if fixed_delta is not None:J=np.delete(J,1+names.index('Delta'),axis=2)
        A=whiten(W,J).reshape(-1,J.shape[-1]);scale=np.linalg.norm(A,axis=0);sv=np.linalg.svd(A/np.maximum(scale,1e-300),compute_uv=False)
        res.update(jacobian_rank=int(np.sum(sv>sv[0]*1e-10)),jacobian_condition=float(sv[0]/max(sv[-1],1e-300)),jacobian_singular_values=sv,profile_grid_Q=qs,profile_grid_chi2=scores,minima=minima)
    return res

def raw_prediction(L,q,params,model='leading',nu=NU,omega=OMEGA,background=True):
    Q=params[0];X,dX,names=design(L,q,Q,model,nu,omega);p=len(names);coef=params[1:1+p]
    z=np.einsum('nip,p->ni',X,coef);dz=np.einsum('nip,p->ni',dX,coef)
    singular=np.exp(np.clip(z[:,1],-300,300));raw_basis=.033*(np.asarray(L)/LREF)**(-3)
    pred=z.copy();pred[:,1]=singular+(params[-1]*raw_basis if background else 0)
    J=np.zeros((len(L),2,len(params)));J[:,:,0]=dz;J[:,:,1:1+p]=X;J[:,1,:1+p]*=(singular*((z[:,1]>-300)&(z[:,1]<300)))[:,None]
    if background:J[:,1,-1]=raw_basis
    return pred,J

def solve_raw(L,q,y,W,initial,model='leading',nu=NU,omega=OMEGA,background=True):
    names=initial['names'];p0=np.array(initial['raw_parameters']) if 'raw_parameters' in initial else np.r_[initial['Q'],initial['coefficients'],[0.] if background else []]
    low=np.r_[QBOUNDS[0],np.full(len(names),-100.),[-100.] if background else []];high=np.r_[QBOUNDS[1],np.full(len(names),100.),[100.] if background else []]
    def fun(p):return whiten(W,raw_prediction(L,q,p,model,nu,omega,background)[0]-y).ravel()
    def jac(p):return whiten(W,raw_prediction(L,q,p,model,nu,omega,background)[1]).reshape(-1,len(p))
    fit=least_squares(fun,np.clip(p0,low+1e-8,high-1e-8),jac=jac,bounds=(low,high),x_scale='jac',ftol=1e-10,xtol=1e-10,gtol=1e-8,max_nfev=400)
    pred,J=raw_prediction(L,q,fit.x,model,nu,omega,background);A=whiten(W,J).reshape(-1,len(fit.x));scale=np.linalg.norm(A,axis=0);sv=np.linalg.svd(A/np.maximum(scale,1e-300),compute_uv=False)
    return dict(Q=float(fit.x[0]),qc=QREF+fit.x[0]/QSCALE,Delta=float(fit.x[1+names.index('Delta')]),coefficients=fit.x[1:1+len(names)],names=names,raw_parameters=fit.x,background_coefficient=float(fit.x[-1]) if background else 0.,prediction=pred,chi2=float(fit.fun@fit.fun),success=bool(fit.success),nfev=int(fit.nfev),boundary=bool(np.any(abs(fit.x-low)<1e-6)|np.any(abs(fit.x-high)<1e-6)),rank=int(np.sum(sv>sv[0]*1e-10)),condition=float(sv[0]/max(sv[-1],1e-300)),multiple_near_minima=None,exponent_clipped=bool(np.any(abs(np.einsum('nip,p->ni',design(L,q,fit.x[0],model,nu,omega)[0],fit.x[1:1+len(names)])[:,1])>=300)))

def load_bank(analysis):
    a=json.loads(Path(analysis).read_text());point={s['id']:s['point'] for s in a['source_jobs']};rows=[];boot=[];cache_records=[]
    for g in a['groups']:
        if g['N']<12 or not .80049<=g['q']<=.801501 or 'mean' not in g['summary']:continue
        points=[point[i] for i in g['summary']['ids']]
        # Require an already-published joint cache. No new raw-data reads or writes.
        key=hashlib.sha256(json.dumps([sorted(p['id'] for p in points),a['settings']['cut'],a['settings']['factor'],cq.VERSION,'common_xi']).encode()).hexdigest()
        meta=OUT/'xi-cache'/(key+'.json');data=OUT/'xi-cache'/(key+'.npz')
        if not meta.exists() or not data.exists():raise RuntimeError('Required existing joint cache absent: '+key)
        s=json.loads(meta.read_text());B=np.load(data)['samples'];mean=np.array(s['mean']);R=float(cq.ratio(mean,g['N']));RB=cq.ratio(B,g['N'])
        values=np.array([R,np.log(mean[0])]);samples=np.column_stack([RB,np.log(B[:,0])])
        if not np.isfinite(samples).all():raise RuntimeError('Nonfinite joint cache: '+key)
        rows.append(dict(L=g['N'],q=g['q'],mean=values,ids=s['ids'],flags=s.get('flags',[]),diagnostics=s['diagnostics'],individual=s.get('individual')));boot.append(samples)
        cache_records.append(dict(metadata=str(meta),metadata_sha256=sha(meta),samples=str(data),samples_sha256=sha(data)))
    return rows,np.stack(boot,axis=1),dict(analysis=str(analysis),sha256=sha(analysis),settings=a['settings'],source_jobs=a['source_jobs'],cache_records=cache_records)


def configs():
    # Declared before fitting; no model is chosen by its p value.
    out=[]
    for window in ['inner','broad']:
        for minimum in [12,16,24]:
            for model in ['leading','mixed']:
                out.append(dict(name=f'{window}-min{minimum}-{model}',window=window,minimum=minimum,maximum=64,model=model,nu=NU,omega=OMEGA,holdout=None))
    for model in ['analytic','second','cubic']:
        out.append(dict(name='inner-min12-'+model,window='inner',minimum=12,maximum=64,model=model,nu=NU,omega=OMEGA,holdout=None))
    for omega in [.6,1.]:out.append(dict(name=f'inner-min12-omega{omega}',window='inner',minimum=12,maximum=64,model='leading',nu=NU,omega=omega,holdout=None))
    for nu in [.63,.71]:out.append(dict(name=f'inner-min12-nu{nu}',window='inner',minimum=12,maximum=64,model='leading',nu=nu,omega=OMEGA,holdout=None))
    for window in ['inner','broad']:
        for maximum in [32,48]:
            for model in ['leading','mixed']:
                out.append(dict(name=f'{window}-train12to{maximum}-{model}',window=window,minimum=12,maximum=maximum,model=model,nu=NU,omega=OMEGA,holdout='larger-than-training'))
    for background in [False,True]:out.append(dict(name='inner-min12-raw-'+('background' if background else 'zero-background'),window='inner',minimum=12,maximum=64,model='leading',nu=NU,omega=OMEGA,holdout=None,raw=True,background=background))
    return out


def select(rows,c,heldout=False):
    lo,hi=WINDOWS[c['window']]
    return np.array([i for i,r in enumerate(rows) if lo-1e-14<=r['q']<=hi+1e-14 and (r['L']>c['maximum'] if heldout else c['minimum']<=r['L']<=c['maximum'])])


def predictive_summary(y,pred,residual_draws):
    r=(y-pred).ravel();C=np.cov(residual_draws.reshape(len(residual_draws),-1),rowvar=False,ddof=1);ev=np.linalg.eigvalsh(C)
    if ev[0]<=0:return dict(status='singular_predictive_covariance')
    stat=float(r@np.linalg.solve(C,r));dof=len(r);sd=np.sqrt(np.diag(C)).reshape(y.shape)
    return dict(status='descriptive_heldout_comparison',observed_minus_predicted=y-pred,predictive_SD=sd,chi2=stat,dof=dof,nominal_p=float(chi2.sf(stat,dof)),max_abs_R=float(abs((y-pred)[:,0]).max()),max_abs_log_m2=float(abs((y-pred)[:,1]).max()),max_abs_fractional_m2=float(abs(np.expm1((y-pred)[:,1])).max()),residual_covariance=C,covariance_eigenvalues=ev,scope='Heldout measurement and training-fit uncertainty, shared bootstrap across all prediction coordinates; nominal Gaussian diagnostic, not confirmatory model selection.')


def fit_one(c,rows,B,draws=2000,central_only=False):
    ix=select(rows,c);L=np.array([rows[i]['L'] for i in ix]);q=np.array([rows[i]['q'] for i in ix]);y=np.array([rows[i]['mean'] for i in ix]);bs=B[:draws,ix,:].copy();W,cov=blocks_whitener(bs)
    raw=c.get('raw',False)
    if raw and c['holdout']:raise ValueError('Raw heldout configuration not supported')
    args=dict(model=c['model'],nu=c['nu'],omega=c['omega'])
    init=solve_profile(L,q,y,W,**args,grid=81,full=True)
    if raw:
        y[:,1]=np.exp(y[:,1]);bs[:,:,1]=np.exp(bs[:,:,1]);W,cov=blocks_whitener(bs)
        starts=[]
        for Q0 in [init['Q'],-.25,.25,.75]:
            for delta0 in [.45,.55,.65]:
                trial=dict(init);trial['Q']=Q0;trial['coefficients']=init['coefficients'].copy();trial['coefficients'][init['names'].index('Delta')]=delta0
                starts.append(solve_raw(L,q,y,W,trial,**args,background=c['background']))
        center=min(starts,key=lambda v:v['chi2'] if v['success'] else np.inf)
        center['multistart_objectives']=[v['chi2'] for v in starts]
        center['multistart_Delta']=[v['Delta'] for v in starts]
    else:center=init
    count=len(center['names'])+1+int(raw and c['background']);dof=2*len(ix)-count
    center.update(dof=dof,p=float(chi2.sf(center['chi2'],dof)) if dof>0 else None)
    result=dict(config=c,groups=len(ix),sizes=sorted(set(L.tolist())),coordinate_indices=ix,source_ids=sorted({s for i in ix for s in rows[i]['ids']}),input_flags=sorted({f for i in ix for f in rows[i]['flags']}),central=center,covariance_blocks=cov,bootstrap_draws=0)
    if central_only:return result,None
    parameters=[];prediction_draws=[];boundary=0;failed=0;multiple=0;failure_reasons={}
    hi=select(rows,c,True) if c['holdout'] else np.array([],dtype=int)
    LH=np.array([rows[i]['L'] for i in hi]);qH=np.array([rows[i]['q'] for i in hi]);yh=np.array([rows[i]['mean'] for i in hi]) if len(hi) else None
    for j,v in enumerate(bs):
        try:
            fit=solve_raw(L,q,v,W,center,**args,background=c['background']) if raw else solve_profile(L,q,v,W,**args,grid=11,full=True)
            joint_count=len(fit['names'])+1+int(raw and c['background'])
            reasons=[]
            if not np.isfinite([fit['Q'],fit['Delta']]).all():reasons.append('nonfinite_parameters')
            if not fit.get('success',True):reasons.append('optimizer_failure')
            if fit.get('exponent_clipped',False):reasons.append('exponent_clipped')
            if fit.get('jacobian_rank',fit['rank'])!=joint_count:reasons.append('joint_rank_failure')
            if fit.get('jacobian_condition',fit['condition'])>1e8:reasons.append('joint_condition_failure')
            if reasons:
                failed+=1;parameters.append(np.full(3,np.nan));prediction_draws.append(None)
                for reason in reasons:failure_reasons[reason]=failure_reasons.get(reason,0)+1
                continue
            parameters.append([fit['qc'],fit['Delta'],fit.get('background_coefficient',0.)]);boundary+=fit['boundary'];multiple+=bool(fit.get('multiple_near_minima'))
            if len(hi):prediction_draws.append(np.einsum('nip,p->ni',design(LH,qH,fit['Q'],**args)[0],fit['coefficients']))
        except (ValueError,np.linalg.LinAlgError,RuntimeError):failed+=1;parameters.append(np.full(3,np.nan));prediction_draws.append(None)
    params=np.array(parameters);valid=np.isfinite(params).all(1);result.update(bootstrap_draws=len(params),successful_draws=int(valid.sum()),failed_draws=failed,bootstrap_failure_reasons=failure_reasons,boundary_draws=boundary,multiple_minima_draws=multiple,bootstrap_valid_fraction=float(valid.mean()))
    if valid.sum()>1:
        result.update(parameter_SD=np.std(params[valid],axis=0,ddof=1),parameter_MC95=np.quantile(params[valid],[.025,.975],axis=0).T,parameter_covariance=np.cov(params[valid],rowvar=False),q_Delta_correlation=float(np.corrcoef(params[valid,:2].T)[0,1]))
    flags=[]
    if not center.get('success',True):flags.append('central_optimizer_failure')
    if center.get('exponent_clipped',False):flags.append('exponent_clipped')
    if dof<=0:flags.append('zero_dof_diagnostic')
    if center['p'] is not None and center['p']<.01:flags.append('rejected_joint_response')
    if center['boundary'] or boundary/max(1,len(params))>.01:flags.append('profile_boundary_unresolved')
    if valid.mean()<.99:flags.append('bootstrap_failure_unresolved')
    if center.get('jacobian_rank',center['rank'])<count:flags.append('joint_parameter_rank_failure')
    if center.get('jacobian_condition',center['condition'])>1e8:flags.append('weak_parameter_identification')
    if center.get('multiple_near_minima') or multiple/max(1,len(params))>.01:flags.append('multiple_profile_minima')
    qc=center['qc'];lo,upper=WINDOWS[c['window']]
    if not lo<=qc<=upper:flags.append('critical_coupling_outside_fitted_window')
    result['flags']=flags;result['status']='diagnostic_not_accepted' if flags else 'conditional_joint_fit'
    if len(hi) and valid.mean()>=.99:
        ph=np.einsum('nip,p->ni',design(LH,qH,center['Q'],**args)[0],center['coefficients']);pd=np.stack([prediction_draws[i] for i in np.flatnonzero(valid)]);res=B[:draws,hi,:][valid]-pd
        result['heldout']=predictive_summary(yh,ph,res);result['heldout'].update(coordinate_indices=hi,sizes=sorted(set(LH.tolist())),predicted=ph,observed=yh)
    return result,params


def main(args):
    primary=OUT/'common-xi-refined-primary-20261015-213346/analysis.json'
    rows,B,provenance=load_bank(primary)
    plans=configs();folder=Path(args.output) if args.output else OUT/('joint-fss-existing-'+time.strftime('%Y%m%d-%H%M%S'));folder.mkdir(exist_ok=False)
    write(folder/'design.json',dict(created=cq.m.run.now(),implementation_sha256=sha(__file__),source_analysis_sha256=sha(primary),scope='Existing data only; no simulation, archive, kernel or old-analysis modifications.',models=plans,windows=WINDOWS,q_profile_bounds=QBOUNDS,q_reference=QREF,q_scale=QSCALE,L_reference=LREF,draws=args.draws,primary_assumptions=dict(nu=NU,omega=OMEGA,WF_informed=True),observable_covariance='Independent physical coordinates, full 2x2 joint covariance within each coordinate; fixed central whitening; no cross-coordinate empirical covariance.',correction_models='Leading L^-omega; mixed xL^-omega in both observables; analytic thermal amplitude drift in logm2; second correction L^-2omega; cubic thermal; rawm2 with and without additive B L^-3. Stress exponents are not confidence ranges.',heldout='Predeclared train max48->test64, and train max32->test48/64. Heldout data not used in each fit; models are not chosen by best predictive p.',uncertainty='Conditional block-bootstrap MC intervals. No fit spread or successful prediction alone calibrates sigma_sys.'))
    write(folder/'sources.json',provenance);np.savez_compressed(folder/'input-bank.npz',L=np.array([r['L'] for r in rows]),q=np.array([r['q'] for r in rows]),means=np.array([r['mean'] for r in rows]),samples=B);write(folder/'input-groups.json',rows)
    shutil.copyfile(__file__,folder/'joint_fss_existing.py');results=[]
    for c in plans:
        path=folder/(c['name']+'.json')
        start=time.monotonic();result,params=fit_one(c,rows,B,args.draws,args.central_only);result['runtime_seconds']=time.monotonic()-start;write(path,result)
        if params is not None:np.savez_compressed(folder/(c['name']+'-draws.npz'),parameters=params)
        results.append(result);write(folder/'status.json',dict(state='running',completed=len(results),total=len(plans),latest=c['name']))
        print(c['name'],result['central']['qc'],result['central']['Delta'],result['central']['p'],result.get('parameter_SD'),result.get('flags'),result.get('heldout',{}).get('nominal_p'),flush=True)
    write(folder/'analysis.json',dict(results=results,scope='Predeclared existing-data joint FSS model/prediction study, conditional on covariance, scaling laws and source history diagnostics. Not calibrated systematic error.',provenance_sha256=sha(folder/'sources.json')))
    write(folder/'status.json',dict(state='central_diagnostics_complete' if args.central_only else 'fits_complete_pending_independent_review',completed=len(results),total=len(plans)))
    print(folder,flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output');ap.add_argument('--draws',type=int,default=2000);ap.add_argument('--central-only',action='store_true');main(ap.parse_args())
