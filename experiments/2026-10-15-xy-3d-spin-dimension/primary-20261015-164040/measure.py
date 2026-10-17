"""Measure the integrated spin scaling dimension at fixed Binder U, from MC."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl');os.environ['MPLBACKEND']='Agg'
from pathlib import Path
import sys,json,hashlib,tempfile,argparse,shutil,csv
import numpy as np,h5py
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.stats import chi2
import matplotlib.pyplot as plt
OUT=Path(__file__).resolve().parent;CAMPAIGN=OUT.parent/'2026-10-13-xy-infinite-volume';sys.path.insert(0,str(CAMPAIGN))
import run
from analyze_calibration import helper
from threshold_calibration import monotone_values
DRAWS=2000

def binder(m):return 1-m[...,1]/(2*m[...,0]**2)

def summarize(arrays,ids,cut,factor):
    means=[];blocks=[];diagnostics=[];flags=[];individual=[]
    seed=int(hashlib.sha256(json.dumps([ids,cut,factor,1]).encode()).hexdigest()[:16],16);rng=np.random.default_rng(seed)
    for a,i in zip(arrays,ids):
        a=a[cut:];m2,m4,ma=a.T;influence=m4.mean()/m2.mean()**3*(m2-m2.mean())-(m4-m4.mean())/(2*m2.mean()**2)
        taus={k:helper.tau_int(v) for k,v in dict(m2=m2,m4=m4,mabs=ma,binder=influence).items()}
        width=factor*max(20,int(np.ceil(10*max([1.]+[float(t) for t in taus.values() if t is not None and np.isfinite(t)]))))
        nb=len(a)//width;diagnostics.append(dict(id=i,rows=len(a),width=width,blocks=nb,taus=taus))
        if nb<2:return dict(status='insufficient_blocks',diagnostics=diagnostics),None
        if nb<50:flags.append('fewer_than_50_blocks_per_start')
        edges=np.linspace(0,len(a),nb+1,dtype=int);sums=np.add.reduceat(a,edges[:-1],axis=0);length=np.diff(edges);ix=rng.integers(nb,size=(DRAWS,nb))
        bs=sums[ix].sum(1);bn=length[ix].sum(1);means.append((bs,bn));blocks.append(len(a));individual.append(dict(values=a.mean(0).tolist(),errors=(bs/bn[:,None]).std(0,ddof=1).tolist()))
        halves=[]
        for half in np.array_split(a,2):
            nbh=len(half)//width
            if nbh<2:flags.append('insufficient_half_blocks');halves.append(None);continue
            edgesh=np.linspace(0,len(half),nbh+1,dtype=int);bm=np.array([v.mean(0) for v in np.split(half,edgesh[1:-1])]);halves.append((half.mean(0),bm.std(0,ddof=1)/np.sqrt(nbh)))
        if all(x is not None for x in halves):
            z=(halves[1][0]-halves[0][0])/np.maximum(np.hypot(halves[1][1],halves[0][1]),1e-15)
            diagnostics[-1]['half_z']=z.tolist()
            if max(abs(z[[0,2]]))>4:flags.append('magnetization_half_shift')
    mean=np.concatenate([a[cut:] for a in arrays]).mean(0);B=sum(v[0] for v in means)/sum(v[1] for v in means)[:,None]
    vals=np.array([v['values'] for v in individual]);errs=np.array([v['errors'] for v in individual])
    for col,name in [(0,'m2'),(2,'mabs')]:
        w=1/np.maximum(errs[:,col],1e-15)**2;center=np.average(vals[:,col],weights=w)
        if len(ids)>1 and chi2.sf(np.sum(w*(vals[:,col]-center)**2),len(ids)-1)<.001:flags.append('start_disagreement_'+name)
    return dict(status='conditional_block_summary',mean=mean.tolist(),flags=sorted(set(flags)),diagnostics=diagnostics,individual=individual,binder=float(binder(mean))),B

def group_summary(points,cut,factor):
    ids=sorted(p['id'] for p in points);key=hashlib.sha256(json.dumps([ids,cut,factor,1]).encode()).hexdigest();cache=OUT/'cache';cache.mkdir(exist_ok=True)
    meta=cache/(key+'.json');npz=cache/(key+'.npz')
    if meta.exists():return json.loads(meta.read_text()),np.load(npz)['samples'] if npz.exists() else None
    arrays=[]
    for i in ids:
        with h5py.File(CAMPAIGN/'chains'/f'{i}.h5') as h:arrays.append(np.column_stack([h[k][:] for k in ['magnetization_squared','magnetization_fourth','magnetization_abs']]))
    summary,B=summarize(arrays,ids,cut,factor);summary['ids']=ids
    if B is not None:
        with tempfile.NamedTemporaryFile(dir=cache,suffix='.npz',delete=False) as f:np.savez_compressed(f,samples=B);f.flush();os.fsync(f.fileno());tmp=Path(f.name)
        os.replace(tmp,npz)
    run.atomic_json(meta,summary);return summary,B

def interpolate_record(x,means,B,target=.25,method='linear',window='local'):
    U=binder(means);UB=binder(B);weights=1/np.maximum(UB.std(0,ddof=1),1e-12)**2;project=monotone_values(U,weights)
    hits=np.flatnonzero((project[:-1]<=target)&(project[1:]>target))
    if len(hits)!=1:return dict(status='unbracketed',valid_fraction=0.),None
    k=int(hits[0]);indices=np.arange(len(x)) if window=='full' else np.arange(max(0,k-1),min(len(x),k+3))
    xx=x[indices];ww=weights[indices]
    def one(m):
        raw=binder(m);u=monotone_values(raw,ww);hits=np.flatnonzero((u[:-1]<=target)&(u[1:]>target))
        if len(hits)!=1:return [np.nan]*3
        j=int(hits[0]);t=(target-u[j])/(u[j+1]-u[j]);root=xx[j]+t*(xx[j+1]-xx[j])
        if method=='pchip':root=brentq(lambda b:float(PchipInterpolator(xx,u)(b)-target),xx[j],xx[j+1])
        if method=='moment_linear':
            if not (raw[j]<=target<raw[j+1]):return [np.nan]*3
            root=brentq(lambda b:float(binder(m[j]+(b-xx[j])/(xx[j+1]-xx[j])*(m[j+1]-m[j]))-target),xx[j],xx[j+1])
        if method=='pchip':values=PchipInterpolator(xx,m,axis=0)(root)
        else:values=m[j]+(root-xx[j])/(xx[j+1]-xx[j])*(m[j+1]-m[j])
        return [root,values[0],values[2]]
    center=np.array(one(means[indices]));draws=np.array([one(b[:,indices,:][0]) for b in B[:,None,:,:]])
    valid=np.isfinite(draws).all(1)&(draws[:,1:] >0).all(1);fraction=float(valid.mean());centralvalid=np.isfinite(center).all() and (center[1:]>0).all()
    record=dict(status='conditional_contour' if centralvalid and fraction>=.99 else 'unresolved',location=center.tolist() if centralvalid else None,MC95=np.quantile(draws[valid],[.025,.975],axis=0).T.tolist() if valid.sum()>1 else None,valid_fraction=fraction,selected_indices=indices.tolist(),beta_window=xx.tolist(),binder_values=U.tolist(),monotone_max_shift=float(np.max(abs(project-U))),node_count=len(indices),moment_ratio_at_interpolated_root=None)
    return record,draws

def fit_dimension(N,y,B,power,minimum=8,omega=None):
    keep=np.array(N)>=minimum;N=np.asarray(N)[keep];y=np.asarray(y)[keep];B=np.asarray(B)[:,keep];valid=np.isfinite(B).all(1)&(B>0).all(1)
    X=np.column_stack([np.ones(len(N)),-power*np.log(N)]+([N.astype(float)**(-omega)] if omega is not None else []))
    if len(N)<X.shape[1]:return dict(status='insufficient_sizes',sizes=N.tolist())
    if valid.mean()<.99:return dict(status='unresolved_joint_bootstrap',sizes=N.tolist(),valid_fraction=float(valid.mean()))
    L=np.log(B[valid]);C=np.atleast_2d(np.cov(L,rowvar=False));e,Q=np.linalg.eigh(C)
    if min(e)<=max(e)*1e-10:return dict(status='rank_unresolved',sizes=N.tolist())
    W=(Q/np.sqrt(e))@Q.T;Z=W@X;scale=np.linalg.norm(Z,axis=0)
    if np.linalg.matrix_rank(Z/scale)<X.shape[1]:return dict(status='design_rank_unresolved',sizes=N.tolist())
    T=np.linalg.pinv(Z/scale)@W;coef=T@np.log(y)/scale;bs=L@T.T/scale;r=W@(np.log(y)-X@coef);stat=float(r@r);dof=len(N)-X.shape[1];prob=float(chi2.sf(stat,dof)) if dof else None
    return dict(status='zero_dof_diagnostic' if dof==0 else 'rejected_power_law' if prob<.01 else 'conditional_dimension_fit',sizes=N.tolist(),Delta=float(coef[1]),MC95=np.quantile(bs[:,1],[.025,.975]).tolist(),chi2=stat,dof=dof,p=prob,valid_fraction=float(valid.mean()),coefficients=coef.tolist(),minimum=minimum,omega=omega,delta_samples=bs[:,1])

def analyze(tag='primary',cut=1001,factor=1,target=.25,starts='all',window='local'):
    folder=OUT/(tag+'-'+run.now()[:19].replace('-','').replace(':','').replace('T','-'));folder.mkdir()
    shutil.copyfile(__file__,folder/'measure.py');groups={};source=[]
    for f in (CAMPAIGN/'jobs').glob('*.json'):
        j=json.loads(f.read_text());p=j['point']
        if p['D']!=3 or p['epsilon_over_pi'] not in [.82,.84,.86,.88,.90,.92,.94,.96,.98,1.] or p.get('campaign_role')=='zero_beta_endpoint' or j['state']!='complete' or not j.get('audit',{}).get('passed'):continue
        if starts!='all' and p['replica']%2!=int(starts=='odd'):continue
        groups.setdefault((p['epsilon_over_pi'],p['N'],p['beta']),[]).append(p);source.append(dict(id=p['id'],audit=j['audit']))
    data=[];boots={};curves=[]
    for key,points in sorted(groups.items()):
        summary,B=group_summary(points,cut,factor);data.append(dict(q=key[0],N=key[1],beta=key[2],summary=summary))
        if B is not None:boots[key]=B
    for q,N in sorted({k[:2] for k in groups}):
        selected=sorted([g for g in data if (g['q'],g['N'])==(q,N) and (q,N,g['beta']) in boots],key=lambda g:g['beta'])
        if len(selected)<2:continue
        x=np.array([g['beta'] for g in selected]);means=np.array([g['summary']['mean'] for g in selected]);B=np.stack([boots[q,N,g['beta']] for g in selected],axis=1)
        for method in ['linear','pchip','moment_linear']:
            r,rb=interpolate_record(x,means,B,target,method,window);r.update(q=q,N=N,method=method,flags=sorted({f for i in r.get('selected_indices',[]) for f in selected[i]['summary']['flags']}));curves.append(r)
            if rb is not None:np.savez_compressed(folder/f'root_q{q:.2f}_N{N}_{method}.npz',samples=rb)
    fits=[]
    for q in sorted({r['q'] for r in curves}):
        for method in ['linear','pchip','moment_linear']:
            rs=sorted([r for r in curves if r['q']==q and r['method']==method and r['status']=='conditional_contour'],key=lambda r:r['N'])
            if not rs:continue
            B=np.stack([np.load(folder/f'root_q{q:.2f}_N{r["N"]}_{method}.npz')['samples'] for r in rs],axis=1)
            for column,moment,power in [(1,'m2',2),(2,'mabs',1)]:
                for minimum,omega in [(6,None),(8,None),(12,None),(6,.5),(6,.8),(6,1.)]:
                    f=fit_dimension([r['N'] for r in rs],[r['location'][column] for r in rs],B[:,:,column],power,minimum,omega)
                    if 'delta_samples' in f:np.savez_compressed(folder/f'fit{len(fits)}.npz',samples=f.pop('delta_samples'))
                    fits.append(dict(q=q,method=method,moment=moment,minimum=minimum,omega=omega,omitted_sizes=[r['N'] for r in curves if r['q']==q and r['method']==method and r['N']>=minimum and r['status']!='conditional_contour'],result=f))
    run.atomic_json(folder/'analysis.json',dict(created=run.now(),settings=dict(tag=tag,cut=cut,factor=factor,target=target,starts=starts,window=window),source_jobs=source,groups=data,contours=curves,fits=fits,scope='Existing MC only. Fixed Binder contour measures an effective integrated spin scaling exponent if it approaches the critical scaling regime and 2Delta<3. No nu or Delta fixed. Conditional block MC uncertainty; interpolation, corrections/contact backgrounds and finite sizes are not certified total confidence coverage.'))
    run.atomic_json(OUT/f'latest-{tag}.json',dict(folder=folder.name));print(folder,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--tag',default='primary');p.add_argument('--cut',type=int,default=1001);p.add_argument('--factor',type=int,default=1);p.add_argument('--target',type=float,default=.25);p.add_argument('--starts',choices=['all','even','odd'],default='all');p.add_argument('--window',choices=['local','full'],default='local');a=p.parse_args();analyze(a.tag,a.cut,a.factor,a.target,a.starts,a.window)
