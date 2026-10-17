"""Joint moment bootstrap for eta_eff and correlation-length inference."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import h5py,numpy as np
from scipy.stats import chi2
from analyze_calibration import helper
import run

def read(point):
    with h5py.File(run.OUT/'chains'/f"{point['id']}.h5") as f:
        n=len(f['magnetization_squared']);w=np.empty(n)
        for first in range(0,n,100):
            stop=min(n,first+100);w[first:stop]=(f['winding'][first:stop].astype(float)**2).mean((1,2))
        return np.column_stack([f['magnetization_squared'][:],f['magnetization_fourth'][:],f['structure_zero'][:],f['structure_min'][:].mean(1),f['mean_link_cosine'][:],f['vortex_density'][:],w])

def quantities(m,N):
    m=np.asarray(m);a,b,f0,fm,cosine,vortex,w=np.moveaxis(m,-1,0)
    with np.errstate(invalid='ignore',divide='ignore'):
        R=np.sqrt(f0/fm-1)/(2*N*np.sin(np.pi/N));U=1-b/(2*a*a)
    return np.stack([U,R,a,b,cosine,vortex,w],axis=-1)

def summarize(points,cut=2001,draws=2000,factor=1):
    points=sorted(points,key=lambda p:p['id']);N=points[0]['N'];arrays=[read(p)[cut:] for p in points]
    assert all(p['N']==N and p['D']==2 and p['beta']==points[0]['beta'] and p['epsilon']==points[0]['epsilon'] for p in points)
    rng=np.random.default_rng([2026101402,N,round(points[0]['epsilon_over_pi']*1000),round(points[0]['beta']*1000),cut,factor])
    total=np.zeros((draws,7));length=np.zeros(draws);diagnostics=[];individual=[];individual_samples=[]
    for p,a in zip(points,arrays):
        m2,m4,f0,fm,cosine,vortex,w=a.T
        influence=(f0-f0.mean())/fm.mean()-f0.mean()/fm.mean()**2*(fm-fm.mean())
        taus={k:helper.tau_int(v) for k,v in dict(m2=m2,xi_ratio=influence,cosine=cosine,vortex=vortex,winding2=w).items()}
        width=factor*max(20,int(np.ceil(10*max(t for t in taus.values() if t is not None))))
        nb=len(a)//width;diagnostics.append(dict(id=p['id'],retained=len(a),block_width=width,blocks=nb,taus=taus))
        if nb<2:return dict(status='insufficient_blocks',diagnostics=diagnostics),None
        edges=np.linspace(0,len(a),nb+1,dtype=int);sums=np.add.reduceat(a,edges[:-1],axis=0);sizes=np.diff(edges)
        ss=np.zeros((draws,7));sn=np.zeros(draws)
        for first in range(0,draws,100):
            stop=min(draws,first+100);ix=rng.integers(nb,size=(stop-first,nb));ss[first:stop]=sums[ix].sum(1);sn[first:stop]=sizes[ix].sum(1)
        total+=ss;length+=sn;individual.append(quantities(a.mean(0),N));individual_samples.append(quantities(ss/sn[:,None],N))
    central=quantities(np.concatenate(arrays).mean(0),N);samples=quantities(total/length[:,None],N);flags=[]
    if any(d['blocks']<50 for d in diagnostics):flags.append('fewer_than_50_blocks_per_start')
    if not np.isfinite(samples[:,1]).all():flags.append('invalid_xi_bootstrap_draws_retained')
    for name,index in [('xi_over_L',1),('m2',2)]:
        v=np.array(individual)[:,index];se=np.nanstd(np.array(individual_samples)[:,:,index],axis=1,ddof=1);eligible=np.isfinite(v)&np.isfinite(se)&(se>0)
        if eligible.sum()>1:
            mean=np.average(v[eligible],weights=1/se[eligible]**2);pv=chi2.sf(np.sum(((v[eligible]-mean)/se[eligible])**2),eligible.sum()-1)
            if pv<.001:flags.append('start_disagreement_'+name)
    keys=['binder','xi_over_L','m2','m4','link_cosine','vortex_density','winding_squared']
    return dict(status='conditional_joint_bootstrap',values=dict(zip(keys,central.tolist())),errors=dict(zip(keys,np.nanstd(samples,axis=0,ddof=1).tolist())),flags=flags,diagnostics=diagnostics,individual_values=np.array(individual).tolist(),retained=sum(len(a) for a in arrays),cut=cut,block_factor=factor,draws=draws,seed_components=[2026101402,N,round(points[0]['epsilon_over_pi']*1000),round(points[0]['beta']*1000),cut,factor]),samples
