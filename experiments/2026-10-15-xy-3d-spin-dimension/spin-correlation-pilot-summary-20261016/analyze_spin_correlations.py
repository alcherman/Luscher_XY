"""Time-aware joint bootstrap for long-distance spin correlation measurements."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ['MPLBACKEND']='Agg';os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-spin-correlation-mpl')
from pathlib import Path
import argparse,json,hashlib,shutil,time
import numpy as np
from scipy.stats import chi2
import common_xi_quotients as cq
import measure_spin_correlations as ex
OUT=Path(__file__).resolve().parent
DRAWS=2000

def clean(x):
    if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [clean(v) for v in x]
    if isinstance(x,np.ndarray):return clean(x.tolist())
    if isinstance(x,(np.floating,float)):return float(x) if np.isfinite(x) else None
    if isinstance(x,np.integer):return int(x)
    if isinstance(x,np.bool_):return bool(x)
    return x

def write(p,a):ex.atomic(p,clean(a))
def tau(v):
    if len(v)<8:return None
    t=cq.m.helper.tau_int(v)
    return float(t) if t is not None and np.isfinite(t) else None

def derived(a,L):
    a=np.asarray(a)
    with np.errstate(invalid='ignore',divide='ignore'):
        R=np.sqrt(a[...,3]/a[...,4]-1)/(2*L*np.sin(np.pi/L));lgq=np.log(a[...,0]);lgh=np.log(a[...,1]);lm=np.log(a[...,2])
    return np.stack([R,lgq,lgh,lm],axis=-1)

def block_width_cycles(checkpoint_tau,dense_tau,cadence=100,factor=2):
    # tau_CP describes the sampled process. cadence*tau_CP is a blocking scale,
    # not a reconstruction of the unobserved sub-cadence ACF.
    values=[1.]+[cadence*v for v in checkpoint_tau.values() if v is not None and np.isfinite(v)]+[v for v in dense_tau.values() if v is not None and np.isfinite(v)]
    requested=factor*max(20,int(np.ceil(10*max(values))))
    cp=max(1,int(np.ceil(requested/cadence)))
    return cp,cp*cadence

def bootstrap_blocks(a,width,rng,draws=DRAWS):
    nb=len(a)//width
    if nb<2:return None,dict(blocks=nb,rows=len(a),width_checkpoints=width)
    edges=np.linspace(0,len(a),nb+1,dtype=int);sums=np.add.reduceat(a,edges[:-1],axis=0);counts=np.diff(edges);ix=rng.integers(nb,size=(draws,nb));ss=sums[ix].sum(1);nn=counts[ix].sum(1)
    return (ss,nn),dict(blocks=nb,rows=len(a),width_checkpoints=width,actual_block_lengths_checkpoints=counts.tolist())

def summarize(entries,cut=1000,factor=2):
    L=entries[0][0]['point']['N'];ids=sorted(e[0]['id'] for e in entries);seed=int(hashlib.sha256(json.dumps([ids,cut,factor,'spin-correlation-v1']).encode()).hexdigest()[:16],16);rng=np.random.default_rng(seed)
    total=np.zeros((DRAWS,6));n=np.zeros(DRAWS);arrays=[];diags=[];individual=[];flags=[];insufficient=False
    axis_means=[];sin_means=[]
    for meta,file in entries:
        if ex.sha(file)!=meta['derived_data_sha256']:raise ValueError('Derived source changed')
        with np.load(file) as d:
            mask=d['sweep']>cut;times=d['sweep'][mask];a=d['series'][mask];axis_means.append(d['axis_real'][mask].mean(0));sin_means.append(d['axis_imag'][mask].mean(0));arrays.append(a)
            if len(a)<4 or not np.all(np.diff(times)==100):raise ValueError('Missing/nonregular measurement time series')
            mu=a.mean(0);influence=(a[:,3]-mu[3])/mu[4]-mu[3]*(a[:,4]-mu[4])/mu[4]**2
            monitor=dict(G_quarter=a[:,0],G_half=a[:,1],m2=a[:,2],R_influence=influence,winding2=a[:,5]);cp_tau={}
            for name,v in monitor.items():
                cp_tau[name+'_postcut']=tau(v);cp_tau[name+'_latehalf']=tau(v[len(v)//2:])
            dense_tau={}
            dm=d['dense_sweep']>cut;m=d['dense_m2'][dm];z=d['dense_S0'][dm];f=d['dense_F'][dm];ri=(z-z.mean())/f.mean()-z.mean()*(f-f.mean())/f.mean()**2
            for name,v in [('m2',m),('R_influence',ri)]:
                dense_tau[name+'_postcut']=tau(v);dense_tau[name+'_latehalf']=tau(v[len(v)//2:])
            for k,v in meta.get('original_retention_tau_estimates_cycles',{}).items():
                if k.startswith('winding2'):dense_tau['archived_retention_'+k]=v
            width,physical_width=block_width_cycles(cp_tau,dense_tau,100,factor);bs,info=bootstrap_blocks(a,width,rng)
            diagnostic=dict(id=meta['id'],samples=len(a),first_cycle=int(times[0]),last_cycle=int(times[-1]),cadence_cycles=100,checkpoint_tau=cp_tau,dense_tau_cycles=dense_tau,block_width_cycles=physical_width,blocking=info,negative_configuration_G_counts={'quarter':int(np.sum(a[:,0]<=0)),'half':int(np.sum(a[:,1]<=0))},constant_winding_axes=meta.get('constant_winding_axes_all_checkpoints'))
            if info['blocks']<50:flags.append('fewer_than_50_blocks_per_start')
            if bs is None:insufficient=True;flags.append('insufficient_blocks_per_start')
            else:
                total+=bs[0];n+=bs[1];v=derived(bs[0]/bs[1][:,None],L);valid=np.isfinite(v).all(1)
                individual.append(dict(id=meta['id'],raw_mean=mu,derived=derived(mu,L),joint_valid_fraction=float(valid.mean()),stat_SD=np.std(v[valid],axis=0,ddof=1) if valid.sum()>1 else None))
            halves=[]
            for half in np.array_split(a,2):
                nb=len(half)//width
                if nb<2:halves.append(None);continue
                bm=np.array([x.mean(0) for x in np.array_split(half,nb)]);values=np.column_stack([bm[:,0],bm[:,1],bm[:,2]])
                halves.append((half[:,[0,1,2]].mean(0),values.std(0,ddof=1)/np.sqrt(nb)))
            if all(v is not None for v in halves):
                zz=(halves[1][0]-halves[0][0])/np.maximum(np.hypot(halves[0][1],halves[1][1]),1e-15);diagnostic['half_z_Gq_Gh_m2']=zz
                if np.max(abs(zz))>4:flags.append('checkpoint_half_shift')
            else:flags.append('insufficient_half_blocks')
            diags.append(diagnostic)
    raw=np.concatenate(arrays).mean(0);center=derived(raw,L);base=dict(L=L,ids=ids,raw_mean=raw,derived=center,flags=sorted(set(flags)),diagnostics=diags,individual=individual,axis_real_mean=np.average(axis_means,axis=0,weights=[len(a) for a in arrays]),axis_imag_mean=np.average(sin_means,axis=0,weights=[len(a) for a in arrays]),samples_total=sum(len(a) for a in arrays))
    if insufficient:return dict(base,status='insufficient_blocks',joint_valid_fraction=0.),None
    B=total/n[:,None];Z=derived(B,L);valid=np.isfinite(Z).all(1)
    if len(individual)==2 and all(v['stat_SD'] is not None for v in individual):
        z=(individual[1]['derived']-individual[0]['derived'])/np.maximum(np.hypot(individual[1]['stat_SD'],individual[0]['stat_SD']),1e-15);base['start_z_R_logGq_logGh_logm2']=z
        if np.max(abs(z))>4:base['flags']=sorted(set(base['flags']+['checkpoint_start_shift']))
    base.update(status='conditional_checkpoint_block_summary' if valid.mean()>=.99 and np.isfinite(center).all() else 'unresolved_joint_tails',joint_valid_fraction=float(valid.mean()),invalid_draws=dict(R=int(np.sum(~np.isfinite(Z[:,0]))),G_quarter=int(np.sum(B[:,0]<=0)),G_half=int(np.sum(B[:,1]<=0)),m2=int(np.sum(B[:,2]<=0))),MC95=np.quantile(Z[valid],[.025,.975],axis=0).T if valid.sum()>1 else None,SD=np.std(Z[valid],axis=0,ddof=1) if valid.sum()>1 else None,raw_MC95=np.quantile(B,[.025,.975],axis=0).T)
    return base,Z


def fixed_q_slope(groups,boot):
    N=np.array([g['L'] for g in groups]);y=np.array([g['derived'][1:3] for g in groups]);B=boot[:,:,1:3];valid=np.isfinite(B).all((1,2))
    if valid.mean()<.99:return dict(status='unresolved_joint_tails',joint_valid_fraction=float(valid.mean()))
    C=np.stack([np.cov(B[valid,i,:].T) for i in range(len(N))]);W=np.linalg.inv(np.linalg.cholesky(C));X=np.zeros((len(N),2,3));X[:,0,0]=1;X[:,1,1]=1;X[:,:,2]=-2*np.log(N[:,None]/32.)
    A=np.einsum('nij,njk->nik',W,X).reshape(-1,3);target=np.einsum('nij,nj->ni',W,y).ravel();coef=np.linalg.lstsq(A,target,rcond=1e-11)[0];T=np.linalg.pinv(A);samples=np.einsum('nij,dnj->dni',W,B[valid]).reshape(valid.sum(),-1)@T.T;res=A@coef-target;stat=float(res@res);dof=len(target)-3
    return dict(status='conditional_fixed_q_effective_dimension',q=groups[0]['q'],sizes=N,Delta=float(coef[-1]),SD=float(samples[:,-1].std(ddof=1)),MC95=np.quantile(samples[:,-1],[.025,.975]),chi2=stat,dof=dof,nominal_p=float(chi2.sf(stat,dof)),joint_valid_fraction=float(valid.mean()),scope='Fixed representative q, two separated correlators with independent amplitudes and common slope. No correction term or critical-q propagation; not the asymptotic CFT dimension or a systematic bound.')


def main(args):
    folder=Path(args.output);folder.mkdir(exist_ok=False);shutil.copyfile(__file__,folder/Path(__file__).name);groups={};provenance=[]
    for source in args.inputs:
        source=Path(source)
        for p in sorted(source.glob('*.json')):
            if p.stem in ['design','status']:continue
            meta=json.loads(p.read_text())
            if 'id' not in meta or 'derived_data_sha256' not in meta:continue
            key=(meta['point']['N'],meta['group_q']);groups.setdefault(key,[]).append((meta,p.with_suffix('.npz')));provenance.append(dict(path=str(p),sha256=ex.sha(p),data_path=str(p.with_suffix('.npz')),data_sha256=ex.sha(p.with_suffix('.npz'))))
    if len({m['id'] for g in groups.values() for m,_ in g})!=sum(len(g) for g in groups.values()):raise ValueError('Duplicated input chains')
    results=[];draws=[]
    for (L,q),entries in sorted(groups.items()):
        r,B=summarize(entries,args.cut,args.factor);r['q']=q;r['source_metadata']=[str(p.with_suffix('.json')) for _,p in entries];results.append(r)
        if B is not None:draws.append((L,q,B));np.savez_compressed(folder/f'L{L}-q{q:.8f}-draws.npz',samples=B)
        print(L,q,r['status'],'G',r['raw_mean'][:2],'valid',r['joint_valid_fraction'],'blocks',[d['blocking']['blocks'] for d in r['diagnostics']],'flags',r['flags'],flush=True)
    selected=[r for r in results if abs(r['q']-.80075)<1e-12];selected_B=[next(b for n,q,b in draws if n==r['L'] and abs(q-r['q'])<1e-12) for r in selected if any(n==r['L'] and abs(q-r['q'])<1e-12 for n,q,b in draws)]
    slope=fixed_q_slope(selected,np.stack(selected_B,axis=1)) if len(selected)==4 and len(selected_B)==4 else dict(status='incomplete_sizes')
    report=dict(created=time.time(),settings=dict(cut_cycles=args.cut,block_factor=args.factor,draws=DRAWS,cadence_cycles=100,log_columns=['R','logG_quarter','logG_half','logm2']),groups=results,fixed_q_slope=slope,provenance=provenance,scope='Existing uniformly sampled checkpoints. Spatial sites/axes not independent MC events; block widths in physical cycles. Sparse-chain autocorrelation and bootstrap intervals remain conditional; not a certified systematic error or CFT precision.')
    write(folder/'analysis.json',report);write(folder/'status.json',dict(state='summaries_complete_pending_review',groups=len(results),fixed_q_slope=slope));print('Fixed q slope',clean(slope),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',nargs='+',required=True);p.add_argument('--output',required=True);p.add_argument('--cut',type=int,default=1000);p.add_argument('--factor',type=int,default=2);main(p.parse_args())
