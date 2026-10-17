"""Exploratory BKT logarithmic extrapolations with correlated eta/R roots."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl');os.environ['MPLBACKEND']='Agg'
import argparse,json,shutil
from pathlib import Path
import numpy as np
from scipy.stats import chi2
import run
from scout import json_finite
OUT=Path(__file__).resolve().parent

def fit_roots(rows,B,minimum=32,offset_R=0,offset_eta=0,mode='joint'):
    keep=np.array([r['N']>=minimum and (mode=='joint' or r['observable']==mode) for r in rows]);chosen=[r for r,k in zip(rows,keep) if k];B=B[:,keep]
    if len(chosen)<(4 if mode=='joint' else 3):return dict(status='insufficient_sizes')
    valid=np.isfinite(B).all(1);fraction=float(valid.mean())
    if fraction<.99:return dict(status='bootstrap_coverage_unresolved',joint_valid_fraction=fraction)
    B=B[valid];y=np.array([r['beta_location'] for r in chosen]);x=[]
    for r in chosen:
        N=r['N']*(np.sqrt(2) if r['observable']=='eta_eff' else 1);C=offset_eta if r['observable']=='eta_eff' else offset_R
        if np.log(N)+C<=0:return dict(status='invalid_log_offset')
        x.append((np.log(N)+C)**(-2))
    X=np.column_stack([np.ones(len(chosen)),*[np.array(x)*np.array([r['observable']==obs for r in chosen]) for obs in (['xi_over_L','eta_eff'] if mode=='joint' else [mode])]])
    covariance=np.cov(B,rowvar=False);eigen,U=np.linalg.eigh(covariance);floor=max(eigen.max()*1e-10,1e-20);rank=int((eigen>floor).sum())
    if rank<len(chosen):return dict(status='covariance_rank_unresolved',covariance_rank=rank,observations=len(chosen))
    W=(U/np.sqrt(eigen))@U.T;Z=W@X;norm=np.linalg.norm(Z,axis=0);scaled=Z/norm
    if np.linalg.matrix_rank(scaled)<X.shape[1]:return dict(status='unidentifiable')
    transform=np.linalg.pinv(scaled)@W;coef=transform@y/norm;draws=(B@transform.T)/norm
    residual=y-X@coef;statistic=float(np.dot(W@residual,W@residual));df=len(chosen)-len(coef)
    return dict(status='exploratory_log_fit',beta_infinite=float(coef[0]),MC95=np.quantile(draws[:,0],[.025,.975]).tolist(),coefficients=coef.tolist(),chi2=statistic,dof=df,p=float(chi2.sf(statistic,df)) if df>0 else None,joint_valid_fraction=fraction,covariance_rank=rank,condition_number=float(np.linalg.cond(scaled)),covariance=covariance.tolist(),root_coordinates=[dict(observable=r['observable'],N=r['N'],pair=r['pair']) for r in chosen],flags=sorted({f for r in chosen for f in r.get('flags',[])}),individual_flags={k:v for r in chosen for k,v in r.get('individual_flags',{}).items()},samples=draws[:,0])

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('snapshot',type=Path);args=parser.parse_args();a=json.loads((args.snapshot/'analysis.json').read_text())
    folder=OUT/('d2-eta-R-log-extrapolation-'+run.now()[:19].replace('-','').replace(':','').replace('T','-'));folder.mkdir();shutil.copyfile(__file__,folder/Path(__file__).name)
    fits=[]
    for q in sorted({r['epsilon_over_pi'] for r in a['locations']}):
        rows=[r for r in a['locations'] if r['epsilon_over_pi']==q and r['status']=='conditional_crossing'];B=np.column_stack([np.load(args.snapshot/f"e{round(q*1000):04d}_{r['observable']}_N{r['N']:03d}.npz")['samples'] for r in rows]) if rows else np.empty((2000,0))
        for minimum in [16,32,64]:
            for CR,CE in [(-1,-1),(0,0),(1,1),(2,2),(4,4),(0,2),(2,0)]:
                for mode in ['joint','xi_over_L','eta_eff']:
                    r=fit_roots(rows,B,minimum,CR,CE,mode);index=len(fits)
                    if 'samples' in r:np.savez_compressed(folder/f'fit{index:04d}.npz',samples=r.pop('samples'))
                    fits.append(dict(epsilon_over_pi=q,minimum_N=minimum,offset_R=CR,offset_eta=CE,mode=mode,result=r))
    run.atomic_json(folder/'analysis.json',json_finite(dict(created=run.now(),source=str(args.snapshot/'analysis.json'),source_sha256=run.sha256(args.snapshot/'analysis.json'),fits=fits,scope='Exploratory fixed-offset inverse-log-squared shifts with shared beta infinity and separate observable amplitudes; effective eta size is geometric mean sqrt2*N. Common bootstrap validity mask and full root covariance retained. Offsets and geometry are sensitivity choices, not a calibrated uncertainty family. No accepted infinite-volume result or precision certification. Calibration against endpoint and size cuts required.')))
    print('snapshot',folder,flush=True)
    for f in fits:
        if f['epsilon_over_pi'] in [.60,.80,1.] and f['minimum_N']==32 and f['mode']=='joint' and (f['offset_R'],f['offset_eta'])==(0,0):print(f,flush=True)
