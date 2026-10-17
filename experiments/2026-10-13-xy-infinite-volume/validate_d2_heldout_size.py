"""Predict the new N512 roots using only smaller boxes; retain shared covariance."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[k]='1'
import json,shutil
from pathlib import Path
import numpy as np
from scipy.stats import chi2
import run
OUT=Path(__file__).resolve().parent

def prediction_check(rows,B,minimum=32,offset=0):
    eligible=np.array([minimum<=r['N']<=(512 if r['observable']=='xi_over_L' else 256) for r in rows]);rows=[r for r,k in zip(rows,eligible) if k];B=B[:,eligible]
    train=np.array([r['N']<(512 if r['observable']=='xi_over_L' else 256) for r in rows]);test=~train
    if test.sum()!=2 or train.sum()<4:return dict(status='insufficient_sizes')
    valid=np.isfinite(B).all(1)
    if valid.mean()<.99:return dict(status='bootstrap_coverage_unresolved',joint_valid_fraction=float(valid.mean()))
    B=B[valid];y=np.array([r['beta_location'] for r in rows]);cov=np.cov(B[:,train],rowvar=False)
    ev,U=np.linalg.eigh(cov)
    if ev.min()<=ev.max()*1e-10:return dict(status='covariance_rank_unresolved')
    x=np.array([(np.log(r['N']*(np.sqrt(2) if r['observable']=='eta_eff' else 1))+offset)**-2 for r in rows])
    X=np.column_stack([np.ones(len(rows)),x*np.array([r['observable']=='xi_over_L' for r in rows]),x*np.array([r['observable']=='eta_eff' for r in rows])])
    W=(U/np.sqrt(ev))@U.T;Z=W@X[train];norm=np.linalg.norm(Z,axis=0);T=(np.linalg.pinv(Z/norm)@W)/norm[:,None]
    coef=T@y[train];drawcoef=B[:,train]@T.T
    predicted=X[test]@coef;predboot=drawcoef@X[test].T;residual=y[test]-predicted;resboot=B[:,test]-predboot
    predictive_cov=np.cov(resboot,rowvar=False)
    statistic=float(residual@np.linalg.solve(predictive_cov,residual));train_res=y[train]-X[train]@coef;trainstat=float(np.sum((W@train_res)**2))
    return dict(status='descriptive_heldout_check',minimum_N=minimum,log_offset=offset,joint_valid_fraction=float(valid.mean()),training_beta_infinite=float(coef[0]),training_chi2=trainstat,training_dof=int(train.sum()-3),heldout_chi2=statistic,heldout_dof=2,heldout_p=float(chi2.sf(statistic,2)),prediction_covariance=predictive_cov.tolist(),rows=[dict(observable=r['observable'],N=r['N'],measured_beta=r['beta_location'],predicted_beta=float(p),difference=float(d),difference_MC95=np.quantile(draw,[.025,.975]).tolist()) for r,p,d,draw in zip([r for r,k in zip(rows,test) if k],predicted,residual,resboot.T)])

def main():
    source=OUT/json.loads((OUT/'latest-d2-eta-correlation.json').read_text())['folder'];a=json.loads((source/'analysis.json').read_text())
    folder=OUT/('d2-heldout-N512-validation-'+run.now()[:19].replace('-','').replace(':','').replace('T','-'));folder.mkdir();shutil.copyfile(__file__,folder/Path(__file__).name)
    checks=[]
    for q in [.60,.80,1.]:
        rows=[r for r in a['locations'] if r['epsilon_over_pi']==q and r['status']=='conditional_crossing']
        B=np.column_stack([np.load(source/f"e{round(q*1000):04d}_{r['observable']}_N{r['N']:03d}.npz")['samples'] for r in rows])
        for minimum in [32,64]:
            for C in [0,1,2,4]:checks.append(dict(epsilon_over_pi=q,**prediction_check(rows,B,minimum,C)))
    run.atomic_json(folder/'analysis.json',dict(created=run.now(),source=str(source),source_sha256=run.sha256(source/'analysis.json'),checks=checks,scope='Smaller-box inverse-log-squared models predict R512 and eta256:512 roots. Shared raw-source bootstrap covariance is retained in training and prediction residuals. This is a descriptive heldout-size check of previously specified offsets/minimum sizes, not a precision guarantee, independent replication of all smaller data, or model selection by best p value. No known benchmark is fitted or imposed.'))
    run.atomic_json(OUT/'latest-d2-heldout-validation.json',dict(folder=folder.name,status='descriptive_not_final'))
    print(folder,flush=True)
    for c in checks:
        if c['minimum_N']==32:print(c['epsilon_over_pi'],c['log_offset'],c.get('heldout_p'),[(r['observable'],r['difference']) for r in c.get('rows',[])],flush=True)

if __name__=='__main__':main()
