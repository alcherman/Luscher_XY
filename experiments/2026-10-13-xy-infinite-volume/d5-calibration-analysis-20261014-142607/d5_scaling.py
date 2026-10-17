"""Periodic above-upper-critical-dimension size fits, independent of D4 logs."""
import numpy as np
from scipy.stats import chi2

def fit(N,y,B,minimum=0,correction=False):
    keep=np.asarray(N)>=minimum;N=np.asarray(N,dtype=float)[keep];y=np.asarray(y)[keep];B=np.asarray(B)[:,keep]
    X=np.column_stack([np.ones(len(N)),N**(-2.5)]+([N**(-3)] if correction else []))
    if len(N)<X.shape[1]:return dict(status='insufficient_sizes')
    valid=np.isfinite(B).all(1)
    if valid.mean()<.99:return dict(status='bootstrap_coverage_unresolved',joint_valid_fraction=float(valid.mean()))
    B=B[valid];C=np.cov(B,rowvar=False);e,U=np.linalg.eigh(C)
    if np.min(e)<=np.max(e)*1e-10:return dict(status='covariance_rank_unresolved')
    W=(U/np.sqrt(e))@U.T;Z=W@X;scale=np.linalg.norm(Z,axis=0)
    if np.linalg.matrix_rank(Z/scale)<X.shape[1]:return dict(status='unidentifiable')
    T=np.linalg.pinv(Z/scale)@W;coef=(T@y)/scale;draws=(B@T.T)/scale
    r=W@(y-X@coef);stat=float(r@r);df=len(N)-X.shape[1]
    return dict(status='zero_dof_diagnostic' if df==0 else 'conditional_D5_size_fit',sizes=N.astype(int).tolist(),minimum_N=minimum,correction=correction,location_infinite=float(coef[0]),MC95=np.quantile(draws[:,0],[.025,.975]).tolist(),coefficients=coef.tolist(),chi2=stat,dof=df,p=float(chi2.sf(stat,df)) if df>0 else None,joint_valid_fraction=float(valid.mean()),samples=draws[:,0],scope='Leading N^-5/2; additive N^-3 correction. Conditional hypothesis for periodic above-critical-dimension scaling; zeroDOFfits diagnostics only.')
