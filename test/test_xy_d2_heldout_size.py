"""Independent predictive covariance oracle for shared-source size holdouts."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments/2026-10-13-xy-infinite-volume'))
from validate_d2_heldout_size import prediction_check

def example():
    rows=[];X=[]
    for obs,sizes,A in [('xi_over_L',[32,64,128,256,512],-.5),('eta_eff',[32,64,128,256],-.3)]:
        for N in sizes:
            x=np.log(N*(np.sqrt(2) if obs=='eta_eff' else 1))**-2
            rows.append(dict(observable=obs,N=N,beta_location=1.12+A*x))
            X.append([1,x if obs=='xi_over_L' else 0,x if obs=='eta_eff' else 0])
    rng=np.random.default_rng(617);B=np.array([r['beta_location'] for r in rows])+.002*rng.normal(size=(2000,1))+.001*rng.normal(size=(2000,len(rows)))
    return rows,B,np.array(X)

def test_predictive_covariance_keeps_cross_correlations():
    rows,B,X=example();r=prediction_check(rows,B)
    assert r['training_beta_infinite']==pytest.approx(1.12,abs=1e-11)
    assert r['heldout_chi2']<1e-15
    train=np.array([v['N']<(512 if v['observable']=='xi_over_L' else 256) for v in rows]);test=~train
    C=np.cov(B,rowvar=False);Ctt=C[np.ix_(train,train)];Cht=C[np.ix_(test,train)];Chh=C[np.ix_(test,test)]
    G=np.linalg.solve(X[train].T@np.linalg.solve(Ctt,X[train]),X[train].T@np.linalg.inv(Ctt));P=X[test]@G
    expected=Chh+P@Ctt@P.T-Cht@P.T-P@Cht.T
    np.testing.assert_allclose(r['prediction_covariance'],expected,rtol=1e-9,atol=1e-16)

def test_bad_shared_draw_coverage_stays_unresolved():
    rows,B,X=example();B[:30,-1]=np.nan
    r=prediction_check(rows,B);assert r['status']=='bootstrap_coverage_unresolved';assert r['joint_valid_fraction']==pytest.approx(.985)
