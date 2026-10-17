"""Analytical normalization and correlated BKT-root regression checks."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments/2026-10-13-xy-infinite-volume'))
from d2_joint_statistics import quantities
from extrapolate_d2_eta_correlation import fit_roots


def test_joint_moment_normalization_and_invalid_R():
    N=32;m2=.3;R=.75;f0=N*N*m2;fm=f0/(1+(2*N*np.sin(np.pi/N)*R)**2)
    q=quantities([m2,.1,f0,fm,.8,.01,0],N)
    assert q[1]==pytest.approx(R);assert q[2]==pytest.approx(m2)
    invalid=quantities([m2,.1,.5,1.,.8,.01,0],N)
    assert np.isnan(invalid[1]);assert invalid[2]==pytest.approx(m2)


def example():
    rows=[];truth=[]
    for observable,sizes,A in [('xi_over_L',[16,32,64,128],-.7),('eta_eff',[16,32,64],-.4)]:
        for N in sizes:
            effective=N*(np.sqrt(2) if observable=='eta_eff' else 1)
            beta=1.12+A/(np.log(effective)+.5)**2
            rows.append(dict(observable=observable,N=N,pair=[N,2*N] if observable=='eta_eff' else None,beta_location=beta,flags=[],individual_flags={}))
            truth.append(beta)
    rng=np.random.default_rng(441);common=rng.normal(size=(2000,1));draws=np.array(truth)+.002*common+.001*rng.normal(size=(2000,len(rows)))
    return rows,draws


def test_correlated_log_fit_recovers_intercept_and_covariance():
    rows,B=example();r=fit_roots(rows,B,16,.5,.5,'joint')
    assert r['status']=='exploratory_log_fit';assert r['beta_infinite']==pytest.approx(1.12,abs=1e-11)
    assert r['coefficients'][1:]==pytest.approx([-.7,-.4],abs=1e-10);assert r['dof']==4
    X=np.column_stack([np.ones(len(rows)),[1/(np.log(v['N'])+.5)**2 if v['observable']=='xi_over_L' else 0 for v in rows],[1/(np.log(v['N']*np.sqrt(2))+.5)**2 if v['observable']=='eta_eff' else 0 for v in rows]])
    C=np.cov(B,rowvar=False);expected=np.linalg.inv(X.T@np.linalg.solve(C,X))[0,0]
    assert np.var(r['samples'],ddof=1)==pytest.approx(expected,rel=1e-9)


def test_common_bootstrap_validity_is_not_silently_conditioned():
    rows,B=example();B[:30,0]=np.nan
    r=fit_roots(rows,B,16,.5,.5,'joint')
    assert r['status']=='bootstrap_coverage_unresolved';assert r['joint_valid_fraction']==pytest.approx(.985)
