"""Check nonlinear estimators and invalid correlation-length handling."""
import numpy as np
from analyze import summarize


def sample():
    n=1001;L=16;rng=np.random.default_rng(20260924)
    a=rng.uniform(.2,.6,n)
    d=dict(magnetization_squared=a,magnetization_fourth=a*a,structure_zero=L*L*a,
        structure_min=np.repeat((L*L*a/5)[:,None],2,axis=1),
        mean_link_cosine=rng.uniform(.5,.7,n),vortex_density=np.zeros(n),
        winding=np.zeros((n,2,L),dtype=int),local_accepted=np.full(n,100),
        cluster_sizes=np.full((n,1),80),winding_accepted=np.zeros(n))
    return d,dict(L=L,seed=1234)


def test_joint_bootstrap_recomputes_binder_and_correlation_length():
    d,p=sample();result,s=summarize(d,p,1,draws=200)
    a=d['magnetization_squared'][1:]
    np.testing.assert_allclose(result['values']['binder'],1-np.mean(a*a)/(2*a.mean()**2))
    np.testing.assert_allclose(s[:,0],1-s[:,3]/(2*s[:,2]**2),rtol=1e-13)
    expected=1/(p['L']*np.sin(np.pi/p['L']))
    np.testing.assert_allclose(s[:,1],expected,rtol=1e-13)
    assert result['retained']==1000 and result['blocks']>=20


def test_negative_structure_factor_ratio_is_not_clipped_into_a_length():
    d,p=sample();d['structure_min']*=10
    result,s=summarize(d,p,1,draws=200)
    assert result['values']['xi_over_L'] is None
    assert result['errors']['xi_over_L'] is None
    assert result['invalid_xi_bootstrap_fraction']==1
    assert np.isnan(s[:,1]).all() and np.isfinite(s[:,0]).all()


def test_retained_tail_contributes_to_nonlinear_estimator():
    d,p=sample();d['magnetization_squared'][-1]=.99;d['magnetization_fourth'][-1]=.99**2
    result,_=summarize(d,p,4,draws=200)
    a=d['magnetization_squared'][4:];b=d['magnetization_fourth'][4:]
    assert result['retained']==997
    np.testing.assert_allclose(result['values']['binder'],1-b.mean()/(2*a.mean()**2),rtol=1e-13)
