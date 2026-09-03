import numpy as np
from pool import ratios,moments,combine
from boundary import crossing

def test_pool_means_before_nonlinear_ratios():
    L=16
    raw=np.array([[.25,.08,1.,.6,.01,.03],[.4,.18,2.,.7,.02,.05]])
    q=ratios(raw,L)
    np.testing.assert_allclose(moments(q,L),raw,rtol=1e-13)
    result=combine(q,[1,3],L)
    np.testing.assert_allclose(result,ratios(np.average(raw,weights=[1,3],axis=0),L),rtol=1e-13)
    assert abs(result[0]-np.average(q[:,0],weights=[1,3]))>1e-3

def test_crossing_is_interpolation_not_extrapolation():
    assert crossing([0,.02,.04],[.48,.49,.50],.495)==.03
    assert np.isnan(crossing([0,.02,.04],[.48,.49,.50],.51))
    assert np.isnan(crossing([0,.02,.04],[.48,.5,.48],.49))
    assert np.isnan(crossing([0,.02,.04],[.49,.49,.49],.49))

def test_joint_replicate_bootstrap_preserves_moment_identity():
    rng=np.random.default_rng(12);raw=np.zeros((4,100,6));raw[:,:,0]=rng.uniform(.2,.4,(4,100))
    raw[:,:,1]=raw[:,:,0]**2*1.018;raw[:,:,2]=rng.uniform(1,2,(4,100));raw[:,:,3:]=.1
    q=ratios(raw,32);pooled=combine(q,[10,10,10,10],32)
    np.testing.assert_allclose(pooled,ratios(raw.mean(0),32),rtol=1e-13)

def test_local_linear_root_recovers_known_crossing_and_rejects_extrapolation():
    from boundary import linear_roots
    x=np.array([.1,.12,.14]);y=.49+.02*(x-.12)
    root,coeff=linear_roots(x,y,np.ones(3)*.0001,.49)
    np.testing.assert_allclose(root,.12,atol=1e-12)
    np.testing.assert_allclose(coeff,[.49,.02],atol=1e-12)
    root,_=linear_roots(x,y,np.ones(3)*.0001,.495)
    assert np.isnan(root)

def test_local_window_excludes_distant_curved_response():
    from boundary import local_window,linear_roots
    x=np.array([.1,.2,.3,.32,.34,.36,.38,.5])
    y=.49+.1*(x-.35)+.2*(x-.35)**2
    err=np.full(len(x),.0001)
    ix=local_window(x,y,err,.49)
    assert x[ix[-1]]-x[ix[0]] < .040001
    root,_=linear_roots(x[ix],y[ix],err[ix],.49)
    assert abs(root-.35)<.001
