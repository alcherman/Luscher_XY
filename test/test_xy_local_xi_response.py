"""Analytic response/quotient tests, including explicit failed-root accounting."""
import sys
from pathlib import Path
import numpy as np
import pytest
P=Path(__file__).resolve().parents[1]/'experiments/2026-10-15-xy-3d-spin-dimension'
sys.path.insert(0,str(P))
import local_xi_response as r


def values(N,t,shift=0,amplitude=1,degree=2):
    m2=amplitude*N**(-2*.519)*(1+.05*t+(.01*t*t if degree==2 else 0))
    R=.59-(N/16)*.05*(t-shift)+(.002*(t-shift)**2 if degree==2 else 0)
    S=N**3*m2;c=(2*N*np.sin(np.pi/N))**2
    return np.column_stack([m2,1.4*m2*m2,np.sqrt(m2),S,S/(1+c*R*R)])


@pytest.mark.parametrize('degree',[1,2])
def test_exact_local_response_and_joint_crossing_preserve_dimension(degree):
    q=np.linspace(.8005,.8015,9);t=np.linspace(-1,1,9);a=values(16,t,degree=degree);b=values(32,t,degree=degree)
    shifts=np.linspace(-.1,.1,151)
    A=np.array([values(16,t,s,1+s/20,degree) for s in shifts]);B=np.array([values(32,t,s,1+s/20,degree) for s in shifts])
    record,samples=r.pair_fit(q,a,A,16,b,B,32,.8005,.8015,degree)
    assert record['status']=='conditional_local_quotient'
    assert record['valid_fraction']==1
    assert record['epsilon_over_pi']==pytest.approx(.801,abs=1e-12)
    assert record['Delta_eff']==pytest.approx(.519,abs=1e-12)
    np.testing.assert_allclose(samples[:,0],.801+.0005*shifts,atol=1e-12)
    np.testing.assert_allclose(samples[:,4],.519,atol=1e-12)
    assert np.corrcoef(samples[:,0],samples[:,2])[0,1]>.99
    assert all(d['dof']==len(q)-degree-1 for fit in record['responses'] for d in fit['response_diagnostics'])


def test_root_failure_causes_are_distinct_and_no_extrapolation():
    # Rows are coefficients for R(t) and m2(t).
    a=np.array([[.6,0.,0.],[.1,0.,0.]])
    assert r.quotient(a,np.array([[.6-.25,0.,1.],[.05,0.,0.]]),16,32,.8,.81)[1]=='multiple_roots_in_window'
    assert r.quotient(a,np.array([[.2,.1,0.],[.05,0.,0.]]),16,32,.8,.81)[1]=='no_root_in_window'
    assert r.quotient(a,np.array([[.6,.1,0.],[-.05,0.,0.]]),16,32,.8,.81)[1]=='nonphysical_root'
    assert r.quotient(a,a,16,32,.8,.81)[1]=='coincident_responses'


def test_rejected_draws_are_counted_instead_of_silently_dropped():
    q=np.linspace(.8005,.8015,5);t=np.linspace(-1,1,5);a=values(16,t,degree=1);b=values(32,t,degree=1)
    shifts=np.linspace(-.1,.1,100);shifts[:10]=2
    A=np.array([values(16,t,s,1+s/20,1) for s in shifts]);B=np.array([values(32,t,s,1+s/20,1) for s in shifts])
    record,samples=r.pair_fit(q,a,A,16,b,B,32,.8005,.8015,1)
    assert record['status']=='unresolved_local_quotient' and record['valid_fraction']==.9
    assert record['bootstrap_root_events']['no_root_in_window']==10
    assert np.isnan(samples[:10]).all()
