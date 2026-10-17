"""Analytic oracles for same-coupling crossing and spin-dimension quotients."""
import sys
from pathlib import Path
import numpy as np
import pytest

P=Path(__file__).resolve().parents[1]/'experiments/2026-10-15-xy-3d-spin-dimension'
sys.path.insert(0,str(P))
import common_xi_quotients as x


def curve(N,q,slope,Delta=.519,shift=0):
    # S0/F-1 is linear so raw-moment interpolation has an exact known crossing.
    c=(2*N*np.sin(np.pi/N))**2
    m2=N**(-2*Delta)*(1+.2*(q-.8))
    F=np.ones(len(q))
    S0=1+c*(.6**2+slope*(q-.8-shift))
    return np.column_stack([m2,1.5*m2**2,np.sqrt(m2),S0,F])


def test_ratio_of_means_and_lattice_normalization():
    q=np.array([.79,.8,.81]);a=curve(8,q,2)
    assert x.ratio(a[1],8)==pytest.approx(.6)
    assert not np.isclose(np.mean(x.ratio(a,8)),x.ratio(a.mean(0),8),atol=1e-8)


@pytest.mark.parametrize('method',['linear','pchip'])
def test_common_coupling_crossing_recovers_delta_and_propagates_root(method):
    q=np.array([.79,.8,.81]);a=curve(8,q,2);b=curve(16,q,5)
    result=x.crossing(q,a,8,q,b,16,.79,.81,method)
    assert result[0]==pytest.approx(.8,abs=1e-11)
    assert result[3]==pytest.approx(.6)
    assert result[4]==pytest.approx(.519,abs=1e-11)
    shifts=np.linspace(-.001,.001,51)
    B1=np.array([curve(8,q,2,shift=s) for s in shifts]);B2=np.array([curve(16,q,5,shift=s) for s in shifts])
    record,B=x.paired_record(q,a,B1,8,q,b,B2,16,method)
    assert record['valid_fraction']==1
    np.testing.assert_allclose(B[:,0],.8+shifts,atol=1e-11)
    np.testing.assert_allclose(B[:,4],.519,atol=1e-10)
    assert np.corrcoef(B[:,0],B[:,1])[0,1]>.999999


def test_no_extrapolation_or_silent_invalid_draw_removal():
    q=np.array([.79,.8,.81]);a=curve(8,q,2);b=curve(16,q,5)
    no=curve(16,q,5,shift=.02)
    assert np.isnan(x.crossing(q,a,8,q,no,16,.79,.81)).all()
    B1=np.tile(a,(100,1,1));B2=np.tile(b,(100,1,1));B2[:10]=no
    record,B=x.paired_record(q,a,B1,8,q,b,B2,16)
    assert record['status']=='unresolved_bootstrap_crossing'
    assert record['valid_fraction']==.9
    assert np.isnan(B[:10]).all()


def test_shared_pair_covariance_and_joint_valid_mask():
    rng=np.random.default_rng(940);N=[6,8,12,16];D=.52;omega=.789
    # Adjacent quotients share a size: analytic correlated oracle.
    e=rng.normal(0,.003,(2000,5));noise=e[:,1:]-e[:,:-1]
    y=D+.1*np.array(N,dtype=float)**(-omega)
    B=y+noise
    records=[dict(N=n,status='conditional_quotient',Delta_eff=t) for n,t in zip(N,y)]
    draws=[np.column_stack([np.ones((2000,4)),B[:,i]]) for i in range(4)]
    f=x.correction_fit(records,draws,6,omega)
    C=np.cov(B,rowvar=False);X=np.column_stack([np.ones(4),np.array(N,dtype=float)**(-omega)])
    T=np.linalg.solve(X.T@np.linalg.solve(C,X),X.T@np.linalg.inv(C))
    assert f['Delta']==pytest.approx(D,abs=1e-12)
    assert f['SD']==pytest.approx((B@T.T)[:,0].std(ddof=1))
    draws[0][:100,4]=np.nan
    assert x.correction_fit(records,draws,6)['status']=='unresolved_joint_bootstrap'
    assert x.correction_fit(records,draws,8)['status']=='conditional_correction_fit'
    assert x.correction_fit(records,draws,12)['status']=='zero_dof_diagnostic'


@pytest.mark.parametrize('method',['linear','pchip'])
def test_two_interior_roots_with_same_sign_endpoints_are_not_hidden(method):
    q=np.array([0.,1.]);N1,N2=8,16;c1=(2*N1*np.sin(np.pi/N1))**2;c2=(2*N2*np.sin(np.pi/N2))**2
    F1=1+q;F2=1+2*q
    S1=(1+5*c1)*F1-35*c1/16;S2=(1+4.5*c2)*F2-15*c2/8
    a=np.column_stack([S1/N1**3,np.ones((2,2)),S1,F1]);b=np.column_stack([S2/N2**3,np.ones((2,2)),S2,F2])
    assert np.isnan(x.crossing(q,a,N1,q,b,N2,0.,1.,method)).all()
    assert x.crossing(q,a,N1,q,b,N2,0.,.5,method)[0]==pytest.approx(.25)


def test_roundoff_merge_does_not_round_actual_scan_coordinates():
    q=.806;ulps=np.spacing(q)
    np.testing.assert_equal(x.merge_roundoff([q,q-ulps,q+1e-8]),[q-ulps,q+1e-8])
