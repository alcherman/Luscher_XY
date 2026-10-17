import sys
from pathlib import Path
import numpy as np
import pytest
from scipy.integrate import quad
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments/2026-10-13-xy-infinite-volume'))
from d5_scaling import fit

def test_D5_leading_fit_and_correlated_intercept_oracle():
    N=np.array([4,6,8,10]);X=np.column_stack([np.ones(4),N**(-2.5)]);y=X@np.array([.23,-.1])
    rng=np.random.default_rng(8);B=y+rng.multivariate_normal(np.zeros(4),1e-6*(np.eye(4)+.3*np.ones((4,4))),size=4000)
    f=fit(N,y,B);assert f['location_infinite']==pytest.approx(.23)
    C=np.cov(B,rowvar=False);T=np.linalg.solve(X.T@np.linalg.solve(C,X),X.T@np.linalg.inv(C))
    np.testing.assert_allclose(f['samples'],B@T[0],atol=1e-12)
    assert f['dof']==2

def test_D5_correction_three_sizes_is_diagnostic():
    N=np.array([4.,6.,8.]);y=.23-.1*N**(-2.5)+.05*N**(-3)
    B=y+np.random.default_rng(4).normal(0,.001,(2000,3));f=fit(N,y,B,correction=True)
    assert f['status']=='zero_dof_diagnostic';assert f['p'] is None
    assert f['location_infinite']==pytest.approx(.23)

def test_D5_common_missing_draw_mask():
    N=np.array([4,6,8]);y=.23-.1*N**(-2.5);B=y+np.random.default_rng(4).normal(0,.001,(2000,3));B[:30,1]=np.nan
    assert fit(N,y,B)['status']=='bootstrap_coverage_unresolved'

def test_O2_quartic_zero_mode_Binder():
    integrals=[quad(lambda r:r**(k+1)*np.exp(-r**4),0,np.inf)[0] for k in [0,2,4]]
    U=1-integrals[2]*integrals[0]/(2*integrals[1]**2)
    assert U==pytest.approx(1-np.pi/4)

def test_D5_local_bracket_retains_reliable_anchors_and_epsilon_orientation():
    from analyze_d5_calibration import local_indices,crossing_record
    x=np.array([.86,.88,.89,.90,.92]);y=.25-5*(x-.891)
    B=y+np.random.default_rng(13).normal(0,.001,(2000,5))
    np.testing.assert_array_equal(local_indices(y,B,.25,-1),[2,3])
    r,draws=crossing_record(x,y,B,sign=-1)
    assert r['location']==pytest.approx(.891);assert r['valid_fraction']==1
    assert np.std(draws)>0
    assert r['selected_window']==[.89,.90]

def test_D5_weak_bracket_uses_neighbors_and_does_not_hide_missing_tails():
    from analyze_d5_calibration import local_indices,crossing_record
    x=np.array([.1,.2,.3,.4]);y=np.array([.24,.249,.251,.26])
    B=np.tile(y,(2000,1));B[::2]+=.004;B[1::2]-=.004
    np.testing.assert_array_equal(local_indices(y,B,.25,1),[0,1,2,3])
    r,draws=crossing_record(x[:2],y[:2],B[:,:2])
    assert r['status']=='unresolved';assert r['valid_fraction']==.5
    assert np.isnan(draws).sum()==1000

def test_D5_linear_and_PCHIP_are_distinct_local_interpolation_checks():
    from analyze_d5_calibration import crossing_record
    from scipy.interpolate import PchipInterpolator
    from scipy.optimize import brentq
    x=np.array([0.,.1,.2,.3]);y=np.array([.05,.20,.28,.40])
    B=y+np.random.default_rng(18).normal(0,.0001,(2000,4))
    linear,_=crossing_record(x,y,B,window='full')
    cubic,_=crossing_record(x,y,B,method='pchip',window='full')
    assert linear['location']==pytest.approx(.1625)
    oracle=brentq(lambda a:PchipInterpolator(x,y)(a)-.25,.1,.2)
    assert cubic['location']==pytest.approx(oracle)
    assert abs(cubic['location']-linear['location'])>.0001

def test_D5_half_shift_diagnostic_keeps_absolute_values_and_flags():
    from analyze_d5_calibration import half_diagnostics
    rng=np.random.default_rng(28);a=np.column_stack([np.ones(800),np.r_[rng.normal(1.4,.005,400),rng.normal(1.6,.005,400)]])
    rows=half_diagnostics([a],[{'id':'synthetic_shift'}],0,{'diagnostics':[{'block_width':20}]})
    assert rows[0]['flag']=='Binder_half_shift';assert rows[0]['blocks']==[20,20]
    assert rows[0]['Binder_halves'][0]==pytest.approx(.30,abs=.001)
    assert rows[0]['Binder_halves'][1]==pytest.approx(.20,abs=.001)


def test_D5_inconsistent_size_ansatz_is_diagnostic_not_accepted():
    N=np.array([6.,8.,10.]);y=.23-.1*N**(-2.5);y[-1]+=.01
    B=y+np.random.default_rng(91).normal(0,.0001,(2000,3))
    result=fit(N,y,B)
    assert result['status']=='rejected_size_fit';assert result['p']<.01
    assert result['dof']==1;assert 'MC95' in result


def test_D5_zero_beta_reuse_is_explicit_exact_coordinate_not_rounding():
    from analyze_d5_calibration import exact_zero_beta_fixed_scan_anchor
    p=dict(campaign_role='zero_beta_endpoint',beta=0,epsilon_over_pi=.90,epsilon_sampling_milli=900)
    assert exact_zero_beta_fixed_scan_anchor(p)
    for q,milli in [(.894,894),(.896,896),(.898,898),(.9001,900),(.89,890)]:
        assert not exact_zero_beta_fixed_scan_anchor(dict(p,epsilon_over_pi=q,epsilon_sampling_milli=milli))
    assert not exact_zero_beta_fixed_scan_anchor(dict(p,beta=.001))
    assert not exact_zero_beta_fixed_scan_anchor({k:v for k,v in p.items() if k!='epsilon_sampling_milli'})


def test_D5_atomic_cache_failure_preserves_previous_complete_data(tmp_path,monkeypatch):
    import analyze_d5_calibration as analysis
    path=tmp_path/'bootstrap.npz';old=np.arange(12,dtype=float).reshape(4,3)
    analysis.atomic_bootstrap_cache(path,old)
    replace=analysis.os.replace
    def failure(*args):raise RuntimeError('injected interruption before replace')
    with monkeypatch.context() as m:
        m.setattr(analysis.os,'replace',failure)
        with pytest.raises(RuntimeError):analysis.atomic_bootstrap_cache(path,old+1)
    np.testing.assert_array_equal(np.load(path)['samples'],old)
    assert list(tmp_path.glob('*.tmp.npz'))
    analysis.atomic_bootstrap_cache(path,old+2)
    np.testing.assert_array_equal(np.load(path)['samples'],old+2)
