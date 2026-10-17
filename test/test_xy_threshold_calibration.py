"""Check threshold bracketing and the infinite-size intercept separately."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

PATH=Path(__file__).resolve().parents[1]/'experiments/2026-10-13-xy-infinite-volume/threshold_calibration.py'
SPEC=importlib.util.spec_from_file_location('xy_threshold_calibration',PATH)
analysis=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(analysis)


@pytest.mark.parametrize('method',['linear','pchip'])
def test_threshold_crossing_is_bracketed_and_does_not_extrapolate(method):
    beta=np.array([.1,.2,.3,.4]);values=2*beta-.25
    assert analysis.crossing(beta,values,method=method)==pytest.approx(.25)
    assert np.isnan(analysis.crossing(beta,[.3,.35,.4,.45],method=method))
    assert np.isnan(analysis.crossing([.3],[.4],method=method))
    assert np.isnan(analysis.crossing(beta,[.1,.4,.1,.4],method=method))


@pytest.mark.parametrize('D',[3,4])
def test_size_fit_recovers_known_intercept_and_identifies_two_size_limit(D):
    N=np.array([6,8,12,16]);x=analysis.coordinate(D,N)
    values=.3017-.15*x
    boot=values+np.random.default_rng(1234).normal(0,.0001,size=(1000,4))
    result=analysis.size_fit(D,N,values,boot)
    assert result['beta_infinite']==pytest.approx(.3017,abs=1e-10)
    assert result['bootstrap_95'][0]<.3017<result['bootstrap_95'][1]
    two=analysis.size_fit(D,N,values,boot,minimum=12)
    assert two['df']==0 and two['status']=='zero_dof_provisional'
    assert analysis.size_fit(D,N,values,boot,minimum=12,correction=True)['status']=='unidentifiable'


def test_weighted_monotone_crossing_handles_noise_without_censoring_multiple_roots():
    beta=np.array([.1,.2,.3,.4]);values=[0,.3,.2,.5]
    np.testing.assert_allclose(analysis.monotone_values(values),[0,.25,.25,.5])
    np.testing.assert_allclose(analysis.monotone_values(values,[1,1,3,1]),[0,.225,.225,.5])
    assert analysis.monotone_crossing(beta,values)==pytest.approx(.25)
    assert np.isnan(analysis.monotone_crossing(beta,[.1,.12,.15,.2]))
