from pathlib import Path
import importlib.util
import numpy as np
P=Path(__file__).resolve().parents[1]/'experiments/2026-10-17-xy-3d-negative-kappa/boundary-resolution/resolve_brackets.py'
s=importlib.util.spec_from_file_location('negative_bracket_oracles',P);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def test_unique_missing_multiple_and_exact_root_ledger():
 assert m.roots(np.array([-3.,-2.,-1.]),np.array([-.5,0.,.5]),0.)==[dict(root=-2.,bracket=[-2.,-2.])]
 assert not m.roots(np.array([-2.,-1.]),np.array([.2,.3]),.5)
 assert len(m.roots(np.array([-3.,-2.,-1.]),np.array([1.,-1.,1.])))==2
 assert len(m.roots(np.array([-3.,-2.,-1.]),np.array([1.,np.nan,-1.])))==0

def test_local_linear_and_pchip_agree_on_analytic_linear_curve():
 x=np.array([-3.,-2.,-1.]);y=2*x+4.3
 for kind in ['linear','pchip']:
  r=m.roots(x,y,0.,kind);np.testing.assert_allclose(r[0]['root'],-2.15,atol=1e-12);assert r[0]['bracket']==[-3.,-2.]

def test_moments_pool_before_nonlinear_binder_and_correlation_ratio():
 n=8;a=np.array([[.1,.015,20.48,0.,0.,0.,0.],[.3,.095,20.48,0.,0.,0.,0.]])
 got=m.moments(a,n);np.testing.assert_allclose(got['Binder_U'],.3125);np.testing.assert_allclose(got['R'],1/(n*np.sin(np.pi/n)))
 assert abs(got['Binder_U']-np.mean([m.moments(v[None],n)['Binder_U'] for v in a]))>.04
 bad=a.copy();bad[:,2]=n**3;assert m.moments(bad,n)['R'] is None

def test_hysteretic_switch_counts_do_not_count_intermediate_frames_as_switches():
 a=m.phase_visits(np.array([.01,.08,.2,.1,.18,.09,.02,.03]));assert a['low_frames']==3 and a['high_frames']==2 and a['intermediate_frames']==3
 assert a['complete_switches']==2;assert [d['postcut_index'] for d in a['switch_events']]==[2,6]

def test_nonfinite_serialization_is_explicit():
 assert m.clean([np.float64(np.nan),np.float64(np.inf),1.])==[None,None,1.]
