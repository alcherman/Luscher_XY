"""Oracles for operator dimension, contour covariance and fit identifiability."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest
P=Path(__file__).resolve().parents[1]/'experiments/2026-10-15-xy-3d-spin-dimension/measure.py'
spec=importlib.util.spec_from_file_location('spin_dimension',P);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def test_exact_power_law_recovers_operator_dimension_not_twice_it():
    rng=np.random.default_rng(489);N=np.array([6,8,12,16,24]);Delta=.519;center=.8*N**(-2*Delta)
    B=center*np.exp(rng.normal(0,.003,(3000,len(N))))
    f=m.fit_dimension(N,center,B,2,8)
    assert f['Delta']==pytest.approx(Delta,abs=1e-12)
    assert f['MC95'][0]<Delta<f['MC95'][1]
    f1=m.fit_dimension(N,np.sqrt(center),np.sqrt(B),1,8)
    assert f1['Delta']==pytest.approx(f['Delta'],abs=1e-12)


def test_contour_root_and_magnetization_uncertainty_are_joint():
    x=np.array([.4,.5,.6]);U=np.array([.2,.25,.3]);m2=np.array([.1,.11,.12]);ma=np.sqrt(m2)*.9
    means=np.column_stack([m2,2*(1-U)*m2**2,ma]);shifts=np.linspace(-.006,.006,1001)
    B=np.broadcast_to(means,(len(shifts),)+means.shape).copy();B[:,:,1]=2*(1-U[None,:]-shifts[:,None])*m2[None,:]**2
    r,b=m.interpolate_record(x,means,B,.25,'linear')
    assert r['valid_fraction']==1
    assert np.corrcoef(b[:,0],b[:,1])[0,1]>.99999
    assert b[:,1].std()>0
    np.testing.assert_allclose(b[:,1],.11-.2*shifts,atol=1e-12)


def test_missing_bracket_is_not_extrapolated():
    x=np.array([.4,.5,.6]);m2=np.array([.1,.11,.12]);means=np.column_stack([m2,2*(1-np.array([.05,.1,.15]))*m2**2,np.sqrt(m2)])
    r,b=m.interpolate_record(x,means,np.tile(means,(200,1,1)))
    assert r['status']=='unbracketed' and b is None


def test_size_cut_recomputes_valid_mask_and_zero_dof_is_diagnostic():
    rng=np.random.default_rng(932);N=np.array([6,8,12,16]);y=N**(-1.04);B=y*np.exp(rng.normal(0,.005,(1000,4)));B[:100,0]=np.nan
    assert m.fit_dimension(N,y,B,2,6)['status']=='unresolved_joint_bootstrap'
    assert m.fit_dimension(N,y,B,2,8)['status']=='conditional_dimension_fit'
    assert m.fit_dimension(N,y,B,2,12)['status']=='zero_dof_diagnostic'
    assert m.fit_dimension(N,y,B,2,8,.8)['status']=='zero_dof_diagnostic'


def test_moment_consistent_contour_and_binder_interpolation_differ():
    x=np.array([.4,.5,.6]);m2=np.array([.05,.15,.35]);U=np.array([.1,.35,.48]);means=np.column_stack([m2,2*(1-U)*m2**2,np.sqrt(m2)*.9]);B=np.tile(means,(100,1,1))
    linear,_=m.interpolate_record(x,means,B,.25,'linear')
    moments,_=m.interpolate_record(x,means,B,.25,'moment_linear')
    assert moments['moment_ratio_at_interpolated_root']==pytest.approx(.25,abs=1e-10)
    assert abs(linear['moment_ratio_at_interpolated_root']-.25)>.01


def test_zero_beta_epsilon_scan_reverses_orientation_and_interval():
    import sys
    sys.path.insert(0,str(P.parent))
    import zero_beta as z
    q=np.array([.79,.80,.81,.82]);U=np.array([.4,.3,.2,.1]);m2=np.array([.15,.12,.09,.06])
    means=np.column_stack([m2,2*(1-U)*m2**2,.9*np.sqrt(m2)])
    B=np.broadcast_to(means,(101,)+means.shape).copy()
    shift=np.linspace(-.005,.005,101)
    B[:,:,1]=2*(1-U[None,:]-shift[:,None])*m2[None,:]**2
    r,b=z.endpoint_contour(q,means,B,.25,'linear','full')
    assert r['epsilon_over_pi']==pytest.approx(.805)
    assert r['epsilon_MC95'][0]<.805<r['epsilon_MC95'][1]
    assert r['valid_fraction']==1
    assert np.corrcoef(-b[:,0],b[:,1])[0,1]<-.999


def test_zero_beta_central_root_without_bootstrap_interval(monkeypatch):
    import sys
    sys.path.insert(0,str(P.parent))
    import zero_beta as z
    monkeypatch.setattr(z.m,'interpolate_record',lambda *args: (dict(location=[-.8,.1,.3],MC95=None,selected_indices=[0,1]),None))
    r,b=z.endpoint_contour(np.array([.79,.81]),np.ones((2,3)),np.ones((5,2,3)))
    assert r['epsilon_over_pi']==.8 and 'epsilon_MC95' not in r


def test_larger_pilot_topology_streaming_and_event_denominators(tmp_path,monkeypatch):
    import sys,json,h5py
    sys.path.insert(0,str(P.parent))
    import diagnose_larger_pilot as d
    out=tmp_path/'study';campaign=tmp_path/'campaign';out.mkdir();(campaign/'chains').mkdir(parents=True);(campaign/'jobs').mkdir()
    monkeypatch.setattr(d.m,'OUT',out);monkeypatch.setattr(d.m,'CAMPAIGN',campaign)
    point=dict(id='toy',N=48,epsilon_over_pi=.8015,replica=30,initial_winding=[0,0,0]);rows=10001;t=np.arange(rows);signal=(t%2).astype(float);accepted=np.where(t%5==0,2,0)
    (out/'larger-zero-beta-pilot-design.json').write_text(json.dumps(dict(points=[point])))
    rng=np.random.default_rng(81);m2=.1+.005*rng.normal(size=rows)
    with h5py.File(campaign/'chains/toy.h5','w') as h:
        h['sweep']=t;h['winding']=np.broadcast_to(signal[:,None,None],(rows,3,4)).astype(int);h['winding_accepted']=accepted;h['vortex_density']=.01*rng.random(rows);h['magnetization_squared']=m2;h['magnetization_fourth']=1.5*m2**2;h['magnetization_abs']=np.sqrt(m2)*.9
    (campaign/'jobs/toy.json').write_text(json.dumps(dict(state='complete',audit=dict(passed=True,bytes=10),attempt_generation_seconds=2,audit_seconds=1)))
    d.main();result=json.loads((out/'larger-pilot-diagnostics.json').read_text())['records'][0]
    for r in result['cut_diagnostics']:
        cut=r['cut'];assert r['rows']==rows-cut
        assert r['cycles_with_accepted_winding_update']==int((accepted[cut:]>0).sum())
        assert r['accepted_winding_updates']==int(accepted[cut:].sum())
        assert r['mean_winding_squared']==pytest.approx(signal[cut:].mean())
    with np.load(out/'larger-pilot-derived-histories/toy.npz') as z:
        np.testing.assert_array_equal(z['directional_mean_winding'][:,0],signal)
