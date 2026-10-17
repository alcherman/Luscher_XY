from pathlib import Path
import importlib.util,numpy as np,h5py,pytest
P=Path(__file__).resolve().parents[1]/'experiments/2026-10-17-xy-3d-negative-kappa/replica-exchange'
spec=importlib.util.spec_from_file_location('negative_pt',P/'pt.py');pt=importlib.util.module_from_spec(spec);spec.loader.exec_module(pt)

def test_joint_replay_and_uncommitted_tail(tmp_path):
 kw=dict(N=4,epsilon=2.1,kappas=[-.5,-.25,0.],warmup=2,seed=741,initial_winding=[1,0,0],checkpoint_every=5)
 a=tmp_path/'a.h5';b=tmp_path/'b.h5';pt.run(a,steps=20,**kw)
 with pytest.raises(RuntimeError,match='Injected crash'):pt.run(b,steps=20,inject_crash_after_flush=13,**kw)
 with h5py.File(b) as h:assert h.attrs['committed_rows']==11 and len(h['sweep'])==14
 pt.run(b,steps=20,resume=True,**kw)
 with h5py.File(a) as aa,h5py.File(b) as bb:
  for key in aa:
   if not key.endswith('_seconds'):np.testing.assert_array_equal(aa[key][:],bb[key][:],err_msg=key)
  assert np.all(aa['exchange_accepted'][:]<=aa['exchange_attempted'][:])

def test_swap_sign_and_ladder():
 assert pt.swap_log_ratio(-.5,0,2,3)==-.5
 assert pt.swap_log_ratio(-.5,0,3,2)==.5
 assert np.array_equal(pt.ladder(-.05),[-.05,-.025,0])

def test_zero_control_and_parameter_guard(tmp_path):
 kw=dict(N=4,epsilon=2.,kappas=[0.],warmup=1,seed=7,initial_winding=[0,0,0],checkpoint_every=2)
 f=tmp_path/'one.h5';pt.run(f,steps=4,**kw)
 with h5py.File(f) as h:assert h['exchange_accepted'].shape==(5,0)
 with pytest.raises(ValueError):pt.run(f,steps=6,resume=True,**(kw|dict(kappas=[-.1])))

def test_full_joint_audit_detects_swap_and_field_corruption(tmp_path):
 spec=importlib.util.spec_from_file_location('pt_audit',P/'audit.py');a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
 kw=dict(N=4,epsilon=2.1,kappas=[-.1,-.05,0.],warmup=1,seed=712,initial_winding=[0,0,0],checkpoint_every=2)
 f=tmp_path/'audit.h5';pt.run(f,steps=6,**kw);s=kw|dict(steps=6)
 assert a.audit(f,s,pt.sources())['passed']
 with h5py.File(f,'r+') as h:h['magnetization_squared'][3,1]+=.1
 with pytest.raises(AssertionError):a.audit(f,s,pt.sources())
