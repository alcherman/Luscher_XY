"""Independent initial-state observables, exact replay, and full-field audits."""
from pathlib import Path
import importlib.util,json,hashlib
import h5py,numpy as np,pytest
P=Path(__file__).resolve().parents[1]/'experiments/2026-10-17-xy-3d-negative-kappa/asymptote'
def load(name,file):
 s=importlib.util.spec_from_file_location(name,P/file);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
pt=load('test_asymptote_pt','pt.py');aud=load('test_asymptote_aud','audit.py')
def kw(kind='checkerboard'):
 return dict(N=6,epsilon=.52*np.pi,kappas=[-8.],warmup=2,seed=718,initial_winding=[0,0,0],initial_kind=kind,checkpoint_every=3)
def spec(k,steps):return dict(k,N=k['N'],steps=steps)
def test_checkerboard_analytic_observables_and_full_audit(tmp_path):
 f=tmp_path/'c.h5';k=kw();pt.run(f,steps=6,**k);a=.98*k['epsilon']
 with h5py.File(f) as h:
  assert h['magnetization_squared'][0,0]==pytest.approx(np.cos(a/2)**2)
  assert h['staggered_squared'][0,0]==pytest.approx(np.sin(a/2)**2)
  assert h['nematic_squared'][0,0]==pytest.approx(np.cos(a)**2)
  assert h['mean_link_cosine'][0,0]==pytest.approx(np.cos(a))
  assert h['vortex_density'][0,0]==0
 assert aud.audit(f,spec(k,6),pt.sources())['passed']
 with h5py.File(f,'r+') as h:h['staggered_squared'][3,0]+=.1
 with pytest.raises(AssertionError):aud.audit(f,spec(k,6),pt.sources())
def test_joint_checkerboard_checkpoint_replay_and_kind_refusal(tmp_path):
 a,b=tmp_path/'a.h5',tmp_path/'b.h5';k=kw();pt.run(a,steps=12,**k)
 with pytest.raises(RuntimeError,match='Injected crash'):pt.run(b,steps=12,inject_crash_after_flush=8,**k)
 pt.run(b,steps=12,resume=True,**k)
 with h5py.File(a) as x,h5py.File(b) as y:
  for name in x:
   if not name.endswith('_seconds'):np.testing.assert_array_equal(x[name][:],y[name][:])
 old=hashlib.sha256(b.read_bytes()).hexdigest();bad=dict(k,initial_kind='winding')
 with pytest.raises(ValueError,match='mismatch'):pt.run(b,steps=15,resume=True,**bad)
 assert hashlib.sha256(b.read_bytes()).hexdigest()==old
@pytest.mark.parametrize('q,w',[(.5,1),(.52,1),(.6,1)])
def test_N6_axial_winding_valid_below_N4_geometry_threshold(tmp_path,q,w):
 k=dict(kw('winding'),epsilon=q*np.pi,initial_winding=[w,w,w]);f=tmp_path/f'q{q}.h5';pt.run(f,steps=3,**k)
 assert aud.audit(f,spec(k,3),pt.sources())['passed']
def test_N8_near_saturation_winding_is_strictly_admissible(tmp_path):
 k=dict(kw('winding'),N=8,initial_winding=[2,2,2]);f=tmp_path/'n8.h5';pt.run(f,steps=3,**k);assert aud.audit(f,spec(k,3),pt.sources())['passed']
