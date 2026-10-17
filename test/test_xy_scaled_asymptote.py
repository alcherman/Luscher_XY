"""Vortex-rich admissible Z4 starts and source-guarded replay."""
from pathlib import Path
import importlib.util,itertools,json,hashlib
import numpy as np,h5py,pytest
P=Path(__file__).resolve().parents[1]/'experiments/2026-10-17-xy-3d-negative-kappa/scaled-asymptote'
def load(name,file):
 s=importlib.util.spec_from_file_location(name,P/file);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
pt=load('scaled_pt_test','pt.py');aud=load('scaled_audit_test','audit.py')
def kw():return dict(N=6,epsilon=np.pi/2+.00625,kappas=[-16.,-8.],warmup=1,seed=1109,initial_winding=[0,0,0],initial_kind='z4_vortex',checkpoint_every=3)
def test_z4_plaquette_exact_vortex_probability():
 charges=[]
 for bits in itertools.product([0,1],repeat=4):
  angles=np.pi*np.array(bits)+np.pi/2*np.array([0,1,0,1]);d=(np.roll(angles,-1)-angles+np.pi)%(2*np.pi)-np.pi;charges.append(round(d.sum()/(2*np.pi)))
 assert charges.count(0)==12 and charges.count(-1)==2 and charges.count(1)==2
 assert np.abs(charges).mean()==.25
@pytest.mark.parametrize('epsilon',[np.pi/2,np.pi/2-.001])
def test_z4_strict_support_refuses_half_and_below_before_writes(tmp_path,epsilon):
 f=tmp_path/'bad.h5';k=dict(kw(),epsilon=epsilon)
 with pytest.raises(ValueError,match='strictly above'):pt.run(f,steps=3,**k)
 assert not f.exists()
def test_z4_initial_admissibility_rng_and_independent_full_audit(tmp_path):
 f=tmp_path/'z.h5';k=kw();pt.run(f,steps=6,**k);spec=dict(k,steps=6)
 with h5py.File(f) as h:
  phi=h['phi'][0,0,0];parity=np.indices(phi.shape).sum(0)%2
  assert h['staggered_nematic_squared'][0,0]==pytest.approx(1.)
  assert h['fourfold_squared'][0,0]==pytest.approx(1.)
  np.testing.assert_allclose(np.cos(2*phi),(-1.)**parity,atol=1e-15)
  for axis in range(3):np.testing.assert_allclose(abs((np.roll(phi,-1,axis)-phi+np.pi)%(2*np.pi)-np.pi),np.pi/2,atol=1e-15)
 assert aud.audit(f,spec,pt.sources())['passed']
def test_z4_two_slot_exact_restart_and_fieldhash_guard(tmp_path):
 f,g=tmp_path/'f.h5',tmp_path/'g.h5';k=kw();pt.run(f,steps=10,**k)
 with pytest.raises(RuntimeError,match='Injected crash'):pt.run(g,steps=10,inject_crash_after_flush=5,**k)
 pt.run(g,steps=10,resume=True,**k)
 with h5py.File(f) as a,h5py.File(g) as b:
  for name in a:
   if not name.endswith('_seconds'):np.testing.assert_array_equal(a[name][:],b[name][:])
 with h5py.File(g,'r+') as h:
  pars=json.loads(h.attrs['parameters']);pars['initial_field_sha256']='wrong';h.attrs['parameters']=json.dumps(pars)
 before=hashlib.sha256(g.read_bytes()).hexdigest()
 with pytest.raises(ValueError,match='mismatch'):pt.run(g,steps=12,resume=True,**k)
 assert hashlib.sha256(g.read_bytes()).hexdigest()==before

def test_derived_harmonics_checkerboard_and_independent_fft_oracles():
 root=P.parent/'asymptote'/'derive_harmonics.py';s=importlib.util.spec_from_file_location('derived_harmonic_test',root);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
 n=6;parity=np.indices((n,)*3).sum(0)%2;angle=1.61
 got=m.harmonics(np.array([angle*parity]));np.testing.assert_allclose(got[0],[np.sin(angle)**2,np.cos(2*angle)**2],atol=2e-15)
 rng=np.random.default_rng(443);phi=rng.uniform(0,2*np.pi,size=(3,n,n,n));got=m.harmonics(phi)
 for i,field in enumerate(phi):
  fft=np.fft.fftn(np.exp(2j*field))/n**3;four=np.fft.fftn(np.exp(4j*field))/n**3
  np.testing.assert_allclose(got[i],[abs(fft[n//2,n//2,n//2])**2,abs(four[0,0,0])**2],rtol=2e-13,atol=1e-15)
