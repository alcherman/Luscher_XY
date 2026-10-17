"""Independent support, quadrature and exact trajectory continuation checks."""
from pathlib import Path
import importlib.util
import numpy as np, h5py, pytest
from scipy.integrate import quad
P=Path(__file__).resolve().parents[1]/'experiments/2026-10-17-xy-3d-negative-kappa/heatbath-resolution'
def load(name,f):
 s=importlib.util.spec_from_file_location(name,P/f);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
hb=load('hb_test','heatbath.py');pt=load('hb_producer_test','pt.py');audit=load('hb_audit_test','audit.py')
def wrap(x):return (x+np.pi)%(2*np.pi)-np.pi

def test_allowed_union_against_independent_grid():
 rng=np.random.default_rng(174)
 for eps,neigh in [(1.7,np.array([0,np.pi]*3)),(.8,np.array([6.2,.1,.2,6.1,0,.05]))]+[(eps,rng.uniform(-eps*.9,eps*.9,6)) for eps in [.3,1.58,2.,3.]]:
  lo,hi=hb.intervals(neigh,eps);grid=np.linspace(0,2*np.pi,20001,endpoint=False)+.000001
  actual=np.any((grid[:,None]>lo)&(grid[:,None]<hi),axis=1)
  expected=np.all(np.abs(wrap(grid[:,None]-neigh))<eps,axis=1)
  np.testing.assert_array_equal(actual,expected)
  assert sum(hi-lo)>0

def test_majorant_and_rng_mutation():
 rng=np.random.default_rng(43)
 for k in [-8,2]:
  neigh=np.array([0,np.pi]*3)+rng.uniform(-.003,.003,6);eps=.505*np.pi;lo,hi=hb.intervals(neigh,eps);A,B,M=hb.coefficients(neigh,k,lo,hi)
  for a,b in zip(lo,hi):assert np.max(A*np.cos(np.linspace(a,b,1001))+B*np.sin(np.linspace(a,b,1001)))<=M
  before=repr(rng.bit_generator.state);hb.conditional(neigh,k,eps,rng);assert repr(rng.bit_generator.state)!=before

@pytest.mark.parametrize('kappa',[-8.,2.])
def test_conditional_moments_against_quadrature(kappa):
 neighbors=np.array([0,np.pi,0,np.pi,0,np.pi])+.003*np.arange(6);eps=.505*np.pi;lo,hi=hb.intervals(neighbors,eps);A,B,M=hb.coefficients(neighbors,kappa,lo,hi)
 weight=lambda x:np.exp(A*np.cos(x)+B*np.sin(x)-M)
 Z=sum(quad(weight,a,b,epsabs=1e-13)[0] for a,b in zip(lo,hi));rng=np.random.default_rng(452);samples=np.array([hb.conditional(neighbors,kappa,eps,rng)[0] for _ in range(12000)])
 assert np.all(abs(wrap(samples[:,None]-neighbors))<eps)
 for fun in [np.cos,np.sin,lambda x:np.cos(2*x)]:
  expected=sum(quad(lambda x:fun(x)*weight(x),a,b,epsabs=1e-13)[0] for a,b in zip(lo,hi))/Z;values=fun(samples)
  assert abs(values.mean()-expected)<6*values.std()/np.sqrt(len(values))+1e-5

def test_full_fields_exact_crash_replay_and_audit(tmp_path):
 k=dict(N=4,epsilon=.505*np.pi,kappas=[-8.,-7.],warmup=2,seed=171414,initial_winding=[1,0,0],initial_kind='winding',checkpoint_every=3)
 a,b=tmp_path/'a.h5',tmp_path/'b.h5';pt.run(a,steps=12,**k)
 with pytest.raises(RuntimeError,match='Injected crash'):pt.run(b,steps=12,inject_crash_after_flush=7,**k)
 pt.run(b,steps=12,resume=True,**k)
 with h5py.File(a) as x,h5py.File(b) as y:
  for name in x:
   if not name.endswith('_seconds'):np.testing.assert_array_equal(x[name][:],y[name][:])
 assert audit.audit(a,dict(k,steps=12),pt.sources())['passed']
 with h5py.File(b,'r+') as h:h['heatbath_draws'][3,0]=0
 with pytest.raises(AssertionError):audit.audit(b,dict(k,steps=12),pt.sources())
