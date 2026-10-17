"""Independent arc-membership and disconnected-conditional distribution oracles."""
from pathlib import Path
import importlib.util
import numpy as np,pytest
from scipy.integrate import quad
P=Path(__file__).resolve().parents[1]/'experiments/2026-10-17-xy-3d-negative-kappa/heatbath-resolution'
s=importlib.util.spec_from_file_location('independent_heatbath',P/'heatbath.py');h=importlib.util.module_from_spec(s);s.loader.exec_module(h)

def test_allowed_union_against_independent_distance_membership():
 rng=np.random.default_rng(4131);grid=(np.arange(10001)+.371)*2*np.pi/10001
 for eps in [.505*np.pi,.8*np.pi]:
  for _ in range(12):
   neighbors=rng.uniform(0,2*np.pi,size=6);lo,hi=h.intervals(neighbors,eps)
   direct=np.all(abs(np.angle(np.exp(1j*(grid[:,None]-neighbors[None,:]))))<eps,axis=1)
   actual=np.array([any(a<v<b for a,b in zip(lo,hi))for v in grid]);np.testing.assert_array_equal(actual,direct)
   if len(lo):
    for k in [-9.6,2.]:
     A,B,M=h.coefficients(neighbors,k,lo,hi);values=np.concatenate([np.linspace(a,b,251)for a,b in zip(lo,hi)]);assert np.max(A*np.cos(values)+B*np.sin(values))<=M

@pytest.mark.parametrize('kappa',[-8.,2.])
def test_disconnected_conditional_against_independent_quadrature(kappa):
 # Two antipodal neighborhoods leave disconnected windows about +/-pi/2.
 neighbors=np.array([0.,0.,.03,np.pi,np.pi+.02,np.pi-.01]);eps=.54*np.pi;lo,hi=h.intervals(neighbors,eps);assert len(lo)==2
 A=kappa*np.cos(neighbors).sum();B=kappa*np.sin(neighbors).sum();weight=lambda x:np.exp(A*np.cos(x)+B*np.sin(x))
 Z=sum(quad(weight,a,b,epsabs=1e-12)[0]for a,b in zip(lo,hi));expected=sum(quad(lambda x:weight(x)*np.cos(2*x),a,b,epsabs=1e-12)[0]for a,b in zip(lo,hi))/Z
 mass0=quad(weight,lo[0],hi[0],epsabs=1e-12)[0]/Z
 rng=np.random.default_rng(71343);draws=np.array([h.conditional(neighbors,kappa,eps,rng)[0]for _ in range(16000)])
 assert abs(np.mean(np.cos(2*draws))-expected)<6*np.std(np.cos(2*draws),ddof=1)/np.sqrt(len(draws))
 observed=np.mean((draws>lo[0])&(draws<hi[0]));assert abs(observed-mass0)<6*np.sqrt(mass0*(1-mass0)/len(draws))
 assert np.all(np.max(abs(np.angle(np.exp(1j*(draws[:,None]-neighbors)))),axis=1)<eps)
