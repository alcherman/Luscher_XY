"""Independent full collective-orbit balance and disabled-extension identity."""
from pathlib import Path
import importlib.util,itertools
import h5py,numpy as np,pytest
from supervillain import XY,Lattice
P=Path(__file__).resolve().parents[1]/'experiments/2026-10-17-xy-3d-negative-kappa/sign-sector'
def load(name,p):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
m=load('independent_sector_moves',P/'sector_moves.py')
@pytest.mark.parametrize('kappa',[-2.,1.])
def test_collective_complete_finite_orbit_stationarity(kappa):
 n=4;x=np.arange(n);parity=x%2;eps=.8*np.pi;action=XY(Lattice(1,n),kappa,eps)
 keys=sorted({tuple(((a*parity+b*x)%4).tolist())for a,b in itertools.product(range(4),repeat=2)})
 keys=[key for key in keys if np.max(abs((np.roll(np.array(key)*np.pi/2,-1)-np.array(key)*np.pi/2+np.pi)%(2*np.pi)-np.pi))<eps]
 assert len(keys)>2
 fields=[np.array(k)[None]*np.pi/2 for k in keys];C=np.array([np.cos(np.roll(v[0],-1)-v[0]).sum()for v in fields]);T=np.zeros((len(keys),len(keys)))
 class Tape:
  def __init__(self,index,u):self.index=index;self.u=u
  def integers(self,*a,**kw):return self.index
  def uniform(self,*a,**kw):return self.u
 choices=m.shift_choices(1)
 for i,field in enumerate(fields):
  for index,(a,w)in enumerate(choices):
   # Independent integer-angle construction of the proposed group action.
   candidate=(field[0]+a*parity-2*np.pi*w[0]*x/n)%(2*np.pi);key=tuple((np.rint(candidate/(np.pi/2)).astype(int)%4).tolist())
   if key not in keys:
    got,info=m.collective(action,{'phi':field},Tape(index,.5));assert not info['accepted'];T[i,i]+=1/len(choices);continue
   j=keys.index(key);accept=min(1.,np.exp(kappa*(C[j]-C[i])))
   for prob,u in [(accept,accept/2),(1-accept,(1+accept)/2)]:
    if prob==0:continue
    got,info=m.collective(action,{'phi':field},Tape(index,u));out=tuple((np.rint(got['phi'][0]/(np.pi/2)).astype(int)%4).tolist());T[i,keys.index(out)]+=prob/len(choices)
 weight=np.exp(kappa*C-np.max(kappa*C));weight/=weight.sum();np.testing.assert_allclose(T.sum(1),1,atol=1e-14);np.testing.assert_allclose(weight[:,None]*T,(weight[:,None]*T).T,atol=1e-14)

def test_disabled_extension_exact_native_trajectory(tmp_path):
 new=load('independent_sector_pt',P/'pt.py');old=load('independent_prior_scaled_pt',P.parent/'scaled-asymptote/pt.py')
 kwargs=dict(N=4,epsilon=.54*np.pi,kappas=[-2.],steps=9,warmup=2,seed=928374,initial_winding=[1,0,0],initial_kind='winding',checkpoint_every=3)
 a,b=tmp_path/'old.h5',tmp_path/'new.h5';old.run(a,**kwargs);new.run(b,sector_updates=False,**kwargs)
 with h5py.File(a)as x,h5py.File(b)as y:
  for name in x:
   if name.endswith('_seconds')or name=='rng_state':continue
   np.testing.assert_array_equal(x[name][:],y[name][:])
