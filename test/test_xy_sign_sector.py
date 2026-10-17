"""Independent balance, supplied randomness, complete audit and replay oracles."""
from pathlib import Path
import importlib.util,itertools,json,hashlib
import numpy as np,h5py,pytest
P=Path(__file__).resolve().parents[1]/'experiments/2026-10-17-xy-3d-negative-kappa/sign-sector'
def load(name,f):
 s=importlib.util.spec_from_file_location(name,P/f);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
m=load('sector_test_moves','sector_moves.py');pt=load('sector_test_producer','pt.py');audit=load('sector_test_audit','audit.py')
def wrapped(x):return (x+np.pi)%(2*np.pi)-np.pi

def test_collective_inverse_pairing_and_unit_escape():
 n=8;x=np.indices((n,)*3);phi=(2*np.pi*x[0]/n)[None];choices=m.shift_choices(3)
 assert len(choices)==14
 for j,(a,w) in enumerate(choices):
  inv=next(i for i,(b,v) in enumerate(choices) if b==-a and np.array_equal(v,-w))
  np.testing.assert_allclose(wrapped(m.shift_candidate(m.shift_candidate(phi,j),inv)-phi),0,atol=4e-15)
 j=next(i for i,(a,w) in enumerate(choices) if a==np.pi/2 and np.array_equal(w,[1,0,0]));y=m.shift_candidate(phi,j)
 assert max(abs(wrapped(np.roll(y,-1,axis=d+1)-y)).max() for d in range(3))<.52*np.pi
 assert -2*(m.cosine_sum(y)-m.cosine_sum(phi))>1000
 for k in [-2.,1.]:
  lx=k*m.cosine_sum(phi);ly=k*m.cosine_sum(y)
  np.testing.assert_allclose(lx+min(0,ly-lx),ly+min(0,lx-ly),atol=1e-12)

@pytest.mark.parametrize('kappa',[-2.,2.])
@pytest.mark.parametrize('hard',[False,True])
def test_pi_cluster_full_16state_transition_matrix_detailed_balance(kappa,hard):
 # Independent exact enumeration of all eight positive bonds, including parallel N=2 edges.
 n=2;coords=np.indices((n,n));flat=np.arange(n*n).reshape(n,n);plus=np.array([np.roll(flat,-1,axis=d).ravel() for d in range(2)]);minus=np.array([np.roll(flat,1,axis=d).ravel() for d in range(2)])
 base=np.array([0.,1.2,1.5,.1]) if hard else 1.45*(coords.sum(0)%2).ravel();eps=1.8 if hard else 1.85;states=[]
 for bits in itertools.product([0,1],repeat=4):
  phi=(base+np.pi*np.array(bits))%(2*np.pi)
  if all(np.all(abs(wrapped(phi[plus[d]]-phi))<eps) for d in range(2)):states.append(bits)
 assert len(states)==(4 if hard else 16)
 T=np.zeros((len(states),len(states)));energies=[]
 for si,bits in enumerate(states):
  phi=(base+np.pi*np.array(bits))%(2*np.pi);energies.append(kappa*sum(np.cos(phi[plus[d]]-phi).sum() for d in range(2)));probs=np.array([1. if abs(wrapped(phi[plus[d,i]]-phi[i]-np.pi))>=eps else 1-np.exp(-2*max(0,kappa*np.cos(phi[plus[d,i]]-phi[i]))) for d in range(2) for i in range(4)])
  for mask in itertools.product([0,1],repeat=8):
   mask=np.array(mask);prob=np.prod(np.where(mask,probs,1-probs))
   if not prob:continue
   tape=np.where(mask,probs/2,(1+probs)/2).reshape(2,4)
   for seed in range(4):
    out=phi.copy();m.pi_cluster_serial(out,plus,minus,seed,tape,kappa,eps);newbits=tuple((abs(wrapped(out-base))>np.pi/2).astype(int));T[si,states.index(newbits)]+=prob/4
 weights=np.exp(np.array(energies)-max(energies));weights/=weights.sum();np.testing.assert_allclose(T.sum(1),1,atol=2e-14);flux=weights[:,None]*T;np.testing.assert_allclose(flux,flux.T,atol=3e-15,rtol=2e-13)

def kw():return dict(N=4,epsilon=.52*np.pi,kappas=[-2.,-1.8],warmup=2,seed=171313,initial_winding=[1,0,0],initial_kind='winding',checkpoint_every=3)
def test_extended_full_field_audit_and_exact_crash_resume(tmp_path):
 a,b=tmp_path/'a.h5',tmp_path/'b.h5';k=kw();pt.run(a,steps=16,**k)
 with pytest.raises(RuntimeError,match='Injected crash'):pt.run(b,steps=16,inject_crash_after_flush=7,**k)
 pt.run(b,steps=16,resume=True,**k)
 with h5py.File(a) as x,h5py.File(b) as y:
  for name in x:
   if not name.endswith('_seconds'):np.testing.assert_array_equal(x[name][:],y[name][:])
 assert audit.audit(a,dict(k,steps=16),pt.sources())['passed']
 with h5py.File(b,'r+') as h:h['collective_log_ratio'][3,0]+=1
 with pytest.raises(AssertionError):audit.audit(b,dict(k,steps=16),pt.sources())

def test_mode_guard_refuses_before_writes(tmp_path):
 f=tmp_path/'a.h5';k=kw();pt.run(f,steps=3,sector_updates=False,**k);before=hashlib.sha256(f.read_bytes()).hexdigest()
 with pytest.raises(ValueError,match='mismatch'):pt.run(f,steps=6,sector_updates=True,resume=True,**k)
 assert hashlib.sha256(f.read_bytes()).hexdigest()==before
 assert audit.audit(f,dict(k,steps=3,sector_updates=False),pt.sources())['passed']
