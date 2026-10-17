"""Exact experiment-local sign-sector and symmetric winding/canting moves."""
import math,hashlib
from pathlib import Path
import numpy as np,numba
from supervillain.generator._xy_kernels import _canonical,_wrap

@numba.njit(cache=False)
def pi_cluster_serial(phi,plus,minus,seed,uniforms,kappa,epsilon):
 v=phi.size;reflected=np.empty(v)
 for i in range(v):reflected[i]=_canonical(phi[i]+math.pi)
 included=np.zeros(v,np.bool_);queue=np.empty(v,np.int64);included[seed]=True;queue[0]=seed;size=1;head=0;tests=0
 while head<size:
  i=queue[head];head+=1
  for mu in range(plus.shape[0]):
   for orientation in range(2):
    j=plus[mu,i] if orientation==0 else minus[mu,i];origin=i if orientation==0 else j
    if included[j]:continue
    tests+=1
    forced=abs(_wrap(phi[j]-reflected[i]))>=epsilon or abs(_wrap(reflected[j]-phi[i]))>=epsilon
    probability=-math.expm1(-2*max(0.,kappa*math.cos(phi[j]-phi[i])))
    if forced or uniforms[mu,origin]<probability:included[j]=True;queue[size]=j;size+=1
 for i in range(v):
  value=reflected[i] if included[i] else phi[i]
  for mu in range(plus.shape[0]):
   j=plus[mu,i];other=reflected[j] if included[j] else phi[j]
   if abs(_wrap(other-value))>=epsilon:raise FloatingPointError('Pi cluster strict boundary; field unchanged')
 for index in range(size):phi[queue[index]]=reflected[queue[index]]
 return size,tests

def shift_choices(D):
 # Full Cartesian set includes every inverse with equal, state-independent probability.
 ws=[np.zeros(D,dtype=int)]
 for d in range(D):
  for sign in [-1,1]:
   w=np.zeros(D,dtype=int);w[d]=sign;ws.append(w)
 return [(sign*np.pi/2,w.copy()) for sign in [-1,1] for w in ws]

def cosine_sum(phi):
 a=np.asarray(phi)[0];return float(sum(np.cos(np.roll(a,-1,axis=d)-a).sum() for d in range(a.ndim)))

def shift_candidate(phi,index):
 a=np.asarray(phi);N=a.shape[-1];D=a.ndim-1;coords=np.indices((N,)*D);angle,w=shift_choices(D)[index]
 g=angle*(coords.sum(0)%2)-sum(2*np.pi*w[d]*coords[d]/N for d in range(D))
 return np.remainder(a+g[None],2*np.pi)

def collective(action,cfg,rng):
 index=int(rng.integers(len(shift_choices(action.Lattice.D))))
 old=cosine_sum(cfg['phi']);candidate=action.canonical(shift_candidate(cfg['phi'],index));new=cosine_sum(candidate);u=float(rng.uniform());valid=bool(action.valid({'phi':candidate}));lr=float(action.kappa*(new-old));accepted=bool(valid and (u==0 or np.log(u)<min(0.,lr)))
 return ({'phi':candidate} if accepted else cfg),dict(index=index,uniform=u,valid=int(valid),accepted=int(accepted),log_ratio=lr,before_C=old,candidate_C=new)

def pi_step(action,cfg,plus,minus,rng):
 phi=np.array(cfg['phi'],copy=True);seed=int(rng.integers(phi.size));uniforms=rng.uniform(size=plus.shape)
 size,tests=pi_cluster_serial(phi.reshape(-1),plus,minus,seed,uniforms,float(action.kappa),float(action.epsilon))
 return {'phi':action.canonical(phi)},size,tests

def definition():
 return dict(sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),pi_activation='forced hard bonds; otherwise 1-exp(-2max(0,kappa*cos(delta)))',collective='uniform14 fixed shifts a*parity-2pi*w.x/N; a=+/-pi/2,w=0 or+/-unit axes; exact inverse pairing; MH exp(kappa*DeltaC)',scope='Isolated exact sector-update experiment; frozen native kernels unchanged')
