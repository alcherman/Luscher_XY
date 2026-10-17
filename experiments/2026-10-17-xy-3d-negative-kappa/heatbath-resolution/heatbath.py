"""Exact restricted-angle conditional heatbath; experiment-local only."""
import math,hashlib
from pathlib import Path
import numpy as np,numba
from supervillain.generator.xy import SiteUpdate as NativeSite
from supervillain.lattice import Form

@numba.njit(cache=False)
def intervals(neighbors,epsilon):
 lo=np.empty(32);hi=np.empty(32);lo[0]=0.;hi[0]=2*math.pi;count=1
 for angle in neighbors:
  begin=(angle-epsilon)%(2*math.pi);end=begin+2*epsilon;left=np.empty(2);right=np.empty(2)
  left[0]=begin;right[0]=min(end,2*math.pi);pieces=1
  if end>2*math.pi:left[1]=0.;right[1]=end-2*math.pi;pieces=2
  nl=np.empty(32);nh=np.empty(32);number=0
  for j in range(count):
   for k in range(pieces):
    a=max(lo[j],left[k]);b=min(hi[j],right[k])
    if b>a:
     if number>=32:raise RuntimeError('Too many allowed components')
     nl[number]=a;nh[number]=b;number+=1
  lo=nl;hi=nh;count=number
 return lo[:count],hi[:count]

@numba.njit(cache=False)
def coefficients(neighbors,kappa,lo,hi):
 A=kappa*np.cos(neighbors).sum();B=kappa*np.sin(neighbors).sum();radius=math.hypot(A,B);maximum=-np.inf;peak=math.atan2(B,A)%(2*math.pi)
 for j in range(len(lo)):
  maximum=max(maximum,A*math.cos(lo[j])+B*math.sin(lo[j]),A*math.cos(hi[j])+B*math.sin(hi[j]))
  if lo[j]<=peak<=hi[j]:maximum=max(maximum,radius)
 # Conservative finite-arithmetic envelope; a constant rejection factor cancels exactly.
 return A,B,maximum+64*np.finfo(np.float64).eps*(1+radius)

@numba.njit(cache=False)
def conditional(neighbors,kappa,epsilon,rng):
 lo,hi=intervals(neighbors,epsilon)
 if not len(lo):raise RuntimeError('Empty conditional support')
 A,B,maximum=coefficients(neighbors,kappa,lo,hi);length=(hi-lo).sum()
 for tries in range(1000000):
  distance=rng.random()*length;value=0.
  for j in range(len(lo)):
   width=hi[j]-lo[j]
   if distance<width or j==len(lo)-1:value=lo[j]+distance;break
   distance-=width
  valid=True
  for angle in neighbors:
   if abs((angle-value+math.pi)%(2*math.pi)-math.pi)>=epsilon:valid=False;break
  u=rng.random()
  if valid and (u==0 or math.log(u)<A*math.cos(value)+B*math.sin(value)-maximum):return value,tries+1
 raise RuntimeError('Conditional rejection budget exceeded')

@numba.njit(cache=False)
def sweep(phi,plus,minus,order,kappa,epsilon,rng):
 draws=0;neighbors=np.empty(2*plus.shape[0])
 for site in order:
  for d in range(plus.shape[0]):neighbors[2*d]=phi[plus[d,site]];neighbors[2*d+1]=phi[minus[d,site]]
  value,tries=conditional(neighbors,kappa,epsilon,rng);phi[site]=value;draws+=tries
 return draws

class HeatbathSite(NativeSite):
 def step(self,cfg):
  phi=cfg['phi'];self._validate_phi(phi);out=Form(np.array(phi,dtype=float,copy=True,order='C'),degree=0,lattice=phi.lattice);order=self.rng.permutation(self.Lattice.sites)
  self.last_draws=int(sweep(np.asarray(out).reshape(-1),self.plus,self.minus,order,self.kappa,self.epsilon,self.rng));self.accepted+=self.Lattice.sites;self.proposed+=self.Lattice.sites;self.sweeps+=1;self.acceptance+=1.
  return dict(cfg,phi=out)

def definition():
 return dict(sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),conditional='exp(kappa*sum_neighbor cos(theta-phi_neighbor)) restricted to intersection of strict hard arcs',sampling='Uniform on exact allowed interval union; rejection using analytic maximum over endpoints/interior stationary peak; random site permutation',scope='Exact conditional heatbath, frozen native kernels unchanged')
