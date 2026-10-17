"""A finite-volume Gibbs probability bound; no mixing or order theorem."""
from pathlib import Path
import json,math
import h5py,numpy as np
P=Path(__file__).resolve().parent

def log_bound_per_site(D,kappa,epsilon,eta,c):
 if not (isinstance(D,int) and D>=1 and kappa<0 and 0<epsilon<math.pi and 0<eta<epsilon/4):raise ValueError('Outside derived domain')
 return min(0.,math.log(math.pi/eta)+D*kappa*(c-math.cos(epsilon-4*eta)))
if __name__=='__main__':
 a=dict(rows=sum([json.loads((P/f).read_text())['rows'] for f in ['scaled-N6-analysis.json','size-validation-analysis.json']],[]));out=[]
 for d in a['rows']:
  p=d['point']
  eta=1/(12*abs(p['target_kappa']));c=math.cos(p['epsilon']-4*eta)+(math.log(math.pi/eta)+1)/(3*abs(p['target_kappa']));per=log_bound_per_site(3,p['target_kappa'],p['epsilon'],eta,c)
  with h5py.File(P/'chains'/f"{p['id']}.h5") as h:x=h['mean_link_cosine'][p['warmup']+1:,0]
  out.append(dict(id=p['id'],N=p['N'],initial_label=p['initial_label'],eta=eta,threshold=c,threshold_rule='cos(epsilon-4eta)+(log(pi/eta)+1)/(3*abs(kappa)); independent of observed trajectory',kappa=p['target_kappa'],epsilon_over_pi=p['epsilon_over_pi'],log_probability_upper_bound=p['N']**3*per,frames_above_threshold=int((x>=c).sum()),frames=len(x),fraction_above_threshold=float((x>=c).mean()),mean_link_cosine=float(x.mean()),min=float(x.min()),max=float(x.max()),full_file_sha256=d['full_file_sha256']))
 result=dict(derivation='For normalized product Haar measure, on an even periodic bipartite box take even angles in (-eta,eta) and odd angles in (epsilon-3eta,epsilon-eta). This admissible product set has measure (eta/pi)^V and C/V <= D*cos(epsilon-4eta). Since kappa<0, Z >= (eta/pi)^V*exp(kappa*DV*cos(epsilon-4eta)). On C/(DV)>=c, exp(kappa*C)<=exp(kappa*DV*c); its Haar measure is at most1. Therefore P(event)<=min(1,exp(V*(log(pi/eta)+D*kappa*(c-cos(epsilon-4eta))))).',scope='Exact instantaneous equilibrium configuration-event bound, not a bound on a trajectory mean, not a relaxation-time estimate and not evidence that the complementary sampled sector has equilibrated or ordered. Observed frame fractions are descriptive correlated trajectory frequencies. Persistent high-event occupancy identifies nonequilibrium sampling at these parameters.',rows=out)
 (P/'energy-sector-check.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(out,indent=2))
