"""Finite-volume Binder contours and common-R diagnostics; not an infinite-volume result."""
import importlib.util
from pathlib import Path
_spec=importlib.util.spec_from_file_location("negative_scout_runtime",Path(__file__).with_name("scout.py"));s=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(s)
import numpy as np,h5py,json,csv
from pathlib import Path
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys=s.sys;sys.path.insert(0,str(s.OLD));from analyze_calibration import helper
P=s.P

def binder(x):return 1-x[...,1]/(2*x[...,0]**2)
def ratio(x,n):
 with np.errstate(invalid='ignore',divide='ignore'):return np.sqrt(x[...,2]/x[...,3]-1)/(2*n*np.sin(np.pi/n))
def root(q,u,target=.25):
 q=np.asarray(q);u=np.asarray(u)
 if not np.isfinite(u).all():return None
 found=[]
 for i in range(len(q)-1):
  a,b=u[i]-target,u[i+1]-target
  if a==0 and b==0:return None
  if a==0:found.append(float(q[i]))
  elif b==0:found.append(float(q[i+1]))
  elif a*b<0:found.append(float(q[i]+(-a)*(q[i+1]-q[i])/(b-a)))
 unique=[]
 for x in found:
  if not any(abs(x-y)<1e-12 for y in unique):unique.append(x)
 return unique[0] if len(unique)==1 else None

def main():
 import time,shutil
 if (P/'analysis.json').exists():
  history=P/'analysis-history'/str(time.time_ns());history.mkdir(parents=True)
  for name in ['analysis.json','finite-volume-contours.csv','negative-kappa-preliminary.png','negative-kappa-preliminary.pdf','negative-kappa-preliminary.svg']:
   if (P/name).exists():shutil.copyfile(P/name,history/name)
 rng=np.random.default_rng(2026101709);groups={};rows=[];B=1000
 points={}
 for manifest in sorted((P/'batches').glob('*.json')):
  for point in s.read(manifest)['points']:points[point['id']]=point
 for p in points.values():
  j=s.read(P/'jobs'/f"{p['id']}.json")
  if j['state']!='complete' or not j['audit']['passed']:raise RuntimeError('Missing full audit')
  with h5py.File(P/'chains'/f"{p['id']}.h5") as h:
   cut=1001;x=np.column_stack([h['magnetization_squared'][cut:],h['magnetization_fourth'][cut:],h['structure_zero'][cut:],h['structure_min'][cut:].mean(1)])
   mean=x.mean(0);influence=-(x[:,1]-mean[1])/(2*mean[0]**2)+mean[1]*(x[:,0]-mean[0])/mean[0]**3
   taus=[]
   ri=(x[:,2]-mean[2])/mean[3]-mean[2]*(x[:,3]-mean[3])/mean[3]**2
   w2=(h['winding'][cut:].astype(float)**2).mean((1,2))
   for v in [x[:,0],influence,ri,w2,h['vortex_density'][cut:],h['mean_link_cosine'][cut:]]:
    for part in [v,v[len(v)//2:]]:
     t=helper.tau_int(part);taus.append(float(t) if t is not None and np.isfinite(t) else 1.)
   width=max(20,int(np.ceil(10*max(taus))));nb=len(x)//width;flags=[]
   if nb<20:flags.append('fewer_than_20_blocks_per_start')
   if nb<2:draws=None;flags.append('insufficient_blocks')
   else:
    blocks=x[:nb*width].reshape(nb,width,4).mean(1);draws=blocks[rng.integers(nb,size=(B,nb))].mean(1)
   spectra=[]
   for i in range(cut,len(h['phi']),100):spectra.append(abs(np.fft.fftn(np.exp(1j*h['phi'][i,0])))**2/p['N']**6)
   spectrum=np.mean(spectra,axis=0);peak=np.unravel_index(np.argmax(spectrum),spectrum.shape)
   if peak!=(0,0,0):flags.append('sampled_nonzero_momentum_peak_ferromagnetic_inference_unresolved')
   halves=[v.mean(0) for v in np.array_split(x,2)];ub=[float(binder(v)) for v in halves]
   sd=float(np.std(binder(draws),ddof=1)) if draws is not None else None
   if sd and abs(ub[0]-ub[1])>10*sd:flags.append('large_half_shift')
   rec=dict(point=p,binder=float(binder(mean)),R=float(ratio(mean,p['N'])),binder_SD=sd,max_tau=max(taus),block_cycles=width,blocks=nb,flags=flags,half_Binder=ub,peak_integer_momentum=list(map(int,peak)),S0=float(spectrum[0,0,0]),Sstaggered=float(spectrum[(p['N']//2,)*3]),Speak=float(spectrum[peak]),FFT_fields=len(spectra),global_twist_accepted=int(h['winding_accepted'][cut:].sum()),attempts=len(x))
  rows.append(rec);groups.setdefault((p['kappa'],p['N'],p['epsilon_over_pi']),[]).append((mean,draws,rec))
 coords={}
 for key,rr in groups.items():
  means=np.mean([v[0] for v in rr],axis=0);draws=np.mean([v[1] for v in rr],axis=0) if all(v[1] is not None for v in rr) else None
  flags=sorted(set(f for _,_,r in rr for f in r['flags']))
  if draws is not None:
   diff=abs(rr[0][2]['binder']-rr[1][2]['binder']);se=np.hypot(rr[0][2]['binder_SD'],rr[1][2]['binder_SD'])
   if diff>5*se:flags.append('large_start_disagreement')
  if not np.isfinite(ratio(means,key[1])):flags.append('undefined_second_moment_ratio')
  coords[key]=dict(U=float(binder(means)),R=float(ratio(means,key[1])),draws=draws,flags=flags)
 roots=[]
 for k in [0.,-.1,-.25,-.5]:
  for n in [8,12]:
   qs=np.array(sorted(q for kk,nn,q in coords if kk==k and nn==n));cc=[coords[k,n,q] for q in qs];u=np.array([c['U'] for c in cc]);central=root(qs,u);sample=[]
   if all(c['draws'] is not None for c in cc):
    for b in range(B):
     v=root(qs,np.array([binder(c['draws'][b]) for c in cc]));
     if v is not None:sample.append(v)
   valid=len(sample)/B;flags=sorted(set(f for c in cc for f in c['flags']))
   if valid<.99:flags.append('unresolved_bootstrap_bracket')
   roots.append(dict(kappa=k,N=n,epsilon_over_pi=central,MC95=np.quantile(sample,[.025,.975]).tolist() if sample else None,valid_fraction=valid,flags=flags,q_knots=qs.tolist(),U_knots=u.tolist()))
 # Common-R crossing is checked independently; no extrapolation is attempted from two boxes.
 crossings=[]
 for k in [0.,-.1,-.25,-.5]:
  qs=np.array(sorted(q for kk,n,q in coords if kk==k and n==8));diff=np.array([coords[k,8,q]['R']-coords[k,12,q]['R'] for q in qs]);crossings.append(dict(kappa=k,epsilon_over_pi=root(qs,diff,0),R8_minus_R12=diff.tolist(),scope='Three coarse knots, diagnostic only'))
 output=dict(created=s.legacy.now(),rows=rows,roots=roots,common_R_crossings=crossings,source_manifests=[str(f) for f in sorted((P/'batches').glob('*.json'))],scope='4000-sweep Metropolis scout on N8/12, not a thermodynamic transition curve; data-selected coarse interpolation; conditional MC intervals do not cover interpolation/mixing/finite-size bias.')
 def clean(v):
  if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
  if isinstance(v,list):return [clean(x) for x in v]
  if isinstance(v,(float,np.floating)) and not np.isfinite(v):return None
  return v
 output=clean(output)
 output['postprocessing_source_sha256']=s.legacy.sha256(__file__)
 s.atomic(P/'analysis.json',output)
 with (P/'finite-volume-contours.csv').open('w') as f:
  w=csv.writer(f);w.writerow(['kappa','N','epsilon_over_pi','MC95_low','MC95_high','valid_fraction','flags'])
  for r in roots:w.writerow([r['kappa'],r['N'],r['epsilon_over_pi'],*(r['MC95'] or [None,None]),r['valid_fraction'],';'.join(r['flags'])])
 fig,axes=plt.subplots(1,2,figsize=(11,4.5))
 for n,color,marker in [(8,'C0','o'),(12,'C1','s')]:
  for r in [v for v in roots if v['N']==n and v['epsilon_over_pi'] is not None]:
   x=r['epsilon_over_pi'];ci=r['MC95']
   axes[0].plot(x,r['kappa'],marker=marker,color=color,mfc='none' if r['flags'] else color,ls='none',label=f'N={n}, U=1/4' if r['kappa']==0 else None)
   if ci:
    axes[0].hlines(r['kappa'],ci[0],ci[1],color=color);axes[0].vlines(ci,[r['kappa']-.004]*2,[r['kappa']+.004]*2,color=color)
 ks=np.linspace(-.5,0,120);qtree=[brentq(lambda q:5*s.v1(k,np.pi*q)-1,.01,1) for k in ks];axes[0].plot(qtree,ks,'k--',label='Tree: 5 v₁=1')
 axes[0].set(xlabel='ε/π',ylabel='κ',title='3D: preliminary finite-box contours');axes[0].legend(fontsize=8)
 for k,color in zip([0.,-.1,-.25,-.5],['C0','C1','C2','C3']):
  for n,marker in [(8,'o'),(12,'s')]:
   qq=sorted(q for kk,nn,q in coords if kk==k and nn==n);axes[1].plot(qq,[coords[k,n,q]['U'] for q in qq],marker=marker,ls='-' if n==8 else '--',color=color,label=f'κ={k:g}, N={n}')
 axes[1].axhline(.25,color='k',lw=.7);axes[1].set(xlabel='ε/π',ylabel='U = 1 − ⟨m⁴⟩/(2⟨m²⟩²)',title='Measured Binder data (pooled starts)');axes[1].legend(fontsize=7,ncol=2)
 for ax in axes:ax.grid(alpha=.2)
 fig.text(.5,.01,'Hollow points carry sampling/bracket flags. No infinite-volume extrapolation; bars are conditional MC intervals.',ha='center',fontsize=8);fig.tight_layout(rect=[0,.04,1,1])
 for ext in ['png','pdf','svg']:fig.savefig(P/f'negative-kappa-preliminary.{ext}',dpi=180)
 print(json.dumps(dict(roots=roots,crossings=crossings)),flush=True)
if __name__=='__main__':main()
