"""Signed-cluster per-start mixing, joint moment bootstrap, local contour checks."""
import runner as r
import numpy as np,h5py,json,time,shutil
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from analyze_calibration import helper
P=r.P

def U(x):return 1-x[...,1]/(2*x[...,0]**2)
def R(x,n):
 with np.errstate(invalid='ignore',divide='ignore'):return np.sqrt(x[...,2]/x[...,3]-1)/(2*n*np.sin(np.pi/n))
def clean(v):
 if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [clean(x) for x in v]
 if isinstance(v,np.ndarray):return clean(v.tolist())
 if isinstance(v,(float,np.floating)) and not np.isfinite(v):return None
 if isinstance(v,np.generic):return v.item()
 return v

def crossing(qs,values,target=.25,method='linear'):
 if not np.isfinite(values).all():return None
 found=[]
 for i in range(len(qs)-1):
  a,b=values[i]-target,values[i+1]-target
  if a==0 and b==0:return None
  if a==0:found.append(qs[i])
  elif b==0:found.append(qs[i+1])
  elif a*b<0:
   if method=='linear':found.append(qs[i]-a*(qs[i+1]-qs[i])/(b-a))
   else:found.append(brentq(lambda q:float(PchipInterpolator(qs,values)(q))-target,qs[i],qs[i+1]))
 unique=[]
 for x in found:
  if not any(abs(x-y)<1e-12 for y in unique):unique.append(float(x))
 return unique[0] if len(unique)==1 else None

def diagnose(job,rng,B=1000,cut=None):
 p=job['point'];n=p['N'];cut=p['warmup']+1 if cut is None else cut
 with h5py.File(P/'chains'/f"{p['id']}.h5") as h:
  x=np.column_stack([h['magnetization_squared'][cut:],h['magnetization_fourth'][cut:],h['structure_zero'][cut:],h['structure_min'][cut:].mean(1)]);m=x.mean(0)
  wi=h['winding'][cut:].astype(float);w2=(wi**2).mean((1,2));wa=wi.mean(2)
  ui=-(x[:,1]-m[1])/(2*m[0]**2)+m[1]*(x[:,0]-m[0])/m[0]**3
  ri=(x[:,2]-m[2])/m[3]-m[2]*(x[:,3]-m[3])/m[3]**2
  taus={}
  for name,v in dict(m2=x[:,0],U_influence=ui,R_influence=ri,winding2=w2,cosine=h['mean_link_cosine'][cut:],vortex=h['vortex_density'][cut:]).items():
   for label,part in [('full',v),('late',v[len(v)//2:])]:
    t=helper.tau_int(part);taus[name+'_'+label]=float(t) if t is not None and np.isfinite(t) else 1.
  width=2*max(20,int(np.ceil(10*max(taus.values()))));nb=len(x)//width;flags=[];draws=None
  if nb<20:flags.append('low_block_count')
  if nb>=2:
   blocks=x[:nb*width].reshape(nb,width,4).mean(1);draws=blocks[rng.integers(nb,size=(B,nb))].mean(1)
  else:flags.append('insufficient_blocks')
  sp=np.mean([abs(np.fft.fftn(np.exp(1j*h['phi'][i,0])))**2/n**6 for i in range(cut,len(h['phi']),100)],axis=0);peak=tuple(map(int,np.unravel_index(np.argmax(sp),sp.shape)))
  if peak!=(0,0,0):flags.append('nonzero_sampled_spectrum_peak')
  hs=[part.mean(0) for part in np.array_split(x,2)];use=float(np.std(U(draws),ddof=1)) if draws is not None else None
  if use and abs(U(hs[0])-U(hs[1]))>10*use:flags.append('large_half_shift')
  sizes=h['cluster_sizes'][cut:];changes=int(np.sum(np.any(np.rint(wa[1:])!=np.rint(wa[:-1]),axis=1)))
  d=dict(point=p,id=p['id'],N=n,cut=cut,m2=float(m[0]),U=float(U(m)),R=float(R(m,n)),U_SD=use,taus=taus,block_cycles=width,blocks=nb,flags=flags,half_U=[float(U(v)) for v in hs],mean_winding_axes=wa.mean(0),rounded_winding_changes=changes,zero_rounded_winding_fraction=float(np.mean(np.all(np.rint(wa)==0,axis=1))),spectrum_peak=peak,S0=float(sp[0,0,0]),Speak=float(sp[peak]),mean_cluster_fraction=float(sizes.mean()/n**3),global_cluster_fraction=float(np.mean(sizes==n**3)),twist_accepted=int(h['winding_accepted'][cut:].sum()),postcut_cycles=len(x),generation_seconds=job['generation_seconds'],audit_seconds=job['audit_seconds'])
 return m,draws,d

def main():
 rng=np.random.default_rng(2026101705);points={}
 for f in sorted((P/'batches').glob('*.json')):
  if f.name.endswith('.review.json'):continue
  for p in r.read(f)['points']:points[p['id']]=p
 groups={};rows=[]
 for p in points.values():
  job=r.read(P/'jobs'/f"{p['id']}.json")
  if job['state']!='complete' or not job['audit']['passed']:raise RuntimeError('Await every full audit before analysis')
  m,draws,d=diagnose(job,rng);rows.append(d);groups.setdefault((p['kappa'],p['N'],p['epsilon_over_pi']),[]).append((m,draws,d))
 coords={};summary=[]
 for key,rr in groups.items():
  m=np.mean([v[0] for v in rr],axis=0);draws=np.mean([v[1] for v in rr],axis=0) if all(v[1] is not None for v in rr) else None;flags=sorted(set(f for v in rr for f in v[2]['flags']))
  for i,a in enumerate(rr):
   for b in rr[i+1:]:
    if a[2]['U_SD'] is not None and b[2]['U_SD'] is not None and abs(a[2]['U']-b[2]['U'])>5*np.hypot(a[2]['U_SD'],b[2]['U_SD']):flags.append('large_start_disagreement')
  c=dict(kappa=key[0],N=key[1],q=key[2],U=float(U(m)),R=float(R(m,key[1])),U_MC95=np.quantile(U(draws),[.025,.975]).tolist() if draws is not None else None,flags=sorted(set(flags)),starts=len(rr));summary.append(c);coords[key]=(c,draws)
 roots=[]
 for k,n in sorted(set((key[0],key[1]) for key in coords)):
  qs=np.array(sorted(key[2] for key in coords if key[:2]==(k,n)));cc=[coords[k,n,q] for q in qs];u=np.array([c['U'] for c,_ in cc])
  for method in ['linear','pchip']:
   central=crossing(qs,u,method=method);samples=[]
   if all(x is not None for _,x in cc):
    for b in range(1000):
     value=crossing(qs,np.array([U(x[b]) for _,x in cc]),method=method)
     if value is not None:samples.append(value)
   roots.append(dict(kappa=k,N=n,method=method,epsilon_over_pi=central,MC95=np.quantile(samples,[.025,.975]).tolist() if samples else None,valid_fraction=len(samples)/1000,flags=sorted(set(f for c,_ in cc for f in c['flags'])),q_knots=qs,U_knots=u))
 result=clean(dict(created=r.legacy.now(),rows=rows,coordinates=summary,roots=roots,source_analysis_sha256=r.legacy.sha256(__file__),scope='Signed-cluster finite-volume U=1/4 contours. Shared-moment blocked bootstrap; no certified total coverage or infinite-volume claim. Mixing diagnostics are empirical scope-limited.'))
 if (P/'analysis.json').exists():
  folder=P/'analysis-history'/str(time.time_ns());folder.mkdir(parents=True);shutil.copyfile(P/'analysis.json',folder/'analysis.json')
 r.atomic(P/'analysis.json',result)
 print(json.dumps(dict(coordinates=result['coordinates'],mixing=[{k:v for k,v in x.items() if k in ['id','U','R','spectrum_peak','rounded_winding_changes','blocks','flags','mean_cluster_fraction','global_cluster_fraction']} for x in result['rows']],roots=result['roots'])),flush=True)
if __name__=='__main__':main()
