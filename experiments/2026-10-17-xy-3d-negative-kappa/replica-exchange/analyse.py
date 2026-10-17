"""Physical slot estimates and shared-time-block replica-exchange diagnostics."""
from pathlib import Path
import importlib.util,json,time,shutil
_spec=importlib.util.spec_from_file_location('joint_pt_runner_analysis',Path(__file__).with_name('runner.py'));r=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(r)
import h5py,numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from analyze_calibration import helper
P=r.P

def U(x):return 1-x[...,1]/(2*x[...,0]**2)
def ratio(x,n):
 with np.errstate(invalid='ignore',divide='ignore'):return np.sqrt(x[...,2]/x[...,3]-1)/(2*n*np.sin(np.pi/n))
def clean(x):
 if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [clean(v) for v in x]
 if isinstance(x,np.ndarray):return clean(x.tolist())
 if isinstance(x,(float,np.floating)) and not np.isfinite(x):return None
 return x.item() if isinstance(x,np.generic) else x

def root(ks,values,target=.25,method='linear',all_roots=False):
 if not np.isfinite(values).all():return [] if all_roots else None
 found=[]
 for i in range(len(ks)-1):
  a,b=values[i]-target,values[i+1]-target
  if a==0 and b==0:return [] if all_roots else None
  if a==0:found.append(ks[i])
  elif b==0:found.append(ks[i+1])
  elif a*b<0:
   if method=='linear':found.append(ks[i]-a*(ks[i+1]-ks[i])/(b-a))
   else:found.append(brentq(lambda k:float(PchipInterpolator(ks,values)(k))-target,ks[i],ks[i+1]))
 unique=[]
 for x in found:
  if not any(abs(x-y)<1e-12 for y in unique):unique.append(float(x))
 # The continuation is the first crossing reached from kappa=0. Preserve all
 # additional roots as a diagnostic rather than silently calling them absent.
 return unique if all_roots else (max(unique) if unique else None)

def trips(walkers):
 T,R=walkers.shape
 if R<2:return dict(round_trips=0,durations=[],target_unique_walkers=1,bridge_unique_walkers=1)
 phase=np.zeros(R,int);began=np.zeros(R,int);count=0;dur=[]
 for t in range(T):
  left,right=int(walkers[t,0]),int(walkers[t,-1])
  if phase[left]==2:count+=1;dur.append(t-int(began[left]));phase[left]=1;began[left]=t
  elif phase[left]==0:phase[left]=1;began[left]=t
  if phase[right]==1:phase[right]=2
 return dict(round_trips=count,durations=dur,target_unique_walkers=len(np.unique(walkers[:,0])),bridge_unique_walkers=len(np.unique(walkers[:,-1])))

def response_changes(series):
 codes=np.where(series>.10,1,np.where(series<.03,-1,0));known=codes[codes!=0]
 return int(np.count_nonzero(np.diff(known))) if len(known)>1 else 0

def diagnose(job,rng,B=1000,cut=None):
 p=job['point'];ks=np.asarray(p['kappas']);S=len(ks);n=p['N'];cut=p['warmup']+1 if cut is None else cut
 with h5py.File(P/'chains'/f"{p['id']}.h5") as h:
  x=np.stack([h['magnetization_squared'][cut:],h['magnetization_fourth'][cut:],h['structure_zero'][cut:],h['structure_min'][cut:].mean(2)],axis=2);means=x.mean(0);T=len(x)
  w2=np.empty((T,S));wa=np.empty((T,S,3))
  for start in range(cut,len(h['sweep']),100):
   stop=min(start+100,len(h['sweep']));w=h['winding'][start:stop].astype(float);w2[start-cut:stop-cut]=(w*w).mean((2,3));wa[start-cut:stop-cut]=w.mean(3)
  taus=[];slot_taus=[]
  for i in range(S):
   m=means[i];ui=-(x[:,i,1]-m[1])/(2*m[0]**2)+m[1]*(x[:,i,0]-m[0])/m[0]**3;ri=(x[:,i,2]-m[2])/m[3]-m[2]*(x[:,i,3]-m[3])/m[3]**2
   values=[]
   for series in [x[:,i,0],ui,ri,w2[:,i],h['mean_link_cosine'][cut:,i]]:
    for part in [series,series[T//2:]]:
     tau=helper.tau_int(part);values.append(float(tau) if tau is not None and np.isfinite(tau) else 1.)
   slot_taus.append(max(values));taus+=values
  width=2*max(20,int(np.ceil(10*max(taus))));nb=T//width;flags=[];draws=None
  if nb<20:flags.append('low_block_count')
  if nb>=2:
   blocks=x[:nb*width].reshape(nb,width,S,4).mean(1);indices=rng.integers(nb,size=(B,nb));draws=blocks[indices].mean(1)
  else:flags.append('insufficient_blocks')
  rates=h['exchange_accepted'][cut:].sum(0)/np.maximum(1,h['exchange_attempted'][cut:].sum(0));walkers=h['walker_id'][cut:];trip=trips(walkers);bywalker=np.take_along_axis(x[:,:,0],np.argsort(walkers,axis=1),axis=1);walker_changes=[response_changes(bywalker[:,i]) for i in range(S)]
  if S>1 and (np.min(rates)<.05 or trip['round_trips']<2):flags.append('inadequate_ladder_traversal')
  initial_target=np.asarray(h['phi'][0,0,0]);spectra=[]
  for start in range(cut,len(h['sweep']),100):
   fields=h['phi'][start,:,0];spectra.append(abs(np.fft.fftn(np.exp(1j*fields),axes=(1,2,3)))**2/n**6)
  spectral_samples=np.asarray(spectra);sample_s0=spectral_samples[:,:,0,0,0].copy();nonzero=spectral_samples.copy();nonzero[:,:,0,0,0]=0;sample_nonzero=nonzero.max(axis=(2,3,4));coherent=(sample_nonzero>.1)&(sample_nonzero>4*sample_s0);coherent_fraction=coherent.mean(0)
  spectrum=np.mean(spectra,axis=0);peaks=[tuple(map(int,np.unravel_index(np.argmax(s),s.shape))) for s in spectrum]
  halves=[part.mean(0) for part in np.array_split(x,2)];u_sd=np.std(U(draws),axis=0,ddof=1) if draws is not None else np.full(S,np.nan);m2_sd=np.std(draws[:,:,0],axis=0,ddof=1) if draws is not None else np.full(S,np.nan)
  if draws is not None and np.any(abs(U(halves[0])-U(halves[1]))>10*u_sd):flags.append('large_half_shift')
  d=dict(id=p['id'],point=p,cut=cut,postcut_cycles=T,bootstrapped_cycles=nb*width,bootstrap_retained_fraction=nb*width/T,bootstrap_remainder_cycles=T-nb*width,slots=S,U=U(means),R=ratio(means,n),m2=means[:,0],U_SD=u_sd,m2_SD=m2_sd,half_U=[U(v) for v in halves],half_m2=[v[:,0] for v in halves],slot_max_tau=slot_taus,max_tau=max(taus),block_cycles=width,blocks=nb,flags=flags,exchange_rates=rates,walker_round_trips=trip,walker_uniform_response_changes=walker_changes,response_change_thresholds='m2>0.10 vs m2<0.03, intermediate rows ignored; descriptive events, not independent phase conversions',spectrum_peaks=peaks,coherent_nonzero_warning_fraction=coherent_fraction,spectrum_sample_S0=sample_s0,spectrum_sample_max_nonzero=sample_nonzero,coherent_threshold='Descriptive only: per-configuration max nonzero power >0.1 and >4*S0; not a phase classifier',S0=spectrum[:,0,0,0],Speak=np.max(spectrum,axis=(1,2,3)),target_zero_rounded_winding_fraction=float(np.mean(np.all(np.rint(wa[:,0])==0,axis=1))),target_winding_changes=int(np.sum(np.any(np.rint(wa[1:,0])!=np.rint(wa[:-1,0]),axis=1))),generation_seconds=job['generation_seconds'],audit_seconds=job['audit_seconds'])
 return means,draws,d

def main():
 rng=np.random.default_rng(2026101711);ps={}
 for manifest in sorted((P/'batches').glob('*.json')):
  if manifest.name.endswith('.review.json'):continue
  for p in r.read(manifest)['points']:ps[p['id']]=p
 rows=[];groups={}
 for p in ps.values():
  jf=P/'jobs'/f"{p['id']}.json"
  if not jf.exists() or r.read(jf)['state']!='complete':continue
  j=r.read(jf)
  if not j['audit']['passed']:raise RuntimeError('Full joint audit required')
  mean,draws,d=diagnose(j,rng);rows.append(d)
  if p.get('role') in ['production','fine_transport_validation']:groups.setdefault((p['N'],p['epsilon_over_pi'],tuple(p['kappas'])),[]).append((mean,draws,d))
 roots=[];coords=[]
 for key,rr in groups.items():
  n,q,kt=key;ks=np.asarray(kt);mean=np.mean([v[0] for v in rr],axis=0);draws=np.mean([v[1] for v in rr],axis=0) if all(v[1] is not None for v in rr) else None;flags=sorted(set(f for v in rr for f in v[2]['flags']))
  for i,a in enumerate(rr):
   for b in rr[i+1:]:
    if draws is not None:
     se=np.hypot(a[2]['m2_SD'],b[2]['m2_SD'])
     if np.any(abs(a[2]['m2']-b[2]['m2'])>5*se):flags.append('large_start_disagreement')
  for i,k in enumerate(ks):coords.append(dict(N=n,q=q,kappa=float(k),U=float(U(mean)[i]),m2=float(mean[i,0]),starts=len(rr)))
  for method in ['linear','pchip']:
   all_central=root(ks,U(mean),method=method,all_roots=True);central=max(all_central) if all_central else None;sample=[];multiple=0
   if draws is not None:
    for b in range(1000):
     br=root(ks,U(draws[b]),method=method,all_roots=True);multiple+=len(br)>1
     if br:sample.append(max(br))
   flags_here=sorted(set(flags+(['unresolved_bootstrap_root'] if len(sample)<990 else [])+(['multiple_central_crossings'] if len(all_central)>1 else [])+(['multiple_bootstrap_crossings'] if multiple>50 else [])))
   local_indices=list(range(max(0,int(np.searchsorted(ks,central))-2),min(len(ks),int(np.searchsorted(ks,central))+2))) if central is not None else []
   local=dict(indices=local_indices,max_tau=max((v[2]['slot_max_tau'][i] for v in rr for i in local_indices),default=None),max_coherent_nonzero_fraction=max((v[2]['coherent_nonzero_warning_fraction'][i] for v in rr for i in local_indices),default=None))
   roots.append(dict(local_diagnostics=local,source_ids=[v[2]['id'] for v in rr],N=n,epsilon_over_pi=q,kappa=central,MC95=np.quantile(sample,[.025,.975]).tolist() if sample else None,valid_fraction=len(sample)/1000,method=method,all_central_roots=all_central,multiple_root_draws=multiple,selected_bracket=([float(ks[max(0,np.searchsorted(ks,central)-1)]),float(ks[min(len(ks)-1,np.searchsorted(ks,central))])] if central is not None else None),flags=flags_here,starts=len(rr),kappa_knots=ks,U_knots=U(mean)))
 result=clean(dict(created=r.pt.legacy.now(),rows=rows,coordinates=coords,roots=roots,analysis_sha256=r.pt.legacy.sha256(__file__),scope='Finite-box U=1/4 contours, first crossing from kappa0. All slots use shared time-block resampling, independent ensemble starts. No certified total or infinite-volume precision.'))
 if (P/'analysis.json').exists():
  folder=P/'analysis-history'/str(time.time_ns());folder.mkdir(parents=True);shutil.copyfile(P/'analysis.json',folder/'analysis.json')
 source_dir=P/'analysis-source-snapshots';source_dir.mkdir(exist_ok=True);shutil.copyfile(__file__,source_dir/(result['analysis_sha256']+'.py'))
 r.atomic(P/'analysis.json',result);print(json.dumps(dict(rows=[{k:v for k,v in d.items() if k in ['id','U','m2','blocks','flags','exchange_rates','walker_round_trips','target_zero_rounded_winding_fraction']} for d in result['rows']],roots=result['roots'])),flush=True)
if __name__=='__main__':main()
