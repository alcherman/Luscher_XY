"""Per-start strong-negative discovery, preserving competing observables."""
from pathlib import Path
import importlib.util,json,hashlib,time,shutil,sys
P=Path(__file__).resolve().parent
_sp=importlib.util.spec_from_file_location('strong_pilot_runner',P/'runner.py');r=importlib.util.module_from_spec(_sp);_sp.loader.exec_module(r)
import h5py,numpy as np
from analyze_calibration import helper

def clean(v):
 if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [clean(x) for x in v]
 if isinstance(v,np.ndarray):return clean(v.tolist())
 if isinstance(v,(float,np.floating)) and not np.isfinite(v):return None
 return v.item() if isinstance(v,np.generic) else v

def diagnose(job,draws=1000):
 p=job['point'];n=p['N'];cut=p['warmup']+1;fn=P/'chains'/f"{p['id']}.h5"
 names=['magnetization_squared','magnetization_fourth','staggered_squared','nematic_squared','mean_link_cosine','vortex_density']
 with h5py.File(fn) as h:
  x=np.column_stack([h[k][cut:,0] for k in names]);m=x.mean(0);T=len(x)
  w=h['winding'][cut:,0].astype(float);w2=(w*w).mean((1,2));wa=w.mean(2)
  influence=-(x[:,1]-m[1])/(2*m[0]**2)+m[1]*(x[:,0]-m[0])/m[0]**3
  taus={}
  for name,series in list(zip(names,x.T))+[('winding_squared',w2),('binder_influence',influence)]:
   for label,part in [('all',series),('late',series[T//2:])]:
    tau=helper.tau_int(part);taus[name+'_'+label]=max(1.,float(tau)) if tau is not None and np.isfinite(tau) else 1.
  width=2*max(20,int(np.ceil(10*max(taus.values()))));nb=T//width;flags=[];ci=None
  if nb<20:flags.append('low_block_count')
  if nb>=2:
   blocks=x[:nb*width].reshape(nb,width,len(names)).mean(1);rng=np.random.default_rng(np.random.SeedSequence([2026101708,p['seed']]));bs=blocks[rng.integers(nb,size=(draws,nb))].mean(1);ci=np.quantile(bs,[.025,.975],axis=0).T
  else:flags.append('insufficient_blocks')
  halves=[v.mean(0) for v in np.array_split(x,2)];sps=[];peakrows=[]
  for i in range(cut,len(h['sweep']),50):
   field=h['phi'][i,0,0];sp=abs(np.fft.fftn(np.exp(1j*field)))**2/n**6;sps.append(sp);peak=np.unravel_index(np.argmax(sp),sp.shape);peakrows.append(dict(cycle=i,momentum_index=peak,power=float(sp[peak]),zero=float(sp[0,0,0]),staggered=float(sp[n//2,n//2,n//2])))
  avg=np.mean(sps,axis=0);peak=np.unravel_index(np.argmax(avg),avg.shape);counts={}
  for d in peakrows:
   key=','.join(map(str,d['momentum_index']));counts[key]=counts.get(key,0)+1
  if peak!=(0,0,0):flags.append('nonzero_average_spectrum_peak')
  mean=dict(zip(names,m));mean['Binder_U']=1-m[1]/(2*m[0]**2);mean['cosine_excess_above_infimum']=m[4]-np.cos(p['epsilon'])
  return dict(point=p,id=p['id'],mean=mean,half_means=[dict(zip(names,v)) for v in halves],conditional_block_MC95=dict(zip(names,ci.tolist())) if ci is not None else None,taus=taus,block_cycles=width,blocks=nb,bootstrap_retained_fraction=nb*width/T,postcut_cycles=T,flags=flags,mean_winding_axes=wa.mean(0),zero_rounded_winding_fraction=float(np.mean(np.all(np.rint(wa)==0,axis=1))),twist_accepted=int(h['winding_accepted'][cut:,0].sum()),mean_cluster_fraction=float(h['cluster_sizes'][cut:,0].mean()/n**3),average_spectrum=avg,average_peak=peak,average_peak_power=float(avg[peak]),spectrum_samples=peakrows,sampled_peak_counts=counts,generation_seconds=job['generation_seconds'],audit_seconds=job['audit_seconds'],full_file_sha256=job['audit']['sha256'],actual_bytes=fn.stat().st_size)

def main(manifest):
 b=r.read(manifest);rows=[]
 for p in b['points']:
  j=r.read(P/'jobs'/f"{p['id']}.json")
  if j['state']!='complete' or not j['audit']['passed']:raise RuntimeError('All new pilot closure audits required')
  rows.append(diagnose(j))
 groups=[]
 for key in sorted({(d['point']['N'],d['point']['epsilon_over_pi'],d['point']['target_kappa']) for d in rows}):
  rr=[d for d in rows if (d['point']['N'],d['point']['epsilon_over_pi'],d['point']['target_kappa'])==key]
  groups.append(dict(N=key[0],q=key[1],kappa=key[2],starts=len(rr),initial_kinds=[d['point']['initial_kind'] for d in rr],means_by_start={d['id']:d['mean'] for d in rr},ranges={name:[min(d['mean'][name] for d in rr),max(d['mean'][name] for d in rr)] for name in rr[0]['mean']},scope='One independently seeded trajectory per initial-condition class, not repeated within-class replication. No pooling or accepted critical roots.'))
 result=clean(dict(created=r.pt.legacy.now(),manifest=str(manifest),manifest_sha256=r.pt.legacy.sha256(manifest),analysis_sha256=r.pt.legacy.sha256(__file__),rows=rows,groups=groups,scope='Strong-negative competing-response discovery only. Separate starts, full spectra and source-preserving histories; no infinite-volume, mixing or asymptote inference from agreement alone.'))
 if (P/'analysis.json').exists():
  target=P/'analysis-history'/str(time.time_ns());target.mkdir(parents=True);shutil.copyfile(P/'analysis.json',target/'analysis.json')
 source=P/'analysis-source-snapshots';source.mkdir(exist_ok=True);shutil.copyfile(__file__,source/(result['analysis_sha256']+'.py'));r.atomic(P/'analysis.json',result)
 print(json.dumps(dict(groups=result['groups'],chains=len(rows))),flush=True)
if __name__=='__main__':main(Path(sys.argv[1]).resolve())
