"""Exploratory measured first-vector-loss locators, with all variant failures retained."""
from pathlib import Path
import json,sys,hashlib,time,argparse
import numpy as np,h5py
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
P=Path(__file__).resolve().parent;S=P.parent/'sign-sector'

def roots(xs,ys,target=0.,kind='linear'):
 if kind=='pchip' and not np.isfinite(ys).all():return []
 pairs=[dict(root=float(x),bracket=[float(x),float(x)]) for x,y in zip(xs,ys) if np.isfinite(y) and y==target]
 for i in range(len(xs)-1):
  if not (np.isfinite(ys[i]) and np.isfinite(ys[i+1])):continue
  if (ys[i]-target)*(ys[i+1]-target)<0:
   x=xs[i]+(target-ys[i])*(xs[i+1]-xs[i])/(ys[i+1]-ys[i])
   if kind=='pchip' and np.isfinite(ys).all():x=brentq(lambda v:float(PchipInterpolator(xs,ys)(v)-target),xs[i],xs[i+1])
   pairs.append(dict(root=float(x),bracket=[float(xs[i]),float(xs[i+1])]))
 return pairs

def moments(x,n):
 m=x.mean(0);u=1-m[1]/(2*m[0]**2);s0=n**3*m[0];r2=(s0/m[2]-1)/(4*n*n*np.sin(np.pi/n)**2) if m[2]>0 else np.nan
 return dict(m2=float(m[0]),m4=float(m[1]),Binder_U=float(u),R=float(np.sqrt(r2)) if r2>=0 else None,C=float(m[3]),vortex=float(m[4]),AF_nematic=float(m[5]),staggered_m2=float(m[6]))

def phase_visits(x):
 states=np.where(x<.05,-1,np.where(x>.15,1,0));last=0;events=[]
 for i,current in enumerate(states):
  if current and current!=last:
   if last:events.append(dict(postcut_index=i,from_state=int(last),to_state=int(current)))
   last=current
 return dict(m2_low_threshold=.05,m2_high_threshold=.15,low_frames=int((states==-1).sum()),high_frames=int((states==1).sum()),intermediate_frames=int((states==0).sum()),frames=len(x),complete_switches=len(events),switch_events=events,scope='Descriptive hysteretic response counts, not independent events or a thermodynamic phase classifier')

def read_chain(folder,point):
 f=folder/'jobs'/f"{point['id']}.json"
 if not f.exists():return None
 job=json.loads(f.read_text())
 if job['state']!='complete' or not job['audit']['passed']:return None
 path=folder/'chains'/f"{point['id']}.h5";cut=point['warmup']+1
 with h5py.File(path) as h:
  x=np.column_stack([h['magnetization_squared'][cut:,0],h['magnetization_fourth'][cut:,0],h['structure_min'][cut:,0].mean(1),h['mean_link_cosine'][cut:,0],h['vortex_density'][cut:,0],h['staggered_nematic_squared'][cut:,0],h['staggered_squared'][cut:,0]])
 return dict(point=point,source=str(path),audit_sha256=job['audit']['sha256'],x=x,means=moments(x,point['N']),halves=[moments(v,point['N']) for v in np.array_split(x,2)],response_visits=phase_visits(x[:,0]))

def clean(v):
 if isinstance(v,dict):return {k:clean(x) for k,x in v.items() if k!='x'}
 if isinstance(v,list):return [clean(x) for x in v]
 if isinstance(v,np.generic):return clean(v.item())
 if isinstance(v,float) and not np.isfinite(v):return None
 return v

def binder_families(coordinates,n,window=None):
 cc=sorted([d for d in coordinates if d['N']==n and (window is None or window[0]<=d['kappa']<=window[1])],key=lambda d:d['kappa']);xs=np.array([d['kappa'] for d in cc]);labels=sorted(set.intersection(*[set(d['variants']) for d in cc])) if cc else []
 result=[]
 for label in labels:
  ys=np.array([d['variants'][label]['Binder_U'] for d in cc])
  for target in [.2,.25,.3]:
   for kind in ['linear','pchip']:
    found=roots(xs,ys,target,kind);result.append(dict(N=n,variant=label,target=target,interpolation=kind,roots=found,status='unique' if len(found)==1 else 'missing' if not found else 'multiple'))
 return result

def main(partial=False,epsilon_over_pi=.52,manifest=None):
 manifest_path=Path(manifest).resolve() if manifest else P/'batches/q052-three-size-local-resolution-v1.json'
 batches=[(manifest_path.parent.parent,manifest_path)]
 if abs(epsilon_over_pi-.52)<1e-12:batches.append((S,S/'batches/q052-sign-sector-mixing-validation-v1.json'))
 chains=[];missing=[];diagnostics={}
 for folder in dict.fromkeys([folder for folder,_ in batches]):
  native=folder/'analysis.json'
  if native.exists():
   content=native.read_bytes();digest=hashlib.sha256(content).hexdigest();archive=folder/'diagnostic-snapshots';archive.mkdir(exist_ok=True);destination=archive/(digest+'.json')
   if not destination.exists():destination.write_bytes(content)
  for f in sorted((folder/'diagnostic-snapshots').glob('*.json')):
   for row in json.loads(f.read_text())['rows']:diagnostics.setdefault(row['id'],(f,row))
 for folder,manifest in batches:
  batch=json.loads(manifest.read_text())
  for point in batch['points']:
   if abs(point['epsilon_over_pi']-epsilon_over_pi)>1e-12:continue
   if folder==S and not(point['sector_updates'] and point['target_kappa'] in [-2.,-2.5]):continue
   d=read_chain(folder,point)
   if d is None:missing.append(point['id'])
   else:
    if point['id'] in diagnostics:
     diagnostic_source,row=diagnostics[point['id']];keys=['flags','taus','blocks','block_cycles','bootstrap_retained_fraction','conditional_block_MC95','mean_winding_axes','twist_accepted','average_peak','sampled_peak_counts'];d['sampling_diagnostics']={k:row[k] for k in keys};d['sampling_diagnostics']['source_analysis']=str(diagnostic_source);d['sampling_diagnostics']['source_analysis_sha256']=hashlib.sha256(diagnostic_source.read_bytes()).hexdigest();d['sampling_diagnostics']['mean_sampled_max_Fourier_power']=float(np.mean([v['power'] for v in row['spectrum_samples']]))
    else:d['sampling_diagnostics']=dict(pending=True)
    chains.append(d)
 if missing and not partial:raise RuntimeError(f'{len(missing)} required trajectories incomplete')
 sizes=sorted({d['point']['N'] for d in chains});coordinates=[]
 for n in sizes:
  for k in sorted({d['point']['target_kappa'] for d in chains if d['point']['N']==n}):
   rr=[d for d in chains if d['point']['N']==n and d['point']['target_kappa']==k]
   if len(rr)!=4:continue
   variants={'pooled':moments(np.concatenate([d['x'] for d in rr]),n)}
   for d in rr:
    label=d['point']['initial_label'];variants[label]=d['means']
    for j,h in enumerate(d['halves']):variants[f'{label}:half{j+1}']=h
   for j in range(2):variants[f'pooled:half{j+1}']=moments(np.concatenate([np.array_split(d['x'],2)[j] for d in rr]),n)
   coordinates.append(dict(N=n,kappa=k,starts=4,variants=variants))
 families=[v for n in sizes for v in binder_families(coordinates,n)]
 # Strict common measured side anchors, ALL individual starts AND halves, all tested thresholds.
 common=[]
 for n in sizes:
  cc=[d for d in coordinates if d['N']==n];lo=[d['kappa'] for d in cc if all(v['Binder_U']<.2 for label,v in d['variants'].items() if label!='pooled')];hi=[d['kappa'] for d in cc if all(v['Binder_U']>.3 for label,v in d['variants'].items() if label!='pooled')]
  first_low=max(lo) if lo else None;higher=[h for h in hi if first_low is not None and h>first_low];best=(first_low,min(higher)) if higher else None
  common.append(dict(N=n,lower_kappa=best[0] if best else None,upper_kappa=best[1] if best else None,status='all_start_and_half_anchor_bracket' if best else 'unresolved',low_phase_candidates=lo,high_phase_candidates=hi))
 pairs=[]
 for i,n in enumerate(sizes):
  for m in sizes[i+1:]:
   a={d['kappa']:d['variants']['pooled']['R'] for d in coordinates if d['N']==n};b={d['kappa']:d['variants']['pooled']['R'] for d in coordinates if d['N']==m};xs=sorted(set(a)&set(b));ys=[a[k]-b[k] if a[k] is not None and b[k] is not None else np.nan for k in xs];rr=roots(np.array(xs),np.array(ys));pairs.append(dict(sizes=[n,m],kappas=xs,difference=clean(ys),roots=rr,status='unique' if len(rr)==1 else 'unresolved_or_multiple'))
 local_families=[]
 for bracket in common:
  if bracket['status']!='unresolved':local_families.extend(binder_families(coordinates,bracket['N'],[bracket['lower_kappa'],bracket['upper_kappa']]))
 primary=next((d for d in local_families if d['N']==max(sizes) and d['variant']=='pooled' and d['target']==.25 and d['interpolation']=='linear'),None) if sizes else None
 result=clean(dict(complete=not missing,missing=missing,manifest_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for _,f in batches},epsilon_over_pi=epsilon_over_pi,trajectories=len(chains),source_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),chains=chains,coordinates=coordinates,families=families,local_families=local_families,local_window_policy='Within the sampled continuation window, take the unanimously suppressed anchor closest to higher kappa, and the nearest higher unanimously high-spin anchor. Keep every measured knot inside. This does not certify absence of transitions in unsampled support toward kappa0. Global scan families and outliers remain diagnostic.',common_anchors=common,xi_pair_crossings=pairs,primary=primary,scope='No root-MC uncertainty is calculated here. Observable conditional-MC intervals and sampling flags are joined from the source analyses; the root bracket/sensitivity ranges are descriptive only. Pooledxi/L centralcrossings are descriptive checks. Exploratory fixed finite-box Binder locators within strong staggered-nematic response. All start/half/threshold/interpolation families retained, including failures. Anchor brackets and variant unions are descriptive, not calibrated confidence intervals or thermodynamic estimates. Pooling never removes incompatible starts. No assumed universality or fitted asymptote.'))
 out=P/'bracket-snapshots'/str(time.time_ns());out.mkdir(parents=True);(out/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
 suffix='' if abs(epsilon_over_pi-.52)<1e-12 else f'-q{round(epsilon_over_pi*1e6):07d}';pointer=P/(('partial-brackets' if partial else 'latest-exploratory-brackets')+suffix+'.json');temp=pointer.with_suffix('.tmp');temp.write_text(json.dumps(dict(analysis=str((out/'analysis.json').resolve()),complete=not missing,trajectories=len(chains)),indent=2)+'\n');temp.replace(pointer)
 print(json.dumps(dict(complete=not missing,trajectories=len(chains),common_anchors=common,primary=primary,xi_pair_crossings=pairs),indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--partial',action='store_true');p.add_argument('--epsilon-over-pi',type=float,default=.52);p.add_argument('--manifest');a=p.parse_args();main(a.partial,a.epsilon_over_pi,a.manifest)
