"""Full-history higher harmonics derived from immutable audited field archives."""
from pathlib import Path
import json,hashlib,os,time,sys
import importlib.util
P=Path(__file__).resolve().parent
_s=importlib.util.spec_from_file_location('harmonic_pilot_runner',P/'runner.py');r=importlib.util.module_from_spec(_s);_s.loader.exec_module(r)
import numpy as np,h5py
from analyze_calibration import helper

def harmonics(phi):
 z=np.exp(1j*phi);parity=(-1.)**np.indices(phi.shape[1:]).sum(0);z2=z*z
 return np.column_stack([abs((z2*parity).mean((1,2,3)))**2,abs((z2*z2).mean((1,2,3)))**2])
def main():
 b=r.read(P/'batches/strong-negative-four-start-discovery-v1.json');folder=P/'derived-harmonic-histories';folder.mkdir(exist_ok=True);rows=[];scriptsha=r.pt.legacy.sha256(__file__)
 names=['magnetization_squared','staggered_squared','mean_link_cosine','vortex_density','staggered_nematic_squared','fourfold_squared']
 for p in b['points']:
  j=r.read(P/'jobs'/f"{p['id']}.json");assert j['state']=='complete' and j['audit']['passed'];source=P/'chains'/f"{p['id']}.h5";sha=r.pt.legacy.sha256(source);assert sha==j['audit']['sha256'];cache=folder/f"{p['id']}.npz"
  if cache.exists():
   with np.load(cache) as z:
    assert str(z['source_sha256'])==sha and str(z['analysis_sha256'])==scriptsha;extra=z['higher_harmonics'].copy()
  else:
   extra=np.empty((p['steps']+1,2))
   with h5py.File(source) as h:
    for first in range(0,p['steps']+1,64):extra[first:first+64]=harmonics(h['phi'][first:first+64,0,0])
   temp=cache.with_suffix('.npz.prepared')
   with temp.open('wb') as out:
    np.savez_compressed(out,higher_harmonics=extra,sweep=np.arange(p['steps']+1),source_sha256=sha,analysis_sha256=scriptsha);out.flush();os.fsync(out.fileno())
   os.replace(temp,cache)
  with h5py.File(source) as h:
   cut=p['warmup']+1;x=np.column_stack([h[k][cut:,0] for k in names[:4]]+[extra[cut:,0],extra[cut:,1]])
  taus=[]
  for v in x.T:
   for part in [v,v[len(v)//2:]]:
    tau=helper.tau_int(part);taus.append(max(1.,float(tau)) if tau is not None and np.isfinite(tau) else 1.)
  width=2*max(20,int(np.ceil(10*max(taus))));nb=len(x)//width;ci=None
  if nb>=2:
   blocks=x[:nb*width].reshape(nb,width,len(names)).mean(1);rng=np.random.default_rng(np.random.SeedSequence([p['seed'],2026101712]));sample=blocks[rng.integers(nb,size=(1000,nb))].mean(1);ci=dict(zip(names,np.quantile(sample,[.025,.975],axis=0).T.tolist()))
  rows.append(dict(point=p,source_sha256=sha,derived_history=str(cache),mean=dict(zip(names,x.mean(0).tolist())),halves=[dict(zip(names,v.mean(0).tolist())) for v in np.array_split(x,2)],blocks=nb,block_cycles=width,conditional_joint_block_MC95=ci,flags=['low_block_count'] if nb<20 else [],scope='One chain per initial class; no pooling or equilibrium inference.'))
  print('derived',p['id'],flush=True)
 result=dict(created=r.pt.legacy.now(),analysis_sha256=scriptsha,rows=rows,scope='ALL4001configuration rows per trajectory derived; original files unchanged. Fourfold spin harmonic is distinct from magnetization_fourth. Joint same-row bootstrap retains covariance. Strong finite-box harmonic responses do not alone establish infinite-volume phase order.')
 r.atomic(P/'harmonic-analysis.json',result);sh=P/'analysis-source-snapshots';sh.mkdir(exist_ok=True);import shutil;shutil.copyfile(__file__,sh/(scriptsha+'.py'))
if __name__=='__main__':main()
