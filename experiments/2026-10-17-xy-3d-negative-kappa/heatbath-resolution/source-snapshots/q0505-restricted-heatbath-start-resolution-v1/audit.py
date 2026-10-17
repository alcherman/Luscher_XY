"""Independent sector histories plus unchanged full-field/native closure audit."""
from pathlib import Path
import importlib.util,json
import numpy as np,h5py
P=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('sector_base_audit',P/'base_audit.py');base=importlib.util.module_from_spec(s);s.loader.exec_module(base)
def audit(path,spec,source_hashes):
 result=base.audit(path,spec,source_hashes);n=spec['N'];enabled=spec.get('sector_updates',True);V=n**3
 with h5py.File(path) as f:
  params=json.loads(f.attrs['parameters']);assert params['sector_updates']==enabled and params['sector_seed_stream']==2026101713
  sizes=f['pi_cluster_sizes'][:];tests=f['pi_cluster_bond_tests'][:];assert not sizes[0].any() and not tests[0].any()
  if enabled:assert np.all((sizes[1:]>=1)&(sizes[1:]<=V)) and np.all((tests[1:]>=0)&(tests[1:]<=6*V))
  else:assert not sizes.any() and not tests.any()
  coords=np.indices((n,)*3);parity=coords.sum(0)%2;vectors=[[0,0,0]]
  for axis in range(3):
   for sign in [-1,1]:w=[0,0,0];w[axis]=sign;vectors.append(w)
  choices=[(sign*np.pi/2,w) for sign in [-1,1] for w in vectors]
  for t in range(len(f['sweep'])):
   state=json.loads(f['rng_state'][t]);assert len(state['sector_rng'])==len(spec['kappas'])
   for g in state['sector_rng']:assert g['bit_generator']=='PCG64'
   if not t:continue
   if not enabled:
    assert np.all(f['collective_index'][t]==-1);continue
   for i,kappa in enumerate(spec['kappas']):
    index=int(f['collective_index'][t,i]);assert 0<=index<14;angle,w=choices[index];old=f['phi'][t-1,i,0];candidate=(old+angle*parity-sum(2*np.pi*w[d]*coords[d]/n for d in range(3)))%(2*np.pi)
    links=np.array([(np.roll(candidate,-1,axis=d)-candidate+np.pi)%(2*np.pi)-np.pi for d in range(3)]);valid=np.all(abs(links)<spec['epsilon']);oldC=sum(np.cos(np.roll(old,-1,axis=d)-old).sum() for d in range(3));newC=np.cos(links).sum();lr=kappa*(newC-oldC);u=float(f['collective_uniform'][t,i])
    assert 0<=u<1;assert bool(f['collective_valid'][t,i])==valid;assert bool(f['collective_accepted'][t,i])==bool(valid and (u==0 or np.log(u)<min(0.,lr)))
    for name,x in [('before_C',oldC),('candidate_C',newC),('log_ratio',lr)]:np.testing.assert_allclose(f['collective_'+name][t,i],x,rtol=2e-12,atol=5e-9)
 
 with h5py.File(path) as f:
  assert json.loads(f.attrs['parameters'])['local_update']=='exact_restricted_heatbath'
  assert np.all(f['local_accepted'][1:]==V) and not f['local_accepted'][0].any()
  assert np.all(f['heatbath_draws'][1:]>=V) and not f['heatbath_draws'][0].any()
 result['sector_history_audit']=True
 result['heatbath_history_audit']=True
 return result
