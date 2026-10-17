"""Fixed-epsilon replica exchange with joint all-slot/RNG checkpoints.
Within-slot updates are frozen Metropolis + four validated signed clusters +
frozen global twist. Fields and walker labels exchange; slot actions and RNGs
remain fixed. Complete histories for every thermodynamic slot are retained.
Reference: Hukushima & Nemoto, JPSJ65,1604(1996), doi10.1143/jpsj.65.1604.
"""
import os
for key in ['NUMBA_NUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import sys,json,time,hashlib
P=Path(__file__).resolve().parent;ROOT=P.parents[2];OLD=ROOT/'experiments/2026-10-13-xy-infinite-volume';sys.path[:0]=[str(P.parent),str(ROOT),str(OLD)]
import numpy as np,h5py
from supervillain import XY,Lattice
from supervillain.xy import measure,_provenance
from supervillain.generator.xy import SiteUpdate,WindingUpdate
from signed_cluster import SignedCluster,definition
import run as legacy
import importlib.util
_se=importlib.util.spec_from_file_location('sector_moves_isolated',P/'sector_moves.py');sector=importlib.util.module_from_spec(_se);_se.loader.exec_module(sector)

def runtime():
    p=_provenance();return {k:p[k] for k in ['python','numpy','numba','h5py','threads']}
def sources():return dict(native=legacy.source_hashes(),signed_cluster=definition(),replica_exchange_sha256=legacy.sha256(__file__),runtime=runtime(),sector_extension=sector.definition(),base_audit_sha256=legacy.sha256(P/'base_audit.py'))
def swap_log_ratio(ka,kb,ca,cb):return float((ka-kb)*(cb-ca))
def cosine_sum(phi):
 a=np.asarray(phi)[0];return float(sum(np.cos(np.roll(a,-1,axis=d)-a).sum() for d in range(a.ndim)))
def ladder(target,spacing=.025):
 if target>0:raise ValueError('Negative continuation only')
 if target==0:return np.array([0.])
 return np.linspace(target,0,max(2,int(np.ceil(abs(target)/spacing))+1))

def run(path,*,N,epsilon,kappas,steps,warmup,seed,initial_winding,checkpoint_every=100,initial_kind="winding",sector_updates=True,resume=False,inject_crash_after_flush=None):
 path=Path(path);ks=np.asarray(kappas,dtype=float);R=len(ks);D=3
 if R<1 or not np.isfinite(ks).all() or np.any(np.diff(ks)<0) or np.any(ks>0):raise ValueError('Finite nondecreasing negative-to-zero ladder required')
 if not 0<=warmup<steps:raise ValueError('Warmup smaller than total steps required')
 if initial_kind not in ['winding','checkerboard','z4_vortex'] or N%2:raise ValueError('Reviewed even boxes and explicit admissible initial kind required')
 params=dict(initial_kind=initial_kind,checkerboard_fraction=.98,z4_seed_stream=2026101711,D=3,N=int(N),epsilon=float(epsilon),kappas=ks.tolist(),warmup=int(warmup),seed=int(seed),initial_winding=list(initial_winding),clusters_per_cycle=4,checkpoint_every=int(checkpoint_every),update='replica_exchange_signed_hybrid',target_slot=0,sector_updates=bool(sector_updates),pi_clusters_per_cycle=4,sector_seed_stream=2026101713)
 source=sources();L=Lattice(3,N);actions=[XY(L,float(k),epsilon) for k in ks];streams=np.random.SeedSequence(seed).spawn(3*R+1)
 locals=[SiteUpdate(a,interval_phi=np.pi,seed=streams[3*i],backend='numba',parallel=False) for i,a in enumerate(actions)]
 clusters=[SignedCluster(a,seed=streams[3*i+1]) for i,a in enumerate(actions)]
 twists=[WindingUpdate(a,seed=streams[3*i+2]) for i,a in enumerate(actions)];exchange=np.random.default_rng(streams[-1])
 sector_rngs=[np.random.default_rng(np.random.SeedSequence([int(seed),2026101713,i])) for i in range(R)]
 coords=np.indices((N,)*3);initial=np.mod(sum(2*np.pi*w*coords[a]/N for a,w in enumerate(initial_winding)),2*np.pi)
 if initial_kind=='checkerboard':initial=.98*epsilon*(coords.sum(axis=0)%2)
 if initial_kind=='z4_vortex':
  if epsilon<=np.pi/2:raise ValueError('Z4 vortex initializer requires epsilon strictly above pi/2')
  initial_rng=np.random.default_rng(np.random.SeedSequence([int(seed),2026101711]));bits=initial_rng.integers(0,2,size=(N,)*3);initial=np.pi*bits+np.pi/2*(coords.sum(axis=0)%2);params['initialization_rng_state']=initial_rng.bit_generator.state
 cfgs=[{'phi':a.canonical(initial)} for a in actions]
 if not all(a.valid(c) for a,c in zip(actions,cfgs)):raise ValueError('Initial configuration is not admissible')
 params['initial_field_sha256']=hashlib.sha256(np.ascontiguousarray(cfgs[0]['phi']).tobytes()).hexdigest()
 walkers=np.arange(R,dtype=np.int64);pending=[]
 if path.exists() and not resume:raise FileExistsError(path)
 path.parent.mkdir(parents=True,exist_ok=True)
 with h5py.File(path,'r+' if resume else 'w') as h:
  if resume:
   if json.loads(h.attrs['parameters'])!=params or json.loads(h.attrs['source_hashes'])!=source:raise ValueError('Joint replica checkpoint parameters/source mismatch')
   committed=int(h.attrs['committed_rows']);assert committed>0
   if committed>steps+1:raise ValueError('Requested target precedes the joint checkpoint')
   for ds in h.values():
    if len(ds)<committed:raise ValueError('Dataset shorter than joint committed epoch')
   for ds in h.values():ds.resize(committed,axis=0)
   start=committed-1;states=json.loads(h['rng_state'][start]);cfgs=[{'phi':actions[i].canonical(h['phi'][start,i,0])} for i in range(R)];walkers=h['walker_id'][start].copy()
   for i in range(R):
    locals[i].rng.bit_generator.state=states['slots'][i]['local'];clusters[i].rng.bit_generator.state=states['slots'][i]['cluster'];twists[i].rng.bit_generator.state=states['slots'][i]['winding']
    clusters[i].clusters=int(h['clusters_total'][start,i]);clusters[i].flipped_sites=int(h['cluster_sites_total'][start,i]);clusters[i].bond_tests=int(h['cluster_bond_tests_total'][start,i])
   
   for i in range(R):sector_rngs[i].bit_generator.state=states['sector_rng'][i]
   exchange.bit_generator.state=states['exchange'];assert states['next_exchange_parity']==start%2
  else:
   h.attrs['runtime_provenance']=json.dumps(_provenance(),sort_keys=True);h.attrs['schema']='xy-replica-exchange-v1';h.attrs['parameters']=json.dumps(params,sort_keys=True);h.attrs['source_hashes']=json.dumps(source,sort_keys=True);h.attrs['committed_rows']=0;start=0
  def state(cycle):return json.dumps(dict(slots=[dict(local=locals[i].rng.bit_generator.state,cluster=clusters[i].rng.bit_generator.state,winding=twists[i].rng.bit_generator.state) for i in range(R)],exchange=exchange.bit_generator.state,next_exchange_parity=cycle%2,sector_rng=[g.bit_generator.state for g in sector_rngs]))
  def record(cycle,accepted,twist_accepted,sizes,tests,pre,attempted,swap_accepted,swap_uniform,logs,elapsed):
   measurements=[measure(a,c['phi']) for a,c in zip(actions,cfgs)]
   row={key:np.stack([m[key] for m in measurements]) for key in measurements[0]}
   spin=np.exp(1j*np.stack([np.asarray(c['phi'])[0] for c in cfgs]));parity=(-1.)**coords.sum(axis=0)
   row.update(staggered_nematic_squared=abs(((spin**2)*parity).mean(axis=(1,2,3)))**2,fourfold_squared=abs((spin**4).mean(axis=(1,2,3)))**2)
   row.update(staggered_squared=abs((spin*parity).mean(axis=(1,2,3)))**2,nematic_squared=abs((spin**2).mean(axis=(1,2,3)))**2)
   row.update(pi_cluster_sizes=np.array(pi_sizes),pi_cluster_bond_tests=np.array(pi_tests),**{'collective_'+key:np.array([d[key] for d in collective_info]) for key in collective_info[0]})
   row.update(phi=np.stack([np.asarray(c['phi']) for c in cfgs]),sweep=np.int64(cycle),rng_state=state(cycle),walker_id=walkers.copy(),local_accepted=np.asarray(accepted,dtype=np.int64),winding_accepted=np.asarray(twist_accepted,dtype=np.int64),cluster_sizes=np.asarray(sizes,dtype=np.int64),cluster_bond_tests=np.asarray(tests,dtype=np.int64),clusters_total=np.array([c.clusters for c in clusters]),cluster_sites_total=np.array([c.flipped_sites for c in clusters]),cluster_bond_tests_total=np.array([c.bond_tests for c in clusters]),pre_exchange_C=pre.copy(),post_exchange_C=np.array([cosine_sum(c['phi']) for c in cfgs]),exchange_attempted=attempted.copy(),exchange_accepted=swap_accepted.copy(),exchange_uniform=swap_uniform.copy(),exchange_log_ratio=logs.copy(),cycle_seconds=float(elapsed))
   pending.append(row)
  def flush(commit=True):
   if not pending:return
   first=len(h['sweep']) if 'sweep' in h else 0;stop=first+len(pending)
   for key in pending[0]:
    value=pending[0][key]
    if key not in h:
     if isinstance(value,str):h.create_dataset(key,shape=(0,),maxshape=(None,),dtype=h5py.string_dtype('utf-8'),chunks=(100,))
     else:
      arr=np.asarray(value);chunks=(1,1,1,N,N,N) if key=='phi' else (min(100,checkpoint_every),)+arr.shape
      # Empty exchange arrays at the one-slot kappa0 control need extendable
      # second dimension with a physical nonzero chunk.
      if 0 in arr.shape:h.create_dataset(key,shape=(0,)+arr.shape,maxshape=(None,)+tuple(None if x==0 else x for x in arr.shape),chunks=(min(100,checkpoint_every),)+tuple(max(1,x) for x in arr.shape),dtype=arr.dtype)
      else:h.create_dataset(key,shape=(0,)+arr.shape,maxshape=(None,)+arr.shape,chunks=chunks,dtype=arr.dtype)
    ds=h[key];ds.resize(stop,axis=0);ds[first:stop]=[x[key] for x in pending] if isinstance(value,str) else np.stack([x[key] for x in pending])
   h.flush()
   if not commit:raise RuntimeError('Injected crash after all-slot flush but before joint commit')
   h.attrs['committed_rows']=stop;h.flush();pending.clear()
  empty=np.zeros(max(0,R-1),dtype=np.int64)
  pi_sizes=np.zeros((R,4),int);pi_tests=pi_sizes.copy();collective_info=[dict(index=-1,uniform=0.,valid=0,accepted=0,log_ratio=0.,before_C=0.,candidate_C=0.) for i in range(R)]
  if not resume:
   pre=np.array([cosine_sum(c['phi']) for c in cfgs]);record(0,np.zeros(R,int),np.zeros(R,int),np.zeros((R,4),int),np.zeros((R,4),int),pre,empty,empty,empty.astype(float),empty.astype(float),0.);flush()
  for cycle in range(start+1,steps+1):
   tick=time.perf_counter();accepted=[];wa=[];ss=[];bt=[];pi_sizes=[];pi_tests=[];collective_info=[]
   for i in range(R):
    if sector_updates:cfgs[i],info=sector.collective(actions[i],cfgs[i],sector_rngs[i])
    else:info=dict(index=-1,uniform=0.,valid=0,accepted=0,log_ratio=0.,before_C=0.,candidate_C=0.)
    collective_info.append(info)
    before=locals[i].accepted;cfgs[i]=locals[i].step(cfgs[i]);accepted.append(locals[i].accepted-before);cs=[];ct=[]
    for _ in range(4):cfgs[i]=clusters[i].step(cfgs[i]);cs.append(clusters[i].last_cluster_size);ct.append(clusters[i].last_bond_tests)
    ps=[];ptests=[]
    for _ in range(4):
     if sector_updates:cfgs[i],p_size,p_tests=sector.pi_step(actions[i],cfgs[i],clusters[i].plus,clusters[i].minus,sector_rngs[i])
     else:p_size,p_tests=0,0
     ps.append(p_size);ptests.append(p_tests)
    pi_sizes.append(ps);pi_tests.append(ptests)
    ss.append(cs);bt.append(ct);before=twists[i].accepted;cfgs[i]=twists[i].step(cfgs[i]);wa.append(twists[i].accepted-before)
   pre=np.array([cosine_sum(c['phi']) for c in cfgs]);cvalues=pre.copy();attempted=empty.copy();swapped=empty.copy();uniforms=empty.astype(float);logs=empty.astype(float)
   for i in range((cycle-1)%2,R-1,2):
    lr=swap_log_ratio(ks[i],ks[i+1],cvalues[i],cvalues[i+1]);u=float(exchange.uniform());attempted[i]=1;uniforms[i]=u;logs[i]=lr
    if u==0 or np.log(u)<min(0.,lr):
     cfgs[i],cfgs[i+1]=cfgs[i+1],cfgs[i];walkers[i],walkers[i+1]=walkers[i+1],walkers[i];cvalues[i],cvalues[i+1]=cvalues[i+1],cvalues[i];swapped[i]=1
   record(cycle,accepted,wa,ss,bt,pre,attempted,swapped,uniforms,logs,time.perf_counter()-tick)
   if inject_crash_after_flush==cycle:flush(commit=False)
   if cycle%checkpoint_every==0 or cycle==steps:
    flush();print(f'joint checkpoint {cycle}/{steps}, slots={R}',flush=True)
 return path
