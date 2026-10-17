"""Experiment-local path-MH tempered transition, with all auxiliary RNG checkpointed.
Each attempted path applies an exactly palindromic sequence of signed reflection
kernels. Path log acceptance is k_target*(C_final-C_initial) plus the sum of
k_leg*(C_before-C_after). C is extensive sum over positive lattice links.
All legs use the same epsilon; no hard-support changes. Rejected paths, C values,
and work are recorded. Numerical cutoff failures propagate before mutation.
"""
from pathlib import Path
from contextlib import contextmanager
import hashlib,json,math,time
import numpy as np
from supervillain import xy
from supervillain.action.xy import XY
from supervillain.generator.xy import WindingUpdate
import signed_cluster as signed
_BASE_CLUSTER=xy.ClusterUpdate;_BASE_WINDING=xy.WindingUpdate;_BASE_MEASURE=xy.measure;_BASE_PROVENANCE=xy._provenance
ATTEMPT_PROBABILITY=.1
LADDER_STEPS=24
_LAST=None

def hot_coupling(target):return 0.0  # Empirically tested algorithmic bridge; not a temperature claim.
def ladder(target):
    up=np.linspace(float(target),hot_coupling(target),LADDER_STEPS+1)[1:]
    return np.concatenate([up,up[-2::-1]])
def cosine_sum(phi):
    field=np.asarray(phi)[0]
    return float(sum(np.cos(np.roll(field,-1,axis=a)-field).sum() for a in range(field.ndim)))
def path_log_ratio(target,levels,values):
    c=np.asarray(values,dtype=float);k=np.asarray(levels,dtype=float)
    if len(c)!=len(k)+1 or not np.array_equal(k,k[::-1]):raise ValueError('Exactly palindromic path with one more C value required')
    return float(target*(c[-1]-c[0])+np.dot(k,c[:-1]-c[1:]))
def empty_record(length):
    return dict(tempered_attempted=0,tempered_accepted=0,tempered_log_ratio=0.,tempered_C_path=np.zeros(length+1),tempered_kappa_path=np.zeros(length),tempered_work_clusters=0,tempered_work_sites=0,tempered_work_bond_tests=0,tempered_seconds=0.,native_twist_accepted=0)
class TemperedWinding(WindingUpdate):
    def __init__(self,action,seed=None):
        super().__init__(action,seed);self.levels=ladder(action.kappa)
        self.legs=[signed.SignedCluster(XY(action.Lattice,float(k),action.epsilon),seed=0) for k in self.levels]
        global _LAST;_LAST=empty_record(len(self.levels))
    def step(self,cfg):
        global _LAST
        before=self.accepted;cfg=super().step(cfg);record=empty_record(len(self.levels));record['native_twist_accepted']=self.accepted-before
        if self.rng.uniform()<ATTEMPT_PROBABILITY:
            start=time.perf_counter();record['tempered_attempted']=1;proposal=cfg['phi'];values=[cosine_sum(proposal)];sites=tests=0
            for update in self.legs:
                seed=int(self.rng.integers(self.Lattice.sites));alpha=self.rng.uniform(0,2*np.pi);tape=self.rng.uniform(size=(self.Lattice.D,self.Lattice.sites))
                proposal=update.cluster(proposal,seed,alpha,tape);values.append(cosine_sum(proposal));sites+=update.last_cluster_size;tests+=update.last_bond_tests
            log_ratio=path_log_ratio(self.Action.kappa,self.levels,values);uniform=float(self.rng.uniform());accepted=(uniform==0 or math.log(uniform)<min(0.,log_ratio))
            record.update(tempered_accepted=int(accepted),tempered_log_ratio=log_ratio,tempered_C_path=np.asarray(values),tempered_kappa_path=self.levels.copy(),tempered_work_clusters=len(self.levels),tempered_work_sites=sites,tempered_work_bond_tests=tests,tempered_seconds=time.perf_counter()-start)
            if accepted:self.accepted+=1;cfg=cfg|{'phi':proposal}
        _LAST=record
        return cfg

def definition():
    return dict(version=1,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),signed_cluster=signed.definition(),attempt_probability=ATTEMPT_PROBABILITY,ladder_steps=LADDER_STEPS,bridge_coupling=0.0,path='linear ladder excluding target, then exact palindrome with hot endpoint once',acceptance='k_target*(C_final-C_initial)+dot(k_legs,C_before-C_after)',C='sum positive-direction cos links, extensive',randomness='Every path choice uses native winding RNG; all rejected and accepted C paths retained',winding_accepted='Native record is native_twist_accepted+tempered_accepted; separate histories retained',checkpoint='Native fields/all RNGs at100-cycle commits; exact replay includes auxiliary RNG')
@contextmanager
def installed():
    if xy.ClusterUpdate is not _BASE_CLUSTER or xy.WindingUpdate is not _BASE_WINDING or xy.measure is not _BASE_MEASURE or xy._provenance is not _BASE_PROVENANCE:raise RuntimeError('Conflicting producer extension')
    spec=definition()
    def measure(action,phi):
        values=_BASE_MEASURE(action,phi)
        if _LAST is None:raise RuntimeError('Tempered diagnostics missing')
        values.update({key:(value.copy() if isinstance(value,np.ndarray) else value) for key,value in _LAST.items()});return values
    def provenance():
        value=_BASE_PROVENANCE();base=value['source_sha256'];value['base_kernel_source_sha256']=base;value['tempered_cluster_extension']=spec;value['source_sha256']=hashlib.sha256(json.dumps([base,spec],sort_keys=True).encode()).hexdigest();return value
    xy.ClusterUpdate=signed.SignedCluster;xy.WindingUpdate=TemperedWinding;xy.measure=measure;xy._provenance=provenance
    try:yield
    finally:xy.ClusterUpdate=_BASE_CLUSTER;xy.WindingUpdate=_BASE_WINDING;xy.measure=_BASE_MEASURE;xy._provenance=_BASE_PROVENANCE

def run_chain(path,**kwargs):
    if not kwargs.get('winding',True):raise ValueError('Tempered extension requires winding=True to construct/reset its checkpointed auxiliary update')
    path=Path(path)
    if kwargs.get('resume',False):
        import h5py
        with h5py.File(path) as h:
            if json.loads(h.attrs['provenance']).get('tempered_cluster_extension')!=definition():raise ValueError('Tempered extension mismatch; refusing resume before writes')
    with installed():return xy.run_chain(path,**kwargs)
