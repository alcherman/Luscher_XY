import copy
import tempfile
from pathlib import Path
import h5py
import numpy as np
from supervillain import Ensemble, Lattice, XY
from supervillain.generator.xy import ClusterUpdate, SiteUpdate, WindingUpdate
from supervillain.generator.combining import Sequentially
from supervillain.xy import run_chain

root = Path(tempfile.mkdtemp(prefix='xy-second-review-', dir='/private/tmp'))
print('artifacts', root)

# A real interrupted append: only a prefix of the datasets reaches row 4.
kwargs = dict(D=2, N=3, kappa=.7, epsilon=2.2, sweeps=7, warmup=1, seed=864,
              winding=True, checkpoint_every=2, update='cluster', clusters_per_step=2)
run_chain(root/'expected.h5', **kwargs)
original = h5py.Dataset.__setitem__
def fail_partial(self, key, value):
    if self.name == '/cluster_sites_total' and key == 4:
        raise RuntimeError('injected mid-row append fault')
    return original(self, key, value)
h5py.Dataset.__setitem__ = fail_partial
try:
    run_chain(root/'interrupted.h5', **kwargs)
except RuntimeError as exc:
    print('scan fault', str(exc))
finally:
    h5py.Dataset.__setitem__ = original
with h5py.File(root/'interrupted.h5') as f:
    print('partial rows', f.attrs['committed_rows'], {k: len(f[k]) for k in ['phi','cluster_sizes','cluster_sites_total','cluster_equivalent_sweeps']})
run_chain(root/'interrupted.h5', resume=True, **kwargs)
with h5py.File(root/'expected.h5') as a, h5py.File(root/'interrupted.h5') as b:
    print('partial resume exact', all(np.array_equal(a[k][...],b[k][...]) for k in a))

# Cross-restored lattice objects with compositions and old local APIs.
for name in ('cluster', 'site', 'combined'):
    action = XY(Lattice(2,3), .7, 2.2)
    g = ClusterUpdate(action, seed=886) if name == 'cluster' else SiteUpdate(action, seed=887, parallel=False)
    if name == 'combined':
        g = Sequentially([g, WindingUpdate(action,seed=889)])
    first = Ensemble(action).generate(3,g)
    path=root/f'{name}.h5'
    with h5py.File(path,'w') as f: first.to_h5(f.create_group('ensemble'))
    with h5py.File(path) as f: restored=Ensemble.from_h5(f['ensemble'])
    # Deepcopy avoids shared generator mutation between independent continuations.
    expected=Ensemble.continue_from(first,2)
    actual=Ensemble.continue_from(restored,2)
    print(name,'snapshot continuation exact',np.array_equal(expected.configuration.phi,actual.configuration.phi))

# Standard public append path writes configurations but may not write RNG state.
action=XY(Lattice(2,3),.7,2.2)
first=Ensemble(action).generate(3,ClusterUpdate(action,seed=63))
path=root/'ensemble-append.h5'
with h5py.File(path,'w') as f: first.to_h5(f.create_group('ensemble'))
second=Ensemble.continue_from(first,2)
with h5py.File(path,'r+') as f: second.extend_h5(f['ensemble'])
expected=Ensemble.continue_from(second,2)
with h5py.File(path) as f: restored=Ensemble.from_h5(f['ensemble'])
print('extended ensemble rows/cluster count',len(restored),restored.generator.clusters)
actual=Ensemble.continue_from(restored,2)
print('extended ensemble continuation exact',np.array_equal(expected.configuration.phi,actual.configuration.phi))
print('extended final counters expected/actual',expected.generator.clusters,actual.generator.clusters)
