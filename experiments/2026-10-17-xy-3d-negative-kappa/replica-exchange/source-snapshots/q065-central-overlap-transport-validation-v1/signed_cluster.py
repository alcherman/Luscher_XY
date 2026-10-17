"""Signed-coupling reflection update, isolated from frozen producer/kernel files.
Boundary bonds activate with 1-exp(-max(0,DeltaS)) or with certainty if a
single-end reflection violates the hard constraint. Internal pair energies
are invariant. Both cut orientations retain the native strict guards.
"""
from pathlib import Path
from contextlib import contextmanager
import hashlib,json
import numpy as np
from supervillain import xy
from supervillain.action.xy import XY
from supervillain.generator.xy import ClusterUpdate as NativeCluster
from supervillain.generator.reference_implementation.xy_cluster import reference_cluster
from supervillain.lattice import Form
from _signed_cluster_kernel import signed_cluster_serial
_BASE_CLASS=xy.ClusterUpdate
_BASE_PROVENANCE=xy._provenance
class SignedCluster(NativeCluster):
    def __init__(self,action,seed=None,backend='numba'):
        if not isinstance(action,XY):raise ValueError('XY action required')
        super().__init__(XY(action.Lattice,abs(action.kappa),action.epsilon),seed=seed,backend=backend)
        self.Action=action;self.kappa=float(action.kappa)
    def cluster(self, phi, seed_site, reflection_angle, uniforms):
        """Replay one supplied-randomness cluster and return a fresh 0-form.

        ``reflection_angle`` specifies the normal ``r=(cos(alpha),sin(alpha))``.
        The independent uniform for the positive bond ``x -> x+mu`` is
        ``uniforms[mu, x]``, with flat site index ``x``.  Neither the field nor
        the supplied random tape is mutated.
        """
        self._validate_phi(phi)
        seed_site, reflection_angle, uniforms = self._validate_randomness(
            seed_site, reflection_angle, uniforms,
        )
        updated = Form(
            np.array(np.asarray(phi), dtype=np.float64, copy=True, order="C"),
            degree=0,
            lattice=phi.lattice,
        )
        kernel = reference_cluster if self.backend == "python" else signed_cluster_serial
        size, tests = kernel(
            np.asarray(updated).reshape(-1), self.plus, self.minus,
            seed_site, reflection_angle, uniforms, self.kappa, self.epsilon,
        )
        self.last_cluster_size = int(size)
        self.last_bond_tests = int(tests)
        self.clusters += 1
        self.flipped_sites += int(size)
        self.bond_tests += int(tests)
        return updated


def definition():
    p=Path(__file__)
    return dict(version=1,signed_cluster_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),kernel_sha256=hashlib.sha256(p.with_name('_signed_cluster_kernel.py').read_bytes()).hexdigest(),activation='forced hard bonds; else 1-exp(-2*max(0,kappa*projection_product))',scope='Experiment-local signed-coupling extension. Frozen original kernels unchanged.')
@contextmanager
def installed():
    if xy.ClusterUpdate is not _BASE_CLASS or xy._provenance is not _BASE_PROVENANCE:raise RuntimeError('Conflicting producer extension')
    spec=definition()
    def provenance():
        value=_BASE_PROVENANCE();base=value['source_sha256'];value['base_kernel_source_sha256']=base;value['signed_cluster_extension']=spec
        value['source_sha256']=hashlib.sha256(json.dumps([base,spec],sort_keys=True).encode()).hexdigest();return value
    xy.ClusterUpdate=SignedCluster;xy._provenance=provenance
    try:yield
    finally:xy.ClusterUpdate=_BASE_CLASS;xy._provenance=_BASE_PROVENANCE

def run_chain(path,**kwargs):
    path=Path(path)
    if kwargs.get('resume',False):
        import h5py
        with h5py.File(path) as h:
            if json.loads(h.attrs['provenance']).get('signed_cluster_extension')!=definition():raise ValueError('Signed cluster extension mismatch; no writes permitted')
    with installed():return xy.run_chain(path,**kwargs)
