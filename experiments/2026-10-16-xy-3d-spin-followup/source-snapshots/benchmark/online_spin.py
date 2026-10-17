"""Pure online measurement extension; frozen numerical kernels are unchanged."""
from pathlib import Path
from contextlib import contextmanager
import hashlib,json,time
import numpy as np
from supervillain import xy

_ORIGINAL_MEASURE=xy.measure
_ORIGINAL_PROVENANCE=xy._provenance

def definition():
    return dict(version=1,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),D=3,separations=['L/4','L/2'],cadence_cycles=1,directions='three positive coordinate axes',observable='mean_x exp(i*(phi(x+r)-phi(x))), complex axes retained; primary real axis mean',resume='Use this extension; native unwrapped resume intentionally fails provenance guard')

def correlations(phi):
    p=np.asarray(phi);p=p[0] if p.ndim==4 and p.shape[0]==1 else p
    if p.ndim!=3 or len(set(p.shape))!=1 or p.shape[0]%4 or not np.isfinite(p).all():raise ValueError('Finite cubic3D field with L divisible by4 required')
    spin=np.exp(1j*p);conj=np.conj(spin);out=np.empty((2,3),complex)
    for j,r in enumerate([p.shape[0]//4,p.shape[0]//2]):
        for axis in range(3):out[j,axis]=np.mean(conj*np.roll(spin,-r,axis=axis))
    return out

@contextmanager
def installed():
    if xy.measure is not _ORIGINAL_MEASURE or xy._provenance is not _ORIGINAL_PROVENANCE:raise RuntimeError('Another measurement/provenance hook is active')
    spec=definition()
    def measured(action,phi):
        values=_ORIGINAL_MEASURE(action,phi);t=time.perf_counter();g=correlations(phi)
        values.update(spin_correlation_real=g.real,spin_correlation_imag=g.imag,spin_extra_measurement_seconds=time.perf_counter()-t)
        return values
    def provenance():
        value=_ORIGINAL_PROVENANCE();base=value['source_sha256'];value['base_kernel_source_sha256']=base;value['spin_measurement_extension']=spec
        value['source_sha256']=hashlib.sha256(json.dumps([base,spec],sort_keys=True).encode()).hexdigest()
        return value
    xy.measure=measured;xy._provenance=provenance
    try:yield
    finally:xy.measure=_ORIGINAL_MEASURE;xy._provenance=_ORIGINAL_PROVENANCE

def run_chain(path,**kwargs):
    path=Path(path)
    if kwargs.get('resume',False):
        import h5py
        with h5py.File(path,'r') as h:
            if json.loads(h.attrs['provenance']).get('spin_measurement_extension')!=definition():raise ValueError('Online spin extension mismatch; refusing resume before writes')
    with installed():return xy.run_chain(path,**kwargs)
