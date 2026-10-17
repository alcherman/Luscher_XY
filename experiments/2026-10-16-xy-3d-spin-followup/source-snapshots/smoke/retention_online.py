"""Conservative extension of the authorized field policy to recorded online G.
All old-selected fields remain; G-based stretches/spacing add a union. All
observable/RNG/checkpoint histories are copied byte-exact by the frozen converter.
"""
import math,json
import numpy as np
import retain_fields as retained
_BASE=retained.selection

def selection(h,evidence=None):
    indices,base=_BASE(h,evidence)
    if 'spin_correlation_real' not in h:return indices,base
    n=len(h['sweep']);cut=min(json.loads(h.attrs['parameters'])['warmup']+1,n//2);G=h['spin_correlation_real'][:].mean(2);taus={}
    for j,name in enumerate(['Gquarter','Ghalf']):
        series=G[cut:,j]
        for label,v in [('postwarmup',series),('latehalf',series[len(series)//2:])]:
            t=retained.helper.tau_int(v) if len(v)>20 else None;taus[name+'_'+label]=float(t) if t is not None and np.isfinite(t) else None
    tau=max([1.,base.get('max_measured_tau',1.)]+[v for v in taus.values() if v is not None]);stride=max(1,math.ceil(10*tau));width=max(64,math.ceil(50*tau));keep=set(map(int,indices));keep.update(range(cut,n,stride));windows=[]
    for first in [0,cut,cut+(n-cut)//2-width//2,n-width]:
        first=max(0,min(first,n-1));last=min(n,first+width);keep.update(range(first,last));windows.append([first,last])
    return np.array(sorted(keep),dtype=np.int64),dict(mode='adaptive_with_online_G_union',base_policy=base,online_G_tau_estimates=taus,max_measured_tau=tau,stride=stride,diagnostic_windows=windows,checkpoint_interval=100,all_base_selected_fields_preserved=True,independence='No independence guarantee; all online G and scalar histories remain complete.')

def compact(path,audit):
    retained.selection=selection
    try:return retained.compact(path,audit,parallel=True)
    finally:retained.selection=_BASE
