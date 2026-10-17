"""Descriptive moment/topology/start and actual cost diagnostics, no mixing proof."""
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np,h5py
import measure as m

def tau(x):
    value=m.helper.tau_int(x)
    return float(value) if value is not None and np.isfinite(value) else None

def main():
    plan=json.loads((m.OUT/'larger-zero-beta-pilot-design.json').read_text());ps=plan['points'];folder=m.OUT/'larger-pilot-derived-histories';folder.mkdir(exist_ok=True)
    def one(p):
        path=m.CAMPAIGN/'chains'/f"{p['id']}.h5";j=json.loads((m.CAMPAIGN/'jobs'/f"{p['id']}.json").read_text());assert j['state']=='complete' and j['audit']['passed']
        with h5py.File(path) as h:
            rows=len(h['sweep']);wm=np.empty((rows,3));w2=np.empty(rows)
            for first in range(0,rows,128):
                stop=min(rows,first+128);w=h['winding'][first:stop];wm[first:stop]=w.mean(axis=-1);w2[first:stop]=np.square(w,dtype=float).mean(axis=(1,2))
            accepted=h['winding_accepted'][:];vortex=h['vortex_density'][:];moments=np.column_stack([h[k][:] for k in ['magnetization_squared','magnetization_fourth','magnetization_abs']])
        np.savez_compressed(folder/f"{p['id']}.npz",directional_mean_winding=wm,mean_winding_squared=w2,vortex_density=vortex,winding_accepted=accepted)
        cuts=[]
        for cut in [1001,4001,8001]:
            a=moments[cut:];m2,m4,ma=a.T;influence=m4.mean()/m2.mean()**3*(m2-m2.mean())-(m4-m4.mean())/(2*m2.mean()**2)
            mt={k:tau(v) for k,v in dict(m2=m2,m4=m4,mabs=ma,binder=influence).items()};tt={**{f'mean_winding_axis{i}':tau(wm[cut:,i]) for i in range(3)},'mean_winding_squared':tau(w2[cut:]),'vortex_density':tau(vortex[cut:])}
            flags=[];moment_max=max([1.]+[v for v in mt.values() if v is not None]);topology_max=max([1.]+[v for v in tt.values() if v is not None])
            if topology_max>2*moment_max:flags.append('topology_tau_exceeds_twice_moment_tau_review_block_uncertainty')
            if np.max(np.ptp(wm[cut:],axis=0))==0:flags.append('no_directional_mean_winding_change_observed')
            cuts.append(dict(cut=cut,rows=len(a),moment_tau_cycles=mt,topology_tau_cycles=tt,directional_mean_winding=wm[cut:].mean(axis=0).tolist(),directional_winding_ranges=np.column_stack([wm[cut:].min(axis=0),wm[cut:].max(axis=0)]).tolist(),constant_directional_winding_axes=[i for i in range(3) if np.ptp(wm[cut:,i])==0],mean_winding_squared=float(w2[cut:].mean()),cycles_with_accepted_winding_update=int((accepted[cut:]>0).sum()),accepted_winding_updates=int(accepted[cut:].sum()),flags=flags))
        return dict(id=p['id'],N=p['N'],epsilon_over_pi=p['epsilon_over_pi'],replica=p['replica'],initial_winding=p['initial_winding'],full_audit=j['audit'],generation_seconds=j.get('attempt_generation_seconds'),audit_seconds=j.get('audit_seconds'),current_bytes=path.stat().st_size,cut_diagnostics=cuts)
    with ThreadPoolExecutor(max_workers=2) as pool:records=list(pool.map(one,ps))
    m.run.atomic_json(m.OUT/'larger-pilot-diagnostics.json',dict(created=m.run.now(),source_point_ids=sorted(p['id'] for p in ps),records=records,scope='Descriptive accepted-event denominators, autocorrelation and start/late-cut diagnostics; audits/moves/tau estimates do not prove mixing, independent sampling or absolute exponent accuracy. Resource samples are RSS/CPU observations, not a measured peak.'))
    print(len(records),'audited trajectory diagnostics',flush=True)

if __name__=='__main__':main()
