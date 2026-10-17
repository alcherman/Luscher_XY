"""Pilot cost, dense observable coverage and autocorrelation diagnostics."""
import campaign as c
import json,time,math
import numpy as np
import h5py
P=c.OUT

def tau(v):
    t=c.retain.helper.tau_int(v) if len(v)>20 else None
    return float(t) if t is not None and np.isfinite(t) else None

def diagnose(job):
    p=job['point'];path=P/'chains'/f"{p['id']}.h5";cut=p['warmup']+1
    with h5py.File(path,'r') as h:
        g=h['spin_correlation_real'][cut:].mean(2);m=h['magnetization_squared'][cut:];z=h['structure_zero'][cut:];f=h['structure_min'][cut:].mean(1);w=np.empty(len(h['sweep'])-cut)
        for first in range(cut,len(h['sweep']),100):
            stop=min(first+100,len(h['sweep']));w[first-cut:stop-cut]=(h['winding'][first:stop].astype(float)**2).mean((1,2))
        ri=(z-z.mean())/f.mean()-z.mean()*(f-f.mean())/f.mean()**2
        series=dict(Gquarter=g[:,0],Ghalf=g[:,1],m2=m,R_influence=ri,vortex=h['vortex_density'][cut:],link_cos=h['mean_link_cosine'][cut:],winding2=w)
        taus={}
        for name,v in series.items():
            taus[name+'_postwarmup']=tau(v);taus[name+'_latehalf']=tau(v[len(v)//2:])
        t=max([1.]+[v for v in taus.values() if v is not None]);width=2*max(20,math.ceil(10*t));nb=len(m)//width
        R=np.sqrt(z.mean()/f.mean()-1)/(2*p['N']*np.sin(np.pi/p['N']))
        flags=[]
        if nb<50:flags.append('fewer_than_50_blocks_pilot_only')
        if nb<2:flags.append('insufficient_blocks_for_statistical_summary')
        halfmeans=[]
        for part in np.array_split(np.column_stack([g,m]),2):halfmeans.append(part.mean(0).tolist())
        stats=dict(global_twist_accepted_postwarmup=int(h['winding_accepted'][cut:].sum()),global_twist_attempts_postwarmup=len(m),constant_winding2=bool(np.ptp(w)==0),G_mean=g.mean(0).tolist(),m2_mean=float(m.mean()),R_mean=float(R),signed_configuration_G_nonpositive_counts=(g<=0).sum(0).tolist(),first_last_half_Gq_Gh_m2=halfmeans)
        if nb>=2:
            blocks=np.array([a.mean(0) for a in np.array_split(np.column_stack([g,m]),nb)]);stats['block_SE_Gq_Gh_m2']=(blocks.std(0,ddof=1)/np.sqrt(nb)).tolist()
        timing={}
        for key in ['metropolis_seconds','cluster_seconds','measurement_seconds','spin_extra_measurement_seconds','cycle_seconds']:
            v=h[key][cut:];timing[key]=dict(sum=float(v.sum()),mean=float(v.mean()),median=float(np.median(v)),p95=float(np.quantile(v,.95)))
    return dict(id=p['id'],N=p['N'],point=p,source_sha256=job['audit']['sha256'],generation_seconds=job['generation_seconds'],audit_seconds=job['audit_seconds'],physical_full_bytes=path.stat().st_size,physical_ceiling_bytes=c.storage.estimated_full_bytes(p),postwarmup_cycles=len(m),tau_cycles=taus,max_monitored_tau_cycles=t,requested_block_cycles=width,available_blocks=nb,flags=flags,statistics=stats,timing=timing)

def main():
    f=P/'batches/pilot-retry-v2.json';batch=c.read(f if f.exists() else P/'batches/pilot.json');jobs=[c.read(P/'jobs'/f"{p['id']}.json") for p in batch['points']]
    if any(j['state']!='complete' or not j['audit'].get('online_spin',{}).get('passed') for j in jobs):raise RuntimeError('Pilot incomplete or not fully audited')
    rows=[diagnose(j) for j in jobs];benchmark=[]
    for p in c.read(P/'batches/benchmark.json')['points']:
        j=c.read(P/'jobs'/f"{p['id']}.json")
        with h5py.File(P/'chains'/f"{p['id']}.h5") as h:
            start=p['warmup']+1;v=h['cycle_seconds'][start:];measurement=h['measurement_seconds'][start:]
            benchmark.append(dict(online=p['online'],generation_seconds=j['generation_seconds'],audit_seconds=j['audit_seconds'],postwarmup_cycle_mean=float(v.mean()),postwarmup_measurement_mean=float(measurement.mean()),scope='One serial same-RNG timing comparison, not a replicated performance bound; startup/JIT and I/O affect full wall time.'))
    profiles=[]
    for N in [24,32,48,64]:
        rr=[r for r in rows if r['N']==N];profiles.append(dict(N=N,generation_seconds_range=[min(r['generation_seconds'] for r in rr),max(r['generation_seconds'] for r in rr)],audit_seconds_range=[min(r['audit_seconds'] for r in rr),max(r['audit_seconds'] for r in rr)],max_tau_cycles=max(r['max_monitored_tau_cycles'] for r in rr),requested_blocks_at_100k_after_10k_cut=int(90000/max(r['requested_block_cycles'] for r in rr)),linear_generation_plus_audit_seconds_100k=max(r['generation_seconds']+r['audit_seconds'] for r in rr)*25))
    master=[dict(D=3,N=N,steps=100000) for N in [24,32,48,64] for _ in range(20)]
    worker_seconds=sum(20*r['linear_generation_plus_audit_seconds_100k'] for r in profiles)
    c.atomic(P/'pilot-analysis.json',dict(created=c.frozen.now(),source_batch=str(f if f.exists() else P/'batches/pilot.json'),rows=rows,benchmark=benchmark,profiles=profiles,worker_seconds=sum(r['generation_seconds']+r['audit_seconds'] for r in rows),full_bytes=sum(r['physical_full_bytes'] for r in rows),scope='Pilot cost/coverage/autocorrelation only; initial warmup and tau estimates do not certify equilibrium or target Delta accuracy.',provisional_100k_scenario=dict(trajectories=80,physical_full_ceiling_bytes=sum(c.storage.estimated_full_bytes(p) for p in master),worker_hours_linear=worker_seconds/3600,perfect14worker_wall_lower_bound_hours=worker_seconds/14/3600,scope='Scenario only. Cost extrapolated from4k at one epsilon; no larger/full main commitment before review, resource-contention/different-epsilon/retention/analysis margins not included.')))
    print(json.dumps(dict(profiles=profiles,worker_hours_linear=worker_seconds/3600,full_G_mean=[r['statistics']['G_mean'] for r in rows])),flush=True)
if __name__=='__main__':main()
