"""Source-pinned finite pilot/benchmark specifications for authorized follow-up."""
import sys,json,shutil,time
from pathlib import Path
import numpy as np
import campaign as c
P=c.OUT

def point(L,q,replica,steps,winding,role='pilot',online=True,seed=None):
    micro=int(round(q*1e6));identifier=f'spinfull3_{role}_q{micro:07d}_N{L:03d}_r{replica:02d}_s{steps:06d}'
    if seed is None:seed=int(np.random.SeedSequence([2026101602,3,micro,L,replica,steps]).generate_state(1,dtype=np.uint64)[0])
    return dict(id=identifier,D=3,N=L,beta=0.,epsilon_over_pi=micro/1e6,epsilon=micro*np.pi/1e6,replica=replica,seed=seed,initial_winding=winding,steps=steps,warmup=min(1000,steps//4),clusters_per_cycle=4,online=online,campaign_role=role)

def freeze(name,points,purpose):
    if len(set(p['id'] for p in points))!=len(points):raise ValueError('Duplicate identifiers')
    if any(p['D']!=3 or p['beta']!=0 or p['N'] not in [8,24,32,48,64] or p['N']%4 for p in points):raise ValueError('Unreviewed geometry/model')
    c.route();c.storage.OUT=P;lock=c.storage.acquire()
    try:budget=c.available_for(points)
    finally:lock.close()
    if not budget['admitted']:raise RuntimeError('Entire finite batch does not fit full physical reserve')
    body=dict(name=name,purpose=purpose,points=points,source_hashes=c.sources(),created=c.frozen.now(),storage_budget=budget,measurement=c.online.definition(),max_total_cpu_workers=int(c.read(P/'control.json')['max_total_cpu_workers']),reserve_bytes=c.RESERVE,new_larger_boxes_authorized=False)
    dest=P/'batches'/f'{name}.json'
    if dest.exists():raise FileExistsError(dest)
    c.atomic(dest,body);folder=P/'source-snapshots'/name;folder.mkdir()
    for f in ['online_spin.py','campaign.py','retention_online.py','design.py']:shutil.copyfile(P/f,folder/f)
    return dest
if __name__=='__main__':
    kind=sys.argv[1]
    if kind=='smoke':ps=[point(8,.80075,98,40,[1,0,0],role='smoke')]
    elif kind=='pilot':ps=[point(L,.80075,90+i,4000,w) for L in [24,32,48,64] for i,w in enumerate([[0,0,0],[1,0,0]])]
    elif kind=='benchmark':
        seed=202610160222;ps=[point(32,.80075,99,600,[0,0,0],role=role,online=flag,seed=seed) for role,flag in [('timing_native',False),('timing_online',True)]]
    else:raise ValueError(kind)
    print(freeze(kind,ps,'User-authorized '+kind+' for online long-distance G; no claim of equilibrium or exponent precision.'))
