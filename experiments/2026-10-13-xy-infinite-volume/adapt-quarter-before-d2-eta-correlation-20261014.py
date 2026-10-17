"""Refine actual Binder-quarter brackets, extend missed brackets conservatively."""
import os
for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS']:
    os.environ[name]='1'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json
from pathlib import Path
import h5py
import numpy as np
import run

OUT=Path(__file__).resolve().parent


def main():
    grouped={}
    deferred_dimensions=json.loads((OUT/'control.json').read_text()).get('deferred_dimensions',[])
    for path in sorted((OUT/'jobs').glob('*.json')):
        j=json.loads(path.read_text());p=j['point']
        if j['state']!='complete' or not j['audit']['passed'] or p['D'] in deferred_dimensions:continue
        with h5py.File(OUT/'chains'/f"{p['id']}.h5") as f:
            m=[f['magnetization_squared'][1001:].mean(),f['magnetization_fourth'][1001:].mean()]
        grouped.setdefault((p['D'],round(p['epsilon_over_pi']*100),p['N'],round(p['beta']*1000)),[]).append(m)
    existing={p['id'] for f in (OUT/'batches').glob('*.json') for p in json.loads(f.read_text())['points']}
    points=[];choices=[]
    for D,q,N in sorted({k[:3] for k in grouped}):
        if D in deferred_dimensions:continue
        rows=[]
        for (d,e,n,b),moments in sorted(grouped.items()):
            if (d,e,n)==(D,q,N) and len(moments)>=2:
                m2,m4=np.mean(moments,axis=0);rows.append((b,1-m4/(2*m2*m2)))
        hits=[(a,b) for a,b in zip(rows[:-1],rows[1:]) if a[1]<=.25<b[1]]
        added=[];status='no_refinement_needed'
        limit=15 if D==2 else (8 if D==3 else (8 if N<=8 else 4))
        if len(hits)==1:
            lo,hi=hits[0][0][0],hits[0][1][0]
            if hi-lo>limit:
                added=sorted({round(lo+(hi-lo)*f) for f in [.25,.5,.75]});status='interior_refinement'
        elif max(r[1] for r in rows)<.25:
            # A missed larger-box upper endpoint; step upward, retain both old anchors.
            last=rows[-1][0];step=max(10,min(25,last-rows[-2][0])) if len(rows)>1 else 20
            added=[last+step,last+2*step];status='extend_upper_bracket'
        elif min(r[1] for r in rows)>.25:
            first=rows[0][0]
            if first>0:
                step=max(10,min(25,rows[1][0]-first)) if len(rows)>1 else 20
                added=sorted({max(0,first-step),max(0,first-2*step)});status='extend_lower_bracket'
            else:status='no_positive_quarter_crossing_yet'
        if added:
            for b in added:
                for r in [0,1]:
                    p=run.point(D,q,b,N,r,4000)
                    if p['id'] not in existing:points.append(p);existing.add(p['id'])
        choices.append(dict(D=D,epsilon_over_pi=q/100,N=N,status=status,new_beta_milli=added,
                            measured_beta_milli=[r[0] for r in rows],measured_U=[r[1] for r in rows]))
    name='binder-quarter-adaptive-followup-'+run.now()[:16].replace('-','').replace(':','').replace('T','-')
    m=run.freeze(name,points,'Refine unique directly measured Binder-quarter brackets by25/50/75percent subdivisions, or extend missed upper/lower brackets at modest beta steps. Zero-beta above-threshold boxes are not declared thermodynamically ordered. Two starts,4000cycles,complete raw retention. No critical-coupling acceptance from these simulation-design chords.')
    run.atomic_json(OUT/(name+'-design.json'),dict(created=run.now(),choices=choices,trajectories=len(points),estimated_GiB=m['estimated_bytes']/2**30,source_sha256=run.sha256(Path(__file__))))
    print(name,'trajectories',len(points),'estimated_GiB',m['estimated_bytes']/2**30,flush=True)


if __name__=='__main__':main()
