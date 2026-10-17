"""Finite-box scouting summaries and adaptive design; no thermodynamic claims."""
import json
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
from pathlib import Path
import numpy as np
import run
from analyze_calibration import read,summarize

OUT=Path(__file__).resolve().parent


def json_finite(value):
    """Keep unavailable ratios explicit rather than emitting nonstandard NaN."""
    if isinstance(value,dict):return {k:json_finite(v) for k,v in value.items()}
    if isinstance(value,list):return [json_finite(v) for v in value]
    if isinstance(value,float) and not np.isfinite(value):return None
    return value


def main():
    folder=OUT/'scout-analysis';folder.mkdir(exist_ok=True)
    groups={}
    for f in (OUT/'jobs').glob('*.json'):
        j=json.loads(f.read_text());p=j['point']
        if j['state']=='complete' and p['epsilon_over_pi']<1 and j['audit']['passed']:
            groups.setdefault((p['D'],p['epsilon_over_pi'],p['beta'],p['N']),[]).append(p)
    rows=[]
    for (D,q,beta,N),points in sorted(groups.items()):
        if len(points)!=2:continue
        points=sorted(points,key=lambda p:p['replica'])
        data=[read(p) for p in points];central,boot=summarize(data,points,1001)
        if boot is None:continue
        rows.append(dict(D=D,epsilon_over_pi=q,beta=beta,N=N,summary=central,source_ids=[p['id'] for p in points]))
        np.savez_compressed(folder/f'd{D}_e{round(q*100):03d}_b{round(beta*1000):04d}_N{N:03d}.npz',samples=boot)
    choices=[]
    # Endpoint pilot values seed design only. They are not fixed critical ratios
    # in an infinite-volume estimator and no error bar is inferred from them.
    for D in [3,4]:
        with __import__('h5py').File(OUT/'chains'/f'd{D}_e100_b{455 if D==3 else 302:04d}_N008_r00_s004000.h5') as f:
            reference=np.sqrt(f['structure_zero'][1001:].mean()/f['structure_min'][1001:].mean()-1)/(16*np.sin(np.pi/8))
        for q in sorted({r['epsilon_over_pi'] for r in rows if r['D']==D}):
            subset=sorted([r for r in rows if r['D']==D and r['epsilon_over_pi']==q and r['N']==8],key=lambda r:r['beta'])
            hits=[(a,b) for a,b in zip(subset[:-1],subset[1:]) if a['summary']['values'][1]<=reference<=b['summary']['values'][1]]
            if len(hits)==1:
                lo,hi=hits[0][0]['beta'],hits[0][1]['beta']
                values=[round((lo+(hi-lo)*fraction)*1000) for fraction in [.2,.4,.6,.8]]
                choices.append(dict(D=D,epsilon_over_pi=q,design_beta_interval=[lo,hi],new_beta_milli=values,method='Finite N=8 R-reference design bracket only; not a critical-coupling estimate.'))
            else:
                choices.append(dict(D=D,epsilon_over_pi=q,status='no_interior_reference_bracket',reference=float(reference),method='Needs endpoint/branch validation; no invented positive-beta location.'))
    run.atomic_json(folder/'summary.json',json_finite(dict(created=run.now(),groups=rows,adaptive_choices=choices,source_sha256=run.sha256(Path(__file__)))))
    points=[run.point(c['D'],round(c['epsilon_over_pi']*100),b,N,r,4000) for c in choices if 'new_beta_milli' in c for b in c['new_beta_milli'] for N in [6,8] for r in [0,1]]
    # The first 3D candidate can already be on the ordered side at beta=0.
    points += [run.point(3,80,0,N,r,4000) for N in [12,16] for r in [0,1]]
    manifest=run.freeze('beta-neighborhood-scout',points,'Adaptive coupling neighborhoods from actual finite N=8 correlation-ratio brackets; N=6,8 four added beta points per bracket, plus zero-beta branch checks at d3 epsilon/pi=.80 N=12,16. These are design/proxy measurements, not infinite-volume results.')
    print('groups',len(rows),'adaptive_choices',json.dumps(choices),'new_trajectories',len(points),'estimated_GiB',manifest['estimated_bytes']/2**30,flush=True)


if __name__=='__main__':main()
