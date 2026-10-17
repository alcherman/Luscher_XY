"""Informed larger-size coverage of the 19 non-pilot epsilon slices."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json
from pathlib import Path
import run
OUT=Path(__file__).resolve().parent
source=OUT/json.loads((OUT/'latest-d2-eta-correlation.json').read_text())['folder']/'analysis.json'
a=json.loads(source.read_text());points=[];choices=[]
for q in range(58,101,2):
    if q in [60,80,100]:continue
    roots=[r['beta_location'] for r in a['locations'] if r['epsilon_over_pi']==q/100 and r['beta_location'] is not None and ((r['observable']=='xi_over_L' and r['N']==64) or (r['observable']=='eta_eff' and r['N']==32))]
    assert len(roots)==2
    desired=[min(roots)-(.05 if q<66 else .025),max(roots)+(.015 if q<66 else .01),max(roots)+(.08 if q<66 else .04)]
    obs=next(r for r in a['observations'] if r['epsilon_over_pi']==q/100 and r['N']==64 and r['observable']=='xi_over_L')
    reusable=[s for s in obs['sources'] if 'scout' not in s['source'] and not s.get('flags') and not any(s.get('individual_flags',{}).values())]
    selected=[]
    for b in desired:
        nearest=min(reusable,key=lambda s:abs(s['b']/1000-b)) if reusable else None
        reuse=nearest if nearest and abs(nearest['b']/1000-b)<=.008 else None
        milli=reuse['b'] if reuse else round(max(0,b)*1000)
        assert milli not in [r['beta_milli'] for r in selected]
        selected.append(dict(beta_milli=milli,reused_N64=reuse))
        for N in [64,128,256]:
            if N==64 and reuse:continue
            for replica in [2,3]:points.append(run.point(2,q,milli,N,replica,12000))
    choices.append(dict(epsilon_over_pi=q/100,existing_roots=roots,desired_beta=desired,selected=selected))
old=[p for f in (OUT/'batches').glob('*.json') for p in json.loads(f.read_text())['points']]
assert not {p['id'] for p in points}&{p['id'] for p in old}
assert not {p['seed'] for p in points}&{p['seed'] for p in old}
name='d2-eta-correlation-remaining-grid-20261014'
design=dict(created=run.now(),source=str(source),source_sha256=run.sha256(source),choices=choices,trajectories=len(points),reused_N64_groups=sum(s['reused_N64'] is not None for c in choices for s in c['selected']),estimated_GiB=sum(run.estimated_bytes(p) for p in points)/2**30,
    rationale='All three larger-box pilot slices now have >=.99 valid root draws at N64/128/256 R and 64:128/128:256 eta, with no pooled sampling flags. Broader low-epsilon brackets reflect observed shallow response. Matching old clean multi-start N64 groups are reused where within .008 of desired beta. Larger boxes use fresh starts. These are local coverage runs, not a guarantee of infinite-volume precision; refine or report unresolved when necessary.')
run.atomic_json(OUT/(name+'-design.json'),design)
if __name__=='__main__':
    import sys
    if '--freeze' in sys.argv:
        m=run.freeze(name,points,design['rationale']+' Preserve full fields, histories, RNG, checkpoints, logs and audits. Do not resume any old D2 Binder-quarter jobs. D5 remains disabled.')
        design['estimated_GiB']=m['estimated_bytes']/2**30;run.atomic_json(OUT/(name+'-design.json'),design)
    print(json.dumps({k:design[k] for k in ['trajectories','reused_N64_groups','estimated_GiB']}),flush=True)
