"""Frozen larger-size 2D eta/R pilot with matching beta values for size pairs."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json
from pathlib import Path
import run
OUT=Path(__file__).resolve().parent
pointer=json.loads((OUT/'latest-d2-eta-correlation.json').read_text());source=OUT/pointer['folder']/'analysis.json';a=json.loads(source.read_text())
points=[];choices=[]
existing=[v for f in (OUT/'batches').glob('*.json') for v in json.loads(f.read_text())['points']]
for q in [60,80,100]:
    roots=[r['beta_location'] for r in a['locations'] if r['epsilon_over_pi']==q/100 and ((r['observable']=='xi_over_L' and r['N']==64) or (r['observable']=='eta_eff' and r['N']==32)) and r['beta_location'] is not None]
    if not roots:raise RuntimeError('No existing-data neighborhood')
    center=round(max(roots)*1000);betas=[center-10,center+10,center+30]
    for N in [64,128,256]:
        for b in betas:
            for replica in [2,3]:points.append(run.point(2,q,b,N,replica,12000))
    choices.append(dict(epsilon_over_pi=q/100,old_crossing_roots=roots,beta_milli=betas,N=[64,128,256],replicas=[2,3],cycles=12000))
assert not {p['id'] for p in points}&{p['id'] for p in existing}
assert not {p['seed'] for p in points}&{p['seed'] for p in existing},'Do not count repeated seeded prefixes as independent runs'
name='d2-eta-correlation-larger-size-pilot-20261014'
m=run.freeze(name,points,'User-directed2D eta_eff and xi/L pilot. Threeepsilon slices .60,.80,1.00; matchingbeta neighborhoods atN64/128/256 create direct64:128 and128:256 magnetic exponents, plus correlation ratios. N64 anchors supplement the larger-box pilot, avoiding interpolation between mismatchedbeta sizes. New independentreplicas2/3,12000cycles. Existing316Binder-quarter trajectories stay deferred. Preserve full fields histories RNG checkpoints logs and audits; provisional warmup is not mixing certification. Measure actual cost plus coverage beforefull22slicecommitment.')
run.atomic_json(OUT/(name+'-design.json'),dict(created=run.now(),source=str(source),source_sha256=run.sha256(source),choices=choices,trajectories=len(points),estimated_GiB=m['estimated_bytes']/2**30,design_sha256=run.sha256(Path(__file__)),primary_observables=['eta_eff','xi_over_L'],scope='Finite-size coverage and cost calibration only; no Binder primary criterion or accepted infinite-volume root.'))
print(name,len(points),m['estimated_bytes']/2**30,flush=True)
