"""Six independent larger-size probes; measure cost before further commitments."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json
from pathlib import Path
import run
OUT=Path(__file__).resolve().parent
coordinates=[(60,281),(80,921),(100,1111)]
points=[run.point(2,q,b,512,r,12000) for q,b in coordinates for r in [6,7]]
old=[p for f in (OUT/'batches').glob('*.json') for p in json.loads(f.read_text())['points']]
assert not {p['id'] for p in points}&{p['id'] for p in old}
assert not {p['seed'] for p in points}&{p['seed'] for p in old}
for q,b in coordinates:
    matches=[p for p in old if p['D']==2 and p['N']==256 and round(p['epsilon_over_pi']*100)==q and round(p['beta']*1000)==b]
    assert len(matches)>=2
    assert all(json.loads((OUT/'jobs'/f"{p['id']}.json").read_text())['state']=='complete' for p in matches)
name='d2-eta-correlation-N512-cost-pilot-20261014'
purpose='Independently reviewed six-run D2 N512 cost/observable pilot, matching completed N256 coordinates q=.60 beta=.281,q=.80 beta=.921,q=1 beta=1.111,12000cycles,freshreplicas6/7. Measure actual generation/audit runtime, stored bytes, R512, eta256:512, autocorrelation and start sensitivity before additional N512 sampling. One beta per slice cannot bracket crossings or validate infinite-volume precision. Motivated by unresolved BKT size dependence, including known q1 benchmark miss under baseline fixedoffset0. Preserve every field, history, RNG, checkpoint, provenance and full log; old316D2Binderjobsremain deferred,D5disabled.'
m=run.freeze(name,points,purpose)
run.atomic_json(OUT/(name+'-design.json'),dict(created=run.now(),coordinates=coordinates,N=512,replicas=[6,7],cycles=12000,trajectories=6,estimated_GiB=m['estimated_bytes']/2**30,purpose=purpose,review='Independent read-only review verified exact N256 matches, unique IDs/seeds, storage reserve and limited inference. No blocker.'))
print(name,len(points),m['estimated_bytes']/2**30,flush=True)
