"""Extend observed endpoint brackets, preserving every pilot and seeded prefix."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json
from pathlib import Path
import numpy as np
import run
OUT=Path(__file__).resolve().parent
source=OUT/json.loads((OUT/'latest-zero-beta-endpoints.json').read_text())['folder']/'analysis.json'
choices={2:dict(epsilon_milli=[546,553],sizes=[64,128,256],steps=12000),
         3:dict(epsilon_milli=[802],sizes=[12,16,24,32],steps=8000),
         4:dict(epsilon_milli=[860,862],sizes=[10,12,16,20,24],steps=4000)}
points=[]
for D,c in choices.items():
    for e in c['epsilon_milli']:
        for N in c['sizes']:
            for replica in [6,7]:
                p=run.point(D,e//10,0,N,replica,c['steps'])
                p.update(id=f'zero_d{D}_e{e:04d}_N{N:03d}_r{replica:02d}_s{c["steps"]:06d}',
                         epsilon_over_pi=e/1000,epsilon=e*np.pi/1000,
                         seed=int(np.random.SeedSequence([2026101403,D,e,N,replica]).generate_state(1,dtype=np.uint64)[0]),
                         campaign_role='zero_beta_endpoint',epsilon_sampling_milli=e,
                         primary_observables=['eta_eff','xi_over_L'] if D==2 else ['binder_quarter'])
                points.append(p)
old=[p for f in (OUT/'batches').glob('*.json') for p in json.loads(f.read_text())['points']]
assert not {p['id'] for p in points}&{p['id'] for p in old}
assert not {p['seed'] for p in points}&{p['seed'] for p in old}
m=run.freeze('zero-beta-endpoint-bracket-refinement-20261014',points,
    'Adaptive beta-zero epsilon refinement from first audited snapshot. D2 root tails lack lower-epsilon coverage; add .546/.553 at matching N64/128/256. D3 add .802 inside coarse .800/.806 brackets. D4 add .860/.862 at uniform N10/12/16/20/24, covering unbracketed N10 and tightening existing roots near .862. Exact epsilon coordinates, endpoint role, fresh replicas6/7, full fields RNG histories logs checkpoints audits retained. Preliminary crossings do not certify thermodynamic precision. Old D2 Binder jobs remain deferred.')
run.atomic_json(OUT/'zero-beta-endpoint-bracket-refinement-20261014-design.json',dict(created=run.now(),source=str(source),source_sha256=run.sha256(source),choices=choices,trajectories=len(points),estimated_GiB=m['estimated_bytes']/2**30))
print(len(points),m['estimated_bytes']/2**30,flush=True)
