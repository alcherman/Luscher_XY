"""Bracket extensions supported by completed low-epsilon grid slices."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json
from pathlib import Path
import run
OUT=Path(__file__).resolve().parent
source=OUT/json.loads((OUT/'latest-d2-eta-correlation.json').read_text())['folder']/'analysis.json'
grid=json.loads((OUT/'batches/d2-eta-correlation-remaining-grid-20261014.json').read_text())
choices={58:[20,280],62:[300],66:[640]}
for p in grid['points']:
    if round(p['epsilon_over_pi']*100) in choices:
        job=json.loads((OUT/'jobs'/f"{p['id']}.json").read_text())
        assert job['state']=='complete' and job['audit']['passed']
points=[run.point(2,q,b,N,r,12000) for q,betas in choices.items() for b in betas for N in [64,128,256] for r in [4,5]]
old=[p for f in (OUT/'batches').glob('*.json') for p in json.loads(f.read_text())['points']]
assert not {p['id'] for p in points}&{p['id'] for p in old}
assert not {p['seed'] for p in points}&{p['seed'] for p in old}
name='d2-eta-correlation-grid-bracket-refinement-20261014'
why='Completed q.58/.62/.66 baseline slices have incomplete bootstrap brackets: eta128 q.58 valid.954; eta64/128 q.62 .9055/.987; R128 and eta64 q.66 .975/.965. Extend q.58 beta to.020/.280 because both shallow eta endpoint margins are small; q.62 lowerbeta.300; q.66 upperbeta.640. Match N64/128/256 for joint eta/R, use fresh independent replicas4/5. Preserve every raw field history RNG checkpoint log and audit. This improves measured bracket coverage only, not a claimed infinite-volume accuracy improvement.'
m=run.freeze(name,points,why)
run.atomic_json(OUT/(name+'-design.json'),dict(created=run.now(),source=str(source),source_sha256=run.sha256(source),choices=choices,trajectories=len(points),estimated_GiB=m['estimated_bytes']/2**30,purpose=why))
print(len(points),m['estimated_bytes']/2**30,flush=True)
