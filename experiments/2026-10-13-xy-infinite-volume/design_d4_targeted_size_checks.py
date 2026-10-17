"""Frozen local interpolation and N24 sensitivity at three flagged d4 slices."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json
from pathlib import Path
import run
OUT=Path(__file__).resolve().parent
existing={p['id'] for f in (OUT/'batches').glob('*.json') for p in json.loads(f.read_text())['points']}
choices={88:{20:[42,44,46],24:[39,42,45]},92:{20:[145,147,149],24:[145,147,149]},96:{20:[230,232,234],24:[230,232,234]}}
points=[run.point(4,q,b,N,r,4000) for q,sizes in choices.items() for N,betas in sizes.items() for b in betas for r in [0,1]]
points=[p for p in points if p['id'] not in existing]
name='binder-quarter-d4-targeted-N20-N24-sensitivity-20261013-2000'
m=run.freeze(name,points,'Local beta spacing and N24 sensitivity at epsilon/pi=.88,.92,.96, where leading all-size shifts fail strict goodness-of-fit despite small absolute deviations. N20 narrower windows distinguish interpolation curvature from size-law sensitivity. New N24 three-beta brackets test extrapolation. Two independent starts,4000cycles,full raw retention; no d2 or d5. Not a prerequisite for arbitrary precision or a global crossover fit.')
run.atomic_json(OUT/(name+'-design.json'),dict(created=run.now(),choices=choices,points=len(points),estimated_GiB=m['estimated_bytes']/2**30,source_snapshot='threshold-curve-snapshot-20261013-195243/analysis.json',source_sha256=run.sha256(Path(__file__)),scope='Local design; N24 beta windows are predictions and may require conservative extension if unbracketed. No accepted critical coupling inferred from design.'))
print(name,len(points),m['estimated_bytes']/2**30,flush=True)
