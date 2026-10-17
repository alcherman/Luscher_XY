"""Narrow the measured D5 pilot brackets and add a modest N10 size check."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json
from pathlib import Path
import numpy as np
import run
OUT=Path(__file__).resolve().parent
points=[]
for q,betas in [(92,[70,90,110]),(100,[210,225,240])]:
 for N in [6,8,10]:
  for b in betas:
   for r in [4,5]:points.append(run.point(5,q,b,N,r,8000))
for e in [884,892,900]:
 for N in [6,8,10]:
  for r in [4,5]:
   p=run.point(5,e//10,0,N,r,8000);p.update(id=f'zero_d5_e{e:04d}_N{N:03d}_r{r:02d}_s008000',epsilon_over_pi=e/1000,epsilon=e*np.pi/1000,seed=int(np.random.SeedSequence([2026101403,5,e,N,r]).generate_state(1,dtype=np.uint64)[0]),campaign_role='zero_beta_endpoint',epsilon_sampling_milli=e);points.append(p)
old=[p for f in (OUT/'batches').glob('*.json') for p in json.loads(f.read_text())['points']]
assert not {p['id'] for p in points}&{p['id'] for p in old}
assert not {p['seed'] for p in points}&{p['seed'] for p in old}
name='d5-local-bracket-and-N10-refinement-20261014'
source=OUT/json.loads((OUT/'latest-d5-calibration.json').read_text())['folder']/'analysis.json'
purpose='After54D5calibrationtrajectoriespassedfullaudits, all9linearrootsbracket butcoarsebeta interpolation isnotfinalprecision. q.92 [.07,.09,.11] narrows measured.06:.12 neighborhood;q1 [.210,.225,.240] narrows.18:.28 neighborhood; beta0epsilon .884/.892/.900 narrows measured.88:.90 region. N6/8/10 addsmodestsizeleveragenoD4logformula. Fresh4/5 with8000cycles address7coordinate lowblockcount flags (notcertify mixing). Reusecompatible4000cyclepilotmoments with proper rawretainedweights, neveraverageratios. Preserveexactepsroles,numericalkernels,histories/checkpoints/logs/provenance,fullfieldaudit beforeauthorizedretention. NoadditionalD2precisionruns.'
estimate=sum(run.estimated_bytes(p) for p in points);s=json.loads((OUT/'status.json').read_text());headroom=s['physical_disk_free_bytes']-s['active_reserved_bytes']-s['retention_reserved_bytes']-estimate-25*2**30
assert headroom>0
r=dict(created=run.now(),source=str(source),source_sha256=run.sha256(source),trajectories=len(points),points=points,purpose=purpose,transient_full_GiB=estimate/2**30,headroom_GiB=headroom/2**30,forecast='Approximately1–2hcalendarallowance frommeasuredN8cost, doubledcycles andvolume-scaledN10 extrapolation plusI/O; noactualN10measurementyet. StagefullD5gridonlyafteractualcost,coverage andstoragecheck.')
run.atomic_json(OUT/(name+'-design.json'),r)
if __name__=='__main__':
 import sys
 if '--freeze' in sys.argv:run.freeze(name,points,purpose)
 print(json.dumps({k:r[k] for k in ['trajectories','transient_full_GiB','headroom_GiB']}))
