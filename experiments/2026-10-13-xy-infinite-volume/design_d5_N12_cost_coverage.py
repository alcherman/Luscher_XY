"""Stage N12 at two measured D5 scans and one missing mid-beta N10 anchor."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json
from pathlib import Path
import numpy as np
import run
import storage_admission
OUT=Path(__file__).resolve().parent
name='d5-N12-q1-endpoint-cost-and-coverage-20261014'
points=[run.point(5,100,b,12,r,4000) for b in [224,228,232] for r in [8,9]]
points += [run.point(5,100,228,10,r,8000) for r in [8,9]]
for e in [894,896,898]:
 for r in [8,9]:
  p=run.point(5,e//10,0,12,r,4000);p.update(id=f'zero_d5_e{e:04d}_N012_r{r:02d}_s004000',epsilon_over_pi=e/1000,epsilon=e*np.pi/1000,seed=int(np.random.SeedSequence([2026101403,5,e,12,r]).generate_state(1,dtype=np.uint64)[0]),campaign_role='zero_beta_endpoint',epsilon_sampling_milli=e);points.append(p)
old=[p for f in (OUT/'batches').glob('*.json') for p in json.loads(f.read_text())['points'] if f.stem!=name]
assert not {p['id'] for p in points}&{p['id'] for p in old}
assert not {p['seed'] for p in points}&{p['seed'] for p in old}
source=OUT/json.loads((OUT/'latest-d5-calibration.json').read_text())['folder']/'analysis.json'
s=json.loads((OUT/'status.json').read_text());estimate=sum(storage_admission.estimated_full_bytes(p) for p in points)
# Reserve all previously pending work in addition to live growth/copy budgets.
pending=[]
for p in old:
 f=OUT/'jobs'/(p['id']+'.json');j=json.loads(f.read_text()) if f.exists() else {'state':'pending'}
 if j['state'] in ['pending','interrupted'] and p['D']==5:pending.append(p)
pending_bytes=sum(storage_admission.estimated_full_bytes(p) for p in pending)
headroom=s['physical_disk_free_bytes']-s['active_reserved_bytes']-s['retention_reserved_bytes']-pending_bytes-estimate-25*2**30
assert headroom>0
purpose='After audited N10 cost/coverage, stage modest N12 only at q1 beta.224/.228/.232 and beta0epsilon/pi.894/.896/.898,4000cycles,fresh8/9. Add N10q1beta.228,8000cycles,twofresh8/9 to test the narrow central response: previous.225:.240 interpolation gives.2316 but response is sharply curved and not finalprecision. N12timing/storage/peakmemory unmeasured; these are cost/coverage probes, not fullgridcommitment. Defer N12q92 until lowerbrackets measured. Dedicatedepsilonroles andexactcoords preserved; frozen kernels,fullclosureaudit,allhistories/logs/checkpoints retained,authorizedfieldretention only.'
report=dict(created=run.now(),storage_estimator_sha256=run.sha256(Path(storage_admission.__file__)),trajectories=len(points),physical_coordinates=len(points)//2,points=points,source_analysis=str(source),source_sha256=run.sha256(source),transient_full_GiB=estimate/2**30,pending_prior_full_GiB=pending_bytes/2**30,headroom_after_all_prior_pending_active_copy_and_25GiB_reserve_GiB=headroom/2**30,status_snapshot=s,purpose=purpose,planning='N12 runtime is extrapolated from actualN10 withvolumeandcyclefactors; not measured. Inspect coverage, cost, currentstorage and mixing before moreN12.')
run.atomic_json(OUT/(name+'-design.json'),report)
if __name__=='__main__':
 import sys
 if '--freeze' in sys.argv:
  batch=json.loads((OUT/'batches/d5-local-bracket-and-N10-refinement-20261014.json').read_text())
  for p in batch['points']:
   j=json.loads((OUT/'jobs'/(p['id']+'.json')).read_text());assert j['state']=='complete' and j['audit']['passed']
  a=json.loads(source.read_text());assert all(any(r['axis']==axis and r['q']==q and r['N']==10 and r['status']=='conditional_crossing' for r in a['roots']) for axis,q in [('beta',1.),('epsilon',0)])
  legacy=run.estimated_bytes
  try:
   run.estimated_bytes=storage_admission.estimated_full_bytes
   run.freeze(name,points,purpose)
  finally:run.estimated_bytes=legacy
 print(json.dumps({k:report[k] for k in ['trajectories','physical_coordinates','transient_full_GiB','pending_prior_full_GiB','headroom_after_all_prior_pending_active_copy_and_25GiB_reserve_GiB']}))
