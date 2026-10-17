"""Minimal direct curvature check inside the existing q1 Binder brackets."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
from pathlib import Path
import json
import run,storage_admission as budget
OUT=Path(__file__).resolve().parent
name='d5-q1-local-curvature-check-20261014'
points=[run.point(5,100,230,N,r,4000) for N in [10,12] for r in [12,13]]
source=OUT/'d5-calibration-analysis-20261014-155318-cut1001-block1-U25-local-all/analysis.json'
a=json.loads(source.read_text());assert all(any(r['axis']=='beta' and r['q']==1 and r['N']==N and r['status']=='conditional_crossing' for r in a['roots']) for N in [10,12])
old=[p for f in (OUT/'batches').glob('*.json') for p in json.loads(f.read_text())['points'] if f.stem!=name]
assert not {p['id'] for p in points}&{p['id'] for p in old}
assert not {p['seed'] for p in points}&{p['seed'] for p in old}
pending=[]
for p in old:
 if p['D']!=5:continue
 f=OUT/'jobs'/(p['id']+'.json');j=json.loads(f.read_text()) if f.exists() else {'state':'pending'}
 if j['state'] in ['pending','interrupted']:pending.append(p)
s=json.loads((OUT/'status.json').read_text());estimate=sum(budget.estimated_full_bytes(p) for p in points);prior=sum(budget.estimated_full_bytes(p) for p in pending)
headroom=s['physical_disk_free_bytes']-s['active_reserved_bytes']-s['retention_reserved_bytes']-prior-estimate-25*2**30;assert headroom>0
purpose='Four minimal q1 beta.230 N10/N12,4000cycles,fresh12/13. ExistingN10quarterroot interpolates .228:.240 (U.107:.499); N12uses .228:.232 (U.092:.496). Midpoint checks localcurvature and maygivecloserupperanchor; notguaranteed. Rejectedsizefit doesnotidentifyinterpolation asitscause; finite-sizecorrections/noise remainpossible. No globalfit requirement,no newD2. Physicalchunkbudgetincludesallpriorpending/live/copyand25GiBreserve; nofuturethinningsavings. Frozenkernels, fullclosureaudit andallhistories/checkpoints/logspreserved.'
report=dict(created=run.now(),trajectories=4,physical_coordinates=2,points=points,source=str(source),source_sha256=run.sha256(source),storage_estimator_sha256=run.sha256(Path(budget.__file__)),full_transient_GiB=estimate/2**30,prior_pending_GiB=prior/2**30,headroom_after_all_reservations_and_25GiB_GiB=headroom/2**30,status_snapshot=s,purpose=purpose)
run.atomic_json(OUT/(name+'-design.json'),report)
if __name__=='__main__':
 import sys
 if '--freeze' in sys.argv:
  legacy=run.estimated_bytes
  try:run.estimated_bytes=budget.estimated_full_bytes;run.freeze(name,points,purpose)
  finally:run.estimated_bytes=legacy
 print(json.dumps({k:report[k] for k in ['trajectories','physical_coordinates','full_transient_GiB','prior_pending_GiB','headroom_after_all_reservations_and_25GiB_GiB']}))
