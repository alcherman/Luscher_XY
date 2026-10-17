"""Replace broad small-box endpoint interpolation with shared local coordinates."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
from pathlib import Path
import json,numpy as np
import run,storage_admission as budget
OUT=Path(__file__).resolve().parent
name='d5-zero-beta-local-window-refinement-20261014'
points=[]
for N in [6,8,10]:
 for e in [893,894]:
  for r in [14,15]:
   p=run.point(5,e//10,0,N,r,4000);p.update(id=f'zero_d5_e{e:04d}_N{N:03d}_r{r:02d}_s004000',epsilon_over_pi=e/1000,epsilon=e*np.pi/1000,seed=int(np.random.SeedSequence([2026101403,5,e,N,r]).generate_state(1,dtype=np.uint64)[0]),campaign_role='zero_beta_endpoint',epsilon_sampling_milli=e);points.append(p)
old=[p for f in (OUT/'batches').glob('*.json') for p in json.loads(f.read_text())['points'] if f.stem!=name]
assert len(points)==12 and len({p['id'] for p in points})==12
assert not {p['id'] for p in points}&{p['id'] for p in old};assert not {p['seed'] for p in points}&{p['seed'] for p in old}
source=OUT/'d5-calibration-analysis-20261014-155318-cut1001-block1-U25-local-all/analysis.json'
a=json.loads(source.read_text());assert all(any(r['axis']=='epsilon' and r['N']==N and r['method']=='linear' and r['selected_window']==[.892,.9] for r in a['roots']) for N in [6,8,10])
s=json.loads((OUT/'status.json').read_text());pending=[]
for p in old:
 if p['D']!=5:continue
 f=OUT/'jobs'/(p['id']+'.json');j=json.loads(f.read_text()) if f.exists() else {'state':'pending'}
 if j['state'] in ['pending','interrupted']:pending.append(p)
estimate=sum(budget.estimated_full_bytes(p) for p in points);prior=sum(budget.estimated_full_bytes(p) for p in pending)
headroom=s['physical_disk_free_bytes']-s['active_reserved_bytes']-s['retention_reserved_bytes']-prior-estimate-25*2**30;assert headroom>0
purpose='Twelve4k endpointtrajectories at exactepsilon/pi.893/.894,N6/8/10,fresh14/15,matchingN12existing.894/queued.893. Currentlower-sizequarterroots interpolate broad.892:.900 (U.892=.363/.432/.477 versusU.900~0). These two sharedcoordinates testlocalcurvature and closerbrackets before infinityinference; notprecisionchasingormixingcertification. N12lower-epsilonanchors alreadyfrozenelsewhere; no duplicateN12. Preservezero_beta_endpointrole, exactradians, allhistories/checkpoints/logs, fullfieldaudit beforeauthorizedretention. No modelkernelchanges,no newD2. Future24N12remainingqs stilluncommitted; refreshstoragebeforethem, nevercountfuturethinningsavings.'
report=dict(created=run.now(),trajectories=12,physical_coordinates=6,points=points,source=str(source),source_sha256=run.sha256(source),storage_estimator_sha256=run.sha256(Path(budget.__file__)),full_transient_GiB=estimate/2**30,prior_pending_GiB=prior/2**30,headroom_after_all_known_reservations_and_25GiB_reserve_GiB=headroom/2**30,status_snapshot=s,purpose=purpose)
run.atomic_json(OUT/(name+'-design.json'),report)
if __name__=='__main__':
 import sys
 if '--freeze' in sys.argv:
  legacy=run.estimated_bytes
  try:run.estimated_bytes=budget.estimated_full_bytes;run.freeze(name,points,purpose)
  finally:run.estimated_bytes=legacy
 print(json.dumps({k:report[k] for k in ['trajectories','physical_coordinates','full_transient_GiB','prior_pending_GiB','headroom_after_all_known_reservations_and_25GiB_reserve_GiB']}))
