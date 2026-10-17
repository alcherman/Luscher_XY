"""Final modest D5 size coverage with matched interpolation resolution."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json
from pathlib import Path
import run,storage_admission as budget
OUT=Path(__file__).resolve().parent
name='d5-final-N12-and-matched-midpoints-20261014'
source=OUT/'d5-calibration-analysis-20261014-165416-cut1001-block1-U25-local-all/analysis.json'
a=json.loads(source.read_text());assert len(a['source_jobs'])==234
windows={};points=[]
for q in [90,94,96,98]:
 r=next(r for r in a['roots'] if r['axis']=='beta' and r['q']==q/100 and r['N']==10 and r['method']=='linear')
 assert r['status']=='conditional_crossing' and r['valid_fraction']>=.99
 lo,hi=[round(v*1000) for v in [min(r['selected_window']),max(r['selected_window'])]];middle=round(r['location']*1000)
 assert lo<middle<hi
 windows[q]=[lo,middle,hi]
 for N in [6,8,10,12]:
  for beta in ([lo,middle,hi] if N==12 else [middle]):
   for replica in [16,17]:points.append(run.point(5,q,beta,N,replica,4000))
old=[p for f in (OUT/'batches').glob('*.json') for p in json.loads(f.read_text())['points'] if f.stem!=name]
assert len(points)==48 and len({p['id'] for p in points})==48
assert not {p['id'] for p in points}&{p['id'] for p in old};assert not {p['seed'] for p in points}&{p['seed'] for p in old}
# Snapshot outstanding commitments before fresh physical free space, under the
# same admission lock as scheduling/copy admission. Future savings are excluded.
import shutil
lock=budget.acquire()
try:
 pending=[];active=[]
 for p in old:
  if p['D']!=5:continue
  f=OUT/'jobs'/(p['id']+'.json');j=json.loads(f.read_text()) if f.exists() else {'state':'pending'}
  if j['state'] in ['pending','interrupted']:pending.append(p)
  if j['state']=='running':active.append(p)
 copy=budget.reserved_remaining();live=sum(max(0,budget.estimated_full_bytes(p)-((OUT/'chains'/(p['id']+'.h5')).stat().st_size if (OUT/'chains'/(p['id']+'.h5')).exists() else 0)) for p in active)
 prior=sum(budget.estimated_full_bytes(p) for p in pending);estimate=sum(budget.estimated_full_bytes(p) for p in points);free=shutil.disk_usage((OUT/'chains').resolve()).free
finally:lock.close()
headroom=free-copy-live-prior-estimate-25*2**30;assert headroom>0
purpose='FinalmodestD5coverage:24N12trajectoriesonq.90/.94/.96/.98 at measuredN10bracketendpoints plusroundedquarter-root midpoint, and24matchingmidpoint trajectories atN6/8/10. All4000cycles,fresh16/17, no biggerboxes. Matchingresolutiondirectlytestslocalcurvature ratherthanconfoundingnewN12midpointswithcoarselower-boxinterpolation; currenttwo-knotPCHIP=linear isnotinterpolationvalidation. Windowsaredesignguesses,noguaranteedbracketsormixing. Afterclosurecomparethreshold/cut/block/start/size/correctionfamilies onpinnedsources, reportattainedwiderbars/unresolvedasneeded ratherthanprecisionchasing. Allprior234trajectoriesaudited; frozennumericalkernels, hist/checkpoint/logprovenanceandfullauditbeforeauthorizedretention. Atomicshared-budgetcheckincludesallpriorpending/live/copyand25GiBreserve; nofuturethinningsavings. No newD2work.'
report=dict(created=run.now(),source=str(source),source_sha256=run.sha256(source),storage_estimator_sha256=run.sha256(Path(budget.__file__)),trajectories=48,physical_coordinates=24,N12_trajectories=24,lower_midpoint_trajectories=24,points=points,beta_milli_windows=windows,physical_free_GiB=free/2**30,copy_remaining_GiB=copy/2**30,live_growth_GiB=live/2**30,prior_pending_full_GiB=prior/2**30,full_transient_GiB=estimate/2**30,headroom_after_all_commitments_copy_and_25GiB_reserve_GiB=headroom/2**30,purpose=purpose,forecast='N12measuredgeneration+auditroughly16–22min per4000cycle trajectory;24N12plus24lower-boxmidpoints imply roughly8–10workerhours and1–2hcalendarplanningallowance under17outsideoffice/15daytimeworkerswhile retentionbusy. Thisisnoboundandno guaranteeofprecision; no furtheruniformproductionstageplanned ifattainedresults suffice.')
run.atomic_json(OUT/(name+'-design.json'),report)
if __name__=='__main__':
 import sys
 if '--freeze' in sys.argv:
  for j in a['source_jobs']:
   state=json.loads((OUT/'jobs'/(j['id']+'.json')).read_text());assert state['state']=='complete' and state['audit']['passed']
  legacy=run.estimated_bytes
  try:run.estimated_bytes=budget.estimated_full_bytes;run.freeze(name,points,purpose)
  finally:run.estimated_bytes=legacy
 print(json.dumps({k:report[k] for k in ['trajectories','physical_coordinates','beta_milli_windows','full_transient_GiB','headroom_after_all_commitments_copy_and_25GiB_reserve_GiB']}))
