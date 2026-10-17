"""Measured-cost staged D5 coverage, with exact beta-zero reuse and endpoint repair."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
from pathlib import Path
import json,numpy as np
import run,storage_admission as budget
OUT=Path(__file__).resolve().parent
name='d5-remaining-low-boxes-and-N12-bracket-coverage-20261014'
windows={90:[13,26],94:[101,113,125],96:[144,156,168],98:[183,196,209]}
points=[run.point(5,q,b,N,r,4000) for q,betas in windows.items() for N in [6,8,10] for b in betas for r in [10,11]]
points += [run.point(5,92,b,12,r,4000) for b in [60,64,68] for r in [10,11]]
for e in [890,892,893]:
 for r in [10,11]:
  p=run.point(5,e//10,0,12,r,4000);p.update(id=f'zero_d5_e{e:04d}_N012_r{r:02d}_s004000',epsilon_over_pi=e/1000,epsilon=e*np.pi/1000,seed=int(np.random.SeedSequence([2026101403,5,e,12,r]).generate_state(1,dtype=np.uint64)[0]),campaign_role='zero_beta_endpoint',epsilon_sampling_milli=e);points.append(p)
old=[p for f in (OUT/'batches').glob('*.json') for p in json.loads(f.read_text())['points'] if f.stem!=name]
assert len(points)==78 and len({p['id'] for p in points})==78
assert not {p['id'] for p in points}&{p['id'] for p in old}
assert not {p['seed'] for p in points}&{p['seed'] for p in old}
reuse=[]
for p in old:
 if p['D']==5 and p['N'] in [6,8,10] and p['beta']==0 and p.get('campaign_role')=='zero_beta_endpoint' and p.get('epsilon_sampling_milli')==900 and p['epsilon_over_pi']==.90:
  j=json.loads((OUT/'jobs'/(p['id']+'.json')).read_text());assert j['state']=='complete' and j['audit']['passed'];reuse.append(p)
assert {p['N'] for p in reuse}=={6,8,10}
source=OUT/json.loads((OUT/'latest-d5-calibration.json').read_text())['folder']/'analysis.json'
a=json.loads(source.read_text())
assert any(r['axis']=='epsilon' and r['N']==12 and r['status']=='unresolved' for r in a['roots'])
assert any(r['axis']=='beta' and r['q']==.92 and r['N']==10 and r['status']=='conditional_crossing' for r in a['roots'])
s=json.loads((OUT/'status.json').read_text());pending=[]
for p in old:
 if p['D']!=5:continue
 f=OUT/'jobs'/(p['id']+'.json');j=json.loads(f.read_text()) if f.exists() else {'state':'pending'}
 if j['state'] in ['pending','interrupted']:pending.append(p)
estimate=sum(budget.estimated_full_bytes(p) for p in points);pending_bytes=sum(budget.estimated_full_bytes(p) for p in pending)
headroom=s['physical_disk_free_bytes']-s['active_reserved_bytes']-s['retention_reserved_bytes']-pending_bytes-estimate-25*2**30
assert headroom>0
purpose='After14N12/midpoint probes passed fullaudits, N12 generation834–1134s/audit119–146s supports practicalstaging. Measurefourremainingq.90/.94/.96/.98 exact.02grid atN6/8/10,4000cycles,fresh10/11. Windows interpolate measured finiteNroots only (design guesses, no rejectedinfinitefituse). q90beta0 reuses exactepsilon_sampling_milli900 endpoint histories at declaredq.90, no rounding or rolemutation; all other endpointcoords stayexcludedfixedeps. 66newpositive trajectories=33coords instead72 byreusing3physicalbeta0coords,10auditedstarttrajectories. Add6N12q92 trajectories beta.060/.064/.068 aftermeasuredN10root.065066. Add6N12beta0lower-epsilonanchors.890/.892/.893 because .894/.896/.898 allbelowUquarter; unresolvedlargestendpoint staysflaggeduntilcoverage. Fullphysicalchunkbudgets,nofuturethinningsavings,25GiBreserve. Frozenkernels; allhistories/checkpoints/logs retained,fullaudit beforeauthorizedfieldretention. No newD2work.'
report=dict(created=run.now(),source_analysis=str(source),source_sha256=run.sha256(source),storage_estimator_sha256=run.sha256(Path(budget.__file__)),trajectories=len(points),physical_coordinates=len(points)//2,positive_new_trajectories=66,positive_new_coordinates=33,N12_q092_trajectories=6,N12_endpoint_trajectories=6,reused_exact_beta_zero_points=reuse,reused_trajectories=len(reuse),reused_physical_coordinates=3,points=points,beta_milli_windows=windows,full_transient_GiB=estimate/2**30,prior_pending_full_GiB=pending_bytes/2**30,headroom_after_prior_pending_live_copy_and_25GiB_reserve_GiB=headroom/2**30,status_snapshot=s,purpose=purpose,forecast='ActualN10/N12cost supports roughly1–2h planning allowance includingqueuedwavesandaudits, notaruntimebound. N12onremainingfourpositiveqnotyetcommitted; assess coverage andactualstorage first. Lowblock/start/half flags persist; no certifiedtotalprecision.')
run.atomic_json(OUT/(name+'-design.json'),report)
if __name__=='__main__':
 import sys
 if '--freeze' in sys.argv:
  legacy=run.estimated_bytes
  try:
   run.estimated_bytes=budget.estimated_full_bytes;run.freeze(name,points,purpose)
  finally:run.estimated_bytes=legacy
 print(json.dumps({k:report[k] for k in ['trajectories','physical_coordinates','reused_trajectories','full_transient_GiB','headroom_after_prior_pending_live_copy_and_25GiB_reserve_GiB']}))
