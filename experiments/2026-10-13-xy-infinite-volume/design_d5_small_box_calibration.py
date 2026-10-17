"""User-prioritized D5 calibration with periodic-box, above-critical-dimension analysis."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json
from pathlib import Path
import numpy as np
import run
OUT=Path(__file__).resolve().parent
points=[]
for q,betas in [(92,[0,60,120]),(100,[180,230,280])]:
 for N in [4,6,8]:
  for beta in betas:
   for r in [0,1]:points.append(run.point(5,q,beta,N,r,4000))
for q in [86,88,90]:
 for N in [4,6,8]:
  for r in [2,3]:
   p=run.point(5,q,0,N,r,4000)
   p.update(id=f'zero_d5_e{q*10:04d}_N{N:03d}_r{r:02d}_s004000',seed=int(np.random.SeedSequence([2026101403,5,q*10,N,r]).generate_state(1,dtype=np.uint64)[0]),campaign_role='zero_beta_endpoint',epsilon_sampling_milli=q*10)
   points.append(p)
old=[p for f in (OUT/'batches').glob('*.json') for p in json.loads(f.read_text())['points']]
assert not {p['id'] for p in points}&{p['id'] for p in old}
assert not {p['seed'] for p in points}&{p['seed'] for p in old}
reason='UserOctober14: stopnewD2above512, extractattained512results, thenmoveD5. SmallD5 cost/coveragecalibration,N4/6/8,4000cycles,twoindependentstarts. FixedBinderU=.25locallybracketed; separateD5FSS usingperiodiczero-mode thermalwindow N^-5/2 withtestedsizecuts andsubleading N^-3 sensitivity. Do notreuseD4logscaling. Pilot isnotfinalprecisionormixingcertification; measureactualgeneration+auditcost,physicalstorage,bracketcoverage andstart/blockdiagnostics before full. Fullfields generatedandaudited beforeauthorizedretention; everyhistory,log,checkpointandprovenance retained. Beta0roleexcludedfromfixedepsilonpositivegrouping. Nmax8initial; noD2precisionchasing.'
name='d5-small-box-cost-and-coverage-calibration-20261014'
s=json.loads((OUT/'status.json').read_text());estimate=sum(run.estimated_bytes(p) for p in points)
headroom=s['physical_disk_free_bytes']-s['active_reserved_bytes']-s['retention_reserved_bytes']-estimate-25*2**30
assert headroom>0
r=dict(created=run.now(),points=points,trajectories=len(points),transient_full_GiB=estimate/2**30,headroom_GiB=headroom/2**30,purpose=reason,forecast='NoD5measuredcost yet; these modest boxes calibrate actualcost and coverage before largerboxes/fullgrid. KeepFridaydeadline and analysis margin.',scaling_references=['https://arxiv.org/abs/cond-mat/9809394','https://arxiv.org/abs/1410.5296'],scaling_scope='Inference from periodic phi4/O(n) zero-mode FSS above d_c4, tested on thisconstrainedXYdataset. Corrections/thresholddriftmustbechecked; notassumed exactatN4/6/8.',zero_mode_O2_U_at_t0=1-np.pi/4)
run.atomic_json(OUT/(name+'-design.json'),r)
if __name__=='__main__':
 import sys
 if '--freeze' in sys.argv:run.freeze(name,points,reason)
 print(json.dumps({k:r[k] for k in ['trajectories','transient_full_GiB','headroom_GiB']}))
