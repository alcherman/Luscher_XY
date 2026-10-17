"""Correct the measured D5 q=.92 window without treating coarse roots as final."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json
from pathlib import Path
import run
import storage_admission
OUT=Path(__file__).resolve().parent
name='d5-q092-lower-bracket-refinement-20261014'
points=[run.point(5,92,b,N,r,8000) for N in [6,8,10] for b in [60,64,68] for r in [6,7]]
source=OUT/'d5-calibration-analysis-20261014-145843-cut1001-block1-U25-local-all'/'analysis.json'
a=json.loads(source.read_text())
assert any(r['axis']=='beta' and r['q']==.92 and r['N']==10 and r['status']=='unresolved' for r in a['roots'])
old=[p for f in (OUT/'batches').glob('*.json') for p in json.loads(f.read_text())['points'] if f.stem!=name]
assert not {p['id'] for p in points}&{p['id'] for p in old}
assert not {p['seed'] for p in points}&{p['seed'] for p in old}
for b in [70,90,110]:
 for r in [4,5]:
  j=json.loads((OUT/'jobs'/f'd5_e092_b{b:04d}_N010_r{r:02d}_s008000.json').read_text())
  assert j['state']=='complete' and j['audit']['passed']
s=json.loads((OUT/'status.json').read_text());estimate=sum(storage_admission.estimated_full_bytes(p) for p in points)
headroom=s['physical_disk_free_bytes']-s['active_reserved_bytes']-s['retention_reserved_bytes']-estimate-25*2**30
assert headroom>0
purpose='Measured q.92 N10 U(.070)=.4933: all .070/.090/.110 are above quarter so no bracket. N6/8 locally cross between .060/.070; coarsepilot .06:.12 roots were not precise and misleading. Sample .060/.064/.068 at N6/8/10,fresh6/7,8000cycles. Coverage pilot, not guaranteed bracket or mixing repair. Preserve N8.070 start-disagreement diagnostics; fresh longer .060 sampling complements old4000cycle starts. No new largerD2 orN12q92. All relevant N10 q92 source runs audited before freeze; remaining independent q1/endpoint refinement may still run. Kernels frozen,fullaudit before authorizedretention,histories/checkpoints/logs retained.'
report=dict(created=run.now(),source_analysis=str(source),source_sha256=run.sha256(source),trajectories=len(points),physical_coordinates=len(points)//2,points=points,transient_full_GiB=estimate/2**30,headroom_after_active_copy_and_25GiB_reserve_GiB=headroom/2**30,status_snapshot=s,purpose=purpose,scope='Adaptive bracket correction within existing refinement allowance; no final precision or Friday timing guarantee.')
run.atomic_json(OUT/(name+'-design.json'),report)
if __name__=='__main__':
 import sys
 if '--freeze' in sys.argv:run.freeze(name,points,purpose)
 print(json.dumps({k:report[k] for k in ['trajectories','physical_coordinates','transient_full_GiB','headroom_after_active_copy_and_25GiB_reserve_GiB']}))
