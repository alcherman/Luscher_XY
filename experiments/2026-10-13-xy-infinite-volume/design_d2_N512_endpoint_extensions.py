"""Extend weak epsilon anchors for audited N512 endpoint threshold coverage."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json
from pathlib import Path
import numpy as np
import run
OUT=Path(__file__).resolve().parent
old=[p for f in (OUT/'batches').glob('*.json') for p in json.loads(f.read_text())['points']]
points=[]
for e in [546,572]:
    for r in [10,11]:
        p=run.point(2,e//10,0,512,r,12000)
        p.update(id=f'zero_d2_e{e:04d}_N512_r{r:02d}_s012000',epsilon_over_pi=e/1000,epsilon=e*np.pi/1000,seed=int(np.random.SeedSequence([2026101403,2,e,512,r]).generate_state(1,dtype=np.uint64)[0]),campaign_role='zero_beta_endpoint',epsilon_sampling_milli=e,primary_observables=['eta_eff','xi_over_L'])
        points.append(p)
assert not {p['id'] for p in points}&{p['id'] for p in old}
assert not {p['seed'] for p in points}&{p['seed'] for p in old}
source=OUT/json.loads((OUT/'latest-zero-beta-endpoints.json').read_text())['folder']/'analysis.json'
reason='Audited N512epsilon .553/.558/.565 gives central R512epsilon .557984 andeta256 .557141 but only .9825/.9815 valid crossing draws after extra-neighbor local-window correction. Extendbothsides at.546/.572,matchingexistingN256, tocheckcoveragebeforeinfinitevolumeacceptance andcut/blockchecks. Allhistories/exactcheckpoints/logsretained; fullfieldsproducedandauditedbeforeauthorizedretention. NoN1024expansion orD5.'
name='d2-eta-correlation-N512-endpoint-bracket-extensions-20261014'
s=json.loads((OUT/'status.json').read_text());estimate=sum(run.estimated_bytes(p) for p in points)
headroom=s['physical_disk_free_bytes']-s['active_reserved_bytes']-s['retention_reserved_bytes']-estimate-25*2**30
assert headroom>0
r=dict(created=run.now(),source=str(source),source_sha256=run.sha256(source),points=points,trajectories=len(points),purpose=reason,transient_full_GiB=estimate/2**30,uncommitted_headroom_GiB=headroom/2**30,forecast='Actual firstsix N512endpointgeneration1073–1083s plusaudit145–147s, so20–30min allowance withI/Ocontention; notguaranteedcrossingprecision.',status_at_design=s)
run.atomic_json(OUT/(name+'-design.json'),r)
if __name__=='__main__':
 import sys
 if '--freeze' in sys.argv:run.freeze(name,points,reason)
 print(json.dumps({k:r[k] for k in ['trajectories','transient_full_GiB','uncommitted_headroom_GiB']}))
