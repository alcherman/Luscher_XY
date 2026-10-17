"""Cost probes after the measured N512 holdout checks; no uniform N1024 grid."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json
from pathlib import Path
import numpy as np
import run
OUT=Path(__file__).resolve().parent
old=[p for f in (OUT/'batches').glob('*.json') for p in json.loads(f.read_text())['points']]
points=[run.point(2,80,921,1024,r,12000) for r in [8,9]]
for e in [553,558,565]:
    for r in [8,9]:
        p=run.point(2,e//10,0,512,r,12000)
        p.update(id=f'zero_d2_e{e:04d}_N512_r{r:02d}_s012000',epsilon_over_pi=e/1000,epsilon=e*np.pi/1000,seed=int(np.random.SeedSequence([2026101403,2,e,512,r]).generate_state(1,dtype=np.uint64)[0]),campaign_role='zero_beta_endpoint',epsilon_sampling_milli=e,primary_observables=['eta_eff','xi_over_L'])
        points.append(p)
assert not {p['id'] for p in points}&{p['id'] for p in old}
assert not {p['seed'] for p in points}&{p['seed'] for p in old}
purpose='Targeted D2 size probes after completed N512 brackets. Smaller-box q.80 C0 predicts R512/eta256 roots lowby .00885/.00863, materialfractionof .02 target; higheroffsetmodels remain plausible, so this isnotproofN1024required. TwoN1024 cost/observable probes atq.80beta.921,matchingN512,12000cycles,fresh8/9. Measure actualruntime/memory/storage/mixingbeforeanymoreN1024; onebeta cannotbracket a crossing. Also sixN512 beta0endpoint trajectories epsilon/pi .553/.558/.565,matchingcompletedN256, fresh8/9,12000cycles, tocheckwidelogarithmicendpoint extrapolation withnewsize leverage. Userauthorizedadaptivefieldretention afterfullgeneration/fullfieldaudit applies. Allhistories/checkpoints/provenance/logs remain complete; shared25GiB diskadmission anddaytimecorereserve. Numericalkernelsfrozen, noD5.'
name='d2-eta-correlation-targeted-N1024-cost-and-N512-endpoint-probes-20261014'
source=OUT/json.loads((OUT/'latest-d2-heldout-validation.json').read_text())['folder']/'analysis.json'
design=dict(created=run.now(),source=str(source),source_sha256=run.sha256(source),points=points,trajectories=len(points),transient_full_GiB=sum(run.estimated_bytes(p) for p in points)/2**30,purpose=purpose,forecast='N1024 maycost4–8times measuredN512generation plusaudit: roughly1.5–3.5h allowance for twoconcurrent probes, not measured bound. N512 endpoints about20–30min underpilotthroughput,withretention/I/O contention allowance. NootherN1024commitment untilactualcost/observables audited.')
run.atomic_json(OUT/(name+'-design.json'),design)
if __name__=='__main__':
    import sys
    if '--freeze' in sys.argv:run.freeze(name,points,purpose)
    print(json.dumps({k:design[k] for k in ['trajectories','transient_full_GiB','forecast']}),flush=True)
