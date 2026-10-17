"""Adaptive simulation design from existing Binder-quarter measurements."""
import os
for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS']:
    os.environ[name]='1'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json
from pathlib import Path
import h5py
import numpy as np
import run
from threshold_calibration import crossing,coordinate

OUT=Path(__file__).resolve().parent


def main():
    grouped={};sources=[]
    for path in sorted((OUT/'jobs').glob('*.json')):
        job=json.loads(path.read_text());p=job['point']
        if job['state']!='complete' or not job['audit']['passed']:continue
        if p['D'] not in [3,4] or p['N'] not in [6,8] or p['epsilon_over_pi']>=1:continue
        with h5py.File(OUT/'chains'/f"{p['id']}.h5") as f:
            a=f['magnetization_squared'][1001:];b=f['magnetization_fourth'][1001:]
            moments=[float(a.mean()),float(b.mean())]
        grouped.setdefault((p['D'],round(100*p['epsilon_over_pi']),p['N'],p['beta']),[]).append(moments)
        sources.append(dict(id=p['id'],audit=job['audit']))
    choices=[];points=[]
    for D,q in sorted({key[:2] for key in grouped}):
        roots=[]
        for N in [6,8]:
            vals=[]
            for (d,e,n,beta),moments in sorted(grouped.items()):
                if (d,e,n)!=(D,q,N) or len(moments)!=2:continue
                m2,m4=np.mean(moments,axis=0);vals.append((beta,1-m4/(2*m2*m2)))
            beta,U=np.array(vals).T
            roots.append(crossing(beta,U))
        if not np.isfinite(roots).all():
            choices.append(dict(D=D,epsilon_over_pi=q/100,status='needs_small_box_bracket',finite_roots=[float(r) if np.isfinite(r) else None for r in roots]));continue
        x=coordinate(D,[6,8]);coeff=np.polyfit(x,roots,1)
        predictions=[]
        for N in [12,16]:
            center=float(np.polyval(coeff,coordinate(D,N)))
            # The extrapolation only designs a bracket; it is not a scientific result.
            betas=sorted({max(0,round(1000*(center-.012))),max(1,round(1000*(center+.012)))})
            predictions.append(dict(N=N,design_center=center,beta_milli=betas))
            points += [run.point(D,q,b,N,r,4000) for b in betas for r in [0,1]]
        choices.append(dict(D=D,epsilon_over_pi=q/100,status='designed_larger_box_brackets',N6_N8_roots=roots,predictions=predictions))
    if len({p['id'] for p in points})!=len(points):raise RuntimeError('Duplicate design points')
    budget=sum(run.estimated_bytes(p) for p in points)
    manifest=run.freeze('binder-quarter-constrained-larger-boxes',points,
        'User-directed Binder-quarter primary estimator. Predict N12,16 threshold neighborhoods from actual N6,8 crossings using leading dimension-specific shift only for simulation design. Two beta anchors +/- .012, two independent starts,4000cycles. Adapt brackets when missed; no thermodynamic estimate is accepted from the two small-box design fit. Preserve full fields, RNG, histories, audits and logs.')
    # Three widely separated 2D slices measure actual cost and threshold coverage.
    # The old near-critical anchors alone do not bracket U=.25.
    pilot=[]
    for q,upper in [(60,300),(80,920),(100,1120)]:
        betas=sorted({0,round(.5*upper),round(.9*upper)})
        pilot += [run.point(2,q,b,N,r,4000) for b in betas for N in [16,32,64] for r in [0,1]]
    pilot_manifest=run.freeze('binder-quarter-d2-cost-and-coverage-pilot',pilot,
        'New d2 threshold pilot at epsilon/pi=.60,.80,1 and N16,32,64: beta0 and half/nine-tenths of historical near-critical design anchors .30,.92,1.12, two starts,4000cycles. Locate disordered-side Binder-quarter brackets and measure actual generation-plus-audit/runtime/storage. No assumption that beta0 is disordered; if no positive bracket, retain the evidence and check size or alternate fixed-threshold sensitivity. No new environment or retention changes.')
    run.atomic_json(OUT/'binder-quarter-adaptive-design.json',dict(created=run.now(),choices=choices,source_jobs=sources,
        design_source_sha256=run.sha256(Path(__file__)),larger_box_trajectories=len(points),larger_box_estimated_GiB=budget/2**30,
        larger_box_cost_scope='120 baseline new trajectories expected if15slice designs bracket. Previous forecast estimates2.42 generation-plus-audit worker-hours from measured endpoint costs; budget2-6h with constrained-slice/refinement overhead. This is a planning allowance, not measured precision.',
        d2_pilot_trajectories=len(pilot),d2_pilot_estimated_GiB=pilot_manifest['estimated_bytes']/2**30,
        d2_cost_scope='Not yet measured in current campaign. Small-box pilot before full d2 sampling or N128 escalation.',
        storage_reserve_GiB=25,interpretation='Finite-box threshold design only; not accepted critical couplings.'))
    print('constrained larger-box trajectories',len(points),'estimated GiB',budget/2**30,flush=True)
    print('d2 pilot trajectories',len(pilot),'estimated GiB',pilot_manifest['estimated_bytes']/2**30,flush=True)
    print('choices',json.dumps(choices),flush=True)


if __name__=='__main__':main()
