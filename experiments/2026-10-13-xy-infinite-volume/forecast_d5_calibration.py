"""Measured pilot throughput and explicit, uncommitted D5 planning scenarios."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json
from pathlib import Path
from datetime import datetime,timedelta
import numpy as np
import run
OUT=Path(__file__).resolve().parent
jobs=[json.loads(f.read_text()) for f in (OUT/'jobs').glob('*d5*.json')]
jobs=[j for j in jobs if j['state']=='complete' and j['point']['steps']==4000]
assert len(jobs)==54 and all(j['audit']['passed'] for j in jobs)
closed=[datetime.fromisoformat(j['updated']) for j in jobs]
start=min(t-timedelta(seconds=j['attempt_generation_seconds']+j['audit_seconds']) for j,t in zip(jobs,closed))
elapsed=(max(closed)-start).total_seconds()
worker=sum(j['attempt_generation_seconds']+j['audit_seconds'] for j in jobs)
rate=worker/elapsed
cost={}
for N in [4,6,8]:
    rows=[j for j in jobs if j['point']['N']==N]
    cost[N]=dict(count=len(rows),generation_median=float(np.median([j['attempt_generation_seconds'] for j in rows])),audit_median=float(np.median([j['audit_seconds'] for j in rows])))
base=cost[8]['generation_median']+cost[8]['audit_median']
remaining_seconds=sum((42 if N==12 else 24)*base*(N/8)**5 for N in [6,8,10,12])
refinement_seconds=sum(18*2*base*(N/8)**5 for N in [6,8,10])
baseline_bytes=sum(42*run.estimated_bytes(run.point(5,100,225,N,0,4000)) for N in [6,8,10,12])
remaining_bytes=sum((42 if N==12 else 24)*run.estimated_bytes(run.point(5,100,225,N,0,4000)) for N in [6,8,10,12])
double_bytes=sum((42 if N==12 else 24)*run.estimated_bytes(run.point(5,100,225,N,0,8000)) for N in [6,8,10,12])
physical=sum((OUT/'chains'/f"{j['point']['id']}.h5").stat().st_size for j in jobs)
report=dict(created=run.now(),pilot_closed=54,audits_passed=54,per_size=cost,worker_hours=worker/3600,elapsed_minutes=elapsed/60,achieved_effective_workers=rate,original_audit_GiB=sum(j['audit']['bytes'] for j in jobs)/2**30,current_physical_GiB=physical/2**30)
report['baseline_plan']=dict(positive_slices=[.90,.92,.94,.96,.98,1],beta_samples_per_slice_size=3,sizes=[6,8,10,12],starts=2,cycles=4000,endpoint_epsilon_samples=3,total_trajectories=168,physical_coordinates=84,conservative_full_GiB=baseline_bytes/2**30)
report['refinement_plan']=dict(trajectories=54,sizes=[6,8,10],cycles=8000,status='frozen_running',conservative_full_GiB=271.1326206754893)
report.update(reused_refinement_trajectories=54,reused_physical_coordinates=27,remaining_new_trajectories=114,remaining_physical_coordinates=57,remaining_conservative_full_GiB=remaining_bytes/2**30,uniform_8000_remaining_GiB=double_bytes/2**30,projected_refinement_worker_hours=refinement_seconds/3600,projected_remaining_worker_hours=remaining_seconds/3600,ideal_simulation_wall_hours=(remaining_seconds+refinement_seconds)/rate/3600,planning_simulation_hours=[2,5],extra_statistics_refinement_hours=[4,12],analysis_validation_delivery_hours=[6,12],planning_total_hours=[12,29],hours_to_deadline=(run.DEADLINE-run.datetime.now(run.CENTRAL)).total_seconds()/3600)
report['scope']='Planning range, not a bound or precision guarantee. Volume scaling above N8 is unmeasured. Stage full grid only after actual N10 throughput, mixing/coverage and fresh reservation-aware free space checks. Baseline remaining space exceeds recorded post-refinement headroom by roughly8–10GiB; do not count hypothetical future thinning. Uniform8k remaining grid does not fit uncompacted budget. Wider bars/unresolved locations permitted within Friday budget.'
run.atomic_json(OUT/'runtime-forecast-d5-calibration-20261014.json',report)
text=f'''# Five-dimensional campaign forecast after audited calibration

All **54 of 54** initial trajectories passed full closure audits. They used N=4,6,8, 4,000 cycles and two independent starts per coordinate. Generation plus audits totalled **{worker/3600:.3f} worker-hours** over **{elapsed/60:.2f} minutes**, equivalent to **{rate:.2f} effective workers** under the 15-worker daytime ceiling with retention active. Original full files occupied **{report['original_audit_GiB']:.2f} GiB**; current physical files occupy **{physical/2**30:.2f} GiB**.

Median generation/audit seconds: N=4 **{cost[4]['generation_median']:.2f}/{cost[4]['audit_median']:.2f}**, N=6 **{cost[6]['generation_median']:.2f}/{cost[6]['audit_median']:.2f}**, N=8 **{cost[8]['generation_median']:.2f}/{cost[8]['audit_median']:.2f}**. Nine coarse linear roots retain all bootstrap draws. Seven coordinate groups have fewer than 50 blocks per start. Coarse beta interpolation and these short runs do not establish final precision or mixing.

The baseline scenario has six positive-beta slices at epsilon/pi=.90,.92,.94,.96,.98,1.00, three local beta samples per size, N=6,8,10,12, two starts and 4,000 cycles. Three epsilon samples at beta=0 on the same four sizes bring the total to **168 trajectories**. The running **54-trajectory (27-coordinate), 8,000-cycle refinement** at N=6,8,10 covers two slices and the endpoint; reusing it leaves **114 new baseline trajectories (57 coordinates)**. Its narrower windows address coarse interpolation, and N=10 adds a modest size check. Reuse credits data, not imagined free space.

Scaling measured N=8 cost by volume predicts **{refinement_seconds/3600:.1f} worker-hours** for the refinement and **{remaining_seconds/3600:.1f}** for the remaining baseline. Ideal simulation wall time at achieved throughput is **{report['ideal_simulation_wall_hours']:.1f} hours**. Planning allows **2–5 hours** for generation/audits, **4–12 hours** for extra statistics/refinement, and **6–12 hours** for analysis, independent review and delivery: **12–29 hours total**, with **{report['hours_to_deadline']:.1f} hours** remaining to Friday October16 at15:00Central. Practical attained estimates appear feasible; this is a planning scenario, not a bound or a guarantee of ±0.02 at every point. Actual N=10/N=12 cost and mixing may change scope. Weekday simulation ceilings are15 while retention is busy and16 otherwise; outside office hours they are17/18.

**Storage remains a staging constraint.** The running refinement reserves **271.13 GiB** conservatively. The baseline full-file estimate is **826.29 GiB**; after coordinate reuse, remaining runs require **690.28 GiB**. This is roughly8–10GiB beyond recorded post-refinement headroom before further actual recovery. Do not commit the full grid yet: check actual N=10 cost, retained bytes, and fresh physical free space minus active/copy reservations and the25GiB reserve. Uniform8,000-cycle expansion of the remaining grid would require **{double_bytes/2**30:.2f} GiB** and does not fit current uncompacted capacity. No future thinning savings are counted.

Anticipated computational constraints are larger-box generation, independent full-field audits and SSD I/O. No productive N=10/N=12 timing or peak memory has yet been measured. Periodic zero-mode scaling motivates an N^-5/2 leading shift and an N^-3 stress correction; three-size correction fits have zero degrees of freedom and remain diagnostic. This is a conditional application to constrained XY data, not the4D logarithmic formula or a claim of small-box universality. See [Chen–Dohm](https://arxiv.org/abs/cond-mat/9809394) and [Wittmann–Young](https://arxiv.org/abs/1410.5296).

Numerical kernels, complete observable histories, checkpoints, provenance and logs remain preserved. Authorized field-only retention occurs after full audits. Detailed calculations and source metadata are in runtime-forecast-d5-calibration-20261014.json and latest-d5-calibration.json.
'''
for old,new in {'October16 at15:00Central':'October 16 at 15:00 Central','are15':'are 15','and16':'and 16','are17/18':'are 17/18','roughly8–10GiB':'roughly 8–10 GiB','the25GiB':'the 25 GiB','Uniform8,000':'Uniform 8,000','the4D':'the 4D'}.items():text=text.replace(old,new)
(OUT/'runtime-forecast-d5-calibration-20261014.md').write_text(text)
print('forecast',report['planning_total_hours'],'available',report['hours_to_deadline'])
