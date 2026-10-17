"""Actual cost and closure accounting for the fixed 54-run D5 refinement."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json,statistics
from datetime import datetime,timedelta
from pathlib import Path
import run
OUT=Path(__file__).resolve().parent
batch=json.loads((OUT/'batches/d5-local-bracket-and-N10-refinement-20261014.json').read_text())
jobs=[json.loads((OUT/'jobs'/(p['id']+'.json')).read_text()) for p in batch['points']]
assert all(j['state']=='complete' and j['audit']['passed'] for j in jobs)
starts=[datetime.fromisoformat(j['updated'])-timedelta(seconds=j['attempt_generation_seconds']+j['audit_seconds']) for j in jobs]
stop=max(datetime.fromisoformat(j['updated']) for j in jobs);elapsed=(stop-min(starts)).total_seconds();worker=sum(j['attempt_generation_seconds']+j['audit_seconds'] for j in jobs)
by_size={}
for N in [6,8,10]:
 jj=[j for j in jobs if j['point']['N']==N];stats={}
 for key in ['attempt_generation_seconds','audit_seconds']:
  vals=[j[key] for j in jj];stats[key]=dict(min=min(vals),median=statistics.median(vals),max=max(vals))
 by_size[N]=dict(trajectories=len(jj),**stats,original_audited_GiB=sum(j['audit']['bytes'] for j in jj)/2**30,current_physical_GiB=sum((OUT/'chains'/(j['point']['id']+'.h5')).stat().st_size for j in jj)/2**30)
r=dict(created=run.now(),batch=batch['name'],trajectories=len(jobs),physical_coordinates=len(jobs)//2,all_full_field_audits_passed=True,numerical_sources_match_frozen=batch['source_hashes']==run.source_hashes(),generation_plus_audit_worker_hours=worker/3600,elapsed_minutes=elapsed/60,effective_workers=worker/elapsed,by_size=by_size,original_audited_GiB=sum(j['audit']['bytes'] for j in jobs)/2**30,current_physical_GiB=sum((OUT/'chains'/(j['point']['id']+'.h5')).stat().st_size for j in jobs)/2**30,status_snapshot=json.loads((OUT/'status.json').read_text()),memory_sample='logs/d5-N10-refinement-memory-sample-20261014.log; RSS 468848–792064KiB observed, not peak memory; near-one-core CPU usage per sampled worker.',timing_scope='Observed wall times include real concurrent simulation/audit/SSD conditions and brief resource pausing; no macOS volume-permission stall in this new cohort. Effective workers is aggregate occupied walltime/elapsed, not CPU utilization. N12 cost remains unmeasured.',interpretation='q92 N10 window missed the Binder-quarter crossing; dedicated lower-bracket stage is already frozen. q1 coarse response needs a midpoint test; dedicated modest N12 q1/endpoint cost/coverage stage must pass source and storage guards. No new D2 precision work. Earlier12–29h D5 planning allowance remains a scenario, not a bound or precision guarantee.')
run.atomic_json(OUT/'d5-N10-refinement-actual-cost-20261014.json',r)
print(json.dumps({k:r[k] for k in ['trajectories','generation_plus_audit_worker_hours','elapsed_minutes','effective_workers','original_audited_GiB','current_physical_GiB','by_size']},indent=1))
