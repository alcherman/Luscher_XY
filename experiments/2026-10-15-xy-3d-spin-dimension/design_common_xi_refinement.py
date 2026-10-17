"""Measured common-coupling crossing refinement, no increase beyond N64."""
import json
import sys
import numpy as np
import measure as m
sys.path.insert(0,str(m.CAMPAIGN))
import storage_admission as storage

NAME='d3-zero-beta-common-xi-crossing-refinement-20261015'
SMALL_GRID=[800250,800500,800625,800750,800875,801000,801125,801250,801500]
LARGE_GRID={48:[800250,800500,800625,800875,801125],64:[800625,800875,801125]}


def make_points(existing):
    points=[];reused=[]
    for N in [6,8,12,16,24,32,48,64]:
        steps=12000 if N<=32 else 20000
        for code in SMALL_GRID if N<=32 else LARGE_GRID[N]:
            q=code/1_000_000
            eligible=[p for p in existing if p['D']==3 and p['beta']==0 and p['N']==N and p['steps']>=steps and abs(float(p['epsilon'])/np.pi-q)<1e-14]
            if len(eligible)>=2 and len({p['initial_winding'][0] for p in eligible})>=2:
                reused.append(dict(N=N,epsilon_over_pi=q,ids=sorted(p['id'] for p in eligible)));continue
            for replica in [40,41]:
                p=m.run.point(3,80,0,N,replica,steps)
                p.update(id=f'spinxi3_e{code:06d}_N{N:03d}_r{replica:02d}_s{steps:06d}',epsilon_over_pi=q,epsilon=q*np.pi,
                         seed=int(np.random.SeedSequence([2026101505,3,code,N,replica]).generate_state(1,dtype=np.uint64)[0]),
                         initial_winding=[replica%2,0,0],campaign_role='zero_beta_endpoint',epsilon_sampling_1e6=code,
                         study_role='spin_dimension_common_coupling_xi_refinement',primary_observables=['magnetization_squared','structure_zero','structure_min','xi_over_L'])
                points.append(p)
    # Start the long measured N64/N48 trajectories first to avoid a long idle tail.
    points.sort(key=lambda p:(-p['N'],p['epsilon_sampling_1e6'],p['replica']))
    return points,reused


def main(freeze=False):
    old={p['id']:p for f in (m.CAMPAIGN/'batches').glob('*.json') for p in json.loads(f.read_text())['points']}
    completed=[];jobs={}
    for i,p in old.items():
        f=m.CAMPAIGN/'jobs'/f'{i}.json';j=json.loads(f.read_text()) if f.exists() else dict(state='pending');jobs[i]=j
        if j['state']=='complete' and j.get('audit',{}).get('passed') and not i.startswith('spinxi3_'):completed.append(p)
    ps,reused=make_points(completed);ids={p['id'] for p in ps}
    assert len(ids)==len(ps) and len({p['seed'] for p in ps})==len(ps)
    assert all(p['id'] not in old or old[p['id']]==p for p in ps)
    assert not {p['seed'] for p in ps}&{p['seed'] for i,p in old.items() if i not in ids}
    assert all(p['D']==3 and p['beta']==0 and p['N']<=64 for p in ps)
    prior={}
    for i,p in old.items():
        if i not in ids and jobs[i]['state'] in ['pending','interrupted','running']:
            path=m.CAMPAIGN/'chains'/f'{i}.h5';prior[i]=max(0,storage.estimated_full_bytes(p)-(path.stat().st_size if path.exists() else 0))
    # Forecast per-size measured median generation+audit/row from preceding two scans.
    timing={}
    for N in sorted({p['N'] for p in ps}):
        records=[(jobs[p['id']].get('attempt_generation_seconds',0)+jobs[p['id']].get('audit_seconds',0))/p['steps'] for p in completed if p['D']==3 and p['beta']==0 and p['N']==N and (p['id'].startswith('spinzero3_') or p['id'].startswith('spinpilot3_')) and jobs[p['id']].get('attempt_generation_seconds')]
        if not records:records=[(jobs[p['id']].get('attempt_generation_seconds',0)+jobs[p['id']].get('audit_seconds',0))/p['steps'] for p in completed if p['D']==3 and p['beta']==0 and p['N']==N and jobs[p['id']].get('attempt_generation_seconds')]
        timing[N]=float(np.median(records)) if records else None
    seconds=sum(p['steps']*timing[p['N']] for p in ps if timing[p['N']] is not None)
    source=m.OUT/'common-xi-delivery-20261015-195403'/'quotient-ledger.json'
    points_hash=m.hashlib.sha256(json.dumps(ps,sort_keys=True).encode()).hexdigest()
    with storage.acquire():
        copy=storage.reserved_remaining();physical=m.run.shutil.disk_usage(m.CAMPAIGN/'chains').free
        ceiling=sum(storage.estimated_full_bytes(p) for p in ps);headroom=physical-copy-sum(prior.values())-ceiling-25*2**30
        if headroom<0:raise RuntimeError('Insufficient full transient storage; no future thinning counted')
        body=dict(name=NAME,created=m.run.now(),points=ps,points_sha256=points_hash,trajectories=len(ps),new_coordinates=len(ps)//2,reused_coordinates=reused,
                  source_ledger=str(source),source_sha256=m.run.sha256(source),full_transient_bytes=ceiling,physical_free_bytes=physical,prior_pending_live_growth_bytes=sum(prior.values()),retention_copy_growth_bytes=copy,headroom_after_25GiB_bytes=headroom,
                  measured_seconds_per_cycle_by_N=timing,baseline_worker_hours=seconds/3600,ideal_14_worker_hours=seconds/3600/14,calendar_allowance_hours=[1,3],forecast_scope='Same D3 beta0 measured median generation+audit/cycle. N6 uses earlier4k runs; startup and resource conditions differ. Baseline is not a bound; allowance includes scheduling, audits and scalar analysis, without guaranteed crossing/exponent precision.',
                  purpose='Direct user authorization: improve coarse epsilon grids and resolve common-coupling xi/L crossing24:48. Fine common knots around measured .8007-.8011 roots at N6/8/12/16/24/32; extend N48 below previous lower edge and add matching interior N48/N64 knots. Reuse clean measured coordinates. Two fresh zero/unit-x starts. No Delta,critical xi/L or nu imposed, no boxes above64 and no automatic high-statistics campaign.',
                  analysis_plan='After full audits: common-coupling xi/L roots and joint m2 quotients, exact root multiplicity checks, linear/PCHIP and block2/block4, late-cut and independent-start-parity diagnostics. Report each native bracket, joint root retention and covariance, missing or ambiguous roots. Refine further only if measured coverage requires it; no systematic sigma inferred from arbitrary model spread.',
                  retention='All observable/winding/RNG/timing/provenance histories, logs and exact checkpoints retained. Full fields through closure audit, then existing adaptive policy. SSD originals held; no relocation/deletion.',
                  resources='At most14 single-thread simulations; preceding finite manager lock must be free. Preserve25GiB reserve and all pending/live/copy reservations without future savings. Friday15Central safety deadline retained.')
        if freeze:
            review=json.loads((m.OUT/'common-xi-refinement-design-review.json').read_text())
            assert review['status']=='accepted_with_qualifications' and review['points_sha256']==points_hash
            previous=m.run.estimated_bytes
            try:m.run.estimated_bytes=storage.estimated_full_bytes;m.run.freeze(NAME,ps,body['purpose'])
            finally:m.run.estimated_bytes=previous
        m.run.atomic_json(m.OUT/'common-xi-refinement-design.json',body)
    print(json.dumps({k:v for k,v in body.items() if k!='points'},indent=2),flush=True)


if __name__=='__main__':main('--freeze' in sys.argv)
