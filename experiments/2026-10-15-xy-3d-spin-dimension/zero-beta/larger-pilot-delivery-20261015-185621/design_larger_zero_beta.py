"""N48/N64 beta-zero pilot: reviewed finite windows, cost and full storage."""
import json,sys
import numpy as np
import measure as m
sys.path.insert(0,str(m.CAMPAIGN))
import storage_admission as storage

NAME='d3-zero-beta-spin-N48-N64-pilot-20261015'
GRID={48:(80075,80225,25),64:(80050,80200,25)}

def points():
    result=[]
    for N,(first,last,stride) in GRID.items():
        for code in range(first,last+1,stride):
            for replica in [30,31]:
                q=code/100000;p=m.run.point(3,80,0,N,replica,20000)
                p.update(id=f'spinpilot3_e{code:06d}_N{N:03d}_r{replica:02d}_s020000',epsilon_over_pi=q,epsilon=q*np.pi,
                         seed=int(np.random.SeedSequence([2026101504,3,code,N,replica]).generate_state(1,dtype=np.uint64)[0]),
                         initial_winding=[replica%2,0,0],campaign_role='zero_beta_endpoint',epsilon_sampling_1e5=code,
                         study_role='spin_scaling_dimension_larger_size_pilot',primary_observables=['magnetization_squared','magnetization_fourth','magnetization_abs','binder_contour'])
                result.append(p)
    return result

def main(freeze=False):
    ps=points();old={p['id']:p for f in (m.CAMPAIGN/'batches').glob('*.json') for p in json.loads(f.read_text())['points']};ids={p['id'] for p in ps}
    assert len(ps)==28 and len(ids)==28 and len({p['seed'] for p in ps})==28
    assert not {p['seed'] for p in ps}&{p['seed'] for p in old.values() if p['id'] not in ids}
    assert all(p['id'] not in old or old[p['id']]==p for p in ps)
    assert all(2*np.pi/p['N']<p['epsilon'] for p in ps)
    commitments={}
    for p in old.values():
        f=m.CAMPAIGN/'jobs'/f"{p['id']}.json";j=json.loads(f.read_text()) if f.exists() else {'state':'pending'}
        if j['state'] in ['pending','interrupted','running'] and p['id'] not in ids:
            chain=m.CAMPAIGN/'chains'/f"{p['id']}.h5";commitments[p['id']]=max(0,storage.estimated_full_bytes(p)-(chain.stat().st_size if chain.exists() else 0))
    prediction=m.OUT/'larger-size-epsilon-window-predictions.json';source=m.OUT/'zero-beta'/json.loads((m.OUT/'zero-beta/latest-fine-primary.json').read_text())['folder']/'analysis.json'
    with storage.acquire():
        copies=storage.reserved_remaining();physical=m.run.shutil.disk_usage(m.CAMPAIGN/'chains').free;estimate=sum(storage.estimated_full_bytes(p) for p in ps);headroom=physical-copies-sum(commitments.values())-estimate-25*2**30
        if headroom<0:raise RuntimeError('Unsafe full-retention budget')
        body=dict(name=NAME,created=m.run.now(),source_analysis=str(source),source_analysis_sha256=m.run.sha256(source),prediction_source=str(prediction),prediction_sha256=m.run.sha256(prediction),points=ps,
                  trajectories=len(ps),physical_coordinate_count=len(ps)//2,full_transient_bytes=estimate,physical_free_bytes=physical,prior_pending_live_growth_bytes=sum(commitments.values()),retention_copy_growth_bytes=copies,headroom_after_25GiB_reserve_bytes=headroom,
                  forecast=dict(worker_hours=17.2,initial_calendar_allowance_hours=[1,3],basis='N32 measured233.323sec/12k generation+audit, volume-only N48/N64 extrapolation plus20k/12k. Actual critical slowing, audit/I/O, storage and memory must be measured; not a bound.'),
                  scientific_scope='User authorized N48^3/N64^3 pilot October15 before longer runs. Seven size-specific epsilon knots per box cover stress predictions rather than blindly filling a large grid; two fresh explicitzero/unit-x winding starts,20kcycles,betafixed0. Assess coverage, interpolation, autocorrelation, starts/cuts/blocks, actual resourcecost and size/correction drift before any longer or larger campaign. No Delta/nu imposed, no claim .005 total precision, no N96 or high-statistics follow-on authorized.',
                  retention='Preserve all scalar/winding/RNG/timing/provenance histories, logs and checkpoints. Full fields through independent closure audit, then existing user-authorized adaptive retention. SSD originals held; no transfer/deletion.',
                  resource_policy='At most14 single-thread simulations; wait for previous finite retention manager to exit before opening gate. Preserve25GiB reserve and all known pending/live/copy commitments. No hypothetical recovery.')
        if freeze:
            review=json.loads((m.OUT/'larger-zero-beta-design-review.json').read_text());assert review['status']=='accepted_with_qualifications'
            assert review['points_sha256']==m.hashlib.sha256(json.dumps(ps,sort_keys=True).encode()).hexdigest()
            previous=m.run.estimated_bytes
            try:
                m.run.estimated_bytes=storage.estimated_full_bytes;m.run.freeze(NAME,ps,body['scientific_scope'])
            finally:m.run.estimated_bytes=previous
        m.run.atomic_json(m.OUT/'larger-zero-beta-pilot-design.json',body)
    print(json.dumps({k:v for k,v in body.items() if k!='points'},indent=2),flush=True)

if __name__=='__main__':main('--freeze' in sys.argv)
