"""Reviewed, exact-coordinate finer beta-zero 3D spin-dimension scan."""
from pathlib import Path
import json,sys
import numpy as np
import measure as m
sys.path.insert(0,str(m.CAMPAIGN))
import storage_admission as storage

NAME='d3-zero-beta-spin-fine-epsilon-20261015-stage1'
GRID={8:(80400,81600,100),12:(80300,80900,100),16:(80250,80500,50),24:(80150,80350,25),32:(80150,80275,25)}

def points():
    result=[]
    for N,(first,last,stride) in GRID.items():
        for code in range(first,last+1,stride):
            for replica in [20,21]:
                q=code/100000
                p=m.run.point(3,80,0,N,replica,12000)
                p.update(id=f'spinzero3_e{code:06d}_N{N:03d}_r{replica:02d}_s012000',epsilon_over_pi=q,epsilon=q*np.pi,
                         seed=int(np.random.SeedSequence([2026101503,3,code,N,replica]).generate_state(1,dtype=np.uint64)[0]),
                         initial_winding=[replica%2,0,0],campaign_role='zero_beta_endpoint',epsilon_sampling_1e5=code,
                         study_role='spin_scaling_dimension_fine_scan',primary_observables=['magnetization_squared','magnetization_fourth','magnetization_abs','binder_contour'])
                result.append(p)
    return result

def main(freeze=False):
    ps=points();old={p['id']:p for f in (m.CAMPAIGN/'batches').glob('*.json') for p in json.loads(f.read_text())['points']}
    assert len(ps)==82 and len({p['id'] for p in ps})==82 and len({p['seed'] for p in ps})==82
    assert not {p['seed'] for p in ps}&{p['seed'] for p in old.values() if p['id'] not in {v['id'] for v in ps}}
    assert all(p['id'] not in old or old[p['id']]==p for p in ps)
    committed={}
    for p in old.values():
        f=m.CAMPAIGN/'jobs'/f"{p['id']}.json"
        j=json.loads(f.read_text()) if f.exists() else {'state':'pending'}
        if j['state'] in ['pending','interrupted','running'] and p['id'] not in {v['id'] for v in ps}:
            chain=m.CAMPAIGN/'chains'/f"{p['id']}.h5"
            committed[p['id']]=max(0,storage.estimated_full_bytes(p)-(chain.stat().st_size if chain.exists() else 0))
    with storage.acquire():
        copies=storage.reserved_remaining();physical=m.run.shutil.disk_usage(m.CAMPAIGN/'chains').free;estimate=sum(storage.estimated_full_bytes(p) for p in ps)
        headroom=physical-copies-sum(committed.values())-estimate-25*2**30
        if headroom<0:raise RuntimeError('No safe full-retention budget')
        source=m.OUT/'zero-beta'/json.loads((m.OUT/'zero-beta/latest-primary.json').read_text())['folder']/'analysis.json'
        body=dict(name=NAME,created=m.run.now(),source_analysis=str(source),source_analysis_sha256=m.run.sha256(source),points=ps,trajectories=len(ps),physical_coordinate_count=len(ps)//2,
                  full_transient_bytes=estimate,physical_free_bytes=physical,prior_pending_live_growth_bytes=sum(committed.values()),retention_copy_growth_bytes=copies,headroom_after_25GiB_reserve_bytes=headroom,
                  forecast=dict(worker_hours=2.45,initial_calendar_allowance_minutes=[20,60],basis='Existing generation+full audit timings, no thinning credit. Actual throughput and statistical precision remain to be measured.'),
                  scientific_scope='First resolve interpolation at fixed beta0 with existing L8/12/16/24/32, two independent seeded zero/unit-x-winding starts,12000cycles,exact fine epsilon. No Delta/nu assumed. Measure coverage/mixing/cut/block/start/threshold/size sensitivity before more statistics or larger boxes. Standard 3D O2 WF comparison only; user clarified the linked paper focuses N>2.',
                  retention='Keep all histories/logs/RNG/provenance/checkpoints. Existing explicitly authorized adaptive field retention applies only after full closure audit. Keep SSD originals; no archive relocation/deletion.',
                  review='Independent scientific review requested before freeze; launch authorized by direct user request October15.')
        if freeze:
            review=json.loads((m.OUT/'fine-zero-beta-design-review.json').read_text());assert review['status']=='accepted_with_qualifications'
            assert review['points_sha256']==m.hashlib.sha256(json.dumps(ps,sort_keys=True).encode()).hexdigest()
            previous=m.run.estimated_bytes
            try:
                m.run.estimated_bytes=storage.estimated_full_bytes
                m.run.freeze(NAME,ps,body['scientific_scope'])
            finally:m.run.estimated_bytes=previous
        m.run.atomic_json(m.OUT/'fine-zero-beta-stage1-design.json',body)
    print(json.dumps({k:v for k,v in body.items() if k!='points'},indent=2),flush=True)

if __name__=='__main__':main('--freeze' in sys.argv)
