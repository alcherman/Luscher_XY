import importlib.util,sys
from pathlib import Path
import pytest
P=Path(__file__).resolve().parents[1]/'experiments/2026-10-16-xy-3d-spin-followup'
sys.path.insert(0,str(P))
import run_production as m

def test_frozen_scope():
    allpoints=[p for _,qs in m.STAGES for p in m.points(qs)]
    assert len(allpoints)==80 and len({p['id'] for p in allpoints})==80
    assert len({p['seed'] for p in allpoints})==80
    assert {p['N'] for p in allpoints}=={24,32,48,64}
    assert {p['epsilon_over_pi'] for p in allpoints}=={.80060,.80068,.80076,.80084,.80092}
    assert all(p['steps']==100000 and p['warmup']==10000 and p['beta']==0 and p['online'] for p in allpoints)

def test_existing_stage_rejects_mutation():
    ps=m.points([.80076]);batch={'points':ps,'source_hashes':m.c.sources()};m.validate_batch(batch,ps)
    changed=[dict(p) for p in ps];changed[0]['steps']=100001
    with pytest.raises(RuntimeError):m.validate_batch(batch,changed)
    with pytest.raises(RuntimeError):m.validate_batch(dict(batch,source_hashes={}),ps)

def test_pending_other_manifest_prevents_new_stage(tmp_path,monkeypatch):
    import json
    (tmp_path/'batches').mkdir();(tmp_path/'jobs').mkdir();monkeypatch.setattr(m,'P',tmp_path)
    p=m.points([.80076])[0]
    (tmp_path/'batches/x.json').write_text(json.dumps({'points':[p]}))
    assert m.remaining_commitments(set())==[p]
    assert m.remaining_commitments({p['id']})==[]
    (tmp_path/'jobs'/f"{p['id']}.json").write_text(json.dumps({'state':'complete'}))
    assert m.remaining_commitments(set())==[]

def test_sampling_stop_screens():
    import copy
    row=dict(id='a',N=24,point={'epsilon_over_pi':.8},available_blocks=100,statistics=dict(G_mean=[.05,.04],m2_mean=.04,R_mean=.6,block_SE_Gq_Gh_m2=[.001]*3,first_last_half_Gq_Gh_m2=[[.05,.04,.04],[.05,.04,.04]]))
    assert not m.sampling_stops([row])
    bad=copy.deepcopy(row);bad['available_blocks']=10
    assert m.sampling_stops([bad])
    bad=copy.deepcopy(row);bad['statistics']['first_last_half_Gq_Gh_m2'][1][0]=.08
    assert m.sampling_stops([bad])
    bad=copy.deepcopy(row);bad['id']='b';bad['statistics']['G_mean'][0]=.08
    assert m.sampling_stops([row,bad])

def test_absent_completed_stage_does_not_skip(tmp_path,monkeypatch):
    monkeypatch.setattr(m,'P',tmp_path)
    assert not m.completed_stage('none',m.points([.80076]))

def test_verified_prior_stage_skips_while_later_stage_is_pending(tmp_path,monkeypatch):
    import json
    monkeypatch.setattr(m,'P',tmp_path)
    for folder in ['batches','jobs','chains','retention-records']:
        (tmp_path/folder).mkdir()
    point=m.points([.80068])[0];later=m.points([.80092])[0]
    def save(path,value):path.write_text(json.dumps(value))
    save(tmp_path/'batches/prior.json',dict(points=[point],source_hashes=m.c.sources()))
    save(tmp_path/'batches/later.json',dict(points=[later],source_hashes=m.c.sources()))
    row=dict(id=point['id'],N=point['N'],point=point,available_blocks=100,statistics=dict(G_mean=[.05,.04],m2_mean=.04,R_mean=.6,block_SE_Gq_Gh_m2=[.001]*3,first_last_half_Gq_Gh_m2=[[.05,.04,.04],[.05,.04,.04]]))
    save(tmp_path/'prior-sampling.json',dict(rows=[row],stop_flags=[]))
    save(tmp_path/'jobs'/f"{point['id']}.json",dict(state='complete',audit=dict(passed=True,sha256='full-source',online_spin=dict(passed=True))))
    data=tmp_path/'chains'/f"{point['id']}.h5";data.write_bytes(b'compact')
    ledger=tmp_path/'retention-records'/f"{point['id']}.json"
    record=dict(committed=True,verified=True,original_sha256='full-source',compact_bytes=7)
    save(ledger,record)
    assert m.remaining_commitments(set())==[later]
    assert m.completed_stage('prior',[point])
    # Re-entering the actual pending stage still excludes its own commitment.
    assert m.freeze_stage('later',[later])==tmp_path/'batches/later.json'
    save(ledger,dict(record,committed=False))
    assert not m.completed_stage('prior',[point])
