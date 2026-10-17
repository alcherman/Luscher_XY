"""Production-window, coordinate reuse and finite-manager admission guards."""
import json
import sys
from pathlib import Path
import numpy as np
import pytest

P=Path(__file__).resolve().parents[1]/'experiments/2026-10-15-xy-3d-spin-dimension'
sys.path.insert(0,str(P))
import design_common_xi_refinement as d
import run_common_xi_refinement as runner


def test_fine_common_knots_and_endpoint_lower_guards():
    ps,reuse=d.make_points([])
    assert len(ps)==124 and not reuse
    assert len({p['id'] for p in ps})==len(ps)==len({p['seed'] for p in ps})
    for N in [6,8,12,16,24,32]:
        codes=sorted({p['epsilon_sampling_1e6'] for p in ps if p['N']==N})
        assert codes==d.SMALL_GRID
        assert max(np.diff([x for x in codes if 800500<=x<=801250]))==125
    assert min(p['epsilon_over_pi'] for p in ps if p['N']==48)==.80025
    assert all(p['D']==3 and p['beta']==0 and p['N']<=64 for p in ps)
    for p in ps:
        assert p['epsilon']==pytest.approx(np.pi*p['epsilon_over_pi'],abs=1e-15)
        assert p['initial_winding']==[p['replica']%2,0,0]
        assert p['steps']==(12000 if p['N']<=32 else 20000)


def test_reuse_requires_sufficient_steps_and_both_winding_starts():
    ps,_=d.make_points([])
    existing=[dict(p) for p in ps if p['N'] in [24,32] and p['epsilon_sampling_1e6']==801500]
    new,reuse=d.make_points(existing)
    assert len(new)==120 and len(reuse)==2
    for p in existing:p['steps']=4000
    assert len(d.make_points(existing)[0])==124
    for p in existing:p['steps']=12000;p['initial_winding']=[0,0,0]
    assert len(d.make_points(existing)[0])==124


def test_manager_refuses_unrelated_pending_but_preserves_deferred(tmp_path,monkeypatch):
    monkeypatch.setattr(runner.m,'CAMPAIGN',tmp_path)
    (tmp_path/'batches').mkdir();(tmp_path/'jobs').mkdir()
    (tmp_path/'batches'/'fake.json').write_text(json.dumps(dict(points=[dict(id='allowed'),dict(id='other')])))
    (tmp_path/'jobs'/'other.json').write_text(json.dumps(dict(state='deferred')))
    runner.verify_no_unrelated([dict(id='allowed')])
    (tmp_path/'jobs'/'other.json').write_text(json.dumps(dict(state='pending')))
    with pytest.raises(RuntimeError,match='Unrelated production commitment'):runner.verify_no_unrelated([dict(id='allowed')])


def test_provenance_accepts_snapshot_folder_or_analysis_file(tmp_path):
    snapshot=tmp_path/'snapshot';snapshot.mkdir();analysis=snapshot/'analysis.json';analysis.write_text('{}')
    assert runner.analysis_file(snapshot)==analysis
    assert runner.analysis_file(analysis)==analysis


def test_running_to_complete_handoff_preserves_provenance_path(tmp_path,monkeypatch):
    """Exercise the live progress loop before the exact-cohort analysis handoff."""
    campaign=tmp_path/'campaign';out=tmp_path/'study'
    (campaign/'batches').mkdir(parents=True);(campaign/'jobs').mkdir();out.mkdir()
    new=dict(id='new',replica=40);old=dict(id='old',replica=20)
    (campaign/'batches'/f'{runner.NAME}.json').write_text(json.dumps(dict(points=[new])))
    for point in [old,new]:
        (campaign/'jobs'/f"{point['id']}.json").write_text(json.dumps(dict(state='complete',point=point,audit=dict(passed=True))))
    prior=out/'prior';prior.mkdir();(prior/'analysis.json').write_text(json.dumps(dict(source_jobs=[dict(id='old')])) )
    ledger=out/'prior-delivery';ledger.mkdir()
    (ledger/'provenance.json').write_text(json.dumps(dict(source_analyses=dict(baseline=dict(path=str(prior),sha256=runner.m.run.sha256(prior/'analysis.json'))))))
    (out/'common-xi-refinement-design.json').write_text(json.dumps(dict(source_ledger=str(ledger/'quotient-ledger.json'))))
    monkeypatch.setattr(runner,'OUT',out);monkeypatch.setattr(runner.m,'CAMPAIGN',campaign)
    # This handoff fixture exercises a live campaign before its historical deadline.
    class BeforeDeadline(runner.m.run.datetime):
        @classmethod
        def now(cls,tz=None):
            return cls(2026,10,15,12,tzinfo=runner.m.run.CENTRAL).astimezone(tz)
    monkeypatch.setattr(runner.m.run,'datetime',BeforeDeadline)
    progress=iter([({'running':1},[]),({'running':1},[]),({'complete':1},[])])
    monkeypatch.setattr(runner,'states',lambda points:next(progress))
    monkeypatch.setattr(runner,'coordinator_pid',lambda:777)
    monkeypatch.setattr(runner.time,'sleep',lambda seconds:None)
    monkeypatch.setattr(runner.os,'kill',lambda pid,signal:None)
    class Awake:
        def poll(self):return None
        def terminate(self):pass
    monkeypatch.setattr(runner.subprocess,'Popen',lambda *args,**kwargs:Awake())
    controls=[];monkeypatch.setattr(runner,'control',lambda **kw:controls.append(kw))
    class HandoffReached(Exception):pass
    def analysis(tag,args):
        assert tag=='primary'
        pinned=json.loads(Path(args[1]).read_text())
        assert {s['id'] for s in pinned['source_jobs']}=={'old','new'}
        assert pinned['previous_analysis']==str(prior/'analysis.json')
        assert pinned['previous_sha256']==runner.m.run.sha256(prior/'analysis.json')
        raise HandoffReached
    monkeypatch.setattr(runner,'analyze',analysis)
    with pytest.raises(HandoffReached):runner.main()
    assert controls[0]['launch_allowed_ids']==['new']
    assert controls[0]['max_workers']==14
    assert any(c.get('launch_new') is False for c in controls[1:])
