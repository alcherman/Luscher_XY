from pathlib import Path
import sys,json
import numpy as np
import h5py,pytest
sys.path.insert(0,str(Path(__file__).parents[1]/'experiments/2026-10-16-xy-3d-spin-followup'))
import online_spin as o
from supervillain import xy

KW=dict(D=3,N=8,kappa=0.,epsilon=.80075*np.pi,warmup=8,seed=20261016001,backend='numba',parallel=False,interval_phi=np.pi,winding=True,initial_winding=[1,0,0],checkpoint_every=5,update='hybrid',clusters_per_step=4,record_timing=True,buffer_rows=4)

def compare(path0,path1):
    with h5py.File(path0) as a,h5py.File(path1) as b:
        for k in a:
            if k.endswith('_seconds'):continue
            np.testing.assert_array_equal(a[k][:],b[k][:],err_msg=k)
        assert a.attrs['committed_rows']==b.attrs['committed_rows']

def test_correlation_uniform_helix_and_nonmutation():
    L=8;p=np.broadcast_to(2*np.pi*np.arange(L)[:,None,None]/L,(L,L,L)).copy();before=p.copy();g=o.correlations(p)
    np.testing.assert_allclose(g.real,[[0,1,1],[-1,1,1]],atol=1e-14);np.testing.assert_allclose(g.imag[1],0,atol=1e-14);np.testing.assert_array_equal(p,before)
    np.testing.assert_allclose(o.correlations(np.zeros_like(p)),1)

def test_every_native_physical_history_and_rng_identical(tmp_path):
    a=tmp_path/'native.h5';b=tmp_path/'online.h5';xy.run_chain(a,sweeps=35,**KW);o.run_chain(b,sweeps=35,**KW);compare(a,b)
    with h5py.File(b) as h:
        assert h['spin_correlation_real'].shape==(36,2,3)
        assert json.loads(h.attrs['provenance'])['spin_measurement_extension']==o.definition()
        for i in range(36):
            g=o.correlations(h['phi'][i]);np.testing.assert_array_equal(h['spin_correlation_real'][i],g.real);np.testing.assert_array_equal(h['spin_correlation_imag'][i],g.imag)
    assert xy.measure is o._ORIGINAL_MEASURE and xy._provenance is o._ORIGINAL_PROVENANCE

def test_extended_resume_all_physical_histories_exact(tmp_path):
    a=tmp_path/'full.h5';b=tmp_path/'resumed.h5';o.run_chain(a,sweeps=40,**KW);o.run_chain(b,sweeps=14,**KW);o.run_chain(b,sweeps=40,resume=True,**KW);compare(a,b)

def test_native_without_extension_cannot_resume_extended_file(tmp_path):
    p=tmp_path/'extended.h5';o.run_chain(p,sweeps=14,**KW)
    with pytest.raises(ValueError,match='provenance mismatch'):xy.run_chain(p,sweeps=40,resume=True,**KW)
    with h5py.File(p) as h:assert h.attrs['committed_rows']==15

def test_crash_before_first_return_has_pinned_provenance_and_replays(tmp_path,monkeypatch):
    a=tmp_path/'full.h5';b=tmp_path/'crashed.h5';o.run_chain(a,sweeps=40,**KW);original=o.correlations;calls=[0]
    def crash(phi):
        calls[0]+=1
        if calls[0]==13:raise RuntimeError('injected measurement failure')
        return original(phi)
    monkeypatch.setattr(o,'correlations',crash)
    with pytest.raises(RuntimeError,match='injected'):o.run_chain(b,sweeps=40,**KW)
    with h5py.File(b) as h:
        assert json.loads(h.attrs['provenance'])['spin_measurement_extension']==o.definition();assert h.attrs['committed_rows']==11
    monkeypatch.setattr(o,'correlations',original);o.run_chain(b,sweeps=40,resume=True,**KW);compare(a,b)

def test_extension_mismatch_refused_before_writes(tmp_path,monkeypatch):
    p=tmp_path/'extended.h5';o.run_chain(p,sweeps=14,**KW);before=p.read_bytes();original=o.definition
    monkeypatch.setattr(o,'definition',lambda:dict(original(),version=999))
    with pytest.raises(ValueError,match='mismatch'):o.run_chain(p,sweeps=40,resume=True,**KW)
    assert p.read_bytes()==before

import campaign as c

def test_independent_extra_audit_detects_corrupted_measurement(tmp_path):
    path=tmp_path/'online.h5';kw=dict(KW,buffer_rows=100);o.run_chain(path,sweeps=35,**kw)
    spec=dict(N=8,steps=35)
    assert c.extra_audit(path,spec)['passed']
    with h5py.File(path,'r+') as h:h['spin_correlation_real'][17,0,1]+=.001
    with pytest.raises(AssertionError):c.extra_audit(path,spec)

def test_admission_reserves_growth_and_copy_before_fresh_free(tmp_path,monkeypatch):
    import types
    (tmp_path/'jobs').mkdir();(tmp_path/'chains').mkdir();(tmp_path/'chains/a.h5').write_bytes(b'x'*300)
    (tmp_path/'jobs/a.json').write_text(json.dumps(dict(state='running',point=dict(id='a'))))
    monkeypatch.setattr(c,'OUT',tmp_path);monkeypatch.setattr(c,'TARGET',tmp_path/'chains');monkeypatch.setattr(c,'RESERVE',100)
    monkeypatch.setattr(c.storage,'estimated_full_bytes',lambda p:1000);monkeypatch.setattr(c.storage,'reserved_remaining',lambda:100)
    monkeypatch.setattr(c.shutil,'disk_usage',lambda p:types.SimpleNamespace(free=2000))
    x=c.available_for([dict(id='b')]);assert x['active_remaining_bytes']==700 and x['admitted']
    monkeypatch.setattr(c.shutil,'disk_usage',lambda p:types.SimpleNamespace(free=1700));assert not c.available_for([dict(id='b')])['admitted']

def test_launch_claim_survives_parent_crash_and_obsolete_child_cannot_write(tmp_path,monkeypatch):
    import types
    for name in ['jobs','logs','chains']:(tmp_path/name).mkdir()
    monkeypatch.setattr(c,'OUT',tmp_path);monkeypatch.setattr(c.storage,'acquire',lambda:(tmp_path/'admission.lock').open('a+'))
    monkeypatch.setattr(c,'available_for',lambda ps:dict(admitted=True));point=dict(id='job1')
    original_atomic=c.atomic
    def fail_after_spawn(path,value):
        if value['state']=='running':raise KeyboardInterrupt('injected parent loss after Popen')
        original_atomic(path,value)
    def fake_popen(*args,**kwargs):
        j=c.read(tmp_path/'jobs/job1.json');assert j['state']=='starting' and j['launch_token'];return types.SimpleNamespace(pid=777)
    monkeypatch.setattr(c.subprocess,'Popen',fake_popen);monkeypatch.setattr(c,'atomic',fail_after_spawn)
    with pytest.raises(KeyboardInterrupt):c.admit_and_spawn(point,tmp_path/'batch.json')
    stale=c.read(tmp_path/'jobs/job1.json');assert stale['state']=='starting';oldtoken=stale['launch_token']
    monkeypatch.setattr(c,'atomic',original_atomic);recovered=c.recover_starting(tmp_path/'jobs/job1.json')
    assert recovered['state']=='interrupted' and recovered['launch_token']!=oldtoken
    manifest=tmp_path/'batch.json';manifest.write_text(json.dumps(dict(source_hashes={'pinned':True},points=[point])))
    monkeypatch.setattr(c,'sources',lambda:{'pinned':True});monkeypatch.setattr(c,'route',lambda:None)
    c.worker(manifest,'job1',oldtoken)
    assert not (tmp_path/'chains/job1.h5').exists();assert c.read(tmp_path/'jobs/job1.json')==recovered

def test_recovery_does_not_overwrite_completed_owner(tmp_path,monkeypatch):
    (tmp_path/'jobs').mkdir();j=tmp_path/'jobs/done.json';state=dict(state='complete',point=dict(id='done'),launch_token='old')
    j.write_text(json.dumps(state));monkeypatch.setattr(c,'OUT',tmp_path)
    assert c.recover_starting(j)==state and c.read(j)==state

def test_native_targets_strictly_exceed_warmup():
    assert c.generation_targets(dict(steps=4000,warmup=1000))==[2000,3000,4000]
    assert c.generation_targets(dict(steps=600,warmup=150))==[600]
    assert c.generation_targets(dict(steps=40,warmup=10))==[40]
    assert min(c.generation_targets(dict(steps=100000,warmup=10000)))==11000
