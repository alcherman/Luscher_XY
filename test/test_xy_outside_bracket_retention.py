"""Outside-bracket policy protects inference neighborhoods and exact restart state."""
import json,shutil,sys
from pathlib import Path
import numpy as np,h5py,pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments/2026-10-13-xy-infinite-volume'))
import retain_fields as retention
import design_outside_bracket_retention as design
from supervillain.xy import run_chain

def test_classifier_union_pairs_unknown_and_endpoints():
    p=dict(D=2,N=128,beta=.9,epsilon_over_pi=.8)
    k=design.key(2,.8,128,False)
    a=dict(resolved=True,grid=[.1,.2,.3,.9],protected=[.1,.3])
    b=dict(resolved=True,grid=[.1,.2,.3,.9],protected=[.2,.9])
    assert design.classify(p,{k:[a]})[0]
    assert not design.classify(p,{k:[a,b]})[0]
    assert not design.classify(p,{k:[dict(a,resolved=False)]})[0]
    assert not design.classify(dict(p,beta=.8),{k:[a]})[0]
    assert not design.classify(dict(p,beta=0),{k:[a]})[0]
    assert design.protected_interval([0,.1,.2,.3,.4,.8],[.2,.3])==pytest.approx([.1,.4])

@pytest.mark.parametrize('previous',[False,True])
def test_checkpoints_only_all_histories_and_native_replay(tmp_path,monkeypatch,previous):
    monkeypatch.setattr(retention,'OUT',tmp_path);monkeypatch.setattr(retention.storage_admission,'OUT',tmp_path)
    monkeypatch.setattr(retention.helper,'tau_int',lambda a:1.)
    full=tmp_path/'full.h5';archive=tmp_path/'archive.h5'
    kwargs=dict(N=4,D=2,kappa=0.,epsilon=np.pi,warmup=10,seed=412,parallel=False,update='hybrid',clusters_per_step=4,checkpoint_every=100,buffer_rows=100)
    run_chain(full,sweeps=600,**kwargs);shutil.copyfile(full,archive)
    audit=dict(passed=True,sha256=retention.run.sha256(archive))
    old=retention.compact(archive,audit) if previous else None
    (tmp_path/'outside-bracket-retention-design-20261014.json').write_text(json.dumps(dict(eligible={'archive':dict(audit_sha256=audit['sha256'],reason='test resolved far bracket')})))
    r=retention.compact(archive,audit)
    assert r['policy']['mode']=='checkpoints_only_outside_brackets'
    assert r['original_sha256']==audit['sha256']
    if previous:assert r['previous_retention']==old
    indices=np.arange(0,601,100)
    with h5py.File(full) as a,h5py.File(archive) as b:
        np.testing.assert_array_equal(b['phi'].attrs['retained_indices'],indices)
        assert np.isnan(b['phi'][1]).all()
        for name in a:
            assert retention.equal(a[name][indices] if name=='phi' else a[name][:],b[name][indices] if name=='phi' else b[name][:])
    # Recovery from replacement before event commit, then idempotence.
    ledger=tmp_path/'retention-records'/'archive.json';prepared=json.loads(ledger.read_text());prepared['committed']=False;retention.run.atomic_json(ledger,prepared)
    assert retention.compact(archive,audit)['committed']
    assert retention.compact(archive,audit)['compact_sha256']==r['compact_sha256']
    working=tmp_path/'working.h5';shutil.copyfile(archive,working)
    for path in [full,working]:
        with h5py.File(path,'r+') as h:h.attrs.modify('committed_rows',101)
        run_chain(path,sweeps=105,resume=True,**kwargs)
    with h5py.File(full) as a,h5py.File(working) as b:
        for name in a:assert retention.equal(a[name][100:],b[name][100:])

def test_second_reduction_hash_mismatch_preserves_archive(tmp_path,monkeypatch):
    monkeypatch.setattr(retention,'OUT',tmp_path);monkeypatch.setattr(retention.storage_admission,'OUT',tmp_path);monkeypatch.setattr(retention.helper,'tau_int',lambda a:1.)
    archive=tmp_path/'archive.h5'
    run_chain(archive,sweeps=200,N=4,D=2,kappa=0.,epsilon=np.pi,warmup=10,seed=99,parallel=False,update='hybrid',clusters_per_step=4,checkpoint_every=100)
    audit=dict(passed=True,sha256=retention.run.sha256(archive));retention.compact(archive,audit)
    (tmp_path/'outside-bracket-retention-design-20261014.json').write_text(json.dumps(dict(eligible={'archive':dict(audit_sha256=audit['sha256'])})))
    with h5py.File(archive,'r+') as h:h['magnetization_squared'][1]=99.
    original=archive.read_bytes()
    with pytest.raises(AssertionError,match='hash mismatch'):retention.compact(archive,audit)
    assert archive.read_bytes()==original

def test_second_reduction_crash_before_replace_recovers(tmp_path,monkeypatch):
    monkeypatch.setattr(retention,'OUT',tmp_path);monkeypatch.setattr(retention.storage_admission,'OUT',tmp_path);monkeypatch.setattr(retention.helper,'tau_int',lambda a:1.)
    archive=tmp_path/'archive.h5'
    run_chain(archive,sweeps=200,N=4,D=2,kappa=0.,epsilon=np.pi,warmup=10,seed=99,parallel=False,update='hybrid',clusters_per_step=4,checkpoint_every=100)
    audit=dict(passed=True,sha256=retention.run.sha256(archive));old=retention.compact(archive,audit)
    (tmp_path/'outside-bracket-retention-design-20261014.json').write_text(json.dumps(dict(eligible={'archive':dict(audit_sha256=audit['sha256'])})))
    replace=retention.os.replace
    def fail(src,dst):
        if str(src).endswith('.retention.tmp'):raise RuntimeError('before replace')
        return replace(src,dst)
    monkeypatch.setattr(retention.os,'replace',fail)
    with pytest.raises(RuntimeError,match='before replace'):retention.compact(archive,audit)
    assert retention.run.sha256(archive)==old['compact_sha256']
    monkeypatch.setattr(retention.os,'replace',replace)
    result=retention.compact(archive,audit)
    assert result['policy']['mode']=='checkpoints_only_outside_brackets'
    assert len(list((tmp_path/'retention-attempts').glob('*.json')))==1
