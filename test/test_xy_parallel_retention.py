"""Concurrent conversion keeps per-file exactness and a shared physical budget."""
import json,multiprocessing,queue,sys
from pathlib import Path
import pytest,numpy as np,h5py
CAMPAIGN=Path(__file__).resolve().parents[1]/'experiments/2026-10-13-xy-infinite-volume'
sys.path.insert(0,str(CAMPAIGN))
import retain_fields as retention
import storage_admission as admission
from supervillain.xy import run_chain


def copy_worker(root,identifier,audit,messages,gate,crash=False):
    import retain_fields as r
    r.OUT=Path(root);r.storage_admission.OUT=r.OUT;r.helper.tau_int=lambda a:1.
    original_hash=r.run.sha256;original_replace=r.os.replace
    def sha(path):
        if Path(path)==r.OUT/'chains'/f'{identifier}.h5':
            messages.put(('hash_entered',identifier));assert gate.wait(20)
        return original_hash(path)
    def replace(src,dst):
        original_replace(src,dst)
        if crash and str(src).endswith('.retention.tmp'):raise RuntimeError('after_replace_injection')
    r.run.sha256=sha;r.os.replace=replace
    try:
        result=r.compact(r.OUT/'chains'/f'{identifier}.h5',audit,parallel=True)
        messages.put(('done',identifier,result['compact_sha256']))
    except BaseException as error:messages.put(('error',identifier,repr(error)))


def maintenance_worker(root,messages):
    import fcntl
    root=Path(root);messages.put('maintenance_attempting')
    with (root/'retention-copy.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX);messages.put('maintenance_acquired')


def create_source(root,identifier,seed):
    p=root/'chains';p.mkdir(exist_ok=True);path=p/f'{identifier}.h5'
    run_chain(path,sweeps=220,N=4,D=2,kappa=0.,epsilon=np.pi,warmup=10,seed=seed,parallel=False,update='hybrid',clusters_per_step=4,checkpoint_every=100)
    return dict(passed=True,sha256=retention.run.sha256(path))


@pytest.mark.parametrize('crash',[False,True])
def test_two_copies_share_lock_and_exactness_with_recovery(tmp_path,monkeypatch,crash):
    audits={i:create_source(tmp_path,i,seed) for i,seed in [('a',432),('b',543)]}
    originals={i:(tmp_path/'chains'/f'{i}.h5').read_bytes() for i in audits}
    ctx=multiprocessing.get_context('spawn');messages=ctx.Queue();gate=ctx.Event()
    procs=[ctx.Process(target=copy_worker,args=(str(tmp_path),i,audits[i],messages,gate,crash and i=='a')) for i in audits]
    for proc in procs:proc.start()
    # Both must pass the SH lock before either hash/copy proceeds.
    entered=[messages.get(timeout=20),messages.get(timeout=20)]
    assert {tuple(x) for x in entered}=={('hash_entered','a'),('hash_entered','b')}
    maintenance_messages=ctx.Queue();maintenance=ctx.Process(target=maintenance_worker,args=(str(tmp_path),maintenance_messages));maintenance.start()
    assert maintenance_messages.get(timeout=10)=='maintenance_attempting'
    with pytest.raises(queue.Empty):maintenance_messages.get(timeout=.2)
    gate.set();results=[messages.get(timeout=20),messages.get(timeout=20)]
    for proc in procs:proc.join(timeout=20);assert proc.exitcode==0
    assert maintenance_messages.get(timeout=20)=='maintenance_acquired';maintenance.join(timeout=10)
    if crash:assert any(x[0]=='error' and 'after_replace_injection' in x[2] for x in results)
    else:assert all(x[0]=='done' for x in results)
    monkeypatch.setattr(retention,'OUT',tmp_path);monkeypatch.setattr(admission,'OUT',tmp_path)
    for identifier in audits:
        path=tmp_path/'chains'/f'{identifier}.h5'
        record=retention.compact(path,audits[identifier],parallel=True)
        assert record['committed'] and record['original_sha256']==audits[identifier]['sha256']
        full=tmp_path/f'original-{identifier}.h5';full.write_bytes(originals[identifier])
        with h5py.File(full) as source,h5py.File(path) as reduced:
            indices=reduced['phi'].attrs['retained_indices']
            for name in source:
                assert retention.equal(source[name][indices] if name=='phi' else source[name][:],reduced[name][indices] if name=='phi' else reduced[name][:])
    events=[json.loads(line) for line in (tmp_path/'retention-events.jsonl').read_text().splitlines()]
    assert {e['id'] for e in events}=={'a','b'}
    assert admission.reserved_remaining()==0


def test_all_temporary_leases_count_and_release_only_owner(tmp_path,monkeypatch):
    monkeypatch.setattr(admission,'OUT',tmp_path)
    a=tmp_path/'a.tmp';b=tmp_path/'b.tmp'
    lock=admission.acquire()
    try:
        admission.reserve('a',a,1000);admission.reserve('b',b,2000)
        a.write_bytes(b'x'*100);b.write_bytes(b'y'*300)
        assert admission.reserved_remaining()==2600
        admission.release('a');assert admission.reserved_remaining()==1700
        # A legacy lease remains additive during the safe upgrade.
        (tmp_path/'retention-space-reservation.json').write_text(json.dumps(dict(id='old',temporary_path=str(tmp_path/'old.tmp'),estimated_total_bytes=500)))
        assert admission.reserved_remaining()==2200
        admission.release('unknown');assert admission.reserved_remaining()==2200
        admission.clear_retention_reservations();assert admission.reserved_remaining()==0
    finally:lock.close()


def budget_worker(root,identifier,estimate,messages):
    import storage_admission as a
    from types import SimpleNamespace
    a.OUT=Path(root)
    a.run.shutil.disk_usage=lambda path:SimpleNamespace(free=150-sum(f.stat().st_size for f in Path(root).glob('*.tmp')))
    lock=a.acquire()
    try:messages.put((identifier,a.admit_retention_copy(identifier,Path(root)/f'{identifier}.tmp',estimate,reserve_bytes=0)))
    finally:lock.close()


def test_cross_process_admission_counts_live_copy_remaining(tmp_path,monkeypatch):
    monkeypatch.setattr(admission,'OUT',tmp_path)
    lock=admission.acquire()
    try:admission.reserve('a',tmp_path/'a.tmp',100)
    finally:lock.close()
    (tmp_path/'a.tmp').write_bytes(b'x'*40)
    ctx=multiprocessing.get_context('spawn');messages=ctx.Queue()
    # Physical free110, existing copy still needs60:80 must fail,40 can fit.
    for identifier,estimate,expected in [('b',80,False),('c',40,True)]:
        child=ctx.Process(target=budget_worker,args=(str(tmp_path),identifier,estimate,messages));child.start()
        assert messages.get(timeout=20)==(identifier,expected);child.join(timeout=10);assert child.exitcode==0
    assert admission.reserved_remaining()==100


def test_copy_growth_between_snapshot_and_free_is_conservative(tmp_path,monkeypatch):
    from types import SimpleNamespace
    monkeypatch.setattr(admission,'OUT',tmp_path);source=tmp_path/'a.tmp'
    admission.reserve('a',source,100)
    original=admission.reserved_remaining
    def snapshot_then_grow():
        result=original();source.write_bytes(b'x'*30);return result
    monkeypatch.setattr(admission,'reserved_remaining',snapshot_then_grow)
    monkeypatch.setattr(admission.run.shutil,'disk_usage',lambda p:SimpleNamespace(free=120))
    assert not admission.admit_retention_copy('b',tmp_path/'b.tmp',30,reserve_bytes=0)


def test_two_workers_cannot_convert_same_archive_twice(tmp_path):
    audit=create_source(tmp_path,'same',654);ctx=multiprocessing.get_context('spawn');messages=ctx.Queue();gate=ctx.Event()
    one=ctx.Process(target=copy_worker,args=(str(tmp_path),'same',audit,messages,gate));one.start()
    assert messages.get(timeout=20)==('hash_entered','same')
    two=ctx.Process(target=copy_worker,args=(str(tmp_path),'same',audit,messages,gate));two.start()
    with pytest.raises(queue.Empty):messages.get(timeout=.3)
    gate.set();results=[messages.get(timeout=20),messages.get(timeout=20)]
    one.join(timeout=10);two.join(timeout=10)
    assert one.exitcode==two.exitcode==0
    assert all(r[0]=='done' for r in results) and results[0][2]==results[1][2]
    assert len((tmp_path/'retention-events.jsonl').read_text().splitlines())==1


def test_space_rejection_creates_no_temporary_archive(tmp_path,monkeypatch):
    audit=create_source(tmp_path,'reserve',765)
    monkeypatch.setattr(retention,'OUT',tmp_path);monkeypatch.setattr(admission,'OUT',tmp_path)
    monkeypatch.setattr(admission,'admit_retention_copy',lambda *args,**kwargs:False)
    path=tmp_path/'chains/reserve.h5';before=path.read_bytes()
    with pytest.raises(retention.StorageAdmissionDeferred):retention.compact(path,audit,parallel=True)
    assert path.read_bytes()==before
    assert not list((tmp_path/'chains').glob('.*.retention.tmp'))


def test_completed_temp_rename_during_admission_is_conservative(tmp_path,monkeypatch):
    monkeypatch.setattr(admission,'OUT',tmp_path)
    temp=tmp_path/'finished.tmp';temp.write_bytes(b'x'*80)
    admission.reserve('a',temp,100)
    # A peer can finish its atomic replacement outside the admission lock.
    # Reproduce a stale positive exists() observation followed by missing stat.
    temp.replace(tmp_path/'completed.h5')
    original_exists=Path.exists
    monkeypatch.setattr(Path,'exists',lambda path:True if path==temp else original_exists(path))
    with pytest.raises(FileNotFoundError):
        _=temp.stat().st_size if temp.exists() else 0
    assert admission.reserved_remaining()==100
    # Full reservation until owner release overcounts; it cannot invent space.
    admission.release('a');assert admission.reserved_remaining()==0
