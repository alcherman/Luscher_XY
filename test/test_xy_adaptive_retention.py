"""Retention preserves histories and exact native checkpoint continuation."""
import sys,shutil
from pathlib import Path
import h5py,numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments/2026-10-13-xy-infinite-volume'))
import retain_fields as retention
from supervillain.xy import run_chain

@pytest.mark.parametrize('checkpoint',[100,600])
def test_sparse_archive_exact_continuation(tmp_path,monkeypatch,checkpoint):
    monkeypatch.setattr(retention,'OUT',tmp_path);(tmp_path/'retention-records').mkdir()
    monkeypatch.setattr(retention.storage_admission,'OUT',tmp_path)
    monkeypatch.setattr(retention.helper,'tau_int',lambda a:1.)
    full=tmp_path/'full.h5';archive=tmp_path/'archive.h5'
    kwargs=dict(N=4,D=2,kappa=0.,epsilon=np.pi,warmup=10,seed=412,parallel=False,update='hybrid',clusters_per_step=4,checkpoint_every=100,buffer_rows=100)
    run_chain(full,sweeps=600,**kwargs);shutil.copyfile(full,archive)
    audit=dict(passed=True,sha256=retention.run.sha256(archive))
    result=retention.compact(archive,audit);assert result['retained_fields']<601
    import json
    ledger=tmp_path/'retention-records'/'archive.json';record=json.loads(ledger.read_text());record['committed']=False;retention.run.atomic_json(ledger,record)
    assert retention.compact(archive,audit)['committed']
    with h5py.File(full) as a,h5py.File(archive) as b:
        indices=b['phi'].attrs['retained_indices']
        for name in a:
            if name!='phi':assert retention.equal(a[name][:],b[name][:])
        np.testing.assert_array_equal(a['phi'][indices],b['phi'][indices])
        missing=next(i for i in range(601) if i not in indices);assert np.isnan(b['phi'][missing]).all()
    working=tmp_path/'working.h5';shutil.copyfile(archive,working)
    for path in [full,working]:
        with h5py.File(path,'r+') as h:h.attrs.modify('committed_rows',checkpoint+1)
        run_chain(path,sweeps=checkpoint+5,resume=True,**kwargs)
    with h5py.File(full) as a,h5py.File(working) as b:
        for name in a:assert retention.equal(a[name][checkpoint:],b[name][checkpoint:]),name

def test_hash_failure_preserves_original(tmp_path,monkeypatch):
    monkeypatch.setattr(retention,'OUT',tmp_path)
    monkeypatch.setattr(retention.storage_admission,'OUT',tmp_path)
    p=tmp_path/'source.h5'
    with h5py.File(p,'w') as h:h.create_dataset('phi',data=np.zeros((2,1,2,2)));h.create_dataset('sweep',data=[0,1])
    old=p.read_bytes()
    import pytest
    with pytest.raises(AssertionError,match='hash mismatch'):retention.compact(p,dict(passed=True,sha256='wrong'))
    assert p.read_bytes()==old
