import importlib.util
from pathlib import Path
import pytest
SPEC=importlib.util.spec_from_file_location('archive',Path(__file__).parents[1]/'experiments/2026-10-16-mc-data-hdd-archive/archive_campaign.py')
a=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(a)

def setup(monkeypatch,tmp_path):
    job=tmp_path/'job';job.mkdir();dst=tmp_path/'destination';dst.mkdir()
    monkeypatch.setattr(a,'JOB',job);monkeypatch.setattr(a,'DEST',dst)
    src=tmp_path/'sample.h5';src.write_bytes(b'full observable histories and retained fields'*100)
    item=dict(source=str(src),relative='raw/sample.h5',signature=a.signature(src),remove=True)
    return src,item

def test_exact_copy_and_full_hash(monkeypatch,tmp_path):
    src,i=setup(monkeypatch,tmp_path)
    records=a.transfer([i]);verified=a.verify([i],records)
    assert (a.DEST/i['relative']).read_bytes()==src.read_bytes()
    assert i['relative'] in verified
    assert src.exists()

def test_corrupt_destination_stops_before_delete(monkeypatch,tmp_path):
    src,i=setup(monkeypatch,tmp_path);records=a.transfer([i])
    (a.DEST/i['relative']).write_bytes(b'corruption')
    with pytest.raises(RuntimeError,match='Checksum mismatch'):a.verify([i],records)
    assert src.exists()

def test_changed_source_rejected(monkeypatch,tmp_path):
    src,i=setup(monkeypatch,tmp_path);src.write_bytes(b'changed')
    with pytest.raises(RuntimeError,match='Source changed'):a.transfer([i])
    assert src.exists()

def test_existing_unjournaled_destination_not_overwritten(monkeypatch,tmp_path):
    src,i=setup(monkeypatch,tmp_path);dst=a.DEST/i['relative'];dst.parent.mkdir();dst.write_bytes(b'other')
    with pytest.raises(RuntimeError,match='Unjournaled'):a.transfer([i])
    assert dst.read_bytes()==b'other'

def test_resume_does_not_skip_verification(monkeypatch,tmp_path):
    src,i=setup(monkeypatch,tmp_path);a.transfer([i]);dst=a.DEST/i['relative']
    dst.write_bytes(b'x'*dst.stat().st_size)
    records=a.transfer([i])
    with pytest.raises(RuntimeError,match='Checksum mismatch'):a.verify([i],records)
    assert src.exists()

def test_verified_migration_relocates_link_then_removes_only_source(monkeypatch,tmp_path):
    src,i=setup(monkeypatch,tmp_path);records=a.transfer([i]);verified=a.verify([i],records)
    cam=tmp_path/'campaign';cam.mkdir();source=tmp_path/'ssd';source.mkdir()
    monkeypatch.setattr(a,'CAM',cam);monkeypatch.setattr(a,'SOURCE',source)
    monkeypatch.setattr(a,'volumes',lambda:None);monkeypatch.setattr(a,'gate',lambda:None)
    (cam/'chains').symlink_to(source/'2026-10-13-xy-infinite-volume/chains')
    untouched=tmp_path/'unrelated';untouched.write_text('preserve')
    assert a.delete_verified([i],records,verified)==i['signature'][0]
    assert not src.exists()
    assert (a.DEST/i['relative']).read_bytes()
    assert untouched.read_text()=='preserve'
    assert str((cam/'chains').readlink()).startswith(str(a.DEST))

def test_changed_archive_blocks_ssd_deletion(monkeypatch,tmp_path):
    src,i=setup(monkeypatch,tmp_path);records=a.transfer([i]);verified=a.verify([i],records)
    monkeypatch.setattr(a,'volumes',lambda:None);monkeypatch.setattr(a,'gate',lambda:None)
    cam=tmp_path/'campaign';cam.mkdir();monkeypatch.setattr(a,'CAM',cam)
    (cam/'chains').symlink_to(a.SOURCE/'2026-10-13-xy-infinite-volume/chains')
    (a.DEST/i['relative']).write_bytes(b'bad')
    with pytest.raises(RuntimeError,match='verified archive changed'):a.delete_verified([i],records,verified)
    assert src.exists()
