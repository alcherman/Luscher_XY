"""Ensure the detached coordinator checks the intended raw-data filesystem."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import sys

import pytest

HERE=Path(__file__).resolve().parents[1]/'experiments/2026-10-13-xy-infinite-volume'
sys.path.insert(0,str(HERE))
SPEC=importlib.util.spec_from_file_location('xy_ssd_coordinator',HERE/'coordinator_ssd.py')
route=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(route)


def test_ssd_usage_never_falls_back_to_internal_disk(tmp_path,monkeypatch):
    volume=tmp_path/'Office_SSD'; raw=volume/'campaign/chains'; raw.mkdir(parents=True)
    internal=tmp_path/'repo'; internal.mkdir();(internal/'chains').symlink_to(raw,target_is_directory=True)
    monkeypatch.setattr(route,'SSD',volume);monkeypatch.setattr(route,'TARGET',raw)
    monkeypatch.setattr(route.run,'OUT',internal)
    monkeypatch.setattr(route.storage_admission,'OUT',internal)
    monkeypatch.setattr(Path,'is_mount',lambda p:p==volume)
    seen=[]
    def usage(path):
        seen.append(path);return SimpleNamespace(free=4*1024**4)
    monkeypatch.setattr(route,'original_disk_usage',usage)
    assert route.guarded_usage(internal).free==4*1024**4
    assert seen==[raw]
    monkeypatch.setattr(Path,'is_mount',lambda p:False)
    with pytest.raises(RuntimeError,match='not mounted'):
        route.guarded_usage(internal)
    assert seen==[raw]


def test_mismatched_raw_path_blocks_launch(tmp_path,monkeypatch):
    volume=tmp_path/'Office_SSD';raw=volume/'campaign/chains';raw.mkdir(parents=True)
    internal=tmp_path/'repo';(internal/'chains').mkdir(parents=True)
    monkeypatch.setattr(route,'SSD',volume);monkeypatch.setattr(route,'TARGET',raw)
    monkeypatch.setattr(route.run,'OUT',internal)
    monkeypatch.setattr(route.storage_admission,'OUT',internal)
    monkeypatch.setattr(Path,'is_mount',lambda p:p==volume)
    with pytest.raises(RuntimeError,match='does not match'):
        route.guarded_usage(internal)


def test_scheduler_counts_unwritten_temporary_reservation(tmp_path,monkeypatch):
    import shutil
    volume=tmp_path/'Office_SSD';raw=volume/'campaign/chains';raw.mkdir(parents=True)
    internal=tmp_path/'repo';internal.mkdir();(internal/'chains').symlink_to(raw,target_is_directory=True)
    monkeypatch.setattr(route,'SSD',volume);monkeypatch.setattr(route,'TARGET',raw);monkeypatch.setattr(route.run,'OUT',internal)
    monkeypatch.setattr(route.storage_admission,'OUT',internal);monkeypatch.setattr(Path,'is_mount',lambda p:p==volume)
    monkeypatch.setattr(route,'original_disk_usage',lambda p:shutil._ntuple_diskusage(20000,10000,10000))
    temp=raw/'copy.tmp';temp.write_bytes(b'x'*300)
    route.storage_admission.reserve('a',temp,1000)
    monkeypatch.setattr(route,'coordinate_mode',True)
    assert route.guarded_usage(internal).free==9300
    assert route.admission_lock is not None
    route.publish_status(internal/'status.json',dict(test=True))
    assert route.admission_lock is None


def test_active_growth_during_snapshot_is_conservative(tmp_path,monkeypatch):
    import shutil,json
    volume=tmp_path/'Office_SSD';raw=volume/'campaign/chains';raw.mkdir(parents=True)
    internal=tmp_path/'repo';internal.mkdir();(internal/'chains').symlink_to(raw,target_is_directory=True);(internal/'jobs').mkdir()
    monkeypatch.setattr(route,'SSD',volume);monkeypatch.setattr(route,'TARGET',raw);monkeypatch.setattr(route.run,'OUT',internal)
    monkeypatch.setattr(route.storage_admission,'OUT',internal);monkeypatch.setattr(Path,'is_mount',lambda p:p==volume)
    chain=raw/'a.h5';chain.write_bytes(b'x'*300);temp=raw/'copy.tmp';temp.write_bytes(b'x'*300)
    (internal/'jobs'/'a.json').write_text(json.dumps(dict(state='running',point=dict(id='a'))))
    monkeypatch.setattr(route,'original_estimated_bytes',lambda p:1000)
    def usage(path):
        chain.write_bytes(b'x'*400);temp.write_bytes(b'x'*400)
        return shutil._ntuple_diskusage(20000,10200,9800)
    monkeypatch.setattr(route,'original_disk_usage',usage);monkeypatch.setattr(route,'coordinate_mode',True)
    route.storage_admission.reserve('copy',temp,1000)
    assert route.guarded_usage(internal).free==8400  # conservative; fresh exact margin is8600
    assert route.scheduling_estimate(dict(id='a'))==400
    assert route.scheduling_estimate(dict(id='new'))==1000
    route.publish_status(internal/'status.json',dict(test=True))

def test_user_d2_size_cap_blocks_pending_admission(tmp_path,monkeypatch):
    import json
    monkeypatch.setattr(route.run,'OUT',tmp_path)
    monkeypatch.setattr(route,'reserved_active_ids',set())
    monkeypatch.setattr(route,'original_estimated_bytes',lambda p:123)
    (tmp_path/'control.json').write_text(json.dumps(dict(d2_max_N=512)))
    assert route.scheduling_estimate(dict(id='new-large',D=2,N=1024))>2**100
    assert route.scheduling_estimate(dict(id='allowed',D=2,N=512))==123
    assert route.scheduling_estimate(dict(id='d5',D=5,N=8))==123


def test_finite_pilot_allowlist_blocks_other_pending_jobs(tmp_path,monkeypatch):
    import json
    monkeypatch.setattr(route.run,'OUT',tmp_path)
    monkeypatch.setattr(route,'reserved_active_ids',set())
    monkeypatch.setattr(route,'original_estimated_bytes',lambda p:123)
    (tmp_path/'control.json').write_text(json.dumps(dict(launch_allowed_ids=['pilot'])))
    assert route.scheduling_estimate(dict(id='pilot',D=3,N=64))==123
    assert route.scheduling_estimate(dict(id='unrelated',D=3,N=32))>2**100
