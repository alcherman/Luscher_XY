"""Both producers count the outstanding temporary copy in shared admission."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments/2026-10-13-xy-infinite-volume'))
import storage_admission as admission

def contender(root,queue):
    admission.OUT=Path(root);queue.put('attempting')
    lock=admission.acquire()
    try:queue.put('acquired')
    finally:lock.close()

def test_admission_serializes_independent_processes(tmp_path,monkeypatch):
    import multiprocessing,queue,pytest
    monkeypatch.setattr(admission,'OUT',tmp_path)
    context=multiprocessing.get_context('spawn');messages=context.Queue();lock=admission.acquire()
    process=context.Process(target=contender,args=(str(tmp_path),messages));process.start()
    try:
        assert messages.get(timeout=10)=='attempting'
        with pytest.raises(queue.Empty):messages.get(timeout=.2)
    finally:lock.close()
    assert messages.get(timeout=10)=='acquired';process.join(timeout=10);assert process.exitcode==0

def test_remaining_reservation_and_owner_release(tmp_path,monkeypatch):
    monkeypatch.setattr(admission,'OUT',tmp_path);temp=tmp_path/'copy.tmp'
    lock=admission.acquire()
    try:
        admission.reserve('a',temp,1000);assert admission.reserved_remaining()==1000
        temp.write_bytes(b'x'*300);assert admission.reserved_remaining()==700
        admission.release('b');assert admission.reserved_remaining()==700
        admission.release('a');assert admission.reserved_remaining()==0
    finally:lock.close()


def test_physical_chunk_budget_includes_N10_padding():
    import storage_admission as budget
    spec=dict(D=5,N=10,steps=8000)
    assert budget.estimated_full_bytes(spec)>16612131848
    # Independent observed producer chunk geometry for phi and winding.
    from math import prod,ceil
    phi=prod((128,1,2,3,3,3,3))*8*prod(ceil(a/b) for a,b in zip((8001,1,10,10,10,10,10),(128,1,2,3,3,3,3)))
    winding=prod((32,1,313))*8*prod(ceil(a/b) for a,b in zip((8001,5,10000),(32,1,313)))
    assert budget.estimated_full_bytes(spec)>=1.15*(phi+winding)


def test_chunk_budget_matches_actual_HDF5_chunk_selector_and_allocation(tmp_path):
    import storage_admission as budget
    from h5py._hl.filters import guess_chunk
    import h5py,numpy as np
    from math import prod,ceil
    for tail in [(1,5,5,5),(3,25)]:
        with h5py.File(tmp_path/'layout.h5','w') as f:
            d=f.create_dataset('x',shape=(0,)+tail,maxshape=(None,)+tail,dtype='f8',chunks=True)
            assert d.chunks==guess_chunk((0,)+tail,(None,)+tail,8)
            d.resize(17,axis=0);d[:]=np.ones((17,)+tail)
            expected=8*prod(d.chunks)*prod(ceil(n/c) for n,c in zip(d.shape,d.chunks))
            assert d.id.get_storage_size()==expected
    assert budget.estimated_full_bytes(dict(D=3,N=5,steps=16))>=expected
