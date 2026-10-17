"""Shared admission lock and outstanding retention temporary-copy reservation."""
import fcntl,json,os
from pathlib import Path
import run
_legacy_estimated_bytes=run.estimated_bytes
OUT=Path(__file__).resolve().parent

def acquire():
    lock=(OUT/'storage-admission.lock').open('a+')
    fcntl.flock(lock,fcntl.LOCK_EX)
    return lock

def reservation_paths():
    legacy=OUT/'retention-space-reservation.json'
    return ([legacy] if legacy.exists() else [])+list((OUT/'retention-space-reservations').glob('*.json'))

def reserved_remaining():
    total=0
    for path in reservation_paths():
        record=json.loads(path.read_text());temp=Path(record['temporary_path'])
        # A peer may atomically rename its completed temp onto the archive even
        # while admission is locked. A missing temp conservatively reserves the
        # full estimate until its owner releases the lease. Never exists()+stat().
        try:written=temp.stat().st_size
        except FileNotFoundError:written=0
        total+=max(0,record['estimated_total_bytes']-written)
    return total

def reserve(identifier,temp,estimate):
    assert Path(identifier).name==identifier and identifier not in ('.','..')
    directory=OUT/'retention-space-reservations';directory.mkdir(exist_ok=True)
    run.atomic_json(directory/f'{identifier}.json',dict(id=identifier,temporary_path=str(temp),estimated_total_bytes=estimate,created=run.now(),pid=os.getpid()))
    fd=os.open(directory,os.O_RDONLY)
    try:os.fsync(fd)
    finally:os.close(fd)
    sync_directory()

def sync_directory():
    fd=os.open(OUT,os.O_RDONLY);os.fsync(fd);os.close(fd)

def release(identifier):
    for path in reservation_paths():
        if json.loads(path.read_text())['id']==identifier:
            path.unlink()
            fd=os.open(path.parent,os.O_RDONLY)
            try:os.fsync(fd)
            finally:os.close(fd)

def clear_retention_reservations():
    """Caller holds EX retention-copy maintenance lock AND admission lock."""
    for path in reservation_paths():
        path.unlink()
        fd=os.open(path.parent,os.O_RDONLY)
        try:os.fsync(fd)
        finally:os.close(fd)

def admit_retention_copy(identifier,temp,estimate,active_remaining=0,reserve_bytes=25*2**30):
    """Caller holds admission lock; include peers before reading fresh free."""
    peers=reserved_remaining()
    free=run.shutil.disk_usage(Path(temp).parent).free
    if free-active_remaining-peers-estimate<reserve_bytes:return False
    reserve(identifier,temp,estimate)
    return True


def estimated_full_bytes(spec):
    """Reserve physical HDF5 chunks, including padding of non-power-of-two boxes.

    Matches the frozen producer's chunks=True on initial zero-length datasets.
    No compression or future retention savings are assumed. The environment's
    h5py chunk selector is verified against actual dataset creation in tests.
    """
    from math import prod,ceil
    from h5py._hl.filters import guess_chunk
    D,N,rows=spec['D'],spec['N'],spec['steps']+1
    def allocation(tail):
        shape=(0,)+tuple(tail);chunks=guess_chunk(shape,(None,)+tuple(tail),8)
        final=(rows,)+tuple(tail)
        return prod(chunks)*8*prod(ceil(n/c) for n,c in zip(final,chunks))
    fields=allocation((1,)+tuple([N]*D))
    winding=allocation((D,N**(D-1)))
    estimate=int(1.15*(fields+winding+4096*rows)+16*2**20)
    return max(_legacy_estimated_bytes(spec),estimate)
