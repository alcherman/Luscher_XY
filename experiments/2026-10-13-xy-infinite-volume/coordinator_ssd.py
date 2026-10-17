"""Storage-aware coordinator wrapper; leave the calibration runner immutable."""
from pathlib import Path
import sys,json
import run
import storage_admission

SSD=Path('/Volumes/Office_SSD')
TARGET=SSD/'Luscher_XY/2026-10-13-xy-infinite-volume/chains'
original_disk_usage=run.shutil.disk_usage
original_atomic_json=run.atomic_json
original_estimated_bytes=storage_admission.estimated_full_bytes
coordinate_mode=False
admission_lock=None
reserved_active_ids=set()
budget_snapshot={}


def guarded_usage(path):
    global admission_lock,reserved_active_ids,budget_snapshot
    if coordinate_mode:
        admission_lock=storage_admission.acquire()
    if not SSD.is_mount() or not TARGET.is_dir():
        raise RuntimeError('Campaign SSD is not mounted; refusing new launches')
    if (run.OUT/'chains').resolve()!=TARGET:
        raise RuntimeError('Campaign raw-data link does not match the verified SSD route')
    # Read outstanding growth BEFORE physical free space: intervening writes
    # then reduce free as well, which is conservative rather than phantom space.
    reserved=storage_admission.reserved_remaining()
    active_reserved=0;reserved_active_ids=set()
    if coordinate_mode:
        for job_path in (run.OUT/'jobs').glob('*.json'):
            job=json.loads(job_path.read_text())
            if job['state']=='running':
                point=job['point'];chain=run.OUT/'chains'/f"{point['id']}.h5"
                active_reserved+=max(0,original_estimated_bytes(point)-(chain.stat().st_size if chain.exists() else 0))
                reserved_active_ids.add(point['id'])
    physical=original_disk_usage(TARGET)
    budget_snapshot=dict(physical_disk_free_bytes=physical.free,retention_reserved_bytes=reserved,active_reserved_bytes=active_reserved)
    withheld=reserved+active_reserved
    if withheld==0:return physical
    return type(physical)(physical.total,physical.used+withheld,max(0,physical.free-withheld))


def scheduling_estimate(point):
    # These active jobs were already reserved before the fresh free-space read.
    # Make the frozen scheduler's subsequent remaining-live subtraction zero.
    if point['id'] in reserved_active_ids:
        chain=run.OUT/'chains'/f"{point['id']}.h5"
        return chain.stat().st_size if chain.exists() else 0
    control=run.OUT/'control.json'
    config=json.loads(control.read_text()) if control.exists() else {}
    allowed=config.get('launch_allowed_ids')
    if allowed is not None and point['id'] not in allowed:return 2**120
    cap=config.get('d2_max_N')
    if point.get('D')==2 and cap is not None and point['N']>cap:
        # Frozen scheduler rejects requirements exceeding available storage.
        # This admission sentinel prevents new launches above the user size cap.
        return 2**120
    return original_estimated_bytes(point)


def publish_status(path,value):
    global admission_lock
    try:
        if Path(path)==run.OUT/'status.json':value=dict(value,**budget_snapshot)
        return original_atomic_json(path,value)
    finally:
        if Path(path)==run.OUT/'status.json' and admission_lock is not None:
            admission_lock.close();admission_lock=None


if __name__=='__main__':
    if '--check-only' in sys.argv:
        print('ssd_free_bytes',guarded_usage(run.OUT).free)
    else:
        coordinate_mode=True
        run.atomic_json=publish_status
        run.estimated_bytes=scheduling_estimate
        run.shutil.disk_usage=guarded_usage
        run.coordinate()
