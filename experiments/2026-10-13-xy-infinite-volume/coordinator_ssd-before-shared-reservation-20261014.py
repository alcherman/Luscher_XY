"""Storage-aware coordinator wrapper; leave the calibration runner immutable."""
from pathlib import Path
import sys
import run

SSD=Path('/Volumes/Office_SSD')
TARGET=SSD/'Luscher_XY/2026-10-13-xy-infinite-volume/chains'
original_disk_usage=run.shutil.disk_usage


def guarded_usage(path):
    if not SSD.is_mount() or not TARGET.is_dir():
        raise RuntimeError('Campaign SSD is not mounted; refusing new launches')
    if (run.OUT/'chains').resolve()!=TARGET:
        raise RuntimeError('Campaign raw-data link does not match the verified SSD route')
    return original_disk_usage(TARGET)


if __name__=='__main__':
    if '--check-only' in sys.argv:
        print('ssd_free_bytes',guarded_usage(run.OUT).free)
    else:
        run.shutil.disk_usage=guarded_usage
        run.coordinate()
