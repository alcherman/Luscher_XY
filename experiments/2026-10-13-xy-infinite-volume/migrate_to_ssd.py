"""Verified raw-data copy; preserve the internal originals and frozen runner."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
from time import perf_counter

OUT = Path(__file__).resolve().parent
SSD = Path('/Volumes/Office_SSD')
TARGET = SSD / 'Luscher_XY/2026-10-13-xy-infinite-volume/chains'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(4*1024**2), b''):
            h.update(block)
    return h.hexdigest()


def atomic(path, value):
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, delete=False) as f:
        json.dump(value, f, indent=2); f.write('\n'); f.flush(); os.fsync(f.fileno())
        tmp=Path(f.name)
    tmp.replace(path)


def main():
    assert SSD.is_mount(), 'External SSD is not mounted'
    state=json.loads((OUT/'status.json').read_text())
    assert state['active']==0 and state['paused']==0 and state['pending']==0
    control=json.loads((OUT/'control.json').read_text())
    assert not control['launch_new'], 'Pause new launches before migration'
    jobs=[json.loads(p.read_text()) for p in (OUT/'jobs').glob('*.json')]
    assert len(jobs)==36 and all(j['state']=='complete' and j['audit']['passed'] for j in jobs)
    source=OUT/'chains'
    assert not source.is_symlink(), 'Migration already performed; inspect storage-routing.json'
    TARGET.mkdir(parents=True,exist_ok=True)
    started=perf_counter(); verified=[]
    for job in sorted(jobs,key=lambda j:j['point']['id']):
        name=job['point']['id']+'.h5'; src=source/name; dst=TARGET/name
        expected=job['audit']
        assert src.stat().st_size==expected['bytes']
        if not dst.exists():
            partial=dst.with_suffix('.h5.copying')
            if partial.exists():
                raise RuntimeError('Preserve and inspect previous partial copy: '+str(partial))
            shutil.copyfile(src,partial)
            assert partial.stat().st_size==expected['bytes'] and digest(partial)==expected['sha256']
            partial.replace(dst)
        else:
            assert dst.stat().st_size==expected['bytes'] and digest(dst)==expected['sha256']
        verified.append(dict(name=name,bytes=expected['bytes'],sha256=expected['sha256']))
        print('verified',name,expected['bytes'],flush=True)
    backup=OUT/'chains-internal-copy-2026-10-13'
    assert not backup.exists()
    link=OUT/'chains.external-link.tmp'
    assert not link.exists()
    link.symlink_to(TARGET,target_is_directory=True)
    source.rename(backup)
    link.replace(source)
    record=dict(volume=str(SSD),raw_directory=str(TARGET),internal_copy=str(backup),
                verified_files=verified,bytes=sum(r['bytes'] for r in verified),
                elapsed_seconds=perf_counter()-started,
                originals_preserved=True,source_runner_unchanged=True)
    atomic(OUT/'storage-routing.json',record)
    atomic(TARGET.parent/'migration-verification.json',record)
    print('migration_complete',len(verified),record['bytes'],record['elapsed_seconds'],flush=True)


if __name__=='__main__':main()
