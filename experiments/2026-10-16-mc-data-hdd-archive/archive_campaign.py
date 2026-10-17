"""Finite, fail-closed campaign migration. Never deletes before full SHA verification."""
import fcntl, hashlib, json, os, plistlib, shutil, stat, subprocess, time, traceback
from pathlib import Path
REPO=Path(__file__).resolve().parents[2]
JOB=Path(__file__).resolve().parent
SOURCE=Path('/Volumes/Office_SSD/Luscher_XY')
VOLUME=Path('/Volumes/MC_Data_Backup')
DEST=VOLUME/'Luscher_XY/archive-2026-10-16'
CAM=REPO/'experiments/2026-10-13-xy-infinite-volume'
STUDY=REPO/'experiments/2026-10-15-xy-3d-spin-dimension'
UUIDS={'/Volumes/Office_SSD':'53650819-819D-42D6-BB15-24760A94AE0B',str(VOLUME):'5E567601-B9EE-4158-BCFF-9EF43D269778'}
def atomic(path,data):
    tmp=path.with_suffix(path.suffix+'.tmp')
    with tmp.open('w') as f: json.dump(data,f,indent=2);f.flush();os.fsync(f.fileno())
    os.replace(tmp,path);syncdir(path.parent)
def syncdir(p):
    fd=os.open(p,os.O_RDONLY)
    try: os.fsync(fd)
    finally: os.close(fd)
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        while b:=f.read(8*1024*1024):h.update(b)
    return h.hexdigest()
def signature(p):
    s=p.lstat()
    if not stat.S_ISREG(s.st_mode):raise RuntimeError(f'Not a regular file: {p}')
    return [s.st_size,s.st_mtime_ns,s.st_ino,s.st_dev]
def volumes():
    infos={}
    for path,uid in UUIDS.items():
        info=plistlib.loads(subprocess.check_output(['/usr/sbin/diskutil','info','-plist',path]))
        if info.get('VolumeUUID')!=uid or info.get('MountPoint')!=path:raise RuntimeError('Volume identity changed')
        infos[path]=info
    return infos
def gate():
    c=json.loads((CAM/'control.json').read_text())
    if c.get('launch_new') or c.get('launch_allowed_ids'):raise RuntimeError('Simulation launch gate is open')
def event(data):
    with (JOB/'events.jsonl').open('a') as f: f.write(json.dumps({'time':time.time(),**data})+'\n');f.flush();os.fsync(f.fileno())
def status(phase,**kw):
    atomic(JOB/'status.json',{'phase':phase,'updated':time.time(),**kw})
def inventory():
    items=[]
    roots=[(SOURCE,'raw/Luscher_XY',True),(REPO/'experiments','records/experiments',False),
           (REPO/'supervillain','records/supervillain',False),(REPO/'test','records/test',False),
           (REPO/'docs','records/docs',False),(REPO/'output','records/output',False)]
    for src,prefix,remove in roots:
        for base,dirs,files in os.walk(src,followlinks=False):
            base=Path(base)
            dirs[:]=sorted(d for d in dirs if base/d!=JOB)
            links=[d for d in dirs if (base/d).is_symlink()]
            dirs[:]=[d for d in dirs if d not in links]
            for name in sorted(files+links):
                p=base/name; rel=str(Path(prefix)/p.relative_to(src))
                if p.is_symlink():
                    target=os.readlink(p)
                    # Archive must remain usable after the SSD copy is removed.
                    if p==CAM/'chains':target=str(DEST/'raw/Luscher_XY/2026-10-13-xy-infinite-volume/chains')
                    items.append(dict(source=str(p),relative=rel,link=target,original_link=os.readlink(p),remove=False))
                else:items.append(dict(source=str(p),relative=rel,signature=signature(p),remove=remove))
    for name in ['README.md','AGENTS.md','PYTHON_ENVIRONMENT.md','.python-version','requirements-xy.lock.txt','pyproject.toml','requirements.txt']:
        p=REPO/name
        if p.exists():items.append(dict(source=str(p),relative='records/'+name,signature=signature(p),remove=False))
    return items

def transfer(items):
    records={}
    ledger=JOB/'checksums.jsonl'
    if ledger.exists():
        for line in ledger.read_text().splitlines():
            r=json.loads(line);records[r['relative']]={'sha256':r['sha256'],'bytes':r['bytes']}
    total=sum(i.get('signature',[0])[0] for i in items); done=0
    for n,i in enumerate(items):
        src=Path(i['source']);dst=DEST/i['relative'];dst.parent.mkdir(parents=True,exist_ok=True)
        if 'link' in i:
            if dst.is_symlink():
                if os.readlink(dst)!=i['link']:raise RuntimeError('Destination link mismatch')
            elif dst.exists():raise RuntimeError('Destination collision')
            else:dst.symlink_to(i['link'])
            continue
        if signature(src)!=i['signature']:raise RuntimeError(f'Source changed: {src}')
        previous=records.get(i['relative'])
        if previous and dst.is_file() and dst.stat().st_size==i['signature'][0]:
            done+=dst.stat().st_size;continue # Full destination hash is checked below, including resumed copies.
        if dst.exists():raise RuntimeError(f'Unjournaled destination file: {dst}')
        tmp=dst.with_name(dst.name+'.archive-part');h=hashlib.sha256()
        with src.open('rb') as fi,tmp.open('wb') as fo:
            while b:=fi.read(8*1024*1024):fo.write(b);h.update(b)
            fo.flush();os.fsync(fo.fileno())
        if signature(src)!=i['signature']:raise RuntimeError(f'Source changed while copying: {src}')
        shutil.copystat(src,tmp);os.replace(tmp,dst);syncdir(dst.parent)
        records[i['relative']]={'sha256':h.hexdigest(),'bytes':i['signature'][0]}
        with ledger.open('a') as f:
            f.write(json.dumps({'relative':i['relative'],**records[i['relative']]})+'\n');f.flush();os.fsync(f.fileno())
        done+=i['signature'][0]
        if n%20==0 or i['signature'][0]>100_000_000:
            status('copying',files_done=n+1,files_total=len(items),bytes_done=done,bytes_total=total)
        print(f'Copied {n+1}/{len(items)} {done}/{total} bytes {i["relative"]}',flush=True)
    return records

def verify(items,records):
    verified={}
    for n,i in enumerate(items):
        p=DEST/i['relative']
        if 'link' in i:
            if not p.is_symlink() or os.readlink(p)!=i['link']:raise RuntimeError('Link verification failed')
        else:
            before=signature(p)
            if digest(p)!=records[i['relative']]['sha256'] or signature(p)!=before:raise RuntimeError(f'Checksum mismatch: {p}')
            verified[i['relative']]=before
        if n%20==0:status('verifying',files_done=n+1,files_total=len(items))
    atomic(JOB/'verified.json',verified)
    return verified

def delete_verified(items,records,verified):
    # Whole-copy verification has finished. Revalidate every remaining SSD source
    # BEFORE relocating links or beginning deletion. Never recurse outside inventory.
    removed=set()
    if (JOB/'events.jsonl').exists():
        for line in (JOB/'events.jsonl').read_text().splitlines():
            e=json.loads(line)
            if e.get('event')=='deleted':removed.add(e['relative'])
    raw=[i for i in items if i.get('remove')]
    for i in raw:
        p=Path(i['source'])
        if not p.exists() and i['relative'] in removed:continue
        if signature(p)!=i['signature'] or digest(p)!=records[i['relative']]['sha256']:raise RuntimeError(f'SSD source changed: {p}')
    volumes();gate()
    link=CAM/'chains';new=DEST/'raw/Luscher_XY/2026-10-13-xy-infinite-volume/chains'
    if not link.is_symlink() or os.readlink(link) not in [str(SOURCE/'2026-10-13-xy-infinite-volume/chains'),str(new)]:raise RuntimeError('Unexpected live chains link')
    tmp=link.with_name('chains.archive-new');tmp.unlink(missing_ok=True);tmp.symlink_to(new);os.replace(tmp,link);syncdir(link.parent)
    event({'event':'live_link_relocated','target':str(new)})
    for n,i in enumerate(raw):
        p=Path(i['source']);q=DEST/i['relative']
        if not p.exists() and i['relative'] in removed:continue
        volumes();gate()
        if signature(p)!=i['signature'] or signature(q)!=verified[i['relative']]:raise RuntimeError('Source or verified archive changed before deletion')
        # Persist intent; an interruption after unlink is detected safely on rerun.
        event({'event':'delete_intent','relative':i['relative']})
        p.unlink();syncdir(p.parent)
        event({'event':'deleted','relative':i['relative'],'bytes':i['signature'][0]})
        status('removing_verified_ssd_files',files_done=n+1,files_total=len(raw))
    return sum(i['signature'][0] for i in raw)

def main():
    locks=[]
    for path in [JOB/'archive.lock',STUDY/'fine-scan-manager.lock',STUDY/'common-xi-refinement-manager.lock',CAM/'retention-copy.lock',CAM/'storage-admission.lock']:
        f=path.open('a+');fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB);locks.append(f)
    gate();info=volumes()
    manifest=JOB/'manifest.json'
    items=json.loads(manifest.read_text()) if manifest.exists() else inventory()
    if not manifest.exists():atomic(manifest,items)
    total=sum(i.get('signature',[0])[0] for i in items)
    needed=sum(i.get('signature',[0])[0] for i in items if not (DEST/i['relative']).exists())
    available=min(shutil.disk_usage(VOLUME).free,info[str(VOLUME)]['APFSContainerFree'])
    if available<needed+25*1024**3:raise RuntimeError(f'Insufficient physical container free: {available}, need {needed}+reserve')
    DEST.mkdir(parents=True,exist_ok=True)
    # Preserve source code and Git identity (no ignored data assumed committed).
    if not (JOB/'git-status.txt').exists():
        (JOB/'git-status.txt').write_bytes(subprocess.check_output(['git','status','--short'],cwd=REPO))
        (JOB/'git-head.txt').write_bytes(subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO))
    records=transfer(items)
    verified=verify(items,records)
    atomic(JOB/'all-files-verified.json',{'time':time.time(),'files':len(items),'bytes':total,'destination':str(DEST)})
    shutil.copytree(JOB,DEST/'transfer-records',dirs_exist_ok=True,ignore=shutil.ignore_patterns('*.lock','*.log','*.plist','status.json'))
    freed=delete_verified(items,records,verified)
    status('complete',verified_files=len(items),verified_bytes=total,ssd_bytes_removed=freed,destination=str(DEST))
    shutil.copytree(JOB,DEST/'transfer-records',dirs_exist_ok=True,ignore=shutil.ignore_patterns('*.lock','*.log','*.plist'))
    print('COMPLETE: all files verified, live link relocated, inventoried SSD copies removed.',flush=True)
if __name__=='__main__':
    try:main()
    except BaseException as exc:
        status('failed_safe',error=repr(exc));traceback.print_exc();raise
