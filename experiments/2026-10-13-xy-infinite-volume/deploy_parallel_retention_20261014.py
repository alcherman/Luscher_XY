"""Safely drain the serial converter before replacing its launchd service."""
import fcntl,json,os,plistlib,re,subprocess,time
from pathlib import Path
OUT=Path(__file__).resolve().parent
LABEL='gui/501/org.luscherxy.adaptive-retention-20261014'
control=json.loads((OUT/'control.json').read_text());status=json.loads((OUT/'status.json').read_text())
assert not control['launch_new'] and status['active']==0 and status['pending']==0
text=subprocess.check_output(['launchctl','print',LABEL],text=True);old=int(re.search(r'\bpid = (\d+)',text).group(1))
print('Waiting for existing copy to finish, pid',old,flush=True)
with (OUT/'retention-copy.lock').open('a+') as maintenance:
    fcntl.flock(maintenance,fcntl.LOCK_EX)
    print('Exclusive maintenance lock acquired; no conversion is active',flush=True)
    subprocess.run(['launchctl','bootout',LABEL],check=True)
    deadline=time.monotonic()+30
    while True:
        try:os.kill(old,0)
        except ProcessLookupError:break
        if time.monotonic()>deadline:raise RuntimeError('Old daemon did not exit; no replacement started')
        time.sleep(.2)
    print('Old daemon exited; originals and committed copies preserved',flush=True)
    # Completed simulation coordinator is no longer needed. This also prevents
    # its old single-lease module instance surviving the storage-helper upgrade.
    stopped=subprocess.run(['launchctl','bootout','gui/501/org.luscherxy.phase-boundary-20261013'],capture_output=True,text=True)
    print('Idle simulation coordinator stop:',stopped.returncode,stopped.stderr,flush=True)
    path=OUT/'retention.plist';data=plistlib.loads(path.read_bytes())
    if '--workers' not in data['ProgramArguments']:data['ProgramArguments']+=['--workers','14']
    else:data['ProgramArguments'][data['ProgramArguments'].index('--workers')+1]='14'
    path.write_bytes(plistlib.dumps(data))
    subprocess.run(['launchctl','bootstrap','gui/501',str(path)],check=True)
    print('Replacement service registered with 14 single-thread workers; releasing maintenance lock',flush=True)
print('Safe deployment complete',flush=True)
