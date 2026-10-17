"""Finite controller for already-running, source-frozen worker pool."""
from pathlib import Path
import subprocess,json,os,signal,time,fcntl,datetime
P=Path(__file__).resolve().parent;W=P/'boundary-resolution';prefix='/Users/acherman/research/Luscher_XY/.venv/bin/python -u '+str(W/'runner.py')+' worker '
def snapshot():
 result=[]
 for line in subprocess.check_output(['ps','-axo','pid=,stat=,command='],text=True).splitlines():
  parts=line.strip().split(None,2)
  if len(parts)==3 and parts[2].startswith(prefix):result.append((int(parts[0]),parts[1]))
 return sorted(result)
def record(rows,active,stopped):
 x=dict(updated=datetime.datetime.now().astimezone().isoformat(),cap=4,workers=[p for p,_ in rows],running=active,suspended=stopped,scope='Own reviewed XY workers only, temporary SIGSTOP/SIGCONT concurrency control; other calculations untouched');temp=P/'four-core-controller-status.tmp';temp.write_text(json.dumps(x,indent=2)+'\n');temp.replace(P/'four-core-controller-status.json')
if __name__=='__main__':
 with (P/'four-core-controller.lock').open('a+') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  while True:
   rows=snapshot()
   if not rows:record([],[],[]);break
   wanted={p for p,_ in rows[:4]};active=[];stopped=[]
   for pid,state in rows:
    if pid in wanted:
     if 'T' in state:os.kill(pid,signal.SIGCONT)
     active.append(pid)
    else:
     if 'T' not in state:os.kill(pid,signal.SIGSTOP)
     stopped.append(pid)
   record(rows,active,stopped);time.sleep(2)
