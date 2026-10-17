"""Finite detached discovery batch, then source-preserving analysis and plots."""
from pathlib import Path
import os,sys,subprocess,fcntl,json,traceback
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P));import runner as r
if __name__=='__main__':
 manifest=Path(sys.argv[1]).resolve()
 with (P/'supervisor.lock').open('a+') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  try:
   r.atomic(P/'supervisor-status.json',dict(state='running',pid=os.getpid(),manifest=str(manifest),updated=r.pt.legacy.now()))
   for args in [['runner.py','run',str(manifest)],['analyse.py',str(manifest)],['plot.py']]:subprocess.run([sys.executable,'-u',str(P/args[0]),*args[1:]],cwd=r.ROOT,check=True)
   r.atomic(P/'supervisor-status.json',dict(state='discovery_complete_pending_scientific_review',manifest=str(manifest),updated=r.pt.legacy.now()))
  except BaseException:
   r.atomic(P/'supervisor-status.json',dict(state='failed_needs_attention',manifest=str(manifest),error=traceback.format_exc(),updated=r.pt.legacy.now()));raise
