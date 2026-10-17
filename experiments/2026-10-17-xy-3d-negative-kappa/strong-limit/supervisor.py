"""Finite strong-limit closure; original reconnaissance snapshots preserved."""
from pathlib import Path
import subprocess,sys,json,time,shutil,fcntl,importlib.util
P=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('r',P/'runner.py');r=importlib.util.module_from_spec(s);s.loader.exec_module(r)
if __name__=='__main__':
 with (P/'strong-limit-supervisor.lock').open('a+') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  try:
   r.atomic(P/'strong-limit-supervisor-status.json',dict(state='running',manifest=sys.argv[1],updated=r.pt.legacy.now()))
   subprocess.run([sys.executable,'-u',str(P/'runner.py'),'run',sys.argv[1]],check=True,cwd=r.ROOT)
   # analyse.py itself preserves the previous primary JSON in analysis-history.
   subprocess.run([sys.executable,'-u',str(P/'analyse.py'),sys.argv[1]],check=True,cwd=r.ROOT)
   shutil.copyfile(P/'analysis.json',P/'strong-limit-analysis.json')
   r.atomic(P/'strong-limit-supervisor-status.json',dict(state='complete_pending_scientific_review',updated=r.pt.legacy.now()))
  except BaseException as e:
   r.atomic(P/'strong-limit-supervisor-status.json',dict(state='failed_needs_attention',error=repr(e),updated=r.pt.legacy.now()));raise
