"""Build an unapproved attained delivery after the existing sensitivity service."""
import fcntl,json,os,subprocess,sys,time
from pathlib import Path
OUT=Path(__file__).resolve().parent
# Lock covers both waiting and generation; app restart cannot duplicate this task.
with (OUT/'d5-delivery-builder.lock').open('a+') as lock:
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:
        print('Another delivery builder is active; no duplicate started',flush=True);sys.exit(0)
    while True:
        status=json.loads((OUT/'d5-final-sensitivity-manager-status.json').read_text())
        if status['state']=='snapshots_complete_pending_scientific_review':break
        if status['state'] in ['failed','incomplete_snapshots']:raise RuntimeError(f'Analysis did not close: {status}')
        time.sleep(30)
    print('Sensitivity snapshots ready; creating pending-review delivery',flush=True)
    result=subprocess.run([sys.executable,str(OUT/'build_attained_d5_delivery.py')],cwd=OUT.parents[1])
    if result.returncode:raise RuntimeError(f'Delivery build failed: {result.returncode}')
    print('Draft complete; independent scientific and visual review is still required',flush=True)
