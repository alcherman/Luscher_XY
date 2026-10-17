"""Finite existing-data analysis; waits for the one frozen field-reading job."""
from pathlib import Path
import json,time,subprocess,sys
P=Path(__file__).resolve().parent
status=P/'spin-correlation-analysis-manager-status-20261016.json'
def write(x):
    tmp=status.with_suffix('.tmp');tmp.write_text(json.dumps(x,indent=2));tmp.replace(status)
try:
    while True:
        s=json.loads((P/'spin-correlation-grid-20261016/status.json').read_text())
        if s['state']=='measurements_complete_pending_statistical_audit':break
        write(dict(state='waiting_for_measurements',measurement_status=s));time.sleep(5)
    for tag,options in [('primary',[]),('checkpoint',['--scalar-mode','checkpoint']),('block4',['--factor','4']),('cut4000',['--cut','4000'])]:
        write(dict(state='analyzing',variant=tag))
        summary=P/f'spin-correlation-grid-summary-20261016-{tag}';report=P/f'spin-correlation-delivery-20261016-{tag}'
        subprocess.run([sys.executable,'-u',str(P/'analyze_spin_correlations.py'),'--inputs',str(P/'spin-correlation-pilot-20261016'),str(P/'spin-correlation-grid-20261016'),'--output',str(summary),*options],check=True)
        subprocess.run([sys.executable,'-u',str(P/'report_spin_correlations.py'),'--summary',str(summary),'--output',str(report)],check=True)
    write(dict(state='analyses_complete_pending_independent_review',variants=['primary','checkpoint','block4','cut4000']))
except BaseException as e:
    write(dict(state='failed',error=repr(e)));raise
