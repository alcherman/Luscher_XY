"""Finite existing-data CV test, including prespecified sampling variants."""
from pathlib import Path
import subprocess,sys
P=Path(__file__).resolve().parent
for tag in ['primary','cut4000','block4']:
    out=P/f'spin-control-variate-20261016-v2-{tag}'
    subprocess.run([sys.executable,'-u',str(P/'spin_control_variate.py'),'--variant',tag,'--output',str(out)],check=True)
    for estimator in ['original','CV']:
        subprocess.run([sys.executable,'-u',str(P/'report_spin_correlations.py'),'--summary',str(out/estimator),'--output',str(out/(estimator+'-quotients'))],check=True)
