"""Descriptive independent long-start window checks; no mixing certification."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
import json
from pathlib import Path
import numpy as np
import run
from analyze_calibration import read,summarize
from scout import json_finite
OUT=Path(__file__).resolve().parent
rows=[]
for D,q,N,beta in [(3,92,16,293),(4,98,10,261)]:
    points=[json.loads(f.read_text())['point'] for f in (OUT/'jobs').glob('*.json') if (lambda j:j['state']=='complete' and j['audit']['passed'] and (j['point']['D'],round(j['point']['epsilon_over_pi']*100),j['point']['N'],round(j['point']['beta']*1000),j['point']['replica'] in [2,3])==(D,q,N,beta,True))(json.loads(f.read_text()))]
    for p in sorted(points,key=lambda p:p['replica']):
        a=read(p)
        for first,stop in [(1001,5001),(5001,9001),(9001,13001),(13001,16001),(8001,16001)]:
            summary,boot=summarize([a[first:stop]],[p],0)
            rows.append(dict(id=p['id'],D=D,q=q,N=N,beta=beta/1000,replica=p['replica'],first=first,stop=stop,summary=summary))
folder=OUT/('long-start-window-check-'+run.now()[:19].replace('-','').replace(':','').replace('T','-'));folder.mkdir()
run.atomic_json(folder/'analysis.json',json_finite(dict(created=run.now(),rows=rows,source_sha256=run.sha256(Path(__file__)),scope='Per-start block-bootstrap windows; last-half overlaps earlier windows. Descriptive stationarity check only, not equilibration certification. Raw data retained.')))
print(folder,flush=True)
for r in rows: print(r['D'],r['replica'],r['first'],r['stop'],r['summary'].get('values',[None])[0],r['summary'].get('bootstrap_95'),r['summary'].get('flags'),flush=True)
