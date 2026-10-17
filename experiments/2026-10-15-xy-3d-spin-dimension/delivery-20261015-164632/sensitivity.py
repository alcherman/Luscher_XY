"""Run a finite, source-preserving set of MC analysis sensitivity snapshots."""
import subprocess,sys,json,concurrent.futures
from pathlib import Path
OUT=Path(__file__).resolve().parent
variants={'full':['--window','full'],'U20':['--target','.20'],'U30':['--target','.30'],'cut501':['--cut','501'],'cut2001':['--cut','2001'],'block2':['--factor','2'],'even':['--starts','even'],'odd':['--starts','odd']}
def one(item):
 tag,args=item
 with (OUT/(tag+'.log')).open('w') as log:r=subprocess.run([sys.executable,str(OUT/'measure.py'),'--tag',tag]+args,stdout=log,stderr=subprocess.STDOUT)
 print(tag,r.returncode,flush=True);return tag,r.returncode
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
 results=dict(pool.map(one,variants.items()))
(OUT/'sensitivity-status.json').write_text(json.dumps(results,indent=2)+'\n')
if any(results.values()):raise RuntimeError(results)
