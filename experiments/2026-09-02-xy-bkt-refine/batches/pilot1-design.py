"""Freeze adaptive batch designs before generation; pilot targets select sampling only."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from run import OUT, ROOT, point, atomic_json, sha256, now, _provenance

SCOUT=OUT.parent/'2026-09-02-xy-bkt-scout'

def freeze(name,points,selection):
    (OUT/'batches').mkdir(exist_ok=True)
    value=dict(schema='xy-refinement.v1',created_utc=now(),batch=name,D=2,points=points,
        model_source_sha256=_provenance()['source_sha256'],runner_sha256=sha256(OUT/'run.py'),
        design_source_sha256=sha256(Path(__file__)),selection=selection,
        update='1 full Metropolis sweep + 4 fixed Wolff clusters + 1 winding proposal per cycle',
        retention='Every generated field, observable, RNG state, work count and timing; initial state included.',
        interval_phi=float(np.pi),checkpoint_every=100,buffer_rows=100,backend='serial Numba; four independent CPU workers')
    path=OUT/'batches'/f'{name}.json'
    if path.exists():raise FileExistsError(path)
    atomic_json(path,value)
    print(name,len(points),'chains',sum(p['steps'] for p in points),'cycles',flush=True)


def pilot():
    rows=json.loads((SCOUT/'results.json').read_text());ps=[];selection=[]
    for e in [580]+list(range(600,1001,20)):
        group=sorted([r for r in rows if r['point']['epsilon_milli']==e and r['point']['L']==64],key=lambda r:r['point']['kappa'])
        ks=np.array([r['point']['kappa'] for r in group]);rs=np.array([r['summary']['values']['xi_over_L'] or 0 for r in group])
        interp=PchipInterpolator(ks,np.maximum.accumulate(rs));center=brentq(lambda k:interp(k)-.765,0,1.2)
        k=int(np.clip(20*round(center*1000/20),20,1180))
        ps.append(point(e,k,64,100,10000));selection.append(dict(epsilon_milli=e,kappa_milli=k,initial_R_target=.765))
    freeze('pilot1',ps,dict(method='PCHIP interpolation of old L64 R to .765 selects sampling neighborhoods only. No transition is inferred from this target.',choices=selection,scout_results_sha256=sha256(SCOUT/'results.json')))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['pilot']);args=parser.parse_args()
    pilot()
