"""Explicitly post-hoc existing-data diagnostics; not independent validation."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import importlib.util,json,time,shutil
import numpy as np
OUT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('frozen_joint',OUT/'joint-fss-existing-20261016-full/joint_fss_existing.py');j=importlib.util.module_from_spec(spec);spec.loader.exec_module(j);j.OUT=OUT

def main():
    folder=OUT/'joint-fss-posthoc-diagnostics-20261016';folder.mkdir(exist_ok=False)
    a=json.loads((OUT/'joint-fss-existing-20261016-full/sources.json').read_text());rows,B,source=j.load_bank(a['analysis'])
    if source['sha256']!=a['sha256']:raise RuntimeError('Source changed')
    plans=[dict(name=f'broad-cubic-max{n}',window='broad',minimum=12,maximum=n,model='cubic',nu=j.NU,omega=j.OMEGA,holdout='larger-than-training' if n<64 else None) for n in [64,48,32]]
    plans += [dict(name=f'inner-raw-background-min{n}',window='inner',minimum=n,maximum=64,model='leading',nu=j.NU,omega=j.OMEGA,holdout=None,raw=True,background=True) for n in [16,24]]
    j.write(folder/'design.json',dict(status='post_hoc_exploratory_not_independent_validation',reason='Broad quadratic train<=48 prediction misses heldout64, while raw regular background shifts Delta. Diagnose thermal cubic/window truncation and raw-background size cuts without selecting an error bar by best p.',configs=plans,source_sha256=source['sha256'],implementation_sha256=j.sha(spec.origin)))
    shutil.copyfile(__file__,folder/Path(__file__).name);results=[]
    for c in plans:
        start=time.monotonic();r,p=j.fit_one(c,rows,B,2000);r['runtime_seconds']=time.monotonic()-start;r['post_hoc']=True;j.write(folder/(c['name']+'.json'),r);np.savez_compressed(folder/(c['name']+'-draws.npz'),parameters=p);results.append(r)
        j.write(folder/'status.json',dict(state='running',completed=len(results),total=len(plans)));print(c['name'],r['central']['Delta'],r.get('parameter_SD'),r['flags'],r.get('heldout',{}).get('nominal_p'),flush=True)
    j.write(folder/'analysis.json',dict(results=results,scope='Post-hoc exploratory diagnosis, not independent heldout validation or a calibrated systematic error.'));j.write(folder/'status.json',dict(state='complete_pending_review',completed=5,total=5))
if __name__=='__main__':main()
