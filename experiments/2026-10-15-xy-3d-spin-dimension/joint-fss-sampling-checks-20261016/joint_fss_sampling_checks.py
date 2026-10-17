"""Source-pinned sampling checks using the frozen joint-FSS implementation."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,importlib.util,json,time,shutil
import numpy as np
OUT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('frozen_joint_fss',OUT/'joint-fss-existing-20261016-full/joint_fss_existing.py');j=importlib.util.module_from_spec(spec);spec.loader.exec_module(j)

def main():
    folder=OUT/'joint-fss-sampling-checks-20261016';folder.mkdir(exist_ok=False)
    sources=json.loads((OUT/'local-xi-response-20261015-213702/analysis.json').read_text())['sources']
    plans=[dict(name=f'{model}-min{minimum}',window='inner',minimum=minimum,maximum=64,model=model,nu=j.NU,omega=j.OMEGA,holdout=None) for minimum,model in [(12,'leading'),(12,'mixed'),(16,'leading')]]
    j.write(folder/'design.json',dict(frozen_implementation_sha256=j.sha(spec.origin),configurations=plans,variants=['block4','latecut','even','odd'],draws=2000,scope='Existing joint caches only. Nonpaired comparisons across different bootstrap seeds/cuts/start sets; not a systematic-error calibration.'))
    shutil.copyfile(__file__,folder/Path(__file__).name);results=[]
    for variant in ['block4','latecut','even','odd']:
        source=sources[variant]
        if j.sha(source['path'])!=source['sha256']:raise RuntimeError('Source analysis changed')
        rows,B,provenance=j.load_bank(source['path']);j.write(folder/(variant+'-sources.json'),provenance)
        for plan in plans:
            c=dict(plan);c['name']=variant+'-'+plan['name'];start=time.monotonic();r,draws=j.fit_one(c,rows,B,2000);r['runtime_seconds']=time.monotonic()-start;r['sampling_variant']=variant
            j.write(folder/(c['name']+'.json'),r);np.savez_compressed(folder/(c['name']+'-draws.npz'),parameters=draws);results.append(r)
            j.write(folder/'status.json',dict(state='running',completed=len(results),total=12));print(c['name'],r['central']['Delta'],r.get('parameter_SD'),r['flags'],flush=True)
    j.write(folder/'analysis.json',dict(results=results,scope='Sampling/source-start checks of conditional joint finite-size estimators; no calibrated sigma_sys.'))
    j.write(folder/'status.json',dict(state='complete_pending_review',completed=12,total=12))
if __name__=='__main__':main()
