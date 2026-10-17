"""Compare shared-root GLS fits with archived anchors held fixed."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json
from pathlib import Path
import numpy as np
import run
import argparse,shutil
parser=argparse.ArgumentParser();parser.add_argument('--log-prefix',default='d2-sensitivity');args=parser.parse_args()
from extrapolate_d2_eta_correlation import fit_roots
OUT=Path(__file__).resolve().parent
snapshots={'baseline':OUT/json.loads((OUT/'latest-d2-eta-correlation.json').read_text())['folder']}
for tag in ['cut501','cut4001','block2']:
    log=OUT/'logs'/f'{args.log_prefix}-{tag}-20261014.log'
    snapshots[tag]=next(Path(s.split(' ',1)[1]) for s in log.read_text().splitlines() if s.startswith('snapshot /Users/'))
analyses={k:json.loads((p/'analysis.json').read_text()) for k,p in snapshots.items()}
folder=OUT/('d2-current-chain-sensitivity-'+run.now()[:19].replace('-','').replace(':','').replace('T','-'));folder.mkdir()
shutil.copyfile(__file__,folder/Path(__file__).name)
fits=[];changes=[]
for q in sorted({r['epsilon_over_pi'] for r in analyses['baseline']['locations']}):
    selected={tag:{(r['observable'],r['N']):r for r in a['locations'] if r['epsilon_over_pi']==q and r['N']>=32 and r['status']=='conditional_crossing'} for tag,a in analyses.items()}
    common=set.intersection(*[set(v) for v in selected.values()]);local=[]
    for tag,lookup in selected.items():
        rows=[lookup[k] for k in sorted(common)]
        if rows:
            B=np.column_stack([np.load(snapshots[tag]/f"e{round(q*1000):04d}_{r['observable']}_N{r['N']:03d}.npz")['samples'] for r in rows]);f=fit_roots(rows,B,32,0,0,'joint')
        else:f=dict(status='no_common_roots')
        if 'samples' in f:np.savez_compressed(folder/f'q{round(q*1000)}_{tag}.npz',samples=f.pop('samples'))
        record=dict(epsilon_over_pi=q,variant=tag,common_roots=sorted(common),result=f);fits.append(record);local.append(record)
    base=next(r['result'] for r in local if r['variant']=='baseline')
    if 'beta_infinite' in base:
        successful=[r for r in local if 'beta_infinite' in r['result']]
        changes.append(dict(epsilon_over_pi=q,baseline=base['beta_infinite'],max_absolute_shift=max(abs(r['result']['beta_infinite']-base['beta_infinite']) for r in successful),successful_variants=len(successful),comparison_complete=len(successful)==len(snapshots),unresolved_variants=[r['variant'] for r in local if 'beta_infinite' not in r['result']],maximum_policy='Among successful variants only; unresolved variants do not validate robustness.'))
result=dict(created=run.now(),sources={k:dict(folder=str(p),sha256=run.sha256(p/'analysis.json')) for k,p in snapshots.items()},fits=fits,changes=changes,max_absolute_shift=max(c['max_absolute_shift'] for c in changes),scope='Descriptive common-root minimumN32 C0 joint GLS sensitivity. September archived anchors remain fixed at original cuts/blocks; only current campaign chains change. Pooled moment weights reflect retained rows under each cut. These checks do not certify archived warmup, mixing, finite-size model adequacy or total uncertainty coverage. Correlated analyses of the same simulations are not independent replication. No adjustment to force the known q1 benchmark.')
run.atomic_json(folder/'analysis.json',result);run.atomic_json(OUT/'latest-d2-current-chain-sensitivity.json',dict(folder=folder.name,status='descriptive_not_final'))
print(folder,flush=True);print('max_shift',result['max_absolute_shift'],flush=True)
print(sorted(changes,key=lambda c:c['max_absolute_shift'],reverse=True)[:5],flush=True)
