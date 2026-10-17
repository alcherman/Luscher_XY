"""Descriptive common-root sensitivity checks; no independent-replication claim."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json
from pathlib import Path
import numpy as np
import run
from threshold_calibration import size_fit
from extrapolate_d2_eta_correlation import fit_roots
OUT=Path(__file__).resolve().parent
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--log-prefix');args=parser.parse_args()
tags=['legacy-local','current-local','cut501','cut3001','block2','U20','U30']
snapshots={}
for tag in tags:
    log=OUT/'logs'/(f'{args.log_prefix}-{tag}.log' if args.log_prefix else f'endpoint-sensitivity-{tag}-20261014.log')
    folder=next(Path(s) for s in log.read_text().splitlines() if s.startswith('/Users/') and 'snapshot' in s)
    snapshots[tag]=(folder,json.loads((folder/'analysis.json').read_text()))
folder=OUT/('endpoint-sensitivity-summary-'+run.now()[:19].replace('-','').replace(':','').replace('T','-'));folder.mkdir()
fits=[];changes=[]
def eligible(a,D):
    return {(r['observable'],r['N']):r for r in a['locations'] if r['D']==D and r['interpolation']=='linear' and r['status']=='conditional_crossing' and r['N']>=({2:64,3:12,4:10}[D])}
for D in [2,3,4]:
    selected_tags=['current-local','cut501','cut3001','block2']
    common=set.intersection(*[set(eligible(snapshots[t][1],D)) for t in selected_tags])
    for tag in selected_tags:
        source,a=snapshots[tag];lookup=eligible(a,D);rows=[lookup[k] for k in sorted(common)]
        if not rows:
            fits.append(dict(D=D,variant=tag,common_roots=[],result=dict(status='no_common_roots'),source=str(source)));continue
        B=np.column_stack([np.load(source/f"d{D}_{r['observable']}_N{r['N']}_linear_epsilon_roots.npz")['samples'] for r in rows])
        if D==2:
            adapted=[dict(r,beta_location=r['epsilon_over_pi_location'],pair=[r['N'],2*r['N']] if r['observable']=='eta_eff' else None) for r in rows]
            f=fit_roots(adapted,B,64,0,0,'joint');draws=f.pop('samples',None)
        else:
            valid=np.isfinite(B).all(1)
            f=size_fit(D,[r['N'] for r in rows],[r['epsilon_over_pi_location'] for r in rows],B[valid],{3:12,4:10}[D]);draws=f.pop('intercept_bootstrap',None)
            f['joint_valid_fraction']=float(valid.mean())
        f['input_flags']=sorted({flag for r in rows for flag in r.get('flags',[])})
        f['individual_flags']={k:v for r in rows for k,v in r.get('individual_flags',{}).items()}
        if draws is not None:
            f['epsilon_over_pi_infinite']=f.pop('beta_infinite');np.savez_compressed(folder/f'd{D}_{tag}_epsilon_infinite.npz',samples=draws)
        fits.append(dict(D=D,variant=tag,common_roots=sorted(common),result=f,source=str(source),source_sha256=run.sha256(source/'analysis.json')))
    local=[f for f in fits if f['D']==D and 'epsilon_over_pi_infinite' in f['result']]
    baseline=next((f for f in local if f['variant']=='current-local'),None)
    if baseline:
        base=baseline['result']['epsilon_over_pi_infinite']
        changes.append(dict(D=D,baseline=base,max_absolute_shift=max(abs(f['result']['epsilon_over_pi_infinite']-base) for f in local),successful_variants=len(local),comparison_complete=len(local)==len(selected_tags),scope='Maximum among successful common-current-root variants, same minimum size and C0/leading model; correlated cut/block/window changes, descriptive only.'))
result=dict(created=run.now(),snapshots={k:str(v[0]) for k,v in snapshots.items()},common_root_fits=fits,cut_block_changes=changes,scope='No accepted final endpoint or calibrated total confidence interval. Common-root current-chain comparisons isolate current-data analysis choices without legacy source composition changes. These are correlated analyses of the same simulations, not independent replications. Threshold, interpolation, local/full windows, minimum-size and correction fits remain available in each source snapshot. Failed cut4001 analysis preserved; replaced by cut3001 because D4 runs have4000 cycles. Empty-retention guard was added. Legacy inclusion comparisons are separate from cut/block effects.')
run.atomic_json(folder/'analysis.json',result)
run.atomic_json(OUT/'latest-endpoint-sensitivity.json',dict(folder=folder.name,status='descriptive_checks_not_final_acceptance'))
run.atomic_json(OUT/'latest-zero-beta-endpoints.json',dict(folder=snapshots['legacy-local'][0].name,status='local_brackets_with_legacy_preliminary'))
print(folder,flush=True)
for r in changes:print(r,flush=True)
