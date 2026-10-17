"""Compare cut/block variants on frozen data and the same retained sizes."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json
from pathlib import Path
import numpy as np
import run
from threshold_calibration import size_fit
from scout import json_finite

OUT=Path(__file__).resolve().parent


def main(reference=None,output=None):
    reference=Path(reference) if reference else OUT/'threshold-curve-snapshot-20261013-183856/analysis.json'
    reference=reference.resolve()
    output=Path(output) if output else OUT/'threshold-cut-block-stability-2026-10-13.json'
    base=json.loads(reference.read_text());variants=[]
    for path in sorted(OUT.glob('threshold-curve-snapshot-*-cut*-block*/analysis.json')):
        a=json.loads(path.read_text())
        if a.get('source_snapshot')==str(reference.relative_to(OUT.parents[1])) or a.get('source_snapshot')==str(reference):
            variants.append((path,a))
    chosen={}
    for path,a in variants:chosen[a['cut'],a['block_factor']]=(path,a)
    required=[(501,1),(2001,1),(1001,2)]
    if any(key not in chosen for key in required):
        raise RuntimeError('Cut/block snapshots incomplete; preserve current output and retry later')
    snapshots=[(reference,base)]+[chosen[key] for key in required]
    results=[]
    for D,q in sorted({(g['D'],g['q']) for g in base['groups']}):
        maps=[]
        for path,a in snapshots:
            maps.append({r['N']:r for r in a['locations'] if r['D']==D and r['q']==q and r['threshold']==.25 and r['interpolation']=='linear' and r['status']=='conditional_crossing'})
        common=sorted(set.intersection(*(set(m) for m in maps)))
        item=dict(D=D,epsilon_over_pi=q/100,common_sizes=common,status='insufficient_common_sizes',variants=[])
        if len(common)>=3:
            for (path,a),rows in zip(snapshots,maps):
                B=np.column_stack([np.load(path.parent/f'd{D}_e{q:03d}_N{N:03d}_U25_linear_crossing.npz')['samples'] for N in common])
                valid=np.isfinite(B).all(1)
                fit=size_fit(D,common,[rows[N]['beta_location'] for N in common],B[valid])
                draws=fit.pop('intercept_bootstrap')
                fit['MC_sd']=float(draws.std(ddof=1))
                item['variants'].append(dict(cut=a.get('cut',1001),block_factor=a.get('block_factor',1),fit=fit,joint_valid_fraction=float(valid.mean())))
            values=[v['fit']['beta_infinite'] for v in item['variants']]
            item.update(status='conditional_stability_comparison',beta_range=[min(values),max(values)],max_shift_from_reference=max(abs(v-values[0]) for v in values),
                        tested_MC_interval_union=[min(v['fit']['bootstrap_95'][0] for v in item['variants']),max(v['fit']['bootstrap_95'][1] for v in item['variants'])],
                        doubled_block_MC_sd_ratio=item['variants'][-1]['fit']['MC_sd']/item['variants'][0]['fit']['MC_sd'])
        results.append(item)
    run.atomic_json(output,json_finite(dict(created=run.now(),source=str(reference),snapshots=[str(p) for p,a in snapshots],results=results,
        scope='Same frozen physical-coordinate IDs; comparison fits use the intersection of bracketed sizes. Descriptive sensitivity, not an independent significance test: cuts share data. Conditional MC intervals assume weighted monotonicity; leading size-law failures and low-block/start flags remain relevant. No total-error or mixing certification.')))
    print('common_size_comparisons',sum(r['status']=='conditional_stability_comparison' for r in results),'of',len(results),flush=True)
    for r in results:print(r['D'],r['epsilon_over_pi'],r['status'],r.get('max_shift_from_reference'),flush=True)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--source-snapshot',type=Path)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();main(args.source_snapshot,args.output)
