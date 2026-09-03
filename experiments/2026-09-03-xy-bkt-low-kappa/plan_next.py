"""Read-only adaptive proposals; interpolation selects runs, never phase labels."""
import json
import shutil
import numpy as np
from scipy.interpolate import PchipInterpolator
from run import OUT,atomic_json,now,sha256
from conservative import data,classify
from boundary import REFERENCE
from pathlib import Path


def main():
    lookup=data();pooled=json.loads((OUT/'pooled.json').read_text())
    eps=sorted({e for e,k,L in lookup if e<=660})
    grid=np.arange(0,651,10);profiles={};groups={}
    for (e,k,L),r in lookup.items():groups.setdefault((e,k),{})[L]=r
    labels={key:classify(g)['label'] for key,g in groups.items() if all(L in g for L in [16,32,64])}
    for e in eps:
        ks=sorted(k for ee,k in labels if ee==e)
        fine=[k for k in ks if all(lookup[e,k,L]['source']!='scout' for L in [32,64])]
        if len(fine)>=2:
            # Coarse low-statistics interior points must not distort selection
            # between the new high-statistics anchors. They remain in inference.
            ks=[k for k in ks if k in fine or k<min(fine) or k>max(fine)]
        if len(ks)<2:continue
        profile={}
        for L in [32,64]:
            for key in ['binder','xi_over_L']:
                y=np.array([lookup[e,k,L]['values'][key] for k in ks]);y=np.maximum.accumulate(y)
                profile[L,key]=PchipInterpolator(ks,y,extrapolate=False)(np.clip(grid,ks[0],ks[-1]))
        y=np.array([-np.log2(lookup[e,k,64]['values']['m2']/lookup[e,k,32]['values']['m2']) for k in ks])
        profile['eta']=PchipInterpolator(ks,np.minimum.accumulate(y),extrapolate=False)(np.clip(grid,ks[0],ks[-1]))
        profiles[e]=profile
    se={(L,key):float(np.median([r['errors'][key] for r in pooled if r['L']==L])) for L in [32,64] for key in ['binder','xi_over_L']}
    # Conservative selection scale; inference always uses actual joint bootstrap errors.
    eta_se=.0013
    available=sorted(profiles);target_new=[]
    existing_phase=json.loads((OUT/'conservative-brackets.json').read_text()) if (OUT/'conservative-brackets.json').exists() else None
    if existing_phase:
        for r in existing_phase['regions']:
            if r['x']<=.64 and (r['kappa_upper']-r['kappa_lower'])/2>.05000001:
                a,b=round(1000*r['x_lower']),round(1000*r['x_upper'])
                if b-a>=4:target_new.append((a+b)//2)
    else:target_new=[585,595,610,630]
    targets=sorted(set([550,555,560,570,575]+target_new+[e for e in available if 575<=e<=640]))
    owned={(p['epsilon_milli'],p['kappa_milli']) for f in (OUT/'batches').glob('*.json') for p in json.loads(f.read_text())['points']}
    result=[]
    for e in targets:
        if e in profiles:pr=profiles[e]
        else:
            a=max([ee for ee in available if ee<e],default=None);b=min([ee for ee in available if ee>e],default=None)
            if a is None or b is None:continue
            t=(e-a)/(b-a);pr={key:(1-t)*profiles[a][key]+t*profiles[b][key] for key in profiles[a]}
        m=np.ones(len(grid),bool);q=m.copy()
        for L in [32,64]:
            ref=REFERENCE[L]
            m &= (pr[L,'binder']+2*se[L,'binder']<ref['U'][0])&(pr[L,'xi_over_L']+2*se[L,'xi_over_L']<ref['R'][0])
            q &= (pr[L,'binder']-2*se[L,'binder']>ref['U'][1])&(pr[L,'xi_over_L']-2*se[L,'xi_over_L']>ref['R'][0])
        m &= pr['eta']-2*eta_se>.25;q &= pr['eta']+2*eta_se<.25
        lower=max([k for (ee,k),label in labels.items() if ee==e and label=='massive_like'],default=-1)
        upper=min([k for (ee,k),label in labels.items() if ee==e and label=='quasi_ordered_like'],default=1201)
        candidates=[]
        if e<=575 and ((m|q)[0]):candidates.append((0,'resolve_zero_kappa_endpoint'))
        if np.any(m):
            k=int(grid[m][-1])
            if k>lower:candidates.append((k,'tighten_massive_edge'))
        if np.any(q):
            k=int(grid[q][0])
            if k<upper:candidates.append((k,'tighten_quasi_ordered_edge'))
        for k,why in candidates:
            if (e,k) in owned:continue
            if any(r['epsilon_milli']==e and r['kappa_milli']==k for r in result):continue
            result.append(dict(epsilon_milli=e,epsilon_over_pi=e/1000,kappa_milli=k,kappa=k/1000,reason=why,
                new_epsilon=e not in profiles,steps_per_chain=32000,replicas=4,sizes=[16,32,64],
                raw_bytes_upper_estimate=4*(32000+1)*sum(8*L*L+2600 for L in [16,32,64])))
    for r in result:
        e=r['epsilon_over_pi'];k=r['kappa_milli']
        if r['new_epsilon'] and existing_phase:
            cells=[cell for cell in existing_phase['regions'] if cell['x_lower']<e<cell['x_upper']]
            r['priority_score']=max([(cell['kappa_upper']-cell['kappa_lower'])/(2*.05)-1 for cell in cells],default=0.)
            r['priority_reason']='Split a currently wide epsilon cell; score is its excess width relative to the target, not a predicted reduction.'
        else:
            column=existing_phase['columns'].get(str(r['epsilon_milli'])) if existing_phase else None
            gain=(k/1000-column['lower'] if r['reason']=='tighten_massive_edge' else column['upper']-k/1000) if column else 0.
            r['priority_score']=max(0.,gain/.10)
            r['priority_reason']='Tighten an existing coupling edge; score uses the candidate edge shift relative to the target width.'
        if r['reason']=='resolve_zero_kappa_endpoint':
            r['priority_score']=.1;r['priority_reason']='Separate endpoint check; lower priority for the current vertical-width request.'
    result.sort(key=lambda r:(-r['priority_score'],r['epsilon_milli'],r['kappa_milli']))
    remaining=0
    for f in (OUT/'batches').glob('*.json'):
        for p in json.loads(f.read_text())['points']:
            j=OUT/'jobs'/f"{p['id']}.json";done=json.loads(j.read_text()).get('cycles',0) if j.exists() else 0
            remaining+=max(0,p['steps']-done)*(8*p['L']**2+2600)
    record=dict(created_utc=now(),source_sha256=sha256(Path(__file__)),method=__doc__,proposals=result,
        selection_uncertainty='Monotone PCHIP in kappa, constant endpoint continuation outside observed kappa support, and linear interpolation in epsilon are selection heuristics only. Confirm all brackets with directly simulated observations.',
        free_bytes=shutil.disk_usage(OUT).free,pending_upper_bytes=remaining,reserve_bytes=25*2**30,
        suggested_bytes=sum(r['raw_bytes_upper_estimate'] for r in result),
        note='Proposals are candidates, not an instruction to run all of them. Prioritize newly split wide cells and verify storage before freezing a bounded batch.')
    atomic_json(OUT/'next-proposals.json',record)
    print(json.dumps({k:v for k,v in record.items() if k not in ['source_sha256','created_utc']},indent=2))


if __name__=='__main__':main()
