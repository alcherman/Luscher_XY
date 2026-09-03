"""Repeat the original joint scouting rule on old plus refined data, unchanged."""
import json
import csv
import numpy as np
from run import OUT, EPSILON_MILLI, atomic_json, now, sha256
from boundary import REFERENCE
from pathlib import Path

SCOUT=OUT.parent/'2026-09-02-xy-bkt-scout'


def data():
    lookup={}
    for r in json.loads((SCOUT/'results.json').read_text()):
        p=r['point'];s=r['summary']
        lookup[(p['epsilon_milli'],p['kappa_milli'],p['L'])]=dict(values=s['values'],errors=s['errors'],flags=r['flags'],enough_blocks=s['enough_blocks'],
            samples=SCOUT/'analysis'/f"{p['id']}_bootstrap.npz",source='scout')
    for folder in [OUT.parent/'2026-09-02-xy-bkt-refine',OUT]:
        for r in json.loads((folder/'pooled.json').read_text()):
            lookup[(r['epsilon_milli'],r['kappa_milli'],r['L'])]=dict(values=r['values'],errors=r['errors'],flags=r['flags'],enough_blocks=r['blocks']>=80,
                samples=folder/'pooled'/f"{r['id']}_bootstrap.npz",source=folder.name)
    return lookup


def classify(g):
    if not all(L in g for L in [16,32,64]):return dict(label='unresolved',reason='missing_volume')
    for L in [32,64]:
        r=g[L]
        if r['flags'] or not r['enough_blocks']:return dict(label='unresolved',reason='sampling_flag')
        if any(r['values'][k] is None or r['errors'][k] is None for k in ['binder','xi_over_L']):return dict(label='unresolved',reason='ratio_unresolved')
    q=m=True
    for L in [32,64]:
        r=g[L];v=r['values'];se=r['errors'];ref=REFERENCE[L]
        q &= v['binder']-2*se['binder']>ref['U'][1] and v['xi_over_L']-2*se['xi_over_L']>ref['R'][0]
        m &= v['binder']+2*se['binder']<ref['U'][0] and v['xi_over_L']+2*se['xi_over_L']<ref['R'][0]
    eta=-np.log2(g[64]['values']['m2']/g[32]['values']['m2'])
    bs={L:np.load(g[L]['samples'])['samples'][:,2] for L in [32,64]}
    se=float(np.std(-np.log2(bs[64]/bs[32]),ddof=1))
    q &= eta+2*se<.25;m &= eta-2*se>.25
    assert not(q and m)
    return dict(label='quasi_ordered_like' if q else 'massive_like' if m else 'unresolved',eta_32_64=eta,eta_error=se)


def main():
    lookup=data();groups={}
    for (e,k,L),v in lookup.items():groups.setdefault((e,k),{})[L]=v
    labels={key:classify(g) for key,g in groups.items()};columns={}
    epsilon_grid=sorted({e for e,k,L in lookup})
    for e in epsilon_grid:
        m=[k/1000 for (ee,k),r in labels.items() if ee==e and r['label']=='massive_like']
        q=[k/1000 for (ee,k),r in labels.items() if ee==e and r['label']=='quasi_ordered_like']
        columns[e]=dict(lower=max(m),upper=min(q)) if m and q and max(m)<min(q) else None
    regions=[]
    for a,b in zip(epsilon_grid[:-1],epsilon_grid[1:]):
        ca,cb=columns[a],columns[b]
        if ca and cb:
            lo,hi=min(ca['lower'],cb['lower']),max(ca['upper'],cb['upper'])
            regions.append(dict(x=(a+b)/2000,x_lower=a/1000,x_upper=b/1000,kappa=(lo+hi)/2,kappa_lower=lo,kappa_upper=hi))
    q0=[e/1000 for (e,k),r in labels.items() if k==0 and r['label']=='quasi_ordered_like']
    m0=[e/1000 for (e,k),r in labels.items() if k==0 and r['label']=='massive_like']
    endpoint=dict(lower=max(q0),upper=min(m0)) if q0 and m0 and max(q0)<min(m0) else None
    atomic_json(OUT/'conservative-brackets.json',dict(created_utc=now(),source_sha256=sha256(Path(__file__)),
        method='Same joint U,R,eta 2-SE rule and adjacent-epsilon-cell union as the original scout. Updated by pooling four independent production starts. This is still a heuristic finite-volume bracket.',
        additional_assumption='A single continuous boundary, locally monotone within each epsilon cell.',
        epsilon_grid=epsilon_grid,columns=columns,regions=regions,kappa_zero_endpoint=endpoint,
        labels=[dict(epsilon_milli=e,kappa_milli=k,**r) for (e,k),r in sorted(labels.items())]))
    table=[dict(kind='interior_grid_region',**r,status='finite_volume_scouting_bracket') for r in regions]
    if endpoint:
        lo,hi=endpoint['lower'],endpoint['upper']
        table.insert(0,dict(kind='zero_kappa_endpoint',x=(lo+hi)/2,x_lower=lo,x_upper=hi,
            kappa=0.,kappa_lower=0.,kappa_upper=0.,status='finite_volume_scouting_bracket'))
    if columns[1000]:
        lo,hi=columns[1000]['lower'],columns[1000]['upper']
        table.append(dict(kind='unconstrained_endpoint',x=1.,x_lower=1.,x_upper=1.,
            kappa=(lo+hi)/2,kappa_lower=lo,kappa_upper=hi,status='finite_volume_scouting_bracket'))
    with (OUT/'phase-brackets.csv').open('w',newline='') as f:
        fields=['kind','x','kappa','x_lower','x_upper','kappa_lower','kappa_upper','status']
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(table)
    print('Conservative regions',len(regions),'endpoint',endpoint,flush=True)

if __name__=='__main__':main()
