"""Pool raw moments across independent production starts; joint block bootstrap."""
import json
import os
import tempfile
from pathlib import Path
import numpy as np
from scipy.stats import chi2
from run import OUT, atomic_json, sha256, now
from analyze import finite

KEYS=['binder','xi_over_L','m2','m4','link_cosine','vortex_density','winding_squared']


def moments(q,L):
    """Invert derived ratios to sufficient means; never average Binder or xi/L ratios."""
    q=np.asarray(q);a=q[...,2];b=q[...,3]
    fm=L*L*a/(1+(2*L*np.sin(np.pi/L)*q[...,1])**2)
    return np.stack([a,b,fm,q[...,4],q[...,5],q[...,6]],axis=-1)


def ratios(m,L):
    a,b,fm,c,v,w=np.moveaxis(np.asarray(m),-1,0)
    return np.stack([1-b/(2*a*a),np.sqrt(L*L*a/fm-1)/(2*L*np.sin(np.pi/L)),a,b,c,v,w],axis=-1)


def combine(quantities,weights,L):
    return ratios(np.average(moments(quantities,L),weights=weights,axis=0),L)


def main():
    grouped={};folder=OUT/'pooled';folder.mkdir(exist_ok=True)
    for path in sorted((OUT/'analysis').glob('*.json')):
        r=json.loads(path.read_text());p=r['point']
        if p['replica']>=100:continue
        grouped.setdefault((p['epsilon_milli'],p['kappa_milli'],p['L']),[]).append(r)
    result=[]
    for (e,k,L),rows in sorted(grouped.items()):
        if len(rows)!=4:continue
        rows.sort(key=lambda r:r['point']['replica']);weights=[r['summary']['retained'] for r in rows]
        if any(r['summary']['invalid_xi_bootstrap_fraction'] for r in rows):raise ValueError('Unresolved R bootstrap at refinement point')
        q=np.array([[r['summary']['values'][key] for key in KEYS] for r in rows])
        center=combine(q,weights,L)
        bootstrap=combine(np.array([np.load(OUT/'analysis'/f"{r['point']['id']}_bootstrap.npz")['samples'] for r in rows]),weights,L)
        se=bootstrap.std(0,ddof=1);flags=[];tests={}
        for i,key in enumerate(KEYS):
            errs=np.array([r['summary']['errors'][key] for r in rows]);x=q[:,i]
            valid=errs>0
            if valid.sum()>1:
                mean=np.average(x[valid],weights=1/errs[valid]**2)
                statistic=float(np.sum(((x[valid]-mean)/errs[valid])**2));pv=float(chi2.sf(statistic,valid.sum()-1))
                tests[key]=dict(chi2=statistic,df=int(valid.sum()-1),p=pv)
                if pv<.001:flags.append('start_disagreement_'+key)
        name=f'e{e:04d}_k{k:04d}_L{L:03d}'
        # Readers may select the next batch while pooling runs. Publish a whole
        # archive atomically so they cannot observe an empty or partial NPZ.
        with tempfile.NamedTemporaryFile(dir=folder,suffix='.npz',delete=False) as tmp:
            np.savez_compressed(tmp,samples=bootstrap)
            tmp.flush();os.fsync(tmp.fileno());temporary=Path(tmp.name)
        temporary.replace(folder/f'{name}_bootstrap.npz')
        r=dict(id=name,epsilon_milli=e,epsilon_over_pi=e/1000,kappa_milli=k,kappa=k/1000,L=L,
            values=dict(zip(KEYS,center)),errors=dict(zip(KEYS,se)),retained=sum(weights),
            chains=[r['point']['id'] for r in rows],individual_flags={r['point']['id']:r['flags'] for r in rows if r['flags']},
            flags=flags,start_consistency=tests,blocks=sum(r['summary']['blocks'] for r in rows),
            tau_ranges={key:[min(r['summary']['taus'][key] for r in rows if r['summary']['taus'][key] is not None),
                             max(r['summary']['taus'][key] for r in rows if r['summary']['taus'][key] is not None)]
                        for key in rows[0]['summary']['taus'] if any(r['summary']['taus'][key] is not None for r in rows)})
        result.append(finite(r))
    atomic_json(OUT/'pooled.json',result)
    atomic_json(OUT/'pool-provenance.json',dict(created_utc=now(),source_sha256=sha256(Path(__file__)),groups=len(result),
        method='Pool raw means with retained-row weights. Bootstrap contiguous blocks independently within each start; recompute ratios after pooling. Pilot data excluded.',start_flag_threshold=.001))
    print(len(result),'pooled parameter/volume groups',flush=True)

if __name__=='__main__':main()
