"""Add direct measurements near both edges of the original diagnostic bracket.

Interpolation is used ONLY to select couplings. Final labels still require
actual simulated U, R and eta values and their blocked errors at that coupling.
"""
import argparse
import json
import numpy as np
from scipy.interpolate import PchipInterpolator
from design import freeze
from conservative import data, classify
from boundary import REFERENCE
from run import OUT, point, sha256


def propose():
    lookup=data();points=[];choices=[]
    fine=json.loads((OUT/'pooled.json').read_text())
    eps=sorted({r['epsilon_milli'] for r in fine if r['kappa']>0})
    grid=np.arange(0,1201)/1000
    for e in eps:
        ks=sorted({k for (ee,k,L) in lookup if ee==e and all((e,k,ll) in lookup for ll in (16,32,64))})
        newks=sorted({r['kappa_milli'] for r in fine if r['epsilon_milli']==e and r['L']==64 and r['kappa']>0})
        if len(newks)<3:continue
        labels={k:classify({L:lookup[e,k,L] for L in (16,32,64)})['label'] for k in ks}
        q=np.ones(len(grid),dtype=bool);m=q.copy();estimates={}
        for L in [32,64]:
            rr=[r for r in fine if r['epsilon_milli']==e and r['L']==L and r['kappa']>0]
            if len(rr)<3:break
            for key,lo,hi in [('binder',*REFERENCE[L]['U']),('xi_over_L',*REFERENCE[L]['R'])]:
                y=np.array([lookup[e,k,L]['values'][key] or 0 for k in ks])
                y=np.maximum.accumulate(y)
                # Binder deficit is smoother over coarse gaps; this is selection only.
                if key=='binder':pred=.5-np.exp(PchipInterpolator(np.array(ks)/1000,np.log(np.maximum(.5-y,1e-9)))(grid))
                else:pred=PchipInterpolator(np.array(ks)/1000,y)(grid)
                se=float(np.median([r['errors'][key] for r in rr]))
                q &= pred-2*se>(hi if key=='binder' else lo)
                m &= pred+2*se<lo
                estimates[f'{L}_{key}_se']=se
        else:
            eta=[];eta_err=[]
            for k in ks:
                eta.append(-np.log2(lookup[e,k,64]['values']['m2']/lookup[e,k,32]['values']['m2']))
            for k in newks:
                if (e,k,32) not in lookup:continue
                bs={L:np.load(lookup[e,k,L]['samples'])['samples'][:,2] for L in [32,64]}
                eta_err.append(float(np.std(-np.log2(bs[64]/bs[32]),ddof=1)))
            es=float(np.median(eta_err));pred=PchipInterpolator(np.array(ks)/1000,np.minimum.accumulate(eta))(grid)
            q &= pred+2*es<.25;m &= pred-2*es>.25
            proposed=[];why=[]
            if np.any(m):
                candidate=max(0,int(np.floor(grid[m][-1]*1000/20))*20)
                old_m=max([k for k,v in labels.items() if v=='massive_like'],default=-1)
                while candidate in newks and labels.get(candidate)!='massive_like' and candidate>0:candidate-=20
                if candidate>old_m and candidate not in newks:
                    proposed.append(candidate);why.append('resolve_massive_edge')
            if np.any(q):
                candidate=min(1200,int(np.ceil(grid[q][0]*1000/20))*20)
                old_q=min([k for k,v in labels.items() if v=='quasi_ordered_like'],default=1201)
                while candidate in newks and labels.get(candidate)!='quasi_ordered_like' and candidate<1200:candidate+=20
                if candidate<old_q and candidate not in newks:
                    proposed.append(candidate);why.append('resolve_quasi_ordered_edge')
            proposed=sorted(set(proposed))
            # A fresh high-statistics measurement may replace an old coarse point,
            # but must never overwrite an existing refinement trajectory.
            proposed=[k for k in proposed if not (OUT/'chains'/f'e{e:04d}_k{k:04d}_L016_r000.h5').exists()]
            choices.append(dict(epsilon_milli=e,added_kappa_milli=proposed,reasons=why,eta_se=es,observable_se=estimates))
            for k in proposed:
                for L in [16,32,64]:
                    for rep in range(4):points.append(point(e,k,L,rep,12000))
    return points,choices


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--freeze');args=parser.parse_args()
    points,choices=propose()
    print(json.dumps(choices,indent=2));print('New parameter points',len(points)//12,'chains',len(points))
    if args.freeze:
        if json.loads((OUT/'production-status.json').read_text())['state']!='complete':raise RuntimeError('Production must finish before freezing additions')
        if not points:raise RuntimeError('No additions needed')
        freeze(args.freeze,points,dict(method=__doc__,decisions=choices,pooled_sha256=sha256(OUT/'pooled.json'),selection_source_sha256=sha256(__file__)))

if __name__=='__main__':main()
