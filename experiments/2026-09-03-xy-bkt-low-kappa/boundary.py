"""Resolve finite-volume Binder locations, with grid/MC and model sensitivity separate.

These are conditional finite-size proxy locations, not an infinite-volume fit.
The endpoint reference envelope is a heuristic retained from the scout.
"""
import csv
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.stats import chi2
from run import OUT, EPSILON_MILLI, atomic_json, sha256, now
from analyze import finite

REFERENCE={16:dict(U=(1-1.02025/2,1-1.017811/2),R=(.7445,.8014)),
           32:dict(U=(1-1.01868/2,1-1.018521/2),R=(.7441,.7913)),
           64:dict(U=(1-1.018871/2,1-1.01780/2),R=(.7434,.7861))}


def crossing(x,y,target):
    """Unique increasing piecewise-linear crossing, without extrapolation or clipping."""
    x=np.asarray(x);y=np.asarray(y)
    if len(x)<2 or np.any(np.diff(x)<=0) or not np.isfinite(y).all():return np.nan
    exact=np.where(y==target)[0]
    if len(exact)>1:return np.nan
    hits=[]
    for i in range(len(x)-1):
        if y[i]<=target<=y[i+1] and y[i+1]>y[i]:
            hits.append(x[i]+(target-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]))
        elif y[i]>target>y[i+1]:return np.nan
    hits=np.unique(hits)
    return float(hits[0]) if len(hits)==1 else np.nan


def linear_roots(x,y,errors,target):
    x=np.asarray(x);y=np.asarray(y);errors=np.asarray(errors)
    origin=x.mean();X=np.column_stack([np.ones(len(x)),x-origin])
    projector=np.linalg.pinv(X/errors[:,None])/errors[None,:]
    coefficients=y@projector.T
    with np.errstate(divide='ignore',invalid='ignore'):
        roots=origin+(target-coefficients[...,0])/coefficients[...,1]
    valid=(coefficients[...,1]>0)&(roots>=x.min())&(roots<=x.max())
    return np.where(valid,roots,np.nan),coefficients


def local_window(x,y,errors,target):
    """Smallest consecutive three-point neighborhood with an interior fitted root."""
    if len(x)<=3:return np.arange(len(x))
    eligible=[]
    for first in range(len(x)-2):
        ix=np.arange(first,first+3);root,_=linear_roots(x[ix],y[ix],errors[ix],target)
        if np.isfinite(root):eligible.append((round(float(x[ix[-1]]-x[ix[0]]),10),abs(float(root)-x[ix[1]]),first))
    if not eligible:return np.arange(len(x))
    first=min(eligible)[2]
    return np.arange(first,first+3)


def locate(rows,key,target,index):
    rows=sorted(rows,key=lambda r:r['kappa']);x=np.array([r['kappa'] for r in rows]);y=np.array([r['values'][key] for r in rows])
    errors=np.array([r['errors'][key] for r in rows])
    if len(x)<3:return dict(status='too_few_points',x=x.tolist(),y=y.tolist(),target=target)
    ix=local_window(x,y,errors,target);rows=[rows[i] for i in ix];x=x[ix];y=y[ix];errors=errors[ix]
    center,coef=linear_roots(x,y,errors,target);center=float(center)
    if not np.isfinite(center):return dict(status='outside_sampled_crossing',x=x.tolist(),y=y.tolist(),target=target,linear_coefficients=coef.tolist(),linear_origin=float(x.mean()))
    draws=np.array([np.load(OUT/'pooled'/f"{r['id']}_bootstrap.npz")['samples'][:,index] for r in rows]).T
    bounded,coefficients=linear_roots(x,draws,errors,target)
    fraction=float(np.isfinite(bounded).mean())
    # Keep uncertainty tails rather than silently conditioning them on lying in the grid.
    with np.errstate(divide='ignore',invalid='ignore'):
        roots=x.mean()+(target-coefficients[:,0])/coefficients[:,1]
    valid=np.isfinite(roots)&(coefficients[:,1]>0)
    if valid.sum()<2:return dict(status='unresolved_bootstrap',target=target)
    bracket=next(i for i in range(len(x)-1) if x[i]<=center<=x[i+1])
    grid_width=float(x[bracket+1]-x[bracket])
    q=np.quantile(roots[valid],[.025,.16,.84,.975]);curvature=None
    if np.all(np.diff(y)>0):
        smooth=PchipInterpolator(x,y)
        if y[0]<=target<=y[-1]:
            other=brentq(lambda z:smooth(z)-target,x[0],x[-1]);curvature=abs(other-center)
    statistic=float(np.sum(((y-(coef[0]+coef[1]*(x-x.mean())))/errors)**2))
    pv=float(chi2.sf(statistic,len(x)-2))
    flags=sorted(set(f for r in rows for f in r['flags']))
    if pv<.001:flags.append('local_line_lack_of_fit')
    status='resolved_finite_volume_proxy' if fraction>=.99 and not flags and grid_width<=.02000001 else 'needs_refinement'
    interpolation=curvature or 0.
    return finite(dict(status=status,kappa=center,bootstrap_95=[q[0],q[3]],bootstrap_68=[q[1],q[2]],
        mc_se=float(np.std(roots[valid],ddof=1)),bootstrap_valid_fraction=fraction,
        bootstrap_positive_slope_fraction=float(valid.mean()),
        grid_spacing=grid_width,grid_halfwidth=grid_width/2,
        vertical_lower=q[0]-grid_width/2-interpolation,vertical_upper=q[3]+grid_width/2+interpolation,
        interpolation_sensitivity=curvature,target=target,flags=flags,x=x.tolist(),y=y.tolist(),
        fit_ids=[r['id'] for r in rows],slope=float(coef[1]),linear_fit_chi2=statistic,linear_fit_df=len(x)-2,linear_fit_p=pv))


def main():
    rows=json.loads((OUT/'pooled.json').read_text());grouped={}
    for r in rows:
        if r['kappa']!=0:grouped.setdefault((r['epsilon_milli'],r['L']),[]).append(r)
    result=[]
    for e in sorted({e for e,L in grouped}):
        tests={};locations={}
        for L in (16,32,64):
            group=grouped.get((e,L),[])
            if len(group)<2:continue
            ref=REFERENCE[L]
            locations[str(L)]=locate(group,'binder',float(np.mean(ref['U'])),0)
            tests[str(L)]={}
            for key,index in [('binder',0),('xi_over_L',1)]:
                window=ref['U' if key=='binder' else 'R']
                tests[str(L)][key]=[locate(group,key,target,index) for target in window]
        if '64' not in locations:continue
        ix=EPSILON_MILLI.index(e);x=e/1000
        lower=(EPSILON_MILLI[ix-1]+e)/2000 if ix>0 else x
        upper=(EPSILON_MILLI[ix+1]+e)/2000 if ix<len(EPSILON_MILLI)-1 else x
        result.append(dict(epsilon_milli=e,epsilon_over_pi=x,x_lower=lower,x_upper=upper,
            primary=locations['64'],by_volume=locations,reference_sensitivity=tests))
    atomic_json(OUT/'boundary.json',dict(created_utc=now(),source_sha256=sha256(Path(__file__)),
        reference=REFERENCE,reference_source='https://arxiv.org/pdf/1307.0485 Tables 2 and 3',
        method='Conditional L64 Binder-envelope midpoint location; Local weighted-linear interpolation in the smallest consecutive three-point neighborhood with an interior fitted root, with joint, stratified block-bootstrap ratios; no extrapolated root is accepted. Vertical bars are 95% MC interval expanded by half the actual kappa grid spacing and the shift from monotone PCHIP where available. Horizontal bars are epsilon sampling cells.',
        limitation='This changes the old all-diagnostic grid-region representation. It isolates attained grid/MC resolution. Reference-choice, finite-volume, and other-observable sensitivity are separate and are not a certified thermodynamic uncertainty bound.',
        regions=result))
    with (OUT/'boundary.csv').open('w',newline='') as f:
        fields=['epsilon_over_pi','kappa','x_lower','x_upper','vertical_lower','vertical_upper','mc_se','grid_spacing','bootstrap_valid_fraction','status']
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for r in result:w.writerow({k:r.get(k,r['primary'].get(k)) for k in fields})
    print('Boundary:',{s:sum(r['primary']['status']==s for r in result) for s in {r['primary']['status'] for r in result}},flush=True)

if __name__=='__main__':main()
