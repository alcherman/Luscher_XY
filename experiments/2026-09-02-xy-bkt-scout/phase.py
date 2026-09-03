"""Conservative finite-volume scouting brackets; not a thermodynamic BKT fit.

Use the envelope of measured finite-volume critical U and xi/L in the ordinary
and constraint-only XY models to identify points clearly on either side.
Assuming this envelope is useful for the mixed action is an explicit scouting
heuristic, not a theorem. Keep ambiguous points unresolved. Full precision
requires validated finite-size matching and independent starts.
"""
import csv
import json
from pathlib import Path
from run import OUT, EPSILON_MILLI, KAPPA_MILLI, atomic_json, sha256, now
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# https://arxiv.org/pdf/1307.0485 Tables 2 and 3. Convert their U4 to our 1-U4/2.
# These observable windows do not prescribe either endpoint's parameter value.
REFERENCE={
    32:dict(U=(1-1.01868/2,1-1.018521/2),R=(.7441,.7913)),
    64:dict(U=(1-1.018871/2,1-1.01780/2),R=(.7434,.7861)),
}


def classify(group):
    if not all(L in group for L in (16,32,64)):return dict(label='unresolved',reason='missing_volume')
    for L in (32,64):
        r=group[L];s=r['summary']
        if r['flags'] or not s['enough_blocks']:
            return dict(label='unresolved',reason='sampling_diagnostic_flag')
        if any(s['values'][k] is None or s['errors'][k] is None for k in ('binder','xi_over_L')):
            return dict(label='unresolved',reason='unresolved_moment_ratio')
    quasi=True;massive=True
    for L in (32,64):
        s=group[L]['summary'];u=s['values']['binder'];r=s['values']['xi_over_L']
        du=s['errors']['binder'];dr=s['errors']['xi_over_L'];ref=REFERENCE[L]
        quasi &= u-2*du>ref['U'][1] and r-2*dr>ref['R'][0]
        massive &= u+2*du<ref['U'][0] and r+2*dr<ref['R'][0]
    a,b=group[32],group[64]
    eta=-np.log(b['summary']['values']['m2']/a['summary']['values']['m2'])/np.log(2)
    sa=np.load(OUT/'analysis'/f"{a['point']['id']}_bootstrap.npz")['samples'][:,2]
    sb=np.load(OUT/'analysis'/f"{b['point']['id']}_bootstrap.npz")['samples'][:,2]
    err=float(np.nanstd(-np.log(sb/sa)/np.log(2),ddof=1))
    quasi &= eta+2*err<.25
    massive &= eta-2*err>.25
    assert not (quasi and massive)
    return dict(label='quasi_ordered_like' if quasi else 'massive_like' if massive else 'unresolved',
                eta_32_64=float(eta),eta_error=err,reason='finite_volume_joint_diagnostics')


def bracket_at_x(e,labels):
    m=[k/1000 for k in KAPPA_MILLI if labels.get((e,k),{}).get('label')=='massive_like']
    q=[k/1000 for k in KAPPA_MILLI if labels.get((e,k),{}).get('label')=='quasi_ordered_like']
    if not m or not q:return None
    low,high=max(m),min(q)
    if low>=high:return dict(status='inconsistent_ordering')
    return dict(status='bracketed',lower=low,upper=high)


def main():
    rows=json.loads((OUT/'results.json').read_text())
    grouped={}
    for r in rows:
        p=r['point'];grouped.setdefault((p['epsilon_milli'],p['kappa_milli']),{})[p['L']]=r
    labels={key:classify(g) for key,g in grouped.items()}
    columns={e:bracket_at_x(e,labels) for e in EPSILON_MILLI}
    regions=[]
    # Boundary regions spanning an epsilon grid cell; y limits enclose both
    # neighboring column brackets. Midpoints are plotting representatives.
    for a,b in zip(EPSILON_MILLI[:-1],EPSILON_MILLI[1:]):
        ca,cb=columns[a],columns[b]
        if ca and cb and ca['status']==cb['status']=='bracketed':
            lo,hi=min(ca['lower'],cb['lower']),max(ca['upper'],cb['upper'])
            regions.append(dict(kind='interior_grid_region',x_lower=a/1000,x_upper=b/1000,
                kappa_lower=lo,kappa_upper=hi,x=(a+b)/2000,kappa=(lo+hi)/2,
                status='finite_volume_scouting_bracket'))
    q0=[e/1000 for e in EPSILON_MILLI if labels.get((e,0),{}).get('label')=='quasi_ordered_like']
    m0=[e/1000 for e in EPSILON_MILLI if labels.get((e,0),{}).get('label')=='massive_like']
    if q0 and m0 and max(q0)<min(m0):
        lo,hi=max(q0),min(m0)
        regions.insert(0,dict(kind='zero_kappa_endpoint',x_lower=lo,x_upper=hi,kappa_lower=0.,
            kappa_upper=0.,x=(lo+hi)/2,kappa=0.,status='finite_volume_scouting_bracket'))
    last=columns[1000]
    if last and last['status']=='bracketed':
        regions.append(dict(kind='unconstrained_endpoint',x_lower=1.,x_upper=1.,
            kappa_lower=last['lower'],kappa_upper=last['upper'],x=1.,
            kappa=(last['lower']+last['upper'])/2,status='finite_volume_scouting_bracket'))
    result=dict(created_utc=now(),source_sha256=sha256(Path(__file__)),
        reference_observable_windows=REFERENCE,
        reference_source='https://arxiv.org/pdf/1307.0485, Tables 2 and 3',
        method='Joint U and xi/L critical-observable envelope on L32,L64, with two-SE margins and eta_eff cross-check; L16 retained as a size check.',
        limitation='Heuristic finite-volume brackets, not confidence intervals or an infinite-volume fit. Mixed-action critical values are not proved to lie in this envelope. No parameter endpoint is forced to its literature value.',
        additional_assumption='A single continuous boundary, locally monotone within adjacent epsilon cells; unresolved cells are left blank.',
        regions=regions,columns={str(k):v for k,v in columns.items()},
        labels=[dict(epsilon_over_pi=e/1000,kappa=k/1000,**v) for (e,k),v in sorted(labels.items())])
    atomic_json(OUT/'phase-brackets.json',result)
    with (OUT/'phase-brackets.csv').open('w',newline='') as f:
        fields=['kind','x','kappa','x_lower','x_upper','kappa_lower','kappa_upper','status']
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(regions)
    fig,ax=plt.subplots(figsize=(10,6.5))
    colors={'quasi_ordered_like':'#009E73','massive_like':'#0072B2','unresolved':'#c5c5c5'}
    names={'quasi_ordered_like':'Quasi-ordered-like','massive_like':'Massive-like','unresolved':'Unresolved'}
    for label,color in colors.items():
        group=[(e/1000,k/1000) for (e,k),v in labels.items() if v['label']==label]
        if group:
            x,y=np.array(group).T
            ax.scatter(x,y,c=color,s=20,alpha=.7,label=names[label],zorder=1)
    for i,r in enumerate(regions):
        ax.errorbar(r['x'],r['kappa'],xerr=[[r['x']-r['x_lower']],[r['x_upper']-r['x']]],
            yerr=[[r['kappa']-r['kappa_lower']],[r['kappa_upper']-r['kappa']]],fmt='o',
            color='#222222',ms=4,capsize=3,lw=1.1,zorder=3,
            label='Scouting boundary regions (grid brackets)' if i==0 else None)
    ax.set(xlabel=r'$\epsilon/\pi$',ylabel=r'$\kappa$',xlim=(.49,1.015),ylim=(-.06,1.27),
        title='XY phase-boundary scouting: L = 16, 32, 64')
    ax.grid(alpha=.15);ax.legend(loc='upper left',fontsize=9,framealpha=.9)
    fig.text(.5,.024,'Bars span diagnostic grid brackets; midpoints are representatives. Finite-size bias is not controlled.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.05,1,1))
    for ext in ('png','pdf','svg'):fig.savefig(OUT/f'phase_boundary_scout.{ext}',dpi=160)
    plt.close(fig)
    print(len(regions),'scouting regions;', {label:sum(v['label']==label for v in labels.values()) for label in colors})


if __name__=='__main__': main()
