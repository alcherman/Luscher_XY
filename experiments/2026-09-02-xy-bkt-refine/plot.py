"""Show attained grid/MC precision separately from provisional finite-size sensitivity."""
from run import OUT
import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from boundary import REFERENCE

COLORS={16:'#0072B2',32:'#D55E00',64:'#009E73'}
SCOUT=OUT.parent/'2026-09-02-xy-bkt-scout'


def main():
    data=json.loads((OUT/'boundary.json').read_text());rows=data['regions']
    old=json.loads((SCOUT/'phase-brackets.json').read_text())
    conservative=json.loads((OUT/'conservative-brackets.json').read_text())
    fig,ax=plt.subplots(figsize=(10,6.5))
    first=True
    for r in conservative['regions']:
        ax.fill_between([r['x_lower'],r['x_upper']],r['kappa_lower'],r['kappa_upper'],color='.87',lw=0,
            label='Joint diagnostic bracket (finite-size heuristic)' if first else None);first=False
    selected=[r for r in rows if r['primary']['status']=='resolved_finite_volume_proxy']
    first=True
    for r in selected:
        p=r['primary'];x=r['epsilon_over_pi'];y=p['kappa']
        ax.errorbar(x,y,xerr=[[x-r['x_lower']],[r['x_upper']-x]],
            yerr=[[y-p['vertical_lower']],[p['vertical_upper']-y]],fmt='o',ms=4,color='#111111',capsize=2.5,lw=1,
            label='L=64 Binder proxy: grid + 95% MC + interpolation' if first else None);first=False
    for L in [16,32]:
        good=[r for r in rows if r['by_volume'].get(str(L),{}).get('status')=='resolved_finite_volume_proxy']
        if good:ax.plot([r['epsilon_over_pi'] for r in good],[r['by_volume'][str(L)]['kappa'] for r in good],
            '.',ms=4,color=COLORS[L],label=f'L={L} Binder reference midpoint (size check)')
    endpoint=conservative['kappa_zero_endpoint']
    if endpoint:
        x=(endpoint['lower']+endpoint['upper'])/2
        ax.errorbar(x,0,xerr=(endpoint['upper']-endpoint['lower'])/2,fmt='s',color='#777777',capsize=3,ms=4,label='κ=0 diagnostic bracket')
    ax.set(xlabel=r'$\epsilon/\pi$',ylabel=r'$\kappa$',xlim=(.49,1.015),ylim=(-.04,1.25),
        title='Refined XY boundary scouting: grid precision and finite-size ambiguity')
    ax.grid(alpha=.15);ax.legend(loc='upper left',fontsize=8,framealpha=.95)
    fig.text(.5,.028,'Small bars are conditional finite-volume precision, not total uncertainty in the thermodynamic BKT curve.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.05,1,1))
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'precision_and_finite_size.{ext}',dpi=180)
    plt.close(fig)
    # Faithful before/after comparison using exactly the same conservative representation.
    fig,ax=plt.subplots(figsize=(10,6.5))
    for group,color,shift,label in [(old['regions'],'#b0b0b0',-.001,'Original Δκ=0.2 scout'),(conservative['regions'],'#0072B2',.001,'Refined data; same joint diagnostic rule')]:
        for i,r in enumerate(group):
            x=r['x'];y=r['kappa']
            ax.errorbar(x+shift,y,xerr=[[x-r['x_lower']],[r['x_upper']-x]],yerr=[[y-r['kappa_lower']],[r['kappa_upper']-y]],fmt='o',color=color,ms=3,capsize=2,lw=1,label=label if i==0 else None)
    ax.set(xlabel=r'$\epsilon/\pi$',ylabel=r'$\kappa$',xlim=(.49,1.015),ylim=(-.04,1.27),title='Conservative grid brackets: before and after refinement')
    ax.grid(alpha=.15);ax.legend(loc='upper left',fontsize=9)
    fig.text(.5,.025,'Bars span neighboring epsilon cells and all joint-diagnostic ambiguity; they are not statistical confidence intervals.',ha='center',fontsize=8.5)
    fig.tight_layout(rect=(0,.05,1,1))
    for ext in ['png','pdf']:fig.savefig(OUT/f'conservative_comparison.{ext}',dpi=170)
    plt.close(fig)
    # Main deliverable keeps the same conservative meaning as the original plot.
    fig,ax=plt.subplots(figsize=(10,6.5))
    for i,r in enumerate(conservative['regions']):
        x=r['x'];y=r['kappa']
        ax.errorbar(x,y,xerr=[[x-r['x_lower']],[r['x_upper']-x]],
            yerr=[[y-r['kappa_lower']],[r['kappa_upper']-y]],fmt='o',color='#0072B2',ms=4,capsize=3,lw=1.1,
            label='Refined joint-diagnostic boundary brackets' if i==0 else None)
    if endpoint:
        x=(endpoint['lower']+endpoint['upper'])/2
        ax.errorbar(x,0,xerr=(endpoint['upper']-endpoint['lower'])/2,fmt='o',color='#0072B2',capsize=3,ms=4)
    last=conservative['columns'].get('1000')
    if last:
        y=(last['lower']+last['upper'])/2
        ax.errorbar(1.,y,yerr=(last['upper']-last['lower'])/2,fmt='o',color='#0072B2',capsize=3,ms=4)
    ax.set(xlabel=r'$\epsilon/\pi$',ylabel=r'$\kappa$',xlim=(.49,1.015),ylim=(-.04,1.25),
        title='XY boundary refinement: L = 16, 32, 64')
    ax.grid(alpha=.15);ax.legend(loc='upper left',fontsize=9)
    fig.text(.5,.026,'Bars include grid-cell extent and unresolved joint diagnostics; finite-size scouting brackets, not 1σ errors.',ha='center',fontsize=8.5)
    fig.tight_layout(rect=(0,.05,1,1))
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'phase_boundary_refined.{ext}',dpi=180)
    plt.close(fig)
    # Visible precision budget; x-grid and kappa resolutions are different parameter units.
    fig,ax=plt.subplots(figsize=(9,4.5))
    if selected:
        x=[r['epsilon_over_pi'] for r in selected]
        ax.plot(x,[r['primary']['mc_se'] for r in selected],'o-',label='κ bootstrap standard error')
        ax.plot(x,[max(r['primary']['kappa']-r['primary']['vertical_lower'],r['primary']['vertical_upper']-r['primary']['kappa']) for r in selected],'s-',label='κ plotted half-width (95% MC + grid + interpolation)')
        ax.plot(x,[max(r['epsilon_over_pi']-r['x_lower'],r['x_upper']-r['epsilon_over_pi']) for r in selected],'--',label='ε/π grid half-width (different coordinate)')
    ax.set(xlabel=r'$\epsilon/\pi$',ylabel='Numerical half-width / standard error',title='Attained sampling and grid resolution');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(OUT/'precision_budget.png',dpi=160);plt.close(fig)
    # Direct data supporting each local interpolation, with all three volumes.
    pooled=json.loads((OUT/'pooled.json').read_text());folder=OUT/'slices';folder.mkdir(exist_ok=True)
    for e in sorted({r['epsilon_milli'] for r in pooled if r['kappa']>0}):
        fig,axes=plt.subplots(1,2,figsize=(11,4.3))
        for L in [16,32,64]:
            group=sorted([r for r in pooled if r['epsilon_milli']==e and r['L']==L and r['kappa']>0],key=lambda r:r['kappa'])
            for key,ax,refkey in [('binder',axes[0],'U'),('xi_over_L',axes[1],'R')]:
                ax.errorbar([r['kappa'] for r in group],[r['values'][key] for r in group],yerr=[r['errors'][key] for r in group],fmt='o-',capsize=3,ms=4,color=COLORS[L],label=f'L={L}')
                lo,hi=REFERENCE[L][refkey];ax.axhspan(lo,hi,color=COLORS[L],alpha=.07)
                ax.set_xlabel(r'$\kappa$');ax.grid(alpha=.15)
        axes[0].set_ylabel(r'$U_L$');axes[1].set_ylabel(r'$\xi_{2\mathrm{nd}}/L$');axes[0].legend(fontsize=8)
        fig.suptitle(f'ε/π = {e/1000:g}; four independent starts per point, errors = 1 SE')
        fig.text(.5,.012,'Shading: finite-size critical reference windows from the two endpoint models, not data error bands.',ha='center',fontsize=8)
        fig.tight_layout(rect=(0,.04,1,1));fig.savefig(folder/f'e{e:04d}.png',dpi=150);plt.close(fig)
    print(len(selected),'resolved finite-volume proxy points plotted',flush=True)

if __name__=='__main__':main()
