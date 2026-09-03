"""Compare low-kappa grid brackets without changing their scientific meaning."""
from run import OUT
import json
import matplotlib.pyplot as plt


def draw(ax,rows,color,label,alpha=1):
    for i,r in enumerate(rows):
        x,y=r['x'],r['kappa']
        ax.errorbar(x,y,xerr=[[x-r['x_lower']],[r['x_upper']-x]],
            yerr=[[y-r['kappa_lower']],[r['kappa_upper']-y]],fmt='o',ms=4,capsize=3,
            color=color,alpha=alpha,lw=1.1,label=label if i==0 else None)


def main():
    data=json.loads((OUT/'conservative-brackets.json').read_text())
    old=json.loads((OUT.parent/'2026-09-02-xy-bkt-refine/conservative-brackets.json').read_text())
    rows=list(data['regions']);endpoint=data['kappa_zero_endpoint']
    if endpoint:
        lo,hi=endpoint['lower'],endpoint['upper'];rows.insert(0,dict(x=(lo+hi)/2,x_lower=lo,x_upper=hi,kappa=0,kappa_lower=0,kappa_upper=0))
    last=data['columns'].get('1000')
    if last:rows.append(dict(x=1,x_lower=1,x_upper=1,kappa=(last['lower']+last['upper'])/2,kappa_lower=last['lower'],kappa_upper=last['upper']))
    for zoom in [False,True]:
        fig,ax=plt.subplots(figsize=(10,6.5))
        if zoom:draw(ax,old['regions'],'#b0b0b0','Previous refinement',.8)
        draw(ax,rows,'#0072B2','Longer low-κ runs + existing data')
        ax.set(xlabel=r'$\epsilon/\pi$',ylabel=r'$\kappa$',
            xlim=(.53,.665) if zoom else (.49,1.015),ylim=(-.025,.63) if zoom else (-.04,1.25),
            title='Low-κ grid and statistics refinement: L = 16, 32, 64')
        ax.grid(alpha=.15);ax.legend(loc='upper left',fontsize=9)
        fig.text(.5,.025,'Same U/R/η scouting rule; bars include epsilon-cell extent and unresolved finite-volume diagnostics.',ha='center',fontsize=8.5)
        fig.tight_layout(rect=(0,.05,1,1));name='low_kappa_comparison' if zoom else 'phase_boundary_refined'
        for ext in ['png','pdf','svg']:fig.savefig(OUT/f'{name}.{ext}',dpi=180)
        plt.close(fig)


if __name__=='__main__':main()
