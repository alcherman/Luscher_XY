"""Magnify the kappa=0 region the user wanted to investigate independently."""
import json
from run import OUT, SIZES
import numpy as np
import matplotlib.pyplot as plt
from analyze import COLORS


def main():
    rows=[r for r in json.loads((OUT/'results.json').read_text())
          if r['point']['kappa_milli']==0 and r['point']['epsilon_milli']<=580]
    exponents=[r for r in json.loads((OUT/'kappa_zero_exponents.json').read_text()) if r['epsilon_over_pi']<=.58]
    fig,axes=plt.subplots(1,2,figsize=(11,4.8))
    for L in SIZES:
        group=[r for r in rows if r['point']['L']==L]
        axes[0].errorbar([r['point']['epsilon_over_pi'] for r in group],
            [r['summary']['values']['binder'] for r in group],
            yerr=[r['summary']['errors']['binder'] for r in group],fmt='o-',capsize=3,
            ms=4,color=COLORS[L],label=f'L={L}')
    for pair,color in [([16,32],COLORS[32]),([32,64],COLORS[64])]:
        group=[r for r in exponents if r['sizes']==pair]
        axes[1].errorbar([r['epsilon_over_pi'] for r in group],[r['eta_eff'] for r in group],
            yerr=[r['error'] for r in group],fmt='o-',capsize=3,ms=4,color=color,label=f'{pair[0]}→{pair[1]}')
    axes[1].axhline(.25,color='.35',ls='--',lw=1,label='Asymptotic BKT η=1/4')
    for ax in axes:
        ax.set_xlim(.497,.583);ax.set_xlabel(r'$\epsilon/\pi$');ax.grid(alpha=.2);ax.legend(fontsize=9)
    axes[0].set_ylabel(r'$U_L=1-\langle m^4\rangle/(2\langle m^2\rangle^2)$')
    axes[1].set_ylabel(r'$\eta_{eff}=-\log_2(\langle m^2\rangle_{2L}/\langle m^2\rangle_L)$')
    fig.suptitle('κ=0: scouting just above the vortex-admissibility threshold')
    fig.text(.5,.014,'One cold-start chain per point; 8,000 retained cycles; conditional 1σ block-bootstrap errors. BKT logarithmic corrections remain.',
             ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.055,1,.95))
    for ext in ('png','pdf','svg'):fig.savefig(OUT/f'near_threshold.{ext}',dpi=160)
    plt.close(fig)


if __name__=='__main__':main()
