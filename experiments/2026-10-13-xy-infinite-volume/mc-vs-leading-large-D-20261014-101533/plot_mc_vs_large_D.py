"""Frozen preliminary MC estimates versus the user-supplied leading large-D relation."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl');os.environ['MPLBACKEND']='Agg'
import json,shutil
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import run
OUT=Path(__file__).resolve().parent

def main():
    source=OUT/json.loads((OUT/'latest-preliminary-all-dimensions.json').read_text())['folder']/'analysis.json'
    a=json.loads(source.read_text());rows=a['rows']
    folder=OUT/('mc-vs-leading-large-D-'+run.now()[:19].replace('-','').replace(':','').replace('T','-'));folder.mkdir()
    shutil.copyfile(__file__,folder/Path(__file__).name)
    shutil.copyfile(source.parent/'preliminary-transition-locations.csv',folder/'preliminary-transition-locations.csv')
    colors={2:'#2563a7',3:'#218463',4:'#b25629'}
    def panel(ax,D):
        vals=sorted([r for r in rows if r['d']==D and r['error_axis']=='vertical'],key=lambda r:r['epsilon_over_pi']);color=colors[D]
        x=np.array([r['epsilon_over_pi'] for r in vals]);y=np.array([r['beta'] for r in vals])
        def errors(low,high,center):return np.array([[c-r[low],r[high]-c] for r,c in zip(vals,center)]).T
        ax.errorbar(x,y,yerr=errors('sensitivity_low','sensitivity_high',y),fmt='none',ecolor=color,alpha=.25,elinewidth=4)
        ax.errorbar(x,y,yerr=errors('MC_low','MC_high',y),fmt='o',color=color,markersize=4,capsize=2,elinewidth=1,label='MC extrapolated points')
        ep=next(r for r in rows if r['d']==D and r['error_axis']=='horizontal');q=ep['epsilon_over_pi']
        for lo,hi,lw,alpha in [(ep['sensitivity_low'],ep['sensitivity_high'],4,.25),(ep['MC_low'],ep['MC_high'],1,1)]:
            ax.errorbar([q],[0],xerr=[[q-lo],[hi-q]],fmt='none',ecolor=color,elinewidth=lw,alpha=alpha,capsize=3)
        ax.plot(q,0,marker='D',ms=5,color=color)
        # (1-q)+beta/2=1/(2D): beta=2(q-1)+1/D, positive branch q>=1-1/(2D).
        prediction_q=np.linspace(1-1/(2*D),1,300);prediction_beta=2*(prediction_q-1)+1/D
        assert abs(prediction_beta[0])<1e-12 and abs(prediction_beta[-1]-1/D)<1e-12
        ax.plot(prediction_q,prediction_beta,'--',color='#222222',lw=1.6,label='Leading large-D prediction')
        ax.set(title=f'd = {D}'+('  (exploratory BKT extrapolation)' if D==2 else ''),xlabel=r'$\epsilon/\pi$',ylabel=r'$\beta_c$')
        ax.set_xlim({2:.52,3:.775,4:.845}[D],1.015);ax.set_ylim(-.025 if D==2 else -.009,{2:1.24,3:.49,4:.33}[D]);ax.grid(alpha=.15);ax.legend(loc='upper left',fontsize=8,frameon=False)
    footer='Dark: conditional MC 95% intervals. Light: descriptive extrapolation sensitivity, not a total confidence interval.\nDiamonds at beta=0 have horizontal bars only. MC points are not connected.'
    fig,axes=plt.subplots(1,3,figsize=(15,5.15))
    for D,ax in zip([2,3,4],axes):panel(ax,D)
    fig.suptitle('Preliminary XY phase boundaries and leading large-D prediction',fontsize=14)
    fig.text(.5,.885,r'$\beta_{\mathrm{large}\ D}=2(\epsilon/\pi-1)+1/D$',ha='center',fontsize=12)
    fig.text(.5,.015,footer,ha='center',fontsize=9);fig.tight_layout(rect=[0,.105,1,.86])
    for ext in ['png','pdf','svg']:fig.savefig(folder/f'combined-MC-vs-large-D.{ext}',dpi=180)
    plt.close(fig)
    for D in [2,3,4]:
        fig,ax=plt.subplots(figsize=(7,5.2));panel(ax,D);fig.text(.5,.015,footer,ha='center',fontsize=7.5);fig.tight_layout(rect=[0,.105,1,1])
        for ext in ['png','pdf','svg']:fig.savefig(folder/f'd{D}-MC-vs-large-D.{ext}',dpi=180)
        plt.close(fig)
    run.atomic_json(folder/'analysis.json',dict(source=str(source),source_sha256=run.sha256(source),rows=rows,prediction=dict(relation='(1-epsilon/pi)+beta/2=1/(2D)',beta='2*(epsilon/pi-1)+1/D',zero_beta_epsilon_over_pi={str(D):1-1/(2*D) for D in [2,3,4]},positive_branch_only=True),scope=a['scope'],plot_policy='No interpolation or connecting MC segments. Dashed user-supplied prediction. Original horizontal/vertical error orientations and all values unchanged.'))
    run.atomic_json(OUT/'latest-MC-vs-large-D.json',dict(folder=folder.name,status='preliminary'))
    print(folder,flush=True)

if __name__=='__main__':main()
