"""Plot the final phase curve and the direct kappa-zero diagnostics."""
import json
import csv
import numpy as np
from run import OUT
import matplotlib.pyplot as plt
from conservative import data,REFERENCE


def save(fig,name):
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'{name}.{ext}',dpi=180)
    plt.close(fig)


def main():
    phase=json.loads((OUT/'conservative-brackets.json').read_text())
    rows=list(csv.DictReader((OUT/'phase-brackets.csv').open()))
    fig,ax=plt.subplots(figsize=(10,6.5))
    after=phase['kappa_zero_endpoint']
    for i,r in enumerate(rows):
        x,y=float(r['x']),float(r['kappa'])
        ax.errorbar(x,y,xerr=[[x-float(r['x_lower'])],[float(r['x_upper'])-x]],
            yerr=[[y-float(r['kappa_lower'])],[float(r['kappa_upper'])-y]],fmt='o',color='#0072B2',ms=4,capsize=3,lw=1.2,
            zorder=3)
    ax.set(xlabel=r'$\epsilon/\pi$',ylabel=r'$\kappa$',xlim=(.49,1.015),ylim=(-.04,1.25),title='2D XY phase boundary: L = 16, 32, 64')
    ax.grid(alpha=.15)
    fig.text(.5,.025,'Error bars: finite-volume boundary-location brackets.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.05,1,1));save(fig,'phase_boundary_refined')

    lookup=data();colors={16:'#0072B2',32:'#D55E00',64:'#009E73'}
    fig,axes=plt.subplots(2,2,figsize=(11,8));table=[]
    for L in [16,32,64]:
        selected=sorted((e,r) for (e,k,ll),r in lookup.items() if k==0 and ll==L and 530<=e<=590)
        for j,key in enumerate(['binder','xi_over_L']):
            ax=axes[0,j]
            ax.errorbar([e/1000 for e,r in selected],[r['values'][key] for e,r in selected],
                yerr=[r['errors'][key] for e,r in selected],fmt='o-',ms=4,capsize=3,color=colors[L],lw=.8,label=f'L = {L}')
        for e,r in selected:
            table.append(dict(epsilon_over_pi=e/1000,kappa=0,L=L,source=r['source'],**r['values'],**{key+'_error':v for key,v in r['errors'].items()}))
    for j,(key,label) in enumerate([('U',r'$U_L$'),('R',r'$\xi/L$')]):
        low=min(REFERENCE[L][key][0] for L in [32,64]);high=max(REFERENCE[L][key][1] for L in [32,64])
        axes[0,j].axhspan(low,high,color='.5',alpha=.12,label='Critical-reference envelope')
        axes[0,j].set(ylabel=label,title='Binder cumulant' if j==0 else 'Correlation-length ratio')
        axes[0,j].legend(fontsize=8)
    labels=[r for r in phase['labels'] if r['kappa_milli']==0 and 530<=r['epsilon_milli']<=590 and 'eta_32_64' in r]
    axes[1,0].errorbar([r['epsilon_milli']/1000 for r in labels],[r['eta_32_64'] for r in labels],yerr=[r['eta_error'] for r in labels],fmt='o-',color='#333333',ms=4,capsize=3,lw=.8)
    axes[1,0].axhline(.25,color='.5',ls='--',lw=1)
    axes[1,0].set(ylabel=r'$\eta_{32,64}$',title='Magnetization scaling between L = 32 and 64')
    for r,label,y,color in [(after,'Endpoint',0,'#0072B2')]:
        lo,hi=r['lower'],r['upper']
        axes[1,1].errorbar((lo+hi)/2,y,xerr=(hi-lo)/2,fmt='o',color=color,capsize=6,lw=2)
        axes[1,1].text((lo+hi)/2,y+.15,f'[{lo:.3f}, {hi:.3f}]',ha='center',fontsize=10)
    axes[1,1].set(yticks=[],ylim=(-.5,.5),title=r'$\kappa=0$ boundary-location bracket')
    for ax in axes.flat:
        ax.set(xlabel=r'$\epsilon/\pi$',xlim=(.535,.585));ax.grid(alpha=.15)
    fig.suptitle('κ = 0 endpoint diagnostics',fontsize=15)
    fig.text(.5,.015,'Observable bars: block-bootstrap 1 SE. Endpoint bar: finite-volume boundary-location bracket.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.035,1,.96));save(fig,'endpoint_diagnostics')
    with (OUT/'endpoint-observables.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(table[0]));writer.writeheader();writer.writerows(table)


if __name__=='__main__':main()
