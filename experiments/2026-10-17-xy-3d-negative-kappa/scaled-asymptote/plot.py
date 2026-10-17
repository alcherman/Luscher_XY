"""Test a proposed strong-coupling scaling coordinate; never fit its limit by design."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
P=Path(__file__).resolve().parent;a=json.loads((P/'analysis.json').read_text())
styles={90:('o','C0','Uniform'),91:('^','C1','Unit winding'),92:('s','C2','Canted sublattices'),93:('D','C3','Vortex-rich Z4/parity-locked')}
def rep(d):return int(d['point']['id'].split('_r')[1].split('_')[0])
def u(d):p=d['point'];return abs(p['target_kappa'])*(p['epsilon']-np.pi/2)
def peak(d):return float(np.mean([v['power'] for v in d['spectrum_samples']]))
fig,axes=plt.subplots(3,3,figsize=(12,10),sharex=True,sharey='row')
metrics=['vector_peak','staggered_nematic_squared','vortex_density'];labels=['Sampled max spin Fourier power','Staggered nematic power','Vortex density']
for col,k in enumerate([-4,-8,-16]):
 for row,(metric,label) in enumerate(zip(metrics,labels)):
  ax=axes[row,col]
  for r,(marker,color,name) in styles.items():
   rr=sorted([d for d in a['rows'] if d['point']['target_kappa']==k and rep(d)==r],key=u)
   xs=[u(d) for d in rr];ys=[peak(d) if metric=='vector_peak' else d['mean'][metric] for d in rr];ax.plot(xs,ys,marker=marker,color=color,lw=.9,ms=4)
   if metric!='vector_peak':ax.fill_between(xs,[min(h[metric] for h in d['half_means']) for d in rr],[max(h[metric] for h in d['half_means']) for d in rr],color=color,alpha=.12)
  ax.grid(alpha=.18);ax.set_xlim(0,1.65)
  if row==0:ax.set_title(f'κ = {k}, N = 6')
  if col==0:ax.set_ylabel(label)
  if row==2:ax.set_xlabel(r'$u=|\kappa|(\epsilon-\pi/2)$')
  ax.set_ylim((-.02,.65) if row==0 else ((-.02,1.03) if row==1 else (-.01,.30)))
fig.suptitle('Strong-negative scan: does the response align at fixed u?',fontsize=15)
fig.legend(handles=[Line2D([],[],marker=m,color=c,label=t) for m,c,t in styles.values()],loc='lower center',bbox_to_anchor=(.5,.064),ncol=4,fontsize=8,frameon=False)
fig.text(.5,.015,'Each line is one initial-condition class, with independent seeds at its coordinates; lines guide the eye.\nShading spans two-half means, not confidence intervals. Fourier maxima are sampled every50cycles.\nThe u grid tests a hypothesis; no collapse, equilibrium, transition coupling or asymptote is assumed.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.115,1,.96])
for ext in ['png','pdf','svg']:fig.savefig(P/f'scaled-response-check.{ext}',dpi=170)
fig,axes=plt.subplots(2,2,figsize=(10,8),sharex=True,sharey=True)
for ax,(r,(marker,color,name)) in zip(axes.ravel(),styles.items()):
 ax.axvline(.5,color='.45',ls=':',lw=1)
 for d in a['rows']:
  if rep(d)!=r:continue
  p=d['point'];sc=ax.scatter(p['epsilon_over_pi'],p['target_kappa'],c=[peak(d)],marker=marker,s=65,cmap='viridis',vmin=0,vmax=.6,edgecolor='.2',lw=.4)
 ax.set_title(name);ax.grid(alpha=.15);ax.set(xlim=(.498,.633),ylim=(-17,-3),yticks=[-16,-8,-4]);ax.set_xlabel(r'$\epsilon/\pi$');ax.set_ylabel(r'$\kappa$')
fig.subplots_adjust(left=.09,right=.86,top=.89,bottom=.13,hspace=.3,wspace=.24);cax=fig.add_axes([.89,.23,.022,.55]);fig.colorbar(sc,cax=cax,label='Sampled maximum spin Fourier power')
fig.suptitle('N = 6: strong-negative vector response by initial condition',fontsize=14)
fig.text(.5,.025,'Low first-harmonic power can coexist with strong staggered nematic order.\nThese are finite-box response probes; no critical line or fitted asymptote is drawn.',ha='center',fontsize=9)
for ext in ['png','pdf','svg']:fig.savefig(P/f'strong-negative-asymptote-probes.{ext}',dpi=170)
