"""Finite-box size comparisons, without fitting a transition or asymptote."""
from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
P=Path(__file__).resolve().parent
rows=sum([json.loads((P/f).read_text())['rows'] for f in ['scaled-N6-analysis.json','size-validation-analysis.json']],[])
styles=[('o','C0','Uniform'),('^','C1','Unit winding'),('s','C2','Canted sublattices'),('D','C3','Vortex-rich parity-locked Z4')]
def start(d):return int(d['point']['id'].split('_r')[1].split('_')[0])%10
def peak(d):return np.mean([x['power'] for x in d['spectrum_samples']])
fig,axes=plt.subplots(3,3,figsize=(12,10),sharex=True,sharey='row')
for col,k in enumerate([-4,-8,-16]):
 for row,(metric,label) in enumerate([('peak','Max spin Fourier power'),('staggered_nematic_squared','Staggered nematic power'),('vortex_density','Vortex density')]):
  ax=axes[row,col]
  for n,ls in [(6,'--'),(8,'-')]:
   for j,(m,c,_) in enumerate(styles):
    rr=sorted([d for d in rows if d['point']['N']==n and d['point']['target_kappa']==k and start(d)==j and d['point']['role']!='off_scaling_energy_sector_control'],key=lambda d:d['point']['scaled_u'])
    ax.plot([d['point']['scaled_u'] for d in rr],[peak(d) if metric=='peak' else d['mean'][metric] for d in rr],marker=m,color=c,ls=ls,lw=.9,ms=4,alpha=.8)
  ax.grid(alpha=.18)
  if row==0:ax.set_title(f'κ = {k}')
  if col==0:ax.set_ylabel(label)
  if row==2:ax.set_xlabel(r'$|\kappa|(\epsilon-\pi/2)$')
fig.legend(handles=[Line2D([],[],marker=m,color=c,label=t) for m,c,t in styles]+[Line2D([],[],color='k',ls='--',label='N=6'),Line2D([],[],color='k',label='N=8')],loc='lower center',bbox_to_anchor=(.5,.05),ncol=3,fontsize=8)
fig.suptitle('Strong-negative response: size and initial-condition checks')
fig.text(.5,.012,'Independent trajectories at each coordinate. Lines guide the eye; no accepted critical line or asymptote.\nFourier maxima sampled every50cycles; finite-volume noise floor and trapping prevent precision claims.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.14,1,.96])
for e in ['png','pdf','svg']:fig.savefig(P/f'scaled-size-response-comparison.{e}',dpi=170)
fig,axes=plt.subplots(1,2,figsize=(12,5),sharex=True,sharey=True)
for ax,n in zip(axes,[6,8]):
 ax.axvline(.5,color='.4',ls=':')
 for j,(m,c,t) in enumerate(styles):
  rr=[d for d in rows if d['point']['N']==n and start(d)==j]
  sc=ax.scatter([d['point']['epsilon_over_pi'] for d in rr],[d['point']['target_kappa'] for d in rr],c=[peak(d) for d in rr],vmin=0,vmax=.6,cmap='viridis',marker=m,s=45,edgecolors='.3',lw=.4)
 ax.set(title=f'N={n}',xlabel=r'$\epsilon/\pi$',ylabel=r'$\kappa$',ylim=(-17,-3));ax.grid(alpha=.15)
fig.subplots_adjust(left=.07,right=.86,bottom=.24,top=.88,wspace=.15);fig.colorbar(sc,cax=fig.add_axes([.89,.27,.02,.55]),label='Sampled max spin Fourier power')
fig.legend(handles=[Line2D([],[],ls='',marker=m,color='k',label=t) for m,c,t in styles],loc='lower center',bbox_to_anchor=(.5,.08),ncol=2,fontsize=8)
fig.text(.5,.025,'Response map, not an equilibrium phase boundary. Persistent initial-condition dependence is retained.',ha='center',fontsize=9)
for e in ['png','pdf','svg']:fig.savefig(P/f'negative-kappa-asymptote-size-check.{e}',dpi=170)
with (P/'size-response-table.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['N','kappa','epsilon_over_pi','u','start','vector_max_power','staggered_nematic_power','fourfold_power','vortex_density','mean_link_cosine','blocks','flags'])
 for d in rows:
  p=d['point'];w.writerow([p['N'],p['target_kappa'],p['epsilon_over_pi'],p['scaled_u'],p['initial_label'],peak(d),d['mean']['staggered_nematic_squared'],d['mean']['fourfold_squared'],d['mean']['vortex_density'],d['mean']['mean_link_cosine'],d['blocks'],';'.join(d['flags'])])
