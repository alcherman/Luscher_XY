"""Source-pinned finite-box response maps. No asymptote is fitted."""
from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
P=Path(__file__).resolve().parent;S=P.parent/'scaled-asymptote'
rows=sum([json.loads(f.read_text())['rows'] for f in [S/'size-validation-analysis.json',P/'analysis.json']],[])
rows=[d for d in rows if d['point']['N']==8]
colors=dict(zip([-4,-8,-16,-32,-64],['C0','C1','C2','C3','C4']))
labels=['Uniform','Unit winding','Canted sublattices','Vortex-rich parity-locked Z4'];markers=['o','^','s','D']
def start(d):return int(d['id'].split('_r')[1].split('_')[0])%10
def peak(d):return float(np.mean([x['power'] for x in d['spectrum_samples']]))
def val(d,key):return peak(d) if key=='peak' else d['mean'][key]
fig,axes=plt.subplots(3,2,figsize=(11,10),sharex=True,sharey='row')
for row,(key,label) in enumerate([('peak','Sampled max spin Fourier power'),('staggered_nematic_squared','Staggered nematic power'),('vortex_density','Vortex density')]):
 for k,color in colors.items():
  rr=[d for d in rows if d['point']['target_kappa']==k and d['point']['role']!='off_scaling_energy_sector_control']
  us=sorted({d['point']['scaled_u'] for d in rr});groups=[[val(d,key) for d in rr if d['point']['scaled_u']==u and start(d)!=1] for u in us]
  axes[row,0].plot(us,[np.mean(g) for g in groups],color=color,marker='o',ms=4,label=f'κ={k}')
  axes[row,0].fill_between(us,[min(g) for g in groups],[max(g) for g in groups],color=color,alpha=.14)
  unit=sorted([d for d in rr if start(d)==1],key=lambda d:d['point']['scaled_u']);axes[row,1].plot([d['point']['scaled_u'] for d in unit],[val(d,key) for d in unit],color=color,marker='^',ms=4)
 for ax in axes[row]:ax.grid(alpha=.18)
 axes[row,0].set_ylabel(label)
axes[0,0].set_title('Uniform, canted and Z4 starts: mean and range');axes[0,1].set_title('Unit-winding starts, shown separately')
for ax in axes[-1]:ax.set_xlabel(r'$u=|\kappa|(\epsilon-\pi/2)$')
fig.legend(handles=[Line2D([],[],color=c,marker='o',label=f'κ={k}') for k,c in colors.items()],loc='lower center',bbox_to_anchor=(.5,.055),ncol=5,frameon=False)
fig.suptitle('N=8: stronger-coupling tests of response alignment')
fig.text(.5,.015,'Shading is the range across three initial-state classes, not a confidence interval. Lines guide the eye.\nFinite-box response and persistent trapping; no equilibrium critical curve or fitted asymptote.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.12,1,.96])
for e in ['png','pdf','svg']:fig.savefig(P/f'strong-negative-scaled-response.{e}',dpi=170)
fig,axes=plt.subplots(2,2,figsize=(11,8),sharex=True,sharey=True)
for j,ax in enumerate(axes.ravel()):
 rr=[d for d in rows if start(d)==j]
 ax.axvline(.5,color='.4',ls=':',lw=1)
 sc=ax.scatter([d['point']['epsilon_over_pi'] for d in rr],[d['point']['target_kappa'] for d in rr],c=[peak(d) for d in rr],cmap='viridis',vmin=0,vmax=.6,marker=markers[j],s=55,edgecolors='.3',lw=.4)
 ax.set(title=labels[j],xlabel=r'$\epsilon/\pi$',ylabel=r'$\kappa$',ylim=(-67,-2),yticks=[-64,-32,-16,-8,-4]);ax.grid(alpha=.15)
fig.subplots_adjust(left=.08,right=.86,bottom=.15,top=.88,hspace=.3,wspace=.22);fig.colorbar(sc,cax=fig.add_axes([.89,.27,.02,.5]),label='Sampled max spin Fourier power')
fig.suptitle('3D negative κ: measured N=8 response probes')
fig.text(.5,.025,'The fixed-u scan approaches½ by design; evidence comes from response alignment, not point coordinates.\nFixed ε/π=.60 controls included. Initial-state dependence remains; this is not an accepted phase boundary.',ha='center',fontsize=9)
for e in ['png','pdf','svg']:fig.savefig(P/f'negative-kappa-asymptote-response-map.{e}',dpi=170)
with (P/'N8-response-table.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['N','kappa','epsilon_over_pi','scaled_u','start','vector_max','staggered_nematic','fourfold','vortex','mean_link_cosine','blocks','flags'])
 for d in rows:
  p=d['point'];w.writerow([8,p['target_kappa'],p['epsilon_over_pi'],p['scaled_u'],p['initial_label'],peak(d),*[d['mean'][x] for x in ['staggered_nematic_squared','fourfold_squared','vortex_density','mean_link_cosine']],d['blocks'],';'.join(d['flags'])])
