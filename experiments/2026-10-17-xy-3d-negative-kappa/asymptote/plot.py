"""Strong-negative discovery map: retain all initial-condition differences."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
P=Path(__file__).resolve().parent;a=json.loads((P/'analysis.json').read_text());fig,axes=plt.subplots(2,3,figsize=(12,8),sharex=True,sharey=True)
styles={70:('o',-.0015,'Uniform start'),71:('^',-.0005,'Unit winding'),72:('s',.0005,'Maximal multi-axis winding'),73:('D',.0015,'Canted sublattices')}
for row,n in enumerate([6,8]):
 for col,field in enumerate(['magnetization_squared','staggered_squared','nonzero_peak']):
  ax=axes[row,col];ax.axvline(.5,ls=':',color='.5',lw=.8)
  for d in a['rows']:
   p=d['point']
   if p['N']!=n:continue
   marker,dx,_=styles[int(p['id'].split('_r')[1].split('_')[0])];q=p['epsilon_over_pi'];value=max(np.array(d['average_spectrum']).ravel()[1:]) if field=='nonzero_peak' else d['mean'][field]
   sc=ax.scatter(q+dx,p['target_kappa'],c=[value],marker=marker,s=65,cmap='viridis',vmin=0,vmax=1,edgecolor='.25',linewidth=.5)
  ax.set_title(f'N={n}: '+{'magnetization_squared':r'$\langle m_0^2\rangle$','staggered_squared':r'$\langle m_\pi^2\rangle$','nonzero_peak':'Max nonzero mean Fourier power'}[field]);ax.set(xlim=(.47,.615),ylim=(-8.8,-1.2),yticks=[-8,-4,-2]);ax.grid(alpha=.15)
  if row==1:ax.set_xlabel(r'$\epsilon/\pi$')
  if col==0:ax.set_ylabel(r'$\kappa$')
fig.subplots_adjust(left=.075,right=.90,bottom=.18,top=.89,hspace=.25,wspace=.16);cax=fig.add_axes([.92,.25,.018,.55]);fig.colorbar(sc,cax=cax,label='Normalized power')
fig.legend(handles=[Line2D([],[],marker=m,ls='none',color='.3',label=label) for m,dx,label in styles.values()],loc='lower center',bbox_to_anchor=(.5,.095),ncol=4,frameon=False,fontsize=9)
fig.suptitle('3D negative-κ discovery: uniform, staggered and competing modes',fontsize=14)
fig.text(.5,.025,'Four independently seeded initial conditions at each coordinate; markers offset slightly in ε/π for visibility.\nOne trajectory per initial-condition class. Agreement is not an equilibrium or infinite-volume certificate; no transition curve is fitted.',ha='center',fontsize=9)
for ext in ['png','pdf','svg']:fig.savefig(P/f'strong-negative-discovery-map.{ext}',dpi=170)
