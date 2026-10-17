"""Evidence and an explicitly imposed half-angle asymptote, shown separately."""
from pathlib import Path
import json,hashlib,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
P=Path(__file__).resolve().parent;ROOT=P.parents[1]
source=ROOT/'experiments/2026-10-13-xy-infinite-volume/d2-d3-d4-attained-N512-with-endpoints-20261014/analysis.json'
rows=[r for r in json.loads(source.read_text())['rows'] if r['d']==3]
q0=next(r['epsilon_over_pi'] for r in rows if r['beta']==0)
scouts=[dict(kappa=-.1,low=.748,high=.798,kind='N8/N12 short-scout response region; N12 uncertainty unresolved'),dict(kappa=-.25,low=.65,high=.69,kind='N12 cold-start response-change interval; sector trapping'),dict(kappa=-.5,low=.57,high=.62,kind='N12 cold-start response-change interval; sector trapping')]
def hypothesis(k):
    t=-np.asarray(k)
    return .5+(q0-.5)/(1+.8*t+7.2*t*t)
fig,(ax,bx)=plt.subplots(1,2,figsize=(12,6.6),gridspec_kw={'width_ratios':[1.2,1]})
for r in rows:
    x,y=r['epsilon_over_pi'],r['beta']
    if r['error_axis']=='vertical':
        ax.vlines(x,r['sensitivity_low'],r['sensitivity_high'],color='C0',alpha=.25,lw=5,zorder=2)
        ax.vlines(x,r['MC_low'],r['MC_high'],color='C0',lw=1.5,zorder=3)
    else:
        ax.hlines(y,r['sensitivity_low'],r['sensitivity_high'],color='C0',alpha=.25,lw=5,zorder=2)
        ax.hlines(y,r['MC_low'],r['MC_high'],color='C0',lw=1.5,zorder=3)
    ax.plot(x,y,'o',color='C0',ms=4,zorder=4)
for a in [ax,bx]:
    a.axvline(.5,color='black',lw=1,ls=':',zorder=1)
    a.axhline(0,color='.65',lw=.7,zorder=0)
    for s in scouts:
        y=s['kappa'];a.hlines(y,s['low'],s['high'],color='C1',lw=5,zorder=5)
        a.vlines([s['low'],s['high']],y-.012,y+.012,color='C1',lw=1.3,zorder=5)
    a.set_xlabel(r'$\epsilon/\pi$',fontsize=13);a.set_ylabel(r'$\kappa$',fontsize=13);a.grid(alpha=.14)
# One schematic curve, not fitted: its asymptote is imposed by its formula.
k=np.linspace(-2.3,0,500);bx.plot(hypothesis(k),k,color='.55',ls='--',lw=1.6,zorder=1)
bx.plot(q0,0,'o',color='C0',ms=5)
ax.set(xlim=(.475,1.02),ylim=(-.64,.515),title='Current evidence: MC points and scout regions')
bx.set(xlim=(.475,.845),ylim=(-2.2,.10),title='Illustrative hypothesis — not fitted')
bx.axhspan(-2.2,-.5,color='.94',zorder=-2)
bx.text(.68,-1.35,'NO MC DATA\nbelow κ = −0.5',ha='center',color='.35',fontsize=10)
bx.text(.519,-2.05,r'$\epsilon/\pi=1/2$',ha='left',fontsize=10)
bx.text(.62,-.91,'Dashed curve: imposed hypothesis\n—not a fit or prediction',ha='left',fontsize=9,color='.35')
handles=[Line2D([],[],marker='o',ls='none',color='C0',label='Attained positive-κ MC estimates'),Line2D([],[],lw=5,color='C1',label='Negative-κ scout regions — not equilibrium intervals'),Line2D([],[],ls='--',color='.55',label='Illustrative continuation forced toward ½')]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.085),ncol=1,frameon=False,fontsize=9)
fig.suptitle('3D constrained XY: continuation into negative κ',fontsize=16,y=.98)
fig.text(.5,.025,'Negative-κ scouts use N=8,12 and 4,000 Metropolis sweeps. Sector trapping prevents equilibrium transition estimates.\nPositive-κ dark bars: conditional MC intervals; pale bars: descriptive sensitivity unions, without certified total coverage.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.205,1,.95])
for ext in ['png','pdf','svg']:fig.savefig(P/f'kappa-epsilon-continuation-overview.{ext}',dpi=180)
metadata=dict(positive_source=str(source),positive_source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),negative_source=str(P/'analysis.json'),negative_source_sha256=hashlib.sha256((P/'analysis.json').read_bytes()).hexdigest(),scouts=scouts,hypothesis=dict(formula='q=0.5+(q0-0.5)/(1+0.8*t+7.2*t*t), t=-kappa',q0=q0,parameters='Chosen for a schematic illustration, not statistically fitted',limit='0.5 imposed, not inferred'),conclusion='Data are compatible with a continuation toward smaller epsilon; the limit1/2 is not established by these scouts.',no_new_simulations=True)
(P/'kappa-epsilon-continuation-overview.json').write_text(json.dumps(metadata,indent=2)+'\n')
with (P/'kappa-epsilon-continuation-overview.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['role','epsilon_over_pi','kappa','low','high','error_axis'])
 for r in rows:w.writerow(['attained_MC',r['epsilon_over_pi'],r['beta'],r['MC_low'],r['MC_high'],r['error_axis']])
 for s in scouts:w.writerow(['descriptive_scout_region',None,s['kappa'],s['low'],s['high'],'horizontal_non_CI'])
