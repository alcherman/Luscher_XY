"""Measured fixed-box continuation; never impose the epsilon/pi=1/2 limit."""
from pathlib import Path
import json,hashlib,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
P=Path(__file__).resolve().parent
SOURCE=P.parents[1]/'2026-10-13-xy-infinite-volume/d2-d3-d4-attained-N512-with-endpoints-20261014/analysis.json'
a=json.loads((P/'analysis.json').read_text());positive=[v for v in json.loads(SOURCE.read_text())['rows'] if v['d']==3]
fig,(ax,bx)=plt.subplots(1,2,figsize=(11.5,6.4))
for panel in [ax,bx]:
 panel.axvline(.5,ls=':',color='.4',lw=1,label=r'$\epsilon/\pi=1/2$ reference')
 panel.axhline(0,color='.65',lw=.8);panel.grid(alpha=.15)
 panel.set_xlabel(r'$\epsilon/\pi$');panel.set_ylabel(r'$\kappa$')
for r in positive:
 for panel in [ax,bx]:
  q,k=r['epsilon_over_pi'],r['beta']
  if r['error_axis']=='vertical':
   panel.vlines(q,r['sensitivity_low'],r['sensitivity_high'],color='C0',alpha=.2,lw=5)
   panel.vlines(q,r['MC_low'],r['MC_high'],color='C0',lw=1.3)
  else:
   panel.hlines(k,r['sensitivity_low'],r['sensitivity_high'],color='C0',alpha=.2,lw=5)
   panel.hlines(k,r['MC_low'],r['MC_high'],color='C0',lw=1.3)
  panel.plot(q,k,'o',color='C0',ms=4)
bands=json.loads((P/'response-bands.json').read_text())['bands']
bands=[r for r in bands if r.get('role','production')=='production']
accepted=[];flagged=[]
for r in bands:
 if not r['resolved_response_band']:
  for panel in [ax,bx]:
   q=r['epsilon_over_pi'];panel.fill_betweenx([-2,0],q-.003,q+.003,color='.82',alpha=.6,hatch='///',edgecolor='.6',linewidth=.5)
  continue
 for panel in [ax,bx]:
  color='C1' if r['N']==8 else 'C2';q=r['epsilon_over_pi'];lo,hi=r['kappa_low'],r['kappa_high']
  panel.vlines(q,lo,hi,color=color,lw=5,alpha=.65)
  panel.hlines([lo,hi],q-.003,q+.003,color=color,lw=1.2)
bx.text(.527,-1.75,'Unresolved\nsector trapping',fontsize=8,color='.35',ha='left')
ax.set(xlim=(.47,1.02),ylim=(-2.05,.52),title='3D: attained MC and negative-κ reconnaissance')
bx.set(xlim=(.47,.83),ylim=(-2.05,.07),title='Negative-κ: observed response bands (N = 8)')
fig.legend(handles=[Line2D([],[],color='C0',marker='o',ls='none',label='Attained positive-κ / endpoint'),Line2D([],[],color='C1',marker='s',ls='none',label='N = 8 response band (not a CI)'),Line2D([],[],color='.4',ls=':',label='½ reference; grey scan unresolved')],loc='lower center',bbox_to_anchor=(.5,.115),ncol=3,frameon=False,fontsize=9)
fig.text(.5,.025,'Orange bars span large/small uniform-magnetization responses agreeing between both starts and halves.\nEquilibrium transition locations remain unresolved where ladders do not traverse; bars are NOT confidence intervals.\nBlue bars: earlier extrapolated conditional MC and pale descriptive sensitivity. No asymptote imposed.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.17,1,1])
for ext in ['png','pdf','svg']:fig.savefig(P/f'measured-negative-kappa-continuation.{ext}',dpi=180)
(P/'plot-provenance.json').write_text(json.dumps(dict(analysis_sha256=hashlib.sha256((P/'analysis.json').read_bytes()).hexdigest(),positive_source=str(SOURCE),positive_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),response_bands=bands,response_bands_sha256=hashlib.sha256((P/'response-bands.json').read_bytes()).hexdigest(),interpolated_roots_rejected_as_headline=True,asymptote_imposed=False),indent=2)+'\n')
with (P/'measured-negative-kappa-continuation.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['N','epsilon_over_pi','kappa_low','kappa_high','observed_response_band','equilibrium_location_resolved'])
 for r in bands:w.writerow([r['N'],r['epsilon_over_pi'],r.get('kappa_low'),r.get('kappa_high'),r['resolved_response_band'],False])
