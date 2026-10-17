"""Transition-location summary from existing reviewed sources; no new root inference."""
from pathlib import Path
import json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
P=Path(__file__).resolve().parent
positive=P.parent/'2026-10-13-xy-infinite-volume/d2-d3-d4-attained-N512-with-endpoints-20261014/analysis.json'
negative=P/'replica-exchange/response-bands.json'
rows=[r for r in json.loads(positive.read_text())['rows'] if r['d']==3]
bands=[r for r in json.loads(negative.read_text())['bands'] if r.get('role','production')=='production']
fig,ax=plt.subplots(figsize=(8,6))
ax.axvline(.5,color='.4',ls=':',lw=1.3);ax.axhline(0,color='.6',lw=.8)
for r in rows:
 q,k=r['epsilon_over_pi'],r['beta']
 if r['error_axis']=='vertical':
  ax.vlines(q,r['sensitivity_low'],r['sensitivity_high'],color='C0',lw=6,alpha=.2);ax.vlines(q,r['MC_low'],r['MC_high'],color='C0',lw=1.4)
 else:
  ax.hlines(k,r['sensitivity_low'],r['sensitivity_high'],color='C0',lw=6,alpha=.2);ax.hlines(k,r['MC_low'],r['MC_high'],color='C0',lw=1.4)
 ax.plot(q,k,'o',color='C0',ms=5)
for r in bands:
 q=r['epsilon_over_pi']
 if r['resolved_response_band']:
  lo,hi=r['kappa_low'],r['kappa_high'];ax.vlines(q,lo,hi,color='C1',lw=6,alpha=.7);ax.hlines([lo,hi],q-.004,q+.004,color='C1',lw=1.2)
 else:
  ax.fill_betweenx([-2,0],q-.003,q+.003,color='.8',alpha=.5,hatch='///',edgecolor='.6')
ax.text(.531,-1.65,'Unresolved\n(start trapping)',fontsize=10,color='.4')
ax.text(.69,-1.7,'Deeper scans to κ = −64:\nno reliable transition locations yet',fontsize=10,color='.35')
ax.set(xlim=(.48,1.015),ylim=(-2.03,.54),xlabel=r'$\epsilon/\pi$',ylabel=r'$\kappa$',title='3D: transition curve and exploratory negative-κ continuation')
ax.grid(alpha=.15)
ax.legend(handles=[Line2D([],[],color='C0',marker='o',ls='',label='Extrapolated MC, κ ≥ 0'),Line2D([],[],color='C1',marker='|',ms=13,mew=3,ls='',label='N=8 response intervals, κ < 0'),Line2D([],[],color='.4',ls=':',label='ε/π = ½ reference')],loc='upper left',frameon=True,fontsize=9)
fig.text(.5,.035,'Orange intervals are exploratory response changes, not confidence intervals or validated critical couplings.\nBlue: conditional MC bars and pale descriptive sensitivity ranges. No asymptote is fitted.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.1,1,1])
for ext in ['png','pdf','svg']:fig.savefig(P/f'transition-curve-summary.{ext}',dpi=190)
(P/'transition-curve-summary-provenance.json').write_text(json.dumps(dict(positive_source=str(positive),negative_source=str(negative),source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [positive,negative]},positive_count=len(rows),observed_negative_intervals=sum(r['resolved_response_band'] for r in bands),new_transition_inference=False,deeper_negative_scans_excluded_because_no_reliable_transition_locations=True),indent=2)+'\n')
