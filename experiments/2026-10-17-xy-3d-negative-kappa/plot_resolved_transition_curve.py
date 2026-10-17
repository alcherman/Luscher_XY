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
 if abs(q-.52)<1e-10:continue
 if r['resolved_response_band']:
  lo,hi=r['kappa_low'],r['kappa_high'];ax.vlines(q,lo,hi,color='C1',lw=6,alpha=.7);ax.hlines([lo,hi],q-.004,q+.004,color='C1',lw=1.2)
 else:
  ax.fill_betweenx([-2,0],q-.003,q+.003,color='.8',alpha=.5,hatch='///',edgecolor='.6')
for point in json.loads((P/'exploratory-transition-points.json').read_text())['points']:
 q,k=point['epsilon_over_pi'],point['kappa'];ax.vlines(q,point['kappa_low'],point['kappa_high'],color='C4',lw=6);ax.plot(q,k,'D',color='C4',ms=5)
ax.annotate('ε/π = 0.52: bracket resolved',xy=(.52,-2.05),xytext=(.61,-1.88),arrowprops=dict(arrowstyle='->',color='C4'),color='C4',fontsize=10)

ax.set(xlim=(.48,1.015),ylim=(-2.25,.54),xlabel=r'$\epsilon/\pi$',ylabel=r'$\kappa$',title='3D: transition curve and exploratory negative-κ continuation')
ax.grid(alpha=.15)
ax.legend(handles=[Line2D([],[],color='C0',marker='o',ls='',label='Extrapolated MC, κ ≥ 0'),Line2D([],[],color='C1',marker='|',ms=13,mew=3,ls='',label='N=8 response intervals, κ < 0'),Line2D([],[],color='C4',marker='D',ls='',label='6³ / 8³ / 12³ common bracket'),Line2D([],[],color='.4',ls=':',label='ε/π = ½ reference')],loc='upper left',frameon=True,fontsize=9)
fig.text(.5,.035,'Purple: all four starts and both chain halves bracket the vector-response loss on three sizes.\nNegative-κ bars are descriptive finite-box intervals, not certified confidence intervals. No asymptote is fitted.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.1,1,1])
for ext in ['png','pdf','svg']:fig.savefig(P/f'transition-curve-resolved-q052.{ext}',dpi=190)
(P/'transition-curve-resolved-q052-provenance.json').write_text(json.dumps(dict(positive_source=str(positive),negative_source=str(negative),source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [positive,negative]},positive_count=len(rows),observed_negative_intervals=sum(r['resolved_response_band'] for r in bands),new_transition_inference=False,deeper_negative_scans_excluded_because_no_reliable_transition_locations=True),indent=2)+'\n')
