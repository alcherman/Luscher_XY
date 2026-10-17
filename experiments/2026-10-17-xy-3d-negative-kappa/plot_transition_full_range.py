"""Presentation-only variant of the reviewed data: full 0..1 axis, no inset."""
from pathlib import Path
import json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
P=Path(__file__).resolve().parent
source=P/'exploratory-negative-kappa-delivery.json';delivery=json.loads(source.read_text());points=delivery['points']
positive_source=Path(delivery['positive_source']);negative_source=Path(delivery['older_negative_source'])
assert hashlib.sha256(positive_source.read_bytes()).hexdigest()==delivery['positive_sha256']
assert hashlib.sha256(negative_source.read_bytes()).hexdigest()==delivery['older_negative_sha256']
positive=[r for r in json.loads(positive_source.read_text())['rows'] if r['d']==3]
old=[r for r in json.loads(negative_source.read_text())['bands'] if r.get('role','production')=='production' and r['resolved_response_band'] and all(abs(r['epsilon_over_pi']-v['epsilon_over_pi'])>1e-10 for v in points)]
fig,ax=plt.subplots(figsize=(9,7))
ax.axvline(.5,ls=':',color='.4',lw=1.2);ax.axhline(0,color='.65',lw=.8);ax.grid(alpha=.15)
for r in positive:
 q,k=r['epsilon_over_pi'],r['beta']
 if r['error_axis']=='vertical':ax.vlines(q,r['sensitivity_low'],r['sensitivity_high'],color='C0',lw=5,alpha=.2);ax.vlines(q,r['MC_low'],r['MC_high'],color='C0',lw=1.2)
 else:ax.hlines(k,r['sensitivity_low'],r['sensitivity_high'],color='C0',lw=5,alpha=.2);ax.hlines(k,r['MC_low'],r['MC_high'],color='C0',lw=1.2)
 ax.plot(q,k,'o',ms=4,color='C0',clip_on=False)
for r in old:
 q=r['epsilon_over_pi'];ax.vlines(q,r['kappa_low'],r['kappa_high'],color='C1',lw=5,alpha=.7);ax.hlines([r['kappa_low'],r['kappa_high']],q-.003,q+.003,color='C1',lw=1.)
for r in points:
 q=r['epsilon_over_pi'];ax.vlines(q,r['kappa_low'],r['kappa_high'],color='C4',lw=5);ax.hlines([r['kappa_low'],r['kappa_high']],q-.003,q+.003,color='C4',lw=1.2);ax.plot(q,r['kappa'],'D',color='C4',ms=5)
ax.set(xlim=(0,1),ylim=(min(r['kappa_low'] for r in points)-.6,.65),xlabel=r'$\epsilon/\pi$',ylabel=r'$\kappa$',title='3D: exploratory continuation toward negative κ')
ax.set_xticks([i/10 for i in range(11)])
ax.text(.04,.67,'Dotted line: ε/π = ½ reference\nNo asymptote imposed or fitted',transform=ax.transAxes,color='.35',fontsize=10)
fig.legend(handles=[Line2D([],[],marker='o',color='C0',ls='',label='Earlier extrapolated MC, κ ≥ 0'),Line2D([],[],marker='|',color='C1',ls='',ms=12,mew=3,label='Earlier N=8 response intervals'),Line2D([],[],marker='D',color='C4',ls='',label='New all-start/half, multi-size brackets')],loc='lower center',bbox_to_anchor=(.5,.07),ncol=1,frameon=False,fontsize=9)
fig.text(.5,.015,'Negative-κ bars are descriptive finite-box ranges, not confidence intervals or certified thermodynamic errors.\nStrong staggered-nematic response remains on the spin-suppressed side. No MC connecting curve is drawn.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.19,1,1])
assert tuple(ax.get_xlim())==(0.,1.) and len(fig.axes)==1
base='exploratory-negative-kappa-transition-curve-full-range'
for ext in ['png','pdf','svg']:fig.savefig(P/f'{base}.{ext}',dpi=190)
(P/f'{base}.json').write_text(json.dumps(dict(source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),xlim=[0,1],inset=False,positive_points=len(positive),older_negative_intervals=len(old),new_brackets=len(points)),indent=2)+'\n')
print(base)
