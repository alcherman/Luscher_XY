"""Independent-start and split-history evidence behind the response bands."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent;a=json.loads((P/'analysis.json').read_text());groups={}
for r in a['rows']:
 if r['point'].get('role') in ['production','fine_transport_validation']:groups.setdefault((r['point']['epsilon_over_pi'],r['point']['role']),[]).append(r)
fig,axes=plt.subplots(2,4,figsize=(13,7),squeeze=False)
for ax,((q,role),rr) in zip(axes.ravel(),sorted(groups.items())):
 for r in rr:
  k=r['point']['kappas'];zero=r['point']['initial_winding'][0]==0;c='C0' if zero else 'C1';m=np.asarray(r['half_m2']);ax.plot(k,r['m2'],'.-',color=c,lw=1,ms=3,label='Zero start' if zero else 'Winding start');ax.fill_between(k,m.min(0),m.max(0),color=c,alpha=.2)
 ax.set_title(f'q={q:.2f}'+(' dense ladder' if role=='fine_transport_validation' else ''))
 ax.set_xlabel('κ');ax.set_ylabel(r'$\langle m^2\rangle$');ax.set_ylim(-.01,.73);ax.grid(alpha=.2)
 ax.text(.03,.96,'Round trips: '+', '.join(str(r['walker_round_trips']['round_trips']) for r in rr),transform=ax.transAxes,va='top',fontsize=8)
for ax in axes.ravel()[len(groups):]:ax.axis('off')
axes[0,0].legend(fontsize=8,loc='lower left')
fig.suptitle('N = 8: mean uniform magnetization and two-half ranges',fontsize=14)
fig.text(.5,.015,'Each curve is one independent replica-exchange ensemble. Shading spans its two post-warmup half-means; it is not a confidence band.\nAgreement of these diagnostics does not establish equilibrium; unresolved transport and metastability remain explicit.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.065,1,.96])
for ext in ['png','pdf','svg']:fig.savefig(P/f'independent-start-response-checks.{ext}',dpi=170)
