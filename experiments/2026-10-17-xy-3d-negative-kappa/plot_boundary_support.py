"""Native measured Binder responses and split-half ranges; no fitted phase curve."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
P=Path(__file__).resolve().parent;B=P/'boundary-resolution'
fig,axes=plt.subplots(1,3,figsize=(13,4.7))
for ax,suffix in zip(axes,['','-q0510000','-q0505000']):
 a=json.loads(Path(json.loads((B/f'latest-exploratory-brackets{suffix}.json').read_text())['analysis']).read_text());n=max(c['N'] for c in a['common_anchors']);cc=[c for c in a['coordinates'] if c['N']==n];labels=list(dict.fromkeys(c['point']['initial_label'] for c in a['chains']))
 br=next(c for c in a['common_anchors'] if c['N']==n);ax.axvspan(br['lower_kappa'],br['upper_kappa'],color='C4',alpha=.12)
 for i,label in enumerate(labels):
  xs=np.array([c['kappa'] for c in cc]);ys=np.array([c['variants'][label]['Binder_U'] for c in cc]);halves=np.array([[c['variants'][f'{label}:half{h}']['Binder_U'] for h in [1,2]] for c in cc]);ax.vlines(xs,halves.min(1),halves.max(1),color=f'C{i}',alpha=.65,lw=2);ax.plot(xs,ys,['o','s','^','D'][i],color=f'C{i}',ms=4,label=label)
 ax.axhline(.25,color='.3',ls='--',lw=1);ax.axhspan(.2,.3,color='.5',alpha=.08);ax.grid(alpha=.15);ax.set(xlabel=r'$\kappa$',title=f"ε/π = {a['epsilon_over_pi']}, {n}³",ylabel=r'$U=1-\langle m^4\rangle/(2\langle m^2\rangle^2)$')
axes[0].legend(loc='best',fontsize=7);fig.suptitle('Measured response: four starts; vertical spans show split-half variation',fontsize=12);fig.text(.5,.015,'Spans are descriptive half-to-half variation, not statistical error bars. Purple windows are measured common brackets.\nThe q=0.505 cohort uses exact local heatbath; earlier Metropolis cohorts remain preserved separately.',ha='center',fontsize=9);fig.tight_layout(rect=[0,.11,1,.95])
for ext in ['png','pdf','svg']:fig.savefig(P/f'exploratory-boundary-start-half-support.{ext}',dpi=180)
