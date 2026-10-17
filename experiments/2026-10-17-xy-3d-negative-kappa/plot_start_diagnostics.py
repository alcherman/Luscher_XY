"""Plot independent start sectors without pooling trapped histories."""
import json,hashlib
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent;d=json.loads((P/'analysis.json').read_text());fig,axes=plt.subplots(2,2,figsize=(10,7),sharey=True)
for ax,k in zip(axes.flat,[0.,-.1,-.25,-.5]):
 for n,marker in [(8,'o'),(12,'s')]:
  for w,color in [(0,'C0'),(1,'C1')]:
   rr=sorted([r for r in d['rows'] if r['point']['kappa']==k and r['point']['N']==n and r['point']['initial_winding'][0]==w],key=lambda r:r['point']['epsilon_over_pi'])
   ax.plot([r['point']['epsilon_over_pi'] for r in rr],[r['binder'] for r in rr],marker=marker,color=color,ls='none',mfc='none' if n==8 else color,label=f'N={n}, initial winding {w}')
 ax.axhline(.25,color='k',lw=.7,ls='--');ax.set(title=f'κ = {k:g}',xlabel='ε/π',ylabel='Binder U');ax.grid(alpha=.2);ax.set_ylim(-.3,.54)
axes[0,0].legend(fontsize=8)
fig.suptitle('3D continuation scout: starts shown separately',fontsize=14)
fig.text(.5,.015,'4,000 Metropolis sweeps; first 1,000 excluded. Points are descriptive means, without error bars.\nPersistent winding-sector disagreement prevents an equilibrium phase-boundary claim.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.07,1,.95])
for ext in ['png','pdf','svg']:fig.savefig(P/f'individual-start-binder-diagnostics.{ext}',dpi=170)
summary=dict(scope='Existing ferromagnetic branch continued below kappa=0; no epsilon-near-pi scan.',completed_full_audited_chains=72,N=[8,12],sweeps_per_chain=4000,warmup_label=1000,max_workers=8,old_spin_campaign_remains_paused=True,observations=[dict(kappa=-.1,qualitative_location_epsilon_over_pi='around0.77',qualification='Short finite-volume measurements; N12 bootstrap interval unresolved.'),dict(kappa=-.25,cold_start_loss_of_order_bracket_epsilon_over_pi=[.65,.69],qualification='Descriptive cold-start bracket only; twisted starts trapped in a different sector; not equilibrium uncertainty.'),dict(kappa=-.5,cold_start_loss_of_order_bracket_epsilon_over_pi=[.57,.62],qualification='Descriptive cold-start bracket only; twisted starts trapped in a different sector; not equilibrium uncertainty.')],headline='Short scouts are consistent with continuation toward lower epsilon, but do not determine an equilibrium negative-kappa phase curve.',next_required='Resolve sector mixing or use a validated negative-kappa-compatible nonlocal sampler before extrapolating the negative branch. Do not assume simply longer Metropolis chains settle this.',analysis_sha256=hashlib.sha256((P/'analysis.json').read_bytes()).hexdigest())
(P/'continuation-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
