"""Reviewable exploratory continuation from common measured side anchors."""
from pathlib import Path
import json,hashlib,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
P=Path(__file__).resolve().parent;B=P/'boundary-resolution'
inputs=[B/'latest-exploratory-brackets.json',B/'latest-exploratory-brackets-q0510000.json',B/'latest-exploratory-brackets-q0505000.json']
points=[];analyses=[]
for f in inputs:
 a=json.loads(Path(json.loads(f.read_text())['analysis']).read_text());analyses.append(a)
 if not a['complete'] or any(v['status']!='all_start_and_half_anchor_bracket' for v in a['common_anchors']) or any(v['status']!='unique' for v in a['local_families']):raise RuntimeError('Unresolved start/half/size/root family; continue campaign rather than drop it')
 brackets=[v['bracket'] for fam in a['local_families'] for v in fam['roots']];lo=min([v['lower_kappa'] for v in a['common_anchors']]+[v[0] for v in brackets]);hi=max([v['upper_kappa'] for v in a['common_anchors']]+[v[1] for v in brackets]);q=a['epsilon_over_pi'];k=a['primary']['roots'][0]['root']
 flags={}
 for chain in a['chains']:
  for flag in chain['sampling_diagnostics'].get('flags',[]):flags[flag]=flags.get(flag,0)+1
 points.append(dict(sampling_flag_counts=flags,diagnostic_pending=sum(bool(c['sampling_diagnostics'].get('pending')) for c in a['chains']),epsilon_over_pi=q,kappa=k,kappa_low=lo,kappa_high=hi,sizes=[d['N'] for d in a['common_anchors']],trajectories=a['trajectories'],global_nonunique_families=sum(v['status']!='unique' for v in a['families']),source=json.loads(f.read_text())['analysis'],source_sha256=hashlib.sha256(Path(json.loads(f.read_text())['analysis']).read_bytes()).hexdigest(),scaled_coefficient_range=[(q-.5)*abs(hi),(q-.5)*abs(lo)],scope='Marker=largest-box pooled U=.25 linear interpolation. Bar=observed size/start/half/threshold/interpolation bracket union, not certifiedCI or infinitevolume uncertainty.'))
positive_source=P.parent/'2026-10-13-xy-infinite-volume/d2-d3-d4-attained-N512-with-endpoints-20261014/analysis.json';positive=[r for r in json.loads(positive_source.read_text())['rows'] if r['d']==3]
negative_source=P/'replica-exchange/response-bands.json';old=[r for r in json.loads(negative_source.read_text())['bands'] if r.get('role','production')=='production' and r['resolved_response_band'] and all(abs(r['epsilon_over_pi']-v['epsilon_over_pi'])>1e-10 for v in points)]
fig,ax=plt.subplots(figsize=(9,7));zoom=ax.inset_axes([.42,.15,.55,.43])
for panel in [ax,zoom]:
 panel.axvline(.5,ls=':',color='.4',lw=1.2);panel.axhline(0,color='.65',lw=.8);panel.grid(alpha=.15)
 for r in positive:
  q,k=r['epsilon_over_pi'],r['beta']
  if r['error_axis']=='vertical':panel.vlines(q,r['sensitivity_low'],r['sensitivity_high'],color='C0',lw=5,alpha=.2);panel.vlines(q,r['MC_low'],r['MC_high'],color='C0',lw=1.2)
  else:panel.hlines(k,r['sensitivity_low'],r['sensitivity_high'],color='C0',lw=5,alpha=.2);panel.hlines(k,r['MC_low'],r['MC_high'],color='C0',lw=1.2)
  panel.plot(q,k,'o',ms=4,color='C0')
 for r in old:
  q=r['epsilon_over_pi'];panel.vlines(q,r['kappa_low'],r['kappa_high'],color='C1',lw=5,alpha=.7);panel.hlines([r['kappa_low'],r['kappa_high']],q-.003,q+.003,color='C1',lw=1.)
 for r in points:
  q=r['epsilon_over_pi'];panel.vlines(q,r['kappa_low'],r['kappa_high'],color='C4',lw=5);panel.hlines([r['kappa_low'],r['kappa_high']],q-.003,q+.003,color='C4',lw=1.2);panel.plot(q,r['kappa'],'D',color='C4',ms=5)
ax.set(xlim=(.492,1.015),ylim=(min(r['kappa_low'] for r in points)-.6,.65),xlabel=r'$\epsilon/\pi$',ylabel=r'$\kappa$',title='3D: exploratory continuation toward negative κ')
zoom.set(xlim=(.495,1.01),ylim=(-2.25,.55),title='Zoom near κ = 0');zoom.tick_params(labelsize=8);zoom.title.set_fontsize(10)
ax.text(.54,.65,'Dotted line: ε/π = ½ reference\nNo asymptote imposed or fitted',transform=ax.transAxes,color='.35',fontsize=10)
fig.legend(handles=[Line2D([],[],marker='o',color='C0',ls='',label='Earlier extrapolated MC, κ ≥ 0'),Line2D([],[],marker='|',color='C1',ls='',ms=12,mew=3,label='Earlier N=8 response intervals'),Line2D([],[],marker='D',color='C4',ls='',label='New all-start/half, multi-size brackets')],loc='lower center',bbox_to_anchor=(.5,.07),ncol=1,frameon=False,fontsize=9)
fig.text(.5,.015,'Negative-κ bars are descriptive finite-box ranges, not confidence intervals or certified thermodynamic errors.\nStrong staggered-nematic response remains on the spin-suppressed side. No MC connecting curve is drawn.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.19,1,1])
for ext in ['png','pdf','svg']:fig.savefig(P/f'exploratory-negative-kappa-transition-curve.{ext}',dpi=190)
with (P/'exploratory-negative-kappa-transition-points.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['epsilon_over_pi','kappa','kappa_low','kappa_high','sizes','trajectories']);w.writeheader();w.writerows([{k:r[k] for k in w.fieldnames} for r in points])
old_metropolis=json.loads((B/'q0505-metropolis-cohort-diagnostic.json').read_text())
result=dict(previous_q0505_metropolis_cohort=old_metropolis,points=points,positive_source=str(positive_source),positive_sha256=hashlib.sha256(positive_source.read_bytes()).hexdigest(),older_negative_source=str(negative_source),older_negative_sha256=hashlib.sha256(negative_source.read_bytes()).hexdigest(),asymptote_fitted=False,root_MC_confidence_intervals_computed=False)
(P/'exploratory-negative-kappa-delivery.json').write_text(json.dumps(result,indent=2)+'\n')
lines=['3D negative-kappa continuation: exploratory resolution', '',
'The previously unresolved epsilon/pi=0.52 scan stopped at kappa=-2, just before the spin-response loss. Extending the range and checking four initial states on several sizes resolves a finite-box bracket. Further measurements at epsilon/pi=0.51 and0.505 test the approach toward one-half.', '',
'Measured continuation',
'Markers are the largest-box pooled U=1-<m^4>/(2<m^2>^2)=1/4 linear-interpolation locators. Bars span the native measured bracket endpoints across every size, all four initial states, both halves, U targets0.20/0.25/0.30 and linear/PCHIP interpolation. Every required local-window family must be identified before publication. All global scan families and outliers remain in the ledger; they are not discarded. The bars are descriptive ranges, not confidence intervals or certified infinite-volume uncertainties.', '']
for r,a in zip(points,analyses):
 lines += [f"epsilon/pi={r['epsilon_over_pi']}: locator kappa={r['kappa']:.4f}; measured range [{r['kappa_low']:.3f}, {r['kappa_high']:.3f}]; sizes {r['sizes']}; {r['trajectories']} audited trajectories.",
 f"Sampling flags across {a['trajectories']} source chains: {r['sampling_flag_counts']}; missing diagnostics: {r['diagnostic_pending']}.",
 f"Common-start/half anchors by size: {a['common_anchors']}",
 f"Central pooled xi/L pair crossings: {a['xi_pair_crossings']}",
 f"Compatible abs(kappa)*(epsilon/pi-0.5) range: {r['scaled_coefficient_range']}", '']
lines += ['Interpretation',
'The response-loss windows are compatible with a1/abs(kappa) approach toward epsilon/pi=1/2 over the tested range. The constant-product comparison uses measured classification on both sides, rather than simply plotting a grid designed to approach one-half. No asymptote is imposed or fitted. These finite-size data do not prove the infinite-coupling or infinite-volume limit, exclude additional branches or reentrance, or determine a universal coefficient.', '',
'Sampling and preserved records',
'All new chains contain16,000 update cycles; the first4,000 are excluded from the primary summaries. Four separately seeded starts are used: uniform, unit winding, canted checkerboard, and vortex-rich parity-locked Z4. The last is not fully disordered: it initially has saturated staggered-nematic and fourth-harmonic power.',
'At q=.52/.51 each cycle has one native local Metropolis sweep. At q=.505 the isolated producer uses an exact restricted-angle conditional heatbath sweep, tested independently against quadrature and exact replay, because the earlier local Metropolis cohort had persistent start-dependent transients. Both cohorts are preserved and the new method cohort is explicitly identified by hashed manifests. The conditional sampler proposes over the complete allowed angular union with exact rejection weights. Completed conditional draws are not Metropolis acceptance or evidence of independence. All cycles also have four signed reflection clusters, four exact pi-rotation clusters, one winding attempt and one symmetric collective Metropolis-Hastings shift. All native model kernels and the environment remain frozen. All fields, observable histories, RNG states, checkpoints, provenance and full logs remain available.',
'Finite-state balance, hard-support, baseline-identity and exact checkpoint-replay tests pass. Every production file passes independent field/observable and bookkeeping audits. These checks establish correctness, not mixing. The fresh baseline comparisons do not establish a performance improvement caused by the new moves.', '',
'Analysis and limits',
'Binder observables and xi/L are formed from pooled moments, never by averaging nonlinear per-configuration ratios. Uniform/staggered spin power, Fourier peaks, vorticity, mean link cosine, staggered-nematic power, split halves and hysteretic spin-response switching counts are retained. Sampling diagnostics include estimated autocorrelation times, block counts, conditional observable-MC intervals and all flags, linked to immutable content-hashed analysis snapshots.',
'Two-knot local windows make linear and PCHIP interpolation identical; this agreement is not an independent interpolation robustness test. Native measured bracket widths set the stated exploratory resolution. Many chains retain low or insufficient block-count flags, recorded per source; agreement of start/half central values does not demonstrate full equilibration. The new plot does not calculate root-MC confidence intervals. Size/start/half/threshold/interpolation envelopes are observed sensitivity ranges. Common xi/L crossings use central pooled ratios and provide descriptive finite-pair cross-checks; no universal critical xi/L target or O(2) finite-size exponent is assumed.',
'Strong staggered-nematic response persists after the first spin loss. The suppressed side must not be called fully disordered. Within-band bimodality, finite-volume sensitivity and limited effective samples remain scientific qualifications. No missing or multiple family in the selected local window is omitted. Global families, including outside-window anomalies, remain explicit; the continuation is connected to the sampled high-spin side, with unsampled support toward kappa=0.', '',
'Artifacts and reproduction',
'exploratory-negative-kappa-transition-curve.png/pdf/svg: main continuation plot with zoom near kappa=0.',
'exploratory-negative-kappa-transition-points.csv and exploratory-negative-kappa-delivery.json: points, ranges and source hashes.',
'Each source bracket analysis contains every start/half/threshold/interpolation family, native supporting brackets, source manifest hashes, immutable diagnostic links and all sampling flags. Older positive-kappa results and earlier snapshots remain unchanged.']
report='\n'.join(lines)+'\n'
for before,after in [('and0.505','and 0.505'),('targets0.20','targets 0.20'),('a1/abs','a 1/abs'),('contain16,000','contain 16,000'),('first4,000','first 4,000')]:report=report.replace(before,after)
(P/'exploratory-negative-kappa-transition-report.txt').write_text(report)
print(json.dumps(points,indent=2))
