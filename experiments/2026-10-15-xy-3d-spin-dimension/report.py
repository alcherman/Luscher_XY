"""Summarize conditional scaling-dimension measurements and their limitations."""
from pathlib import Path
import sys,json,shutil,csv
import numpy as np
OUT=Path(__file__).resolve().parent;sys.path.insert(0,str(OUT))
import measure as m
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator
TAGS=['primary','full','U20','U30','cut501','cut2001','block2','even','odd','evenfull','oddfull','U30full']

def actual_binder(a,r):
    if r.get('location') is None:return None
    gs=sorted([g for g in a['groups'] if (g['q'],g['N'])==(r['q'],r['N']) and 'mean' in g['summary']],key=lambda g:g['beta']);ix=r['selected_indices'];x=np.array([gs[i]['beta'] for i in ix]);means=np.array([gs[i]['summary']['mean'] for i in ix]);root=r['location'][0]
    values=PchipInterpolator(x,means,axis=0)(root) if r['method']=='pchip' else np.array([np.interp(root,x,means[:,j]) for j in range(3)])
    return float(m.binder(values))

def main():
    snapshots={}
    for tag in TAGS:
        pointer=OUT/f'latest-{tag}.json'
        if not pointer.exists():raise RuntimeError(f'Analysis still pending: {tag}')
        f=OUT/json.loads(pointer.read_text())['folder']/'analysis.json';snapshots[tag]=dict(source=str(f),sha256=m.run.sha256(f),analysis=json.loads(f.read_text()))
    primary=snapshots['primary']['analysis'];folder=OUT/('delivery-'+m.run.now()[:19].replace('-','').replace(':','').replace('T','-'));folder.mkdir()
    for n in ['report.py','measure.py','sensitivity.py']:shutil.copyfile(OUT/n,folder/n)
    base={};ledger=[];coverage=[];rows=[]
    for r in primary['contours']:
        if r['method']=='linear' and r['status']=='conditional_contour':base.setdefault(r['q'],[]).append(r['N'])
    for tag,s in snapshots.items():
        a=s['analysis']
        for r in a['contours']:
            achieved=actual_binder(a,r)
            coverage.append(dict(variant=tag,q=r['q'],N=r['N'],method=r['method'],status=r['status'],valid_fraction=r['valid_fraction'],achieved_Binder=achieved,target=a['settings']['target'],Binder_discrepancy=achieved-a['settings']['target'] if achieved is not None else None,flags=r.get('flags',[])))
        for f in a['fits']:
            r=f['result'];expected=sorted(n for n in base[f['q']] if n>=f['minimum']);same=r.get('sizes')==expected and not f['omitted_sizes'];flags=sorted({flag for c in a['contours'] if c['q']==f['q'] and c['method']==f['method'] and c['N'] in r.get('sizes',[]) for flag in c.get('flags',[])})
            ledger.append(dict(variant=tag,fit=f,expected_sizes=expected,same_complete_sizes=same,input_flags=flags,identifiable_positive_dof=r.get('dof',0)>0 and r['status'] in ['conditional_dimension_fit','rejected_power_law'],scope='Missing sizes or zeroDOF/rejected fits never establish robustness; no p-value model selection.'))
    for f in primary['fits']:
        if f['method']!='linear' or f['minimum']!=8 or f['omega'] is not None:continue
        r=f['result'];rows.append(dict(epsilon_over_pi=f['q'],moment=f['moment'],Delta=r.get('Delta'),MC_low=r.get('MC95',[None,None])[0],MC_high=r.get('MC95',[None,None])[1],status=r['status'],sizes=r.get('sizes'),omitted=f['omitted_sizes'],p=r.get('p'),valid_fraction=r.get('valid_fraction')))
    fields=list(rows[0]);
    with (folder/'spin-dimension-estimates.csv').open('w',newline='') as stream:w=csv.DictWriter(stream,fieldnames=fields);w.writeheader();w.writerows(rows)
    m.run.atomic_json(folder/'contour-and-moment-consistency.json',coverage);m.run.atomic_json(folder/'full-sensitivity-ledger.json',ledger)
    source_provenance=[]
    for source in primary['source_jobs']:
        r=json.loads((m.CAMPAIGN/'retention-records'/f"{source['id']}.json").read_text());source_provenance.append(dict(id=source['id'],full_audit_sha256=source['audit']['sha256'],current_retained_sha256=r['compact_sha256'],current_retained_bytes=r['compact_bytes']))
    m.run.atomic_json(folder/'source-provenance.json',source_provenance)
    q1=next(r for r in rows if r['epsilon_over_pi']==1 and r['moment']=='m2')
    f1=next(f for f in primary['fits'] if (f['q'],f['method'],f['moment'],f['minimum'],f['omega'])==(1.,'linear','m2',8,None));roots=sorted([r for r in primary['contours'] if r['q']==1 and r['method']=='linear' and r['status']=='conditional_contour'],key=lambda r:r['N'])
    fig,(ax,bx)=plt.subplots(1,2,figsize=(12,5.4))
    for r in roots:
        y=r['location'][1];lo,hi=r['MC95'][1];ax.errorbar(r['N'],y,yerr=[[y-lo],[hi-y]],fmt='o',color='#2563a7' if r['N']>=8 else '#999999',capsize=3)
    L=np.geomspace(6,26,100);coef=f1['result']['coefficients'];ax.plot(L,np.exp(coef[0])*L**(-2*coef[1]),'--',color='#222222',label=rf'$L\geq8$ fit: $\Delta={q1["Delta"]:.3f}$');ax.set(xscale='log',yscale='log',xlabel=r'$L$',ylabel=r'$\langle |m|^2\rangle$ at $U=1/4$',title=r'$\epsilon=\pi$: direct finite-size measurement');ax.set_xticks([6,8,12,16,24],labels=['6','8','12','16','24']);ax.grid(alpha=.15);ax.legend(frameon=False,fontsize=9)
    for moment,color,marker,label in [('m2','#2563a7','o',r'from $\langle |m|^2\rangle$'),('mabs','#ba6a21','^',r'from $\langle |m|\rangle$')]:
        rr=[r for r in rows if r['moment']==moment and r['Delta'] is not None]
        for i,r in enumerate(rr):
            if r['MC_low']<=r['Delta']<=r['MC_high']:
                bx.errorbar(r['epsilon_over_pi'],r['Delta'],yerr=[[r['Delta']-r['MC_low']],[r['MC_high']-r['Delta']]],fmt=marker,color=color,mfc=color if r['status']=='conditional_dimension_fit' and moment=='m2' else 'white',ms=5,capsize=2,label=label if i==0 else None)
        
    bx.axhline(.51905,ls='--',lw=1,color='#333333',label=r'3D XY reference $\Delta\approx0.519$');bx.set(xlabel=r'$\epsilon/\pi$',ylabel=r'effective $\Delta_{e^{i\phi}}$',title=r'Fixed $U=1/4$, leading fits with $L\geq8$',ylim=(.30,.75));bx.grid(alpha=.15);bx.legend(frameon=False,fontsize=8,loc='lower left')
    fig.suptitle('3D spin scaling dimension from existing Monte Carlo histories',fontsize=13);fig.text(.5,.025,'Bars: conditional MC 95%, not total uncertainty. Finite-size/interpolation sensitivity remains; q=.86 has only two resolved sizes (diagnostic).',ha='center',fontsize=8);fig.tight_layout(rect=[0,.08,1,.94])
    for ext in ['png','svg']:fig.savefig(folder/f'3d-spin-scaling-dimension.{ext}',dpi=180)
    plt.close(fig)
    q1_simple=[e for e in ledger if e['fit']['q']==1 and e['fit']['minimum']==8 and e['fit']['omega'] is None and e['same_complete_sizes'] and e['identifiable_positive_dof']]
    centers=[e['fit']['result']['Delta'] for e in q1_simple];intervals=[e['fit']['result']['MC95'] for e in q1_simple]
    q1_correction=[e for e in ledger if e['variant']=='primary' and e['fit']['q']==1 and e['fit']['moment']=='m2' and e['fit']['method']=='linear' and e['fit']['omega'] is not None]
    m.run.atomic_json(folder/'analysis.json',dict(created=m.run.now(),status='pending_independent_review',primary_rows=rows,snapshots={k:{kk:vv for kk,vv in v.items() if kk!='analysis'} for k,v in snapshots.items()},q1_same_size_leading_sensitivity=dict(count=len(q1_simple),central_range=[min(centers),max(centers)],interval_union=[min(v[0] for v in intervals),max(v[1] for v in intervals)]),scope='Effective integrated-order-parameter scaling dimension. Conditional MC intervals, incomplete threshold/start coverage and correction/contact backgrounds prevent a precision/CFT/universality claim.',reference=dict(Delta=.51905,eta=.0381,url='https://arxiv.org/abs/cond-mat/0605083',used_in_fits=False)))
    report=f'''3D SCALING DIMENSION OF exp(i phi): EXISTING-MC EXPLORATORY MEASUREMENT

Main result: epsilon=pi, N=8,12,16,24, fixed Binder U=.25, linear contour interpolation and leading log-size fit give Delta={q1['Delta']:.6f}, conditional block-bootstrap 95% [{q1['MC_low']:.6f},{q1['MC_high']:.6f}]. The independent-in-moment-form check using recorded |m| gives Delta=.515651, MC95 [.461721,.559438]. Both use the same MC histories, so they are correlated checks, not independent replicates. No known nu or Delta is inserted into the fits. The established 3D XY reference Delta=(1+eta)/2 approximately .519 follows eta=.0381 in Campostrini et al., https://arxiv.org/abs/cond-mat/0605083; it is a comparison only.

METHOD AND WHAT IS MEASURED
m is the spatial average of exp(i phi). If a fixed dimensionless Binder contour approaches a continuous critical scaling regime, <|m|^2> scales as L^(-2Delta) and <|m|> as L^(-Delta) at fixed scaling variable. We fit these logarithms with a free Delta and amplitude. The integrated-moment relation assumes IR dominance 2Delta<3; contact/analytic terms proportional to L^-3 and irrelevant corrections can bias finite-box slopes. This is an effective integrated spin-exponent measurement, not direct proof of a local conformal operator dimension or universality across the whole constrained curve.

{len(primary['source_jobs'])} existing audited 3D trajectories contribute {len(primary['groups'])} coordinate groups over ten fixed epsilon/pi slices. No new MC simulations, field data changes or SSD removal. N=24 exists on q=.84,.92,1; other positive slices reach N=16. The current primary has joint m2,m4,recorded|m| block resampling, monitoring all moments and the Binder influence function for autocorrelation; ten-times-largest-tau blocks, minimum20rows. Each draw's own Binder root is used to interpolate its own magnetization before taking logarithms. Source IDs seed distinct coordinate resamples. The log-fit GLS covariance and validity mask are recomputed after every size cut. Fixed central local windows include a bracketing pair and one neighbor on either side; isotonic Binder smoothing is an explicit assumption. Intervals are conditional on these choices, not total coverage. Pooled/individual-start and half-history flags appear in the derived ledger.

SENSITIVITY AND LIMITATIONS
Linear interpolation of Binder and magnetization separately need not give exactly the target Binder when the interpolated m2,m4 ratio is formed. contour-and-moment-consistency.json quantifies that discrepancy; moment_linear solves the Binder ratio from interpolated raw moments as an explicit alternative, not a claim of zero interpolation bias. PCHIP and full measured windows are further shape checks. The primary linear interpolation at N24 actually gives Binder .22463 at q1 and .19458 at q=.92 when reconstructed from its interpolated moments, rather than .25. These are material finite-grid inconsistencies, not extra Monte Carlo noise.
For epsilon=pi, same-complete-size leading fits across available methods, moments and variants have central values spanning {min(centers):.6f} to {max(centers):.6f}; their MC interval union is [{min(v[0] for v in intervals):.6f},{max(v[1] for v in intervals):.6f}]. This union includes rejected positive-DOF fits as diagnostics, is descriptive and correlated, and excludes incomparable/missing-size comparisons; it is NOT total confidence coverage. In particular full-window U=.30 with moment-linear interpolation gives Delta=.32560, MC95 [.28789,.35335], but is rejected by the size-fit diagnostic (p=.000165); full-window linear/PCHIP at that threshold give .45584/.48591. This prevents an all-threshold robustness claim even at q1.
At q1, threshold .20 misses N12 coverage; local .30 and start-parity studies may omit sizes. Full-window follow-ups restore all four sizes for odd starts and U=.30; even starts still miss N12. They are separate source-preserving diagnostics; full-sensitivity-ledger.json records exact included/omitted sizes and every fit status. Successful comparisons cannot establish all-variant robustness. q=.86 primary N8 is unresolved; its N12/16-only fit has zero DOF and is plotted only as a diagnostic.
Adding N6 shifts q1's leading m2 exponent to .4821; using N>=12 gives .4716 with much larger errors. Adding an assumed log-amplitude correction c*N^(-omega), with omega=.5,.8,1, gives central q1 values .6013,.5655,.5538 and broad conditional intervals (roughly .20-.88 across those three fits). These models are weakly constrained and are stress tests, not selected by favorable p-values. The analytic/contact background was not independently fitted or excluded. Thus the present data support a coarse exponent of order one-half, not a high-precision or certified absolute determination.
Constrained-slice interpolation sensitivity is materially larger: at q=.92 the primary N>=8 m2 results are .46381 (linear Binder), .47689 (PCHIP), and .36417 (Binder from linearly interpolated moments). The ~.100 spread exceeds typical conditional MC errors despite acceptable size-fit diagnostics. Those constrained-slice values are not robust exponent determinations. No pooling across epsilon was used to force a universal result.

FILES AND REPRODUCIBILITY
spin-dimension-estimates.csv lists the leading reference values and flags zero-DOF diagnostics. full-sensitivity-ledger.json retains all12 snapshots' methods, cuts, block factors, threshold targets, separate starts, size cuts, correction powers, statuses, required-size comparisons and selected-contour flags. source-provenance.json records original audit/current retained-file hashes. analysis.json pins all source snapshot hashes. Every original snapshot, cache and full execution log is preserved. The field configurations are not needed for this analysis; complete retained scalar MC histories suffice. Five analytic/covariance/interpolation/identifiability oracle tests pass in tests.log. Report and figure remain pending independent scientific review until review.json says accepted_with_qualifications.
'''
    (folder/'spin-dimension-measurement-report.txt').write_text(report);m.run.atomic_json(OUT/'latest-delivery.json',dict(folder=folder.name,status='pending_independent_review'));print(folder,flush=True)

if __name__=='__main__':main()
