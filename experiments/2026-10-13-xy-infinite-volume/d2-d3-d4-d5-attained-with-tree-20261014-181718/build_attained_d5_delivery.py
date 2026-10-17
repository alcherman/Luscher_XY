"""Source-pinned attained D5 and combined figures, with rejected fits explicit."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
os.environ['MPLBACKEND']='Agg'
import argparse,csv,json,shutil
from pathlib import Path
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import run
from tree_level_prediction import endpoint,beta_critical
OUT=Path(__file__).resolve().parent


def eligible_stress(f,expected):
    r=f['result']
    if f['included_root_sizes']!=expected:return False,'missing_or_different_sizes'
    if r.get('dof',0)<=0:return False,'zero_dof_or_unidentified'
    if r.get('status') not in ['conditional_D5_size_fit','rejected_size_fit']:return False,r['status']
    if not np.isfinite([r.get('location_infinite',np.nan),*r.get('MC95',[np.nan,np.nan])]).all():return False,'nonfinite'
    return True,'rejected_fit_stress_only' if r['status']=='rejected_size_fit' else 'conditional_fit_stress_only'


def assemble_d5(snapshots):
    primary=snapshots['primary']['analysis'];rows=[];ledger=[]
    axes=[('beta',q) for q in [.90,.92,.94,.96,.98,1.]]+[('epsilon',0)]
    for axis,q in axes:
        matches=[f for f in primary['size_fits'] if (f['axis'],f['q'],f['method'],f['minimum_N'],f['correction'])==(axis,q,'linear',6,False)]
        if len(matches)!=1:raise ValueError(f'Missing unique primary {axis} {q}')
        f=matches[0];r=f['result']
        if r.get('status') not in ['conditional_D5_size_fit','rejected_size_fit']:raise ValueError(f'Unresolved primary {axis} {q}: {r}')
        center=r['location_infinite'];intervals=[];counts=Counter();central_shifts=[]
        for tag,snapshot in snapshots.items():
            for variant in snapshot['analysis']['size_fits']:
                if (variant['axis'],variant['q'])!=(axis,q):continue
                # Uniform family, excluding N4 and its duplicate aliases. No p-based model ranking.
                if (variant['minimum_N'],variant['correction']) not in [(6,False),(8,False),(6,True)]:continue
                expected=[8,10,12] if variant['minimum_N']==8 else [6,8,10,12]
                use,reason=eligible_stress(variant,expected);counts[reason]+=1
                entry=dict(variant=tag,axis=axis,q=q,method=variant['method'],minimum_N=variant['minimum_N'],correction=variant['correction'],included=use,reason=reason,fit=variant)
                ledger.append(entry)
                if use:
                    intervals.append(variant['result']['MC95']);central_shifts.append(abs(variant['result']['location_infinite']-center))
        for tag,snapshot in snapshots.items():
            present={(v['method'],v['minimum_N'],v['correction']) for v in snapshot['analysis']['size_fits'] if (v['axis'],v['q'])==(axis,q)}
            for method in ['linear','pchip']:
                for minimum,correction in [(6,False),(8,False),(6,True)]:
                    if (method,minimum,correction) not in present:
                        ledger.append(dict(variant=tag,axis=axis,q=q,method=method,minimum_N=minimum,correction=correction,included=False,reason='missing_fit_family',fit=None));counts['missing_fit_family']+=1
        if not intervals:raise ValueError('No identifiable descriptive intervals')
        low=min(center,r['MC95'][0],min(v[0] for v in intervals));high=max(center,r['MC95'][1],max(v[1] for v in intervals))
        accepted=r['status']=='conditional_D5_size_fit'
        rows.append(dict(d=5,epsilon_over_pi=q if axis=='beta' else center,beta=center if axis=='beta' else 0.,MC_low=r['MC95'][0],MC_high=r['MC95'][1],sensitivity_low=low,sensitivity_high=high,error_axis='vertical' if axis=='beta' else 'horizontal',status='attained_conditional_D5_extrapolation' if accepted else 'rejected_size_fit_diagnostic_not_accepted',accepted_extrapolation=accepted,p=r['p'],sizes=f['included_root_sizes'],input_flags=r['input_flags'],sensitivity_counts=dict(counts),max_tested_central_shift=max(central_shifts),baseline_fit=f))
    return rows,ledger


def write_csv(path,rows,fields):
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,extrasaction='ignore');writer.writeheader();writer.writerows(rows)


def draw_figures(folder,rows):
    colors={2:'#2563a7',3:'#218463',4:'#b25629',5:'#7952a3'}
    def panel(ax,D):
        color=colors[D];vals=[r for r in rows if r['d']==D]
        for r in vals:
            horizontal=r['error_axis']=='horizontal';center=r['epsilon_over_pi'] if horizontal else r['beta'];diagnostic=r.get('accepted_extrapolation') is False
            for lo,hi,lw,alpha in [(r['sensitivity_low'],r['sensitivity_high'],4,.24),(r['MC_low'],r['MC_high'],1,1)]:
                if not lo<=center<=hi:raise ValueError(f'Interval does not contain plotted center: {D}, {r}, {lo}, {hi}')
                errors=np.array([[center-lo],[hi-center]])
                ax.errorbar([r['epsilon_over_pi']],[r['beta']],**({'xerr':errors} if horizontal else {'yerr':errors}),fmt='none',ecolor=color,elinewidth=lw,alpha=alpha,capsize=2)
            ax.plot(r['epsilon_over_pi'],r['beta'],marker='s' if diagnostic else ('D' if horizontal else 'o'),ms=5,markerfacecolor='white' if diagnostic else color,markeredgecolor=color,linestyle='none')
        qs=np.linspace(endpoint(D),1,300);ax.plot(qs,[beta_critical(D,q) for q in qs],'--',color='#222222',lw=1.5,label='Tree-level prediction')
        ax.plot([],[],marker='o',color=color,ls='none',ms=4,label='Attained MC extrapolations')
        if D==5:ax.plot([],[],marker='s',mfc='white',mec=color,ls='none',label=r'$q=1$: rejected-fit diagnostic')
        ax.set(title=f'd = {D}'+('  (exploratory BKT)' if D==2 else ''),xlabel=r'$\epsilon/\pi$',ylabel=r'$\beta$')
        lows=[r['sensitivity_low'] for r in vals if r['error_axis']=='horizontal'];ax.set_xlim(min([endpoint(D),*lows])-.015,1.012)
        ymax=max([beta_critical(D,1),*[r['sensitivity_high'] for r in vals if r['error_axis']=='vertical']]);ax.set_ylim(-ymax*.035,ymax*1.08);ax.grid(alpha=.16);ax.legend(loc='upper left',fontsize=8,frameon=False)
    footer='Dark: conditional MC 95%. Pale: descriptive stress envelope, not total confidence coverage.\nD5 open square: rejected size-fit diagnostic. Endpoints: horizontal bars only; MC points unconnected.'
    fig,axes=plt.subplots(2,2,figsize=(11,8.5))
    for D,ax in zip([2,3,4,5],axes.flat):panel(ax,D)
    fig.suptitle('Attained XY phase boundaries and tree-level prediction',fontsize=14);fig.text(.5,.942,r'$(2D-1)v_1(\beta,\epsilon)=1$',ha='center',fontsize=12)
    fig.text(.5,.018,footer,ha='center',fontsize=8);fig.tight_layout(rect=[0,.075,1,.93])
    for ext in ['png','pdf','svg']:fig.savefig(folder/f'combined-MC-vs-tree-level.{ext}',dpi=180)
    plt.close(fig)
    for D in [2,3,4,5]:
        fig,ax=plt.subplots(figsize=(7,5.3));panel(ax,D);fig.text(.5,.018,footer,ha='center',fontsize=7);fig.tight_layout(rect=[0,.10,1,1])
        for ext in ['png','pdf','svg']:fig.savefig(folder/f'd{D}-MC-vs-tree-level.{ext}',dpi=180)
        plt.close(fig)



def draw_size_diagnostics(folder,analysis):
    fig,axes=plt.subplots(2,4,figsize=(14,6.8))
    keys=[('beta',q) for q in [.90,.92,.94,.96,.98,1.]]+[('epsilon',0)]
    for ax,(axis,q) in zip(axes.flat,keys):
        roots=sorted([r for r in analysis['roots'] if (r['axis'],r['q'],r['method'])==(axis,q,'linear') and r['N']>=6],key=lambda r:r['N'])
        f=next(f for f in analysis['size_fits'] if (f['axis'],f['q'],f['method'],f['minimum_N'],f['correction'])==(axis,q,'linear',6,False))['result']
        for r in roots:
            x=r['N']**(-2.5);y=r['location'];lo,hi=r['MC95']
            if not lo<=y<=hi:raise ValueError('Finite-root interval misses its center')
            ax.errorbar(x,y,yerr=[[y-lo],[hi-y]],fmt='o',ms=3,color='#7952a3',capsize=2);ax.annotate(f"N={r['N']}",(x,y),xytext=(2,5),textcoords='offset points',fontsize=6)
        xs=np.linspace(0,6**(-2.5)*1.04,120);color='#b52d35' if f['status']=='rejected_size_fit' else '#333333'
        ax.plot(xs,f['coefficients'][0]+f['coefficients'][1]*xs,'--',color=color,lw=1)
        label='beta=0 endpoint' if axis=='epsilon' else f"q={q:.2f}"
        ax.set(title=f"{label}; p={f['p']:.3g}",xlabel=r'$N^{-5/2}$',ylabel=r'$\epsilon_N/\pi$' if axis=='epsilon' else r'$\beta_N$');ax.grid(alpha=.15)
        if f['status']=='rejected_size_fit':ax.text(.03,.94,'REJECTED FIT',transform=ax.transAxes,va='top',color=color,fontsize=7)
    axes.flat[-1].axis('off');axes.flat[-1].text(0,.9,'Linear Binder U=1/4 roots\nN >= 6, leading size ansatz\n\nBars: conditional MC 95%\nDashed: size fit, not MC interpolation\nRed: rejected fit\n\nLow-block flags remain;\nno certified total accuracy.',va='top',fontsize=9)
    fig.suptitle('5D finite-size locations and leading extrapolation diagnostics');fig.tight_layout(rect=[0,0,1,.95])
    for ext in ['png','pdf','svg']:fig.savefig(folder/f'd5-finite-size-diagnostics.{ext}',dpi=180)
    plt.close(fig)


def main():
    manager=OUT/json.loads((OUT/'latest-d5-sensitivity.json').read_text())['folder'];manifest=json.loads((manager/'manifest.json').read_text())
    if manifest.get('status')!='snapshots_complete_pending_scientific_review':raise RuntimeError('Wait for all source-pinned snapshots')
    snapshots={}
    for tag,info in manifest['snapshots'].items():
        source=OUT/info['folder']/'analysis.json'
        if info['returncode'] or not info.get('source_ids_match_expected') or run.sha256(source)!=info['analysis_sha256']:raise RuntimeError('Snapshot failed or changed')
        snapshots[tag]=dict(source=str(source),sha256=run.sha256(source),analysis=json.loads(source.read_text()))
    d5,ledger=assemble_d5(snapshots)
    attained=OUT/json.loads((OUT/'latest-attained-d2-d3-d4-delivery.json').read_text())['folder'];old=json.loads((attained/'analysis.json').read_text());rows=old['rows']+d5
    stamp=run.now()[:19].replace('-','').replace(':','').replace('T','-');folder=OUT/('d2-d3-d4-d5-attained-with-tree-'+stamp);folder.mkdir()
    for name in [Path(__file__).name,'tree_level_prediction.py','d5_scaling.py']:shutil.copyfile(OUT/name,folder/name)
    shutil.copyfile(attained/'attained-phase-boundary-report.md',folder/'prior-d2-d3-d4-qualifications.md')
    shutil.copyfile(manager/'comparisons.json',folder/'d5-complete-variant-coverage.json')
    shutil.copyfile(manager/'manifest.json',folder/'d5-source-cohort-manifest.json')
    fields=['d','epsilon_over_pi','beta','MC_low','MC_high','sensitivity_low','sensitivity_high','error_axis','status','p','accepted_extrapolation']
    write_csv(folder/'phase-boundary-locations.csv',rows,fields);write_csv(folder/'d5-transition-locations.csv',d5,fields)
    finite=[dict(**{k:r.get(k) for k in ['axis','q','N','method','status','location','MC95','valid_fraction','selected_window','flags']}) for r in snapshots['primary']['analysis']['roots']]
    write_csv(folder/'d5-finite-size-locations.csv',finite,['axis','q','N','method','status','location','MC95','valid_fraction','selected_window','flags'])
    shutil.copyfile(attained/'unresolved-locations.csv',folder/'prior-unresolved-locations.csv')
    unresolved=[dict(d=5,epsilon_over_pi=r['epsilon_over_pi'],reason=r['status'],diagnostic_beta=r['beta'],accepted_infinite_beta=None) for r in d5 if not r['accepted_extrapolation']]
    unresolved.append(dict(d=5,epsilon_over_pi=.88,reason='Not measured as a positive-beta slice; no claim of absence from missing finite-box crossing',diagnostic_beta=None,accepted_infinite_beta=None))
    write_csv(folder/'d5-unresolved-locations.csv',unresolved,['d','epsilon_over_pi','reason','diagnostic_beta','accepted_infinite_beta'])
    draw_figures(folder,rows)
    draw_size_diagnostics(folder,snapshots['primary']['analysis'])
    prediction=dict(relation='(2D-1)*v1(beta,epsilon)=1',v1='normalized integral exp(beta*cos(theta))*cos(theta) on [-epsilon,epsilon]',zero_beta={str(D):endpoint(D) for D in [2,3,4,5]},at_epsilon_pi={str(D):beta_critical(D,1) for D in [2,3,4,5]},document='/Users/acherman/Downloads/constrained_xy_large_d_note.pdf',document_sha256=run.sha256(Path('/Users/acherman/Downloads/constrained_xy_large_d_note.pdf')))
    run.atomic_json(folder/'analysis.json',dict(created=run.now(),rows=rows,prior_source=str(attained/'analysis.json'),prior_sha256=run.sha256(attained/'analysis.json'),snapshots={k:{x:y for x,y in v.items() if x!='analysis'} for k,v in snapshots.items()},manager=str(manager),prediction=prediction,baseline_policy='Uniform minN6 leading N^-2.5, locally bracketed linear U=.25. q1 rejected baseline shown only as open diagnostic, not an accepted coupling.',stress_policy='Union intervals across all9 snapshots, linear/PCHIP, min6/min8 leading and min6 additiveN^-3. Same required sizes, identifiable positiveDOF only; rejected/flagged fits explicitly included as stress diagnostics, never accepted estimates. Not calibrated confidence coverage.',prior_rows_unchanged=rows[:len(old['rows'])]==old['rows']))
    run.atomic_json(folder/'d5-stress-ledger.json',ledger)
    p=snapshots['primary']['analysis'];flags=Counter(flag for g in p['groups'] for flag in g['summary']['flags']);minimum_blocks=min(d['blocks'] for g in p['groups'] for d in g['summary']['diagnostics'])
    table='\n'.join(f"| {'beta=0' if r['error_axis']=='horizontal' else format(r["epsilon_over_pi"],".2f")} | {r['epsilon_over_pi'] if r['error_axis']=='horizontal' else r['beta']:.6f} | [{r['MC_low']:.6f}, {r['MC_high']:.6f}] | [{r['sensitivity_low']:.6f}, {r['sensitivity_high']:.6f}] | {'conditional' if r['accepted_extrapolation'] else 'rejected-fit diagnostic'} |" for r in d5)
    report=f'''# Attained 5D phase boundary and combined comparison

All282 frozen5D trajectories passed full-field closure audits. There are six positive-beta slices q=epsilon/pi=.90,.92,.94,.96,.98,1.00 and a beta=0 scan, with N=6,8,10,12 and independent starts. N4 pilot data are preserved but excluded from the uniform baseline. The primary contains31linear finite-size roots, all100% bootstrap-valid, with62linear/PCHIP roots total. Audits establish recorded-data correctness, not mixing or infinite-volume accuracy.

Five positive-beta slices and the endpoint have conditional baseline size fits; q=1 fails that diagnostic. Its nominal extrapolation is shown as an OPEN SQUARE only and is not an accepted critical coupling. No larger-box campaign is required for this practical attained delivery. Rejected q1 leading/correction fits and the exploratory minN8 alternative remain visible in the complete stress ledger.

## Locations and uncertainty

Positive rows give beta, the beta=0 row gives epsilon/pi. These are descriptive attained estimates under the method below.

| Scan | Baseline location | Conditional MC95 | Descriptive stress envelope | Status |
|---|---:|---|---|---|
{table}

The uniform reference is locally bracketed linear U=.25 followed by beta_N=beta_infinity+a*N^-2.5 (or epsilon_N at beta=0), N>=6. This periodic above-upper-critical-dimension hypothesis is not the4D logarithmic law. O(2) quartic zero-mode Binder U*=1-pi/4=.2146: U=.25 is a shifted ordered-side fixed criterion, not the universal finite-box critical value. An additiveN^-3 correction is a stress hypothesis. Source references are Chen–Dohm cond-mat/9809394 and Wittmann–Young1410.5296, as documented in the campaign design.

Dark bars are conditional2000-draw block-bootstrap intervals, conditional on central-data local-window selection, autocorrelation/block estimates and the size ansatz. For q1 that ansatz is rejected and its interval is only diagnostic. Pale bars unite intervals from9source-pinned snapshots: baseline/fullwindow/cut2001/cut3001/doubledblocks/U.20/U.30/even/oddstarts; each tests linear/PCHIP, N>=6 or8 leading, and N>=6 correction. Only identifiable, positive-DOF fits with all required sizes enter. Rejected fits and flagged inputs remain INCLUDED in this stress envelope and explicitly flagged in d5-stress-ledger.json. Missing sizes, insufficient bootstrap retention, rank failures and zeroDOF fits are separately recorded. The envelope includes amplified Monte Carlo variance, threshold dependence and rejected-model diagnostics: it is neither pure systematic uncertainty nor a certified total95% interval. A small envelope does not demonstrate totalaccuracy. No model is selected by bestp. Even/odd draws share helper seeds and are not paired for difference significance.

## Sampling and interpolation limits

There are{len(p['groups'])}pooled coordinate groups, with group flag counts{dict(flags)}. The smallest retained block count per start is{minimum_blocks}. Selected-root flags and all-measured flags are preserved; absence of start/half flags does not certify equilibration. Long autocorrelations leave few blocks for some N12starts, limiting the reliability of small conditional intervals. The complete first/last-half, per-start, pooled, autocorrelation and block diagnostics remain in the pinned analysis snapshots. Shared beta=0 q=.90 sources used on both scan axes are correlated evidence, not additional independent replication.

The final matched midpoints changed finite-size roots by as much as .003584 relative to the earlier two-knot windows; this comparison also changes measured sources and is not a pure interpolation-bias estimate. Most local positive-beta linear/PCHIP comparisons still use two knots and are identical by construction. The full measured-window PCHIP comparison adds a shape sensitivity check, but does not certify interpolation coverage or exclude curvature between samples. Local3SD windows are data-selected. All unresolved/tail losses and skipped groups in every variant remain explicit; successful common-root comparisons are not all-root validation.

## Comparison and prior dimensions

Dashed theory solves the full supplied (2D-1)v1=1 criterion, v1 the normalized exponential-cosine integral over[-epsilon,epsilon]. It is the regular-tree instability prediction, not an exact hypercubic critical curve. The supplied PDF definition is onpage2 and orangebox onpage5; itsSHA is recorded in analysis.json. MC points are unconnected. Fixed-epsilon points have vertical bars only; beta=0 endpoints horizontal bars only. Open D5square denotes rejected-fit diagnostic, not a validated extrapolation.

The42 attained D2/D3/D4 rows and their errors are copied exactly, with no new2Danalysis. See prior-d2-d3-d4-qualifications.md. In particular D2q1's dark interval excludes1.1199; its wider descriptive envelope includes it, so certified absolute BKT accuracy is not claimed. D2 uses N<=512 in its primary, with512only on representative slices and endpoint. Completed1024cost probes remain preserved/excluded. Prior unresolved neighbors remain in prior-unresolved-locations.csv. The fixed D5q=.88 positive slice was not measured; no absence inference is made from missing finite-box crossings.

## Provenance and preservation

analysis.json pins all9snapshot paths and hashes, the prior delivery hash and theory PDF hash. d5-finite-size-locations.csv contains root windows/status/flags; d5-stress-ledger.json identifies every included or excluded stress fit. All histories, logs, provenance and exact checkpoints remain intact. User-authorized field-only retention runs after full audits. No numerical simulation kernels or running environment were changed. This folder remains pending independent scientific and visual review until delivery-status.json says reviewed_attained_with_qualifications.
'''
    (folder/'attained-phase-boundary-report.md').write_text(report)
    run.atomic_json(folder/'delivery-status.json',dict(status='pending_independent_scientific_and_visual_review',source_count=282,prior_estimates=42,d5_conditional_extrapolations=6,d5_rejected_diagnostic=1,total_plotted=49,report='attained-phase-boundary-report.md'))
    run.atomic_json(OUT/'latest-draft-d5-delivery.json',dict(folder=folder.name,status='pending_independent_scientific_and_visual_review'))
    print(folder,flush=True)

if __name__=='__main__':main()
