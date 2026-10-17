"""Report beta-zero spin dimension without selecting an interpolation by fit p."""
from pathlib import Path
import json,csv,shutil
import numpy as np
import measure as m
import matplotlib.pyplot as plt

ROOT=m.OUT/'zero-beta'
TAGS=['primary','full','cut501','cut2001','block2','U20','U30','even','odd']

def main():
    snapshots={}
    for tag in TAGS:
        path=ROOT/json.loads((ROOT/f'latest-{tag}.json').read_text())['folder']/'analysis.json'
        snapshots[tag]=dict(path=str(path),sha256=m.run.sha256(path),analysis=json.loads(path.read_text()))
    a=snapshots['primary']['analysis'];folder=ROOT/('delivery-'+m.run.now()[:19].replace('-','').replace(':','').replace('T','-'));folder.mkdir()
    for name in ['zero_beta.py','zero_beta_report.py','measure.py']:shutil.copyfile(m.OUT/name,folder/name)
    ledger=[dict(variant=t,**f) for t,s in snapshots.items() for f in s['analysis']['fits']]
    m.run.atomic_json(folder/'sensitivity-ledger.json',ledger)
    rows=[]
    for r in a['contours']:
        rows.append(dict(N=r['N'],method=r['method'],status=r['status'],epsilon_over_pi=r.get('epsilon_over_pi'),m2=r.get('location',[None]*3)[1],mabs=r.get('location',[None]*3)[2],achieved_Binder=r.get('moment_ratio_at_interpolated_root'),joint_valid_fraction=r['valid_fraction'],flags=';'.join(r['flags'])))
    with (folder/'finite-size-contours.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    fig,(ax,bx)=plt.subplots(1,2,figsize=(12,5))
    styles=[('linear','#777777','o','Binder-linear'),('pchip','#2563a7','s','PCHIP'),('moment_linear','#b45e20','^','Moment-consistent linear')]
    for method,color,marker,label in styles:
        rs=sorted([r for r in a['contours'] if r['method']==method and r.get('location')],key=lambda r:r['N'])
        for i,r in enumerate(rs):
            y=r['location'][1];lo,hi=r['MC95'][1]
            ax.errorbar(r['N'],y,yerr=[[y-lo],[hi-y]],fmt=marker,color=color,capsize=2,label=label if i==0 else None)
        fs=[f for f in a['fits'] if f['method']==method and f['moment']=='m2' and f['omega'] is None]
        for f in fs:
            r=f['result']
            if 'Delta' not in r:continue
            y=r['Delta'];lo,hi=r['MC95'];bx.errorbar(f['minimum'],y,yerr=[[y-lo],[hi-y]],fmt=marker,color=color,mfc='white' if r['status']!='conditional_dimension_fit' else color,capsize=2)
    ax.set(xscale='log',yscale='log',xlabel='$L$',ylabel=r'$\langle |m|^2\rangle$ on nominal $U=1/4$ contour',title=r'$\beta=0$: finite-size magnetization');ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.15)
    bx.axhline(.51905,color='#333333',ls='--',lw=1,label='3D XY reference, not fitted');bx.set(xlabel='Minimum included L',ylabel=r'effective $\Delta_{e^{i\phi}}$',title='Leading power fits: interpolation / size-cut sensitivity');bx.set_xticks([6,8,12,16]);bx.grid(alpha=.15);bx.legend(frameon=False,fontsize=8)
    fig.suptitle('3D spin dimension at zero beta: present data do not give a robust precision estimate');fig.text(.5,.015,'Bars: conditional MC 95% only. Open symbols: rejected size fits. Methods are correlated alternatives, not independent measurements.',ha='center',fontsize=8);fig.tight_layout(rect=[0,.055,1,.94])
    for ext in ['png','svg','pdf']:fig.savefig(folder/f'3d-zero-beta-spin-dimension.{ext}',dpi=180)
    plt.close(fig)
    report='''3D SPIN SCALING DIMENSION AT EXACTLY BETA=0

Conclusion: the existing endpoint MC data do not establish a robust precision value of Delta_exp(i phi). Smooth-curve interpolation gives an effective exponent near .54, compatible with the ordinary 3D XY reference .519, but another explicitly tested interpolation gives .62 and rejected size fits. The available data do not establish a discrepancy from XY universality, or verify it precisely. No new MC was launched.

METHOD AND SOURCES
72 audited trajectories at36 exact physical coordinates, L=6,8,12,16,24,32. Each start is block bootstrapped using all retained m2,m4,recorded|m| scalar histories and the Binder influence function; 2000 joint draws. Independently resampled starts are pooled by row counts. Each draw locates its own epsilon crossing and its own magnetization there. U=1-m4/(2*m2^2) decreases with epsilon: the code uses x=-epsilon/pi. Actual epsilon radians are converted directly, never rounded onto a grid. All18 primary method/size contours retain100% draws; this proves bracket retention, not mixing or scaling accuracy. L32 has fewer than50 blocks per start. Histories/checkpoints/fields/SSD originals remain untouched.

At a fixed dimensionless contour approaching continuous critical scaling, m is the spatial average of exp(i phi), <|m|^2>~L^(-2Delta), and recorded <|m|>~L^(-Delta). GLS propagates root/magnetization covariance through shared bootstrap draws and recomputes joint-valid masks after size cuts. No nu, eta or Delta is fixed. This inference assumes infrared dominance 2Delta<3; analytic/contact L^-3 and corrections may matter. This is not a direct local correlator measurement. Universal comparison Delta=(1+eta)/2=.51905 uses eta=.0381 from Campostrini et al., https://arxiv.org/abs/cond-mat/0605083, as a reference only, not a beta-zero universality assumption.

LEADING M2 FITS, L=8,12,16,24,32, NOMINAL U=.25
Binder-linear: Delta=.541603, conditional MC95 [.516621,.565099]; chi2=11.4003 for3DOF, p=.009747, rejected by the declared p<.01 check.
PCHIP: Delta=.537298, MC95 [.507131,.567521]; chi2=1.04116 for3DOF, p=.79129. This is a conditional diagnostic, not a preferred estimator selected for a good p.
Raw-moment-linear, enforcing the Binder ratio exactly: Delta=.624408, MC95 [.611496,.637351]; chi2=26.8490 for3DOF, p=.00000633, rejected. Its tight conditional error is not a trustworthy final exponent error.

WHY A PRECISE CLAIM IS UNSUPPORTED
Independent interpolation of Binder and raw moments is not algebraically consistent: primary Binder-linear interpolation at L8 reconstructs U=.11792 rather than .25, and L24 gives .19061. Moment-linear enforces U=.25 but assumes moments linear across coarse epsilon intervals, which need not hold. Neither procedure is guaranteed accurate; smooth PCHIP and a favorable size-fit p do not settle this. The measured grid is insufficient to control the interpolation contribution for the operator-dimension inference. Finite-size effects/noise remain possible, and are not isolated as the sole cause.

Size cuts alone also matter. At L>=12 the three methods give .52587/.53573/.57823 (linear/PCHIP/moment-linear); the last is rejected, the first is borderline p=.0109. At L>=16 they give .58377/.55849/.60945 with only3sizes and1DOF. Do not select a size cut simply for a larger p. Correction-power fits with omega=.5,.8,1 are retained as stress diagnostics, never precision validation.

Eight additional immutable analyses vary full/local window, cuts501/1001/2001, doubled block widths, U=.20/.25/.30, and individual-start parity. The exact source cohorts and per-fit omitted sizes, flags, statuses, validity, coefficients and intervals are in sensitivity-ledger.json and source-pinned analysis paths. These are correlated descriptive comparisons, not certified total95% intervals. Warmup/block changes do not remove the central interpolation disagreement; threshold variation can magnify it. Even/odd helper draw IDs are not paired for covariance or treated as independent new production. All72 existing sources have two starts per coordinate, but favorable pooled flags do not prove equilibration. No pooling across epsilon or fit to a known exponent was used.

REPRODUCIBILITY AND CHECKS
Finite-size-contours.csv records roots, moments, achieved Binder, validity and flags. Source snapshots retain original audit hashes and point roles. Saved zero_beta.py/measure.py and current complete scalar MC histories are sufficient to rerun this analysis. Seven focused oracle tests pass, including epsilon orientation/interval reversal, exact power normalization, covariance, missing brackets and size-cut masking. The figure was visually inspected. Independent scientific review is recorded separately in review.json.
'''
    (folder/'zero-beta-spin-dimension-report.txt').write_text(report)
    m.run.atomic_json(folder/'analysis.json',dict(status='pending_independent_review',snapshots={t:{k:v for k,v in s.items() if k!='analysis'} for t,s in snapshots.items()},source_count=len(a['source_jobs']),coordinate_count=len(a['groups']),headline='No robust precision estimate at beta=0 from present grid; conditional PCHIP .537[.507,.568], material interpolation sensitivity and rejected alternatives.'))
    m.run.atomic_json(ROOT/'latest-delivery.json',dict(folder=folder.name,status='pending_independent_review'));print(folder,flush=True)

if __name__=='__main__':main()
