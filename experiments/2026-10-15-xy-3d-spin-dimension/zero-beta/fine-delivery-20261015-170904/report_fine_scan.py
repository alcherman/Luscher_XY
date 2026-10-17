"""Source-pinned report for the finer epsilon scan, without fixing an exponent."""
from pathlib import Path
import json,csv,shutil
import numpy as np
import measure as m
import matplotlib.pyplot as plt

TAGS=['fine-primary','fine-full','fine-cut2001','fine-cut4001','fine-block2','fine-U20','fine-U30','fine-even','fine-odd','fine-pooled','fine-evenfull','fine-oddfull','fine-U20full','fine-U30full','fine-pooledfull']

def main():
    root=m.OUT/'zero-beta';snapshots={}
    for tag in TAGS:
        path=root/json.loads((root/f'latest-{tag}.json').read_text())['folder']/'analysis.json'
        snapshots[tag]=dict(path=str(path),sha256=m.run.sha256(path),analysis=json.loads(path.read_text()))
    a=snapshots['fine-primary']['analysis'];plan=json.loads((m.OUT/'fine-zero-beta-stage1-design.json').read_text());expected={p['id'] for p in plan['points']}
    assert {p['id'] for p in a['source_jobs']}==expected
    oldpath=Path(plan['source_analysis']);old=json.loads(oldpath.read_text())
    assert {p['id'] for p in snapshots['fine-pooled']['analysis']['source_jobs']}==expected|{p['id'] for p in old['source_jobs']}
    for tag,s in snapshots.items():
        if tag in ['fine-pooled','fine-pooledfull']:
            assert {p['id'] for p in s['analysis']['source_jobs']}==expected|{p['id'] for p in old['source_jobs']}
            continue
        starts=s['analysis']['settings']['starts']
        subset={p['id'] for p in plan['points'] if starts=='all' or p['replica']%2==int(starts=='odd')}
        assert {p['id'] for p in s['analysis']['source_jobs']}==subset
    folder=root/('fine-delivery-'+m.run.now()[:19].replace('-','').replace(':','').replace('T','-'));folder.mkdir()
    for name in ['report_fine_scan.py','zero_beta.py','measure.py','run_fine_scan.py','design_fine_zero_beta.py']:shutil.copyfile(m.OUT/name,folder/name)
    ledger=[dict(variant=tag,**f) for tag,s in snapshots.items() for f in s['analysis']['fits']]
    roots=[dict(variant=tag,**r) for tag,s in snapshots.items() for r in s['analysis']['contours']]
    m.run.atomic_json(folder/'sensitivity-ledger.json',ledger);m.run.atomic_json(folder/'contour-ledger.json',roots)
    rows=[]
    for r in a['contours']:
        loc=r.get('location');rows.append(dict(N=r['N'],method=r['method'],status=r['status'],epsilon_over_pi=r.get('epsilon_over_pi'),m2=loc[1] if loc else None,mabs=loc[2] if loc else None,valid_fraction=r['valid_fraction'],achieved_Binder=r.get('moment_ratio_at_interpolated_root'),flags=';'.join(r['flags'])))
    with (folder/'finite-size-contours.csv').open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    headline=[f for f in a['fits'] if f['minimum']==8 and f['omega'] is None]
    fig,(ax,bx)=plt.subplots(1,2,figsize=(12,5))
    styles=[('linear','#777777','o','Binder-linear'),('pchip','#2563a7','s','PCHIP'),('moment_linear','#b45e20','^','Moment-consistent linear')]
    for method,color,marker,label in styles:
        for i,r in enumerate(sorted([r for r in a['contours'] if r['method']==method and r.get('location')],key=lambda r:r['N'])):
            y=r['location'][1];lo,hi=r['MC95'][1];ax.errorbar(r['N'],y,yerr=[[max(0,y-lo)],[max(0,hi-y)]],fmt=marker,color=color,capsize=2,label=label if i==0 else None)
        for f in a['fits']:
            if (f['method'],f['moment'],f['omega'])!=(method,'m2',None):continue
            r=f['result']
            if 'Delta' not in r:continue
            y=r['Delta'];lo,hi=r['MC95'];bx.errorbar(f['minimum'],y,yerr=[[max(0,y-lo)],[max(0,hi-y)]],fmt=marker,color=color,mfc=color if r['status']=='conditional_dimension_fit' else 'white',capsize=2)
    ax.set(xscale='log',yscale='log',xlabel='$L$',ylabel=r'$\langle |m|^2\rangle$ on nominal $U=1/4$ contour',title='Finer epsilon scan: beta=0');ax.set_xticks([8,12,16,24,32],labels=['8','12','16','24','32']);ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.15)
    bx.axhline(.51905,color='#333333',ls='--',lw=1,label='3D XY reference, not fitted');bx.set(xlabel='Minimum included L',ylabel=r'effective $\Delta_{e^{i\phi}}$',title='Size-cut and interpolation diagnostics');bx.set_xticks([6,8,12,16]);bx.legend(frameon=False,fontsize=8);bx.grid(alpha=.15)
    fig.suptitle('3D beta=0 spin scaling dimension: new finer-grid Monte Carlo');fig.text(.5,.015,'Conditional MC 95% bars; open symbols mark rejected or zero-DOF fits. Finite-size and interpolation checks are separate from precision.',ha='center',fontsize=8);fig.tight_layout(rect=[0,.055,1,.94])
    for ext in ['png','svg','pdf']:fig.savefig(folder/f'3d-zero-beta-fine-spin-dimension.{ext}',dpi=180)
    plt.close(fig)
    lines=['FINER 3D BETA=0 EPSILON SCAN: SPIN SCALING DIMENSION','', 'Conclusion: interpolation alternatives are now much closer descriptively, but finite-size/correction sensitivity dominates and a controlled two-significant-digit dimension is NOT established. No discrepancy from the standard O2 Wilson-Fisher value is demonstrated.',f"Sources: {len(a['source_jobs'])} fresh audited trajectories on {len(a['groups'])} exact coordinates,12000cycles each, two independent RNG starts with zero/unit-x winding. All original72 trajectories remain separate; pooled variant has154 sources.",'No Delta or nu is fixed. Predeclared primary family is moment-consistent linear,U=.25,L>=8. Conditional and rejected fit statuses are retained; no estimator or size cut is selected by favorable p.','', 'LEADING L>=8 CONDITIONAL DIAGNOSTICS']
    for f in headline:
        r=f['result'];lines.append(f"{f['method']} {f['moment']}: {json.dumps(r)}; flags={f['input_flags']}; omitted={f['omitted_sizes']}")
    lines+=['','SIZE AND CORRECTION LIMITATIONS','For moment-consistent m2 fits, minimum L=8/12/16 gives Delta=.505400/.530758/.550138. The nominal MC95 upper bound at L>=8 is .518262, slightly below the .51905 reference, but the much larger size/correction effects preclude interpreting this as exclusion of Wilson-Fisher behavior. Adding an assumed amplitude correction with omega=.5/.8/1 gives central Delta=.731199/.647079/.619292 and very broad intervals. These high-p correction fits are weakly constrained and are not chosen as measurements; they show that an asymptotic intercept is not controlled by five small sizes. An interval spanning leading-fit alternatives alone is not total uncertainty.']
    lines+=['','INTERPOLATION SUPPORT']
    for method,_,_,_ in styles:
        rs=[r for r in a['contours'] if r['method']==method];dis=[abs(r['moment_ratio_at_interpolated_root']-.25) for r in rs if r.get('moment_ratio_at_interpolated_root') is not None]
        lines.append(f"{method}: {sum(r['status']=='conditional_contour' for r in rs)}/{len(rs)} conditional roots; max actual-Binder discrepancy {max(dis) if dis else None}; min valid fraction {min(r['valid_fraction'] for r in rs)}")
    lines+=['','COVERAGE QUALIFICATIONS','Full-window follow-ups restore all five sizes for separate starts; their moment-linear exponents are .494535 and .512559. Do not infer a significant start difference by pairing shared helper-draw indices. U=.20 remains unbracketed at L16 and has only93% valid L24 draws; U=.30 retains only68.4% at L32. Threshold robustness across all sizes is therefore incomplete, even though the primary U=.25 roots all retain100%. The pooled full-window follow-up restores all six roots including oldL6, but its L>=8 moment-linear size fit is rejected (p=.007616); it is diagnostic. Sources at floating-roundoff-close legacy coordinates remain explicit, which can affect local window support; they are not silently rounded onto the fine grid.','', 'LIMITATIONS','All intervals are conditional block-bootstrap MC95, not certified total errors. m2~L^(-2Delta) and recorded|m|~L^(-Delta) at a fixed dimensionless contour require continuous critical scaling and infrared dominance2Delta<3. Corrections/contact L^-3 terms and critical-contour drift remain possible. GLS uses joint root/moment samples and masks after size cuts. Leading correction powers .5/.8/1 are diagnostic stress tests, not fixed measured universal exponents.','Cut1001/2001/4001, doubled blocks, full/local windows, U=.20/.25/.30, start parity and linear/PCHIP/moment-linear are preserved in the ledger. Missing roots/omitted sizes and start/half/low-block flags are not hidden. Even/odd results are not paired for start-difference covariance. Favorable audits/bracket retention/p-values alone do not prove mixing or precision.','The new fine grid changes source composition and sampling precision; old-versus-new shifts are descriptive, not independent replicated improvement claims. The original paper motivation was clarified by the user: this O2 study compares to the standard WF value, not a claimed second O2CFT. Reference Delta .51905 from eta .0381 (Campostrini et al. https://arxiv.org/abs/cond-mat/0605083) is for comparison only.','All observable/winding/RNG/provenance/log histories and checkpoints retained. User-authorized field thinning follows full independent closure audits; SSD originals are held and no archive relocation/deletion is authorized.','All fifteen snapshot hashes and complete full/parity/pooled source sets are pinned in analysis.json. Figure and substantive interpretation require independent review; numerical output alone is not completion.']
    (folder/'fine-scan-measurement-report.txt').write_text('\n'.join(lines)+'\n')
    m.run.atomic_json(folder/'analysis.json',dict(status='pending_independent_review',snapshots={t:{k:v for k,v in s.items() if k!='analysis'} for t,s in snapshots.items()},old_source=str(oldpath),old_source_sha256=m.run.sha256(oldpath),headline_fits=headline,source_ids=sorted(expected),scope='Finer epsilon scan conditional effective spin dimension; no automated universality or total precision claim.'))
    m.run.atomic_json(root/'latest-fine-delivery.json',dict(folder=folder.name,status='pending_independent_review'));print(folder,flush=True)

if __name__=='__main__':main()
