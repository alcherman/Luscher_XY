"""Source-pinned delivery of common-coupling xi/L quotients, without new MC."""
import argparse
import csv
import json
import shutil
from pathlib import Path
import numpy as np
import common_xi_quotients as x


def main(baseline,block4):
    sources={name:Path(path) for name,path in [('baseline',baseline),('block4',block4)]}
    data={name:json.loads((path/'analysis.json').read_text()) for name,path in sources.items()}
    assert {v['id'] for v in data['baseline']['source_jobs']}=={v['id'] for v in data['block4']['source_jobs']}
    folder=x.m.OUT/('common-xi-delivery-'+x.m.run.now()[:19].replace('-','').replace(':','').replace('T','-'));folder.mkdir()
    shutil.copyfile(__file__,folder/'report_common_xi.py');shutil.copyfile(x.m.OUT/'common_xi_quotients.py',folder/'plotting_source.py')
    ledger=[]
    for name,a in data.items():
        for raw in a['pairs']:
            r=dict(raw,variant=name)
            if 'epsilon_over_pi' in r:
                native=[]
                for N in [r['N'],r['N2']]:
                    q=np.array([g['q'] for g in a['groups'] if g['N']==N and 'mean' in g['summary']]);j=max(0,min(len(q)-2,np.searchsorted(q,r['epsilon_over_pi'],side='right')-1))
                    native.append(dict(N=N,knots=q[j:j+2].tolist(),width=float(q[j+1]-q[j])))
                r['native_interpolation_brackets']=native
            ledger.append(r)
    a=data['baseline'];curves={}
    for g in a['groups']:
        if 'R' in g['summary']:curves.setdefault(g['N'],[]).append((g['q'],g['summary'],None))
    x.plot(folder,curves,[r for r in ledger if r['variant']=='baseline'])
    x.m.run.atomic_json(folder/'quotient-ledger.json',ledger)
    with (folder/'quotients.csv').open('w') as f:
        cols=['variant','N','N2','method','window','status','epsilon_over_pi','xi_over_L','Delta_eff','SD','valid_fraction','input_flags','native_interpolation_brackets']
        w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore');w.writeheader();w.writerows(ledger)
    text=['3D beta=0 spin dimension from common-epsilon xi/L crossings','',
          'METHOD AND CURRENT CONCLUSION',
          'This is the standard phenomenological-renormalization quotient construction: solve xi(L)/L=xi(2L)/(2L) at one epsilon, then Delta_eff=-log(m2(2L)/m2(L))/(2 log 2) at that same epsilon. It does not impose Delta or a critical xi/L value. Available estimates are finite-size effective dimensions. No final asymptotic Delta +/- stat +/- sys is established by this report.',
          'Both sizes use the same epsilon, not separate U=.25 contours. This analysis supersedes treating the earlier arbitrary Binder-fit/correction spread as a systematic standard deviation. The new data do not retrospectively validate that procedure.',
          '', 'ATTAINED QUOTIENTS (one bootstrap SD; conditional on interpolation)',
          'pair      linear Delta_eff +/- SD       PCHIP Delta_eff +/- SD       validity linear/PCHIP']
    for N in sorted({r['N'] for r in ledger}):
        rs={r['method']:r for r in ledger if r['variant']=='baseline' and r['window']=='local' and r['N']==N}
        def value(r):return f"{r['Delta_eff']:.6f} +/- {r['SD']:.6f}" if 'Delta_eff' in r else r['status']
        text.append(f"{N}:{2*N}  {value(rs['linear'])}   {value(rs['pchip'])}   {rs['linear'].get('valid_fraction')}/{rs['pchip'].get('valid_fraction')}")
    text+=['', 'INTERPOLATION AND COVERAGE',
           'Crossing status and retained fractions are reported separately for every pair/method/variant. Unresolved pairs are excluded from fits and shown open in the plot. No missing root is extrapolated. Initial central crossings found by interpolating R itself are not sufficient to establish moment-consistent bootstrap coverage.',
           'These are finite-pair estimates, not extrapolated CFT measurements. The small-pair and interpolation differences remain visible; no fit is selected for proximity to the expected Wilson-Fisher result.',
           'Native interpolation brackets, not just union search windows, are in quotient-ledger.json and quotients.csv. Another size’s finer knots do not improve a coarser size’s native interpolation interval. Inspect those intervals rather than assuming the union-knot spacing is the resolution for both sizes.',
           'Raw S0,F,m2 are interpolated linearly or with PCHIP, then the correlation-length ratio is formed; each draw solves its own root and evaluates m2 at that root. Exact piecewise polynomial-numerator roots detect multiple or tangent roots; absent, nonphysical or ambiguous roots are invalid. No isotonic smoothing or fixed expected exponent is used. Only <=8-ulp representational duplicates are pooled; original source coordinates remain serialized.',
           '', 'STATISTICAL AND FINITE-SIZE QUALIFICATIONS',
           f"{len(a['source_jobs'])} completed/audited trajectories contribute to the scalar-history cohort; per-pair contributing IDs and flags are serialized. Baselinecut={a['settings']['cut']},draws={a['settings']['draws']}. Blocks are factors {a['settings']['factor']} and {data['block4']['settings']['factor']} times max(20,ceil(10 max measured tau)) for m2,m4,|m|,S0,F and Binder/R influences. All joint moments share resampled blocks; independent coordinates have independent draws. Reusing a size across overlapping pairs reuses its draws, preserving their covariance. Comparing block widths is a finite-data diagnostic, not a mixing proof.",
           'Low-block-count flags occur at supporting coordinates, particularly the old shorter anchor runs. Start/half diagnostics and every source audit/provenance record are retained in the analysis snapshots. Statistical SDs and percentile MC95 are conditional; no interpolation/systematic/total coverage is implied.',
           'Diagnostic GLS extrapolations use Delta_eff(L)=Delta+a L^(-.789). The omega choice is Wilson-Fisher-informed, not measured from this study. Fits recompute a common validity mask and covariance after each size cut, record omitted pairs, and mark zeroDOF/rejected cases. These extrapolations do not supply a calibrated sigma_sys or choose a preferred answer.',
           '', 'VALIDATION AND SCOPE',
           'Focused oracle/analysis tests are recorded in common-xi-tests.log, including exact normalization, shared root/magnetization resampling, no extrapolation, two hidden roots with same-sign endpoints, shared-pair covariance, size-cut validity masks and roundoff-only coordinate merging. Tests establish numerical properties, not physical accuracy. Independent scientific review is recorded separately in review.json.',
           'This report analyses completed trajectories only. No numerical kernel, retention policy, original raw archive or previous result was changed by analysis. No additional production is authorized by this report.',
           '', 'PRIMARY METHOD REFERENCE',
           'Ballesteros et al., Finite size effects on measures of critical exponents in d=3 O(N) models, section3: https://arxiv.org/html/cond-mat/9606203',
           '', 'SOURCE SNAPSHOTS']
    for name,path in sources.items():text.append(f"{name}: {path} SHA256(analysis.json)={x.m.run.sha256(path/'analysis.json')}")
    (folder/'common-coupling-spin-dimension-report.txt').write_text('\n'.join(text)+'\n')
    x.m.run.atomic_json(folder/'provenance.json',dict(status='pending_independent_review',source_analyses={k:dict(path=str(p),sha256=x.m.run.sha256(p/'analysis.json')) for k,p in sources.items()},test_log=str(x.m.OUT/'common-xi-tests.log'),sources=len(a['source_jobs'])))
    x.m.run.atomic_json(x.m.OUT/'latest-common-xi-delivery.json',dict(folder=folder.name,status='pending_independent_review'))
    print(folder,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('baseline');p.add_argument('block4');a=p.parse_args();main(a.baseline,a.block4)
