"""Draft N48/N64 pilot output; scientific acceptance remains explicit."""
import json,shutil,statistics
import numpy as np
import measure as m
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

TAGS=['larger-primary','larger-full','larger-cut4001','larger-cut8001','larger-block2','larger-U20full','larger-U30full','larger-evenfull','larger-oddfull','larger-only','larger-only-late']

def main():
    root=m.OUT/'zero-beta';snapshots={}
    for tag in TAGS:
        path=root/json.loads((root/f'latest-{tag}.json').read_text())['folder']/'analysis.json';snapshots[tag]=dict(path=str(path),sha256=m.run.sha256(path),analysis=json.loads(path.read_text()))
    a=snapshots['larger-primary']['analysis'];plan=json.loads((m.OUT/'larger-zero-beta-pilot-design.json').read_text());diagnostics=json.loads((m.OUT/'larger-pilot-diagnostics.json').read_text());cohort=json.loads((m.OUT/'larger-pilot-analysis-cohort.json').read_text());expected={j['id'] for j in cohort['source_jobs']};pilot={p['id'] for p in plan['points']}
    for tag,s in snapshots.items():
        starts=s['analysis']['settings']['starts'];ids=pilot if tag.startswith('larger-only') else {j['id'] for j in cohort['source_jobs'] if starts=='all' or json.loads((m.CAMPAIGN/'jobs'/f"{j['id']}.json").read_text())['point']['replica']%2==int(starts=='odd')}
        assert {j['id'] for j in s['analysis']['source_jobs']}==ids
    folder=root/('larger-pilot-delivery-'+m.run.now()[:19].replace('-','').replace(':','').replace('T','-'));folder.mkdir()
    for name in ['report_larger_pilot.py','diagnose_larger_pilot.py','zero_beta.py','measure.py','run_larger_pilot.py','design_larger_zero_beta.py']:shutil.copyfile(m.OUT/name,folder/name)
    shutil.copyfile(m.OUT/'larger-pilot-diagnostics.json',folder/'trajectory-cost-and-topology-diagnostics.json')
    fits=[dict(variant=t,**f) for t,s in snapshots.items() for f in s['analysis']['fits']];contours=[dict(variant=t,**r) for t,s in snapshots.items() for r in s['analysis']['contours']]
    m.run.atomic_json(folder/'fit-sensitivity-ledger.json',fits);m.run.atomic_json(folder/'contour-ledger.json',contours)
    fig,(ax,bx)=plt.subplots(1,2,figsize=(12,5.2))
    for method,color,marker in [('linear','#777777','o'),('pchip','#2563a7','s'),('moment_linear','#b45e20','^')]:
        rs=sorted([r for r in a['contours'] if r['method']==method and r.get('location') and r.get('MC95')],key=lambda r:r['N'])
        for i,r in enumerate(rs):
            y=r['location'][1];lo,hi=r['MC95'][1];ax.errorbar(r['N'],y,yerr=[[max(0,y-lo)],[max(0,hi-y)]],fmt=marker,color=color,mfc=color if r['status']=='conditional_contour' else 'white',capsize=2,label=method if i==0 else None)
        for f in a['fits']:
            if f['method']!=method or f['moment']!='m2' or f['omega'] is not None or f['minimum']==6:continue
            r=f['result']
            if 'Delta' not in r:continue
            y=r['Delta'];lo,hi=r['MC95'];bx.errorbar(f['minimum'],y,yerr=[[max(0,y-lo)],[max(0,hi-y)]],fmt=marker,color=color,mfc=color if r['status']=='conditional_dimension_fit' else 'white',capsize=2)
    ax.set(xscale='log',yscale='log',xlabel='$L$',ylabel=r'$\langle |m|^2\rangle$ on nominal $U=1/4$ contour',title='3D beta=0: attained finite boxes');ax.set_xticks([8,12,16,24,32,48,64],labels=['8','12','16','24','32','48','64']);ax.xaxis.set_minor_formatter(NullFormatter());ax.grid(alpha=.15);ax.legend(frameon=False,fontsize=8)
    bx.axhline(.51905,color='#333333',ls='--',lw=1,label='WF reference, not imposed');bx.set(xlabel='Minimum included L',ylabel=r'effective $\Delta_{e^{i\phi}}$',title='Leading power: size-cut sensitivity');bx.set_xticks([8,12,16,24,32]);bx.grid(alpha=.15);bx.legend(frameon=False,fontsize=8)
    fig.suptitle('N48/N64 pilot — preliminary, pending independent scientific review');fig.text(.5,.015,'Conditional MC 95% bars. Open markers: unresolved contours or rejected/zero-DOF size fits. No automatic systematic or total-precision claim.',ha='center',fontsize=8);fig.tight_layout(rect=[0,.055,1,.94])
    for ext in ['png','svg','pdf']:fig.savefig(folder/f'N48-N64-pilot-spin-dimension.{ext}',dpi=180)
    plt.close(fig)
    lines=['N48/N64 BETA=0 SPIN-DIMENSION PILOT: DRAFT PENDING SCIENTIFIC REVIEW','',f"New audited trajectories={len(pilot)}; combined fine-grid source count={len(expected)}. No longer runs, N96, Delta or nu imposed. All pre-existing data and earlier deliveries remain preserved.",'','ACTUAL GENERATION AND AUDIT COST']
    cost={}
    for N in [48,64]:
        rr=[r for r in diagnostics['records'] if r['N']==N];tt=[r['generation_seconds']+r['audit_seconds'] for r in rr]
        cost[N]=dict(trajectories=len(rr),total_worker_hours=sum(tt)/3600,seconds_min_median_max=[min(tt),statistics.median(tt),max(tt)],full_original_bytes=sum(r['full_audit']['bytes'] for r in rr));lines.append(f'L{N}: '+json.dumps(cost[N]))
    lines+=['','ATTAINED BRACKETS AND PRIMARY LEADING FITS']
    for r in a['contours']:
        if r['N'] in [48,64]:lines.append(json.dumps(r))
    for f in a['fits']:
        if f['method']=='moment_linear' and f['moment']=='m2':lines.append(json.dumps(f))
    lines+=['','REQUIRED INTERPRETATION','This is a cost, coverage, mixing and size-leverage pilot. Nominal root existence, exact closure audit, accepted updates, small p or good p do not certify precision or mixing. Assess actual Binder reconstructed from moments, linear/PCHIP/moment-consistent differences, omitted sizes, bootstrap tails, start/half/late-cut diagnostics and all tested corrections before drawing a scientific conclusion. The .00025 epsilon spacing may be too coarse at larger L. Do not select a size cut or correction exponent to match the WF value .519.','The U=.20/.30 variations may lack smaller-size or new-size brackets; missing cases remain explicitly incomplete comparisons. Two-size leading fits and three-size correction fits have zero degrees of freedom and are diagnostic only. Corrections omega=.5,.8,1,1.5,2 and size cuts are stress tests, not measured universal correction powers or a certified systematic envelope.','Block-bootstrap histories preserve m2,m4,|m| covariance and each draws own Binder root. Intervals are conditional MC95, not total errors. Check topology tau against the moment block basis in trajectory-cost-and-topology-diagnostics.json; larger blocks or additional sampling may be needed before quoting errors. Event counts have denominators and are not a claim of effective independent sampling. Current source amplitudes/contact corrections and finite-variable drift may remain.','Resource log records sampled per-process RSS/CPU, not a measured peak. The initial planning forecast was 1–3h,17.2workerh using volume scaling only. Report measured time separately; no extrapolation of pilot speed to guaranteed asymptotic precision.','No additional runs are launched by this analysis. Production has closed the launch gate. User-authorized adaptive field retention follows all closure audits and diagnostic extraction, retaining every observable/winding/RNG/provenance/timing history, log and checkpoint; SSD originals remain held. Independent scientific review must inspect headline, absolute effects, noise, scope and coverage before this draft is accepted.']
    (folder/'N48-N64-pilot-report.txt').write_text('\n'.join(lines)+'\n')
    m.run.atomic_json(folder/'analysis.json',dict(status='pending_independent_review',snapshots={t:{k:v for k,v in s.items() if k!='analysis'} for t,s in snapshots.items()},cost=cost,source_ids=sorted(expected),pilot_source_ids=sorted(pilot),scope='Attained pilot diagnostics only; no controlled total/systematic uncertainty claimed without scientific review.'))
    m.run.atomic_json(root/'latest-larger-pilot-delivery.json',dict(folder=folder.name,status='pending_independent_review'));print(folder,flush=True)

if __name__=='__main__':main()
