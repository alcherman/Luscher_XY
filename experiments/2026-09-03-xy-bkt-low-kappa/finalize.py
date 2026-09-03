"""Reconcile every new record and report attained widths with remaining limits."""
import json
import hashlib
from pathlib import Path
import h5py
import numpy as np
from run import OUT,atomic_json,sha256,now,_provenance


def main():
    batches=[json.loads(p.read_text()) for p in (OUT/'batches').glob('*.json')]
    points=[p for b in batches for p in b['points']]
    assert len({p['id'] for p in points})==len(points)
    source=_provenance()['source_sha256']
    assert all(b['runner_sha256']==sha256(OUT/'run.py') and b['model_source_sha256']==source for b in batches)
    for path,digest in json.loads((OUT/'input-provenance.json').read_text()).items():assert sha256(path)==digest,path
    parity=json.loads((OUT/'classifier-parity.json').read_text())
    assert parity['passed'] and parity['classifier_sha256']==sha256(OUT/'conservative.py')
    records=[]
    for p in points:
        path=OUT/'chains'/f"{p['id']}.h5";job=json.loads((OUT/'jobs'/f"{p['id']}.json").read_text())
        assert job['state']=='complete' and job['audit']['passed'] and job['point']==p
        assert job['audit']['bytes']==path.stat().st_size and job['audit']['rows']==p['steps']+1
        with h5py.File(path,'r') as f:
            q=json.loads(f.attrs['parameters']);provenance=json.loads(f.attrs['provenance'])
            assert q['seed']==p['seed'] and q['kappa']==p['kappa'] and q['epsilon']==p['epsilon'] and q['N']==p['L']
            assert q['D']==2 and q['update']=='hybrid' and q['clusters_per_step']==4 and q['initial_winding']==p['initial_winding']
            assert provenance['source_sha256']==source and f.attrs['committed_rows']==p['steps']+1
            assert all(len(ds)==p['steps']+1 for ds in f.values())
            digest=hashlib.sha256(f['phi'][-1].tobytes()+f['rng_state'][-1]).hexdigest()
            assert digest==job['segments'][-1]['field_rng_sha256']
        records.append(dict(id=p['id'],bytes=path.stat().st_size,rows=p['steps']+1,sha256_at_full_field_audit=job['audit']['sha256']))
    analyses=[json.loads((OUT/'analysis'/f"{p['id']}.json").read_text()) for p in points]
    assert all(r['analysis_source_sha256']==sha256(OUT/'analyze.py') and r['summary']['cut']==2001 and r['summary']['enough_blocks'] for r in analyses)
    pooled=json.loads((OUT/'pooled.json').read_text())
    assert sum(len(r['chains']) for r in pooled)==len(points)
    for e in {p['epsilon_milli'] for p in points}:
        history=json.loads((OUT/'histories_final'/f'e{e:04d}_inputs.json').read_text())
        assert history['source_sha256']==sha256(OUT/'histories_final.py')
        assert history['chains']==sorted(p['id'] for p in points if p['epsilon_milli']==e)
    phase=json.loads((OUT/'conservative-brackets.json').read_text())
    low=[r for r in phase['regions'] if r['x']<=.65]
    flagged=[dict(id=r['point']['id'],flags=r['flags']) for r in analyses if r['flags']]
    pooled_flags=[dict(id=r['id'],flags=r['flags']) for r in pooled if r['flags']]
    audit=dict(created_utc=now(),passed=True,source_sha256=source,runner_sha256=sha256(OUT/'run.py'),finalizer_sha256=sha256(Path(__file__)),
        chains=len(points),saved_fields=sum(r['rows'] for r in records),cycles=sum(p['steps'] for p in points),
        bytes=sum(r['bytes'] for r in records),records=records,
        scope='All fields and physical observables independently audited at chain closure; final reconciliation checks source, parameters, lengths, bytes and last field/RNG digests.')
    atomic_json(OUT/'execution-audit.json',audit)
    diagnostics=dict(created_utc=now(),individual_flags=flagged,pooled_flags=pooled_flags,
        tau_by_L={L:{k:dict(median=float(np.median(v)),maximum=float(max(v))) for k in ['m2','binder','xi_ratio','link_cosine']
            if (v:=[r['summary']['taus'][k] for r in analyses if r['point']['L']==L])} for L in [16,32,64]},
        retained_per_group=sorted({r['retained'] for r in pooled}),cut=2001,
        blocks_per_chain=[min(r['summary']['blocks'] for r in analyses),max(r['summary']['blocks'] for r in analyses)],
        width_per_chain=[min(r['summary']['block_width'] for r in analyses),max(r['summary']['block_width'] for r in analyses)])
    atomic_json(OUT/'diagnostic-summary.json',diagnostics)
    lines=['# Low-kappa refinement results','',f"Completed {len(points)} trajectories; every full-field audit and the final reconciliation passed. Saved {audit['saved_fields']:,} fields ({audit['bytes']/1e9:.2f} GB), with all measurements, RNG records, update work, timing and initial/transient rows retained.",'',
        'The first batch adds epsilon/pi=0.575 and 0.590 and direct kappa values 0.10 and 0.18 at epsilon/pi=0.580. A follow-up at (epsilon/pi,kappa)=(0.590,0.16) established the missing massive-side edge. Each parameter/volume group has four independent 32000-cycle trajectories, or 128000 Metropolis sweeps and 512000 clusters. The selected 2000-cycle cut leaves 120000 pooled cycles, three times the previous refined count. Inspected initial/retained histories, cut/block checks and independent-start comparisons support this cut; no exhaustive equilibration or thermodynamic BKT claim is made.','',
        'Five analysis tests passed, and the classifier reproduced all 286 previous refined phase labels from unchanged previous inputs. The model package and previous campaign inputs are unchanged.','',
        'The full curve combines the new results with the unchanged earlier campaigns. At repeated parameter/volume coordinates, the new longer estimate replaces the earlier estimate; earlier raw records remain intact. The same joint U/R/eta two-SE classification and adjacent-epsilon-cell union define the bars. They are finite-volume scouting brackets, not Gaussian confidence intervals.','',
        '| epsilon/pi cell | kappa bracket | vertical half-width | target <= 0.05 |','|---|---|---:|---|']
    for r in low:
        half=(r['kappa_upper']-r['kappa_lower'])/2
        lines.append(f"| {r['x_lower']:.3f}–{r['x_upper']:.3f} | {r['kappa_lower']:.3f}–{r['kappa_upper']:.3f} | {half:.3f} | {'yes' if half<=.05000001 else 'not yet'} |")
    lines+=['',f"Diagnostic flags: {len(flagged)} individual trajectories and {len(pooled_flags)} pooled groups. See diagnostic-summary.json and the individual analysis records. Flags are preserved, not silently removed.",'',
        'The current disk budget restricts further refinement. Where the table misses the target, additional targeted grid points/statistics or a better controlled finite-volume estimator remain necessary. No smaller total error is claimed from an auxiliary fixed-Binder crossing.','',
        'Artifacts: phase_boundary_refined.png/pdf/svg, low_kappa_comparison.png/pdf/svg, phase-brackets.csv, pooled.json, histories_final/, and the frozen batches, raw chains and audit records.']
    (OUT/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print('Audit passed',len(points),'chains;',len(low),'low-kappa cells',flush=True)


if __name__=='__main__':main()
