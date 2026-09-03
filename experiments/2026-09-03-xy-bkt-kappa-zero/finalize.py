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
    old=json.loads((OUT.parent/'2026-09-03-xy-bkt-low-kappa/conservative-brackets.json').read_text())['kappa_zero_endpoint']
    endpoint=phase['kappa_zero_endpoint']
    labels=[r for r in phase['labels'] if r['kappa_milli']==0 and r['epsilon_milli'] in [550,570,575]]
    lines=['# Kappa-zero endpoint refinement','',
        f"Completed {len(points)} trajectories on L=16,32,64 at epsilon/pi=0.550,0.570,0.575, all at kappa=0. Each coordinate/volume uses two independent 7000-cycle starts (cold and winding (1,0)). Every cycle contains one Metropolis sweep, four fixed Wolff clusters and one winding proposal. The 2000-cycle cut leaves 10000 pooled cycles per coordinate/volume.",'',
        f"The previous horizontal bracket was [{old['lower']:.3f}, {old['upper']:.3f}]. The updated bracket is [{endpoint['lower']:.3f}, {endpoint['upper']:.3f}], with plotting midpoint {(endpoint['upper']+endpoint['lower'])/2:.4f} and half-width {(endpoint['upper']-endpoint['lower'])/2:.4f}. Its width is smaller by a factor {(old['upper']-old['lower'])/(endpoint['upper']-endpoint['lower']):.2f}.",'',
        '| New epsilon/pi | Joint U/R/eta label | eta from L32,64 |','|---:|---|---:|']
    for r in labels:
        lines.append(f"| {r['epsilon_milli']/1000:.3f} | {r['label']} | {r.get('eta_32_64',float('nan')):.5f} +/- {r.get('eta_error',float('nan')):.5f} |")
    lines+=['',
        'The same joint U/R/eta two-standard-error rule is used, including unresolved coordinates. These are finite-volume scouting brackets; they are not Gaussian confidence intervals or a thermodynamic extrapolation. Existing higher-statistics estimates at other coordinates remain in use. This small added batch has lower statistics per coordinate than the longer low-kappa campaign.', '',
        f"All {len(points)} full-field audits and the final reconciliation passed: {audit['saved_fields']:,} saved fields, {audit['bytes']/1e9:.3f} GB, with initial/transient fields, physical observables, RNG, work and timing records retained. Five analysis tests passed. All {parity['previous_labels']} prior phase classifications were reproduced from unchanged input files. The model package is unchanged.",'',
        f"Diagnostic flags: {len(flagged)} individual trajectories and {len(pooled_flags)} pooled groups. Blocks per trajectory range from {diagnostics['blocks_per_chain'][0]} to {diagnostics['blocks_per_chain'][1]}. See diagnostic-summary.json and per-chain analysis for autocorrelations, alternate cuts, doubled blocks, retained-half checks and independent-start comparisons. Initial/full/retained histories are saved in histories_final/. Visual inspection of the epsilon/pi=0.570 initial and retained histories supports the cut; doubled-block standard errors range from 0.862 to 1.119 times the primary estimates. See thermalization-review.json.",'',
        'Artifacts: phase_boundary_refined.png/pdf/svg, endpoint_diagnostics.png/pdf/svg, phase-brackets.csv, endpoint-observables.csv, pooled.json, execution-audit.json, and complete frozen manifests and raw chains. The separate target of shrinking all low-kappa vertical bars remains unfinished; this batch targets only the kappa=0 horizontal bracket.']
    (OUT/'REPORT.md').write_text('\n'.join(lines)+'\n')
    print('Audit passed',len(points),'chains; endpoint',endpoint,flush=True)


if __name__=='__main__':main()
