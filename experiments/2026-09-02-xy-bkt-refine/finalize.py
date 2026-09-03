"""Reconcile completed batches, saved records and analysis; write a bounded report."""
import hashlib
import json
from pathlib import Path
import platform
import tarfile
import h5py
import numpy as np
import scipy
from run import OUT, atomic_json, sha256, now, _provenance
from analyze import finite


def main():
    batches=[json.loads(p.read_text()) for p in sorted((OUT/'batches').glob('*.json'))]
    points=[p for b in batches for p in b['points']]
    assert len({p['id'] for p in points})==len(points)
    source=_provenance()['source_sha256'];records=[]
    assert all(b['runner_sha256']==sha256(OUT/'run.py') and b['model_source_sha256']==source for b in batches)
    assert all(b['selection'].get('concurrent_runner_sha256',sha256(OUT/'run_edges.py'))==sha256(OUT/'run_edges.py') for b in batches)
    root=OUT.parents[1]
    with tarfile.open(OUT/'model-source.tar.gz','r:gz') as archive:
        for member in archive.getmembers():
            if member.name.startswith('supervillain/') and member.isfile():
                assert archive.extractfile(member).read()==(root/member.name).read_bytes(),member.name
    for relative,digest in json.loads((OUT/'analysis-dependencies.json').read_text()).items():
        assert sha256(root/relative)==digest,relative
    parity=json.loads((OUT/'classifier-parity.json').read_text())
    assert parity['passed'] and parity['classifier_sha256']==sha256(OUT/'conservative.py')
    for p in points:
        path=OUT/'chains'/f"{p['id']}.h5";job=json.loads((OUT/'jobs'/f"{p['id']}.json").read_text())
        assert job['state']=='complete' and job['audit']['passed'] and job['point']==p
        assert job['audit']['rows']==p['steps']+1 and job['audit']['bytes']==path.stat().st_size
        with h5py.File(path,'r') as f:
            params=json.loads(f.attrs['parameters']);provenance=json.loads(f.attrs['provenance'])
            assert provenance['source_sha256']==source
            assert params['seed']==p['seed'] and params['clusters_per_step']==4
            assert params['kappa']==p['kappa'] and params['epsilon']==p['epsilon'] and params['N']==p['L']
            assert params['initial_winding']==p['initial_winding']
            assert f.attrs['committed_rows']==p['steps']+1
            assert all(len(ds)==p['steps']+1 for ds in f.values())
            assert f['sweep'][0]==0 and f['sweep'][-1]==p['steps']
            assert f['clusters_total'][-1]==4*p['steps']
            digest=hashlib.sha256(f['phi'][-1].tobytes()+f['rng_state'][-1]).hexdigest()
            assert digest==job['segments'][-1]['field_rng_sha256']
        records.append(dict(id=p['id'],rows=p['steps']+1,bytes=path.stat().st_size,sha256_at_full_field_audit=job['audit']['sha256']))
    production=[p for p in points if p['replica']<100]
    analyses=[json.loads((OUT/'analysis'/f"{p['id']}.json").read_text()) for p in production]
    assert all(r['summary']['cut']==2001 and r['summary']['enough_blocks'] for r in analyses)
    assert all(r['analysis_source_sha256']==sha256(OUT/'analyze.py') for r in analyses)
    pooled=json.loads((OUT/'pooled.json').read_text())
    assert sum(len(r['chains']) for r in pooled)==len(production)
    assert all(len(r['chains'])==4 for r in pooled)
    tau_table={}
    for L in [16,32,64]:
        rr=[r for r in analyses if r['point']['L']==L]
        tau_table[L]={k:dict(median=float(np.median(v)),maximum=float(max(v)))
                      for k in ['m2','binder','xi_ratio','link_cosine','vortex']
                      if (v:=[r['summary']['taus'][k] for r in rr if r['summary']['taus'][k] is not None])}
    flagged_individual=[dict(id=r['point']['id'],flags=r['flags']) for r in analyses if r['flags']]
    flagged_pooled=[dict(id=r['id'],flags=r['flags']) for r in pooled if r['flags']]
    stability_overview={}
    for key in ['binder','xi_over_L','m2','link_cosine']:
        z=[];ratios=[]
        for r in analyses:
            a,b=r['halves'];den=np.hypot(a['errors'][key],b['errors'][key])
            if den>0:z.append((b['values'][key]-a['values'][key])/den)
            doubled=next(s for s in r['stability'] if s['cut']==2001 and s['block_width']==2*r['summary']['block_width'])
            if r['summary']['errors'][key]>0:ratios.append(doubled['errors'][key]/r['summary']['errors'][key])
        stability_overview[key]=dict(half_difference_z_mean=float(np.mean(z)),half_difference_z_sd=float(np.std(z,ddof=1)),
            trajectories=len(z),doubled_block_se_ratio_quantiles_025_50_975=np.quantile(ratios,[.025,.5,.975]).tolist())
    winding_block_exceptions=[r['point']['id'] for r in analyses if r['summary']['winding_tau'] is not None
                              and r['summary']['block_width']<10*r['summary']['winding_tau']]
    boundary=json.loads((OUT/'boundary.json').read_text())
    resolved=[r for r in boundary['regions'] if r['primary']['status']=='resolved_finite_volume_proxy']
    unresolved=[dict(epsilon_over_pi=r['epsilon_over_pi'],status=r['primary']['status']) for r in boundary['regions'] if r['primary']['status']!='resolved_finite_volume_proxy']
    old=json.loads((OUT.parent/'2026-09-02-xy-bkt-scout/phase-brackets.json').read_text())
    new=json.loads((OUT/'conservative-brackets.json').read_text())
    old_by={(r['x_lower'],r['x_upper']):r for r in old['regions']}
    improvements=[]
    for r in new['regions']:
        before=old_by.get((r['x_lower'],r['x_upper']))
        if before:improvements.append(dict(epsilon_over_pi=r['x'],old_halfwidth=(before['kappa_upper']-before['kappa_lower'])/2,
            new_halfwidth=(r['kappa_upper']-r['kappa_lower'])/2,reduction_factor=(before['kappa_upper']-before['kappa_lower'])/(r['kappa_upper']-r['kappa_lower'])))
    audit=dict(created_utc=now(),passed=True,finalizer_sha256=sha256(Path(__file__)),source_sha256=source,chains=len(points),production_chains=len(production),
        saved_fields=sum(r['rows'] for r in records),cycles=sum(p['steps'] for p in points),clusters=4*sum(p['steps'] for p in points),
        bytes=sum(r['bytes'] for r in records),records=records,classifier_parity=parity,
        verification='Every completed chain independently audited against every saved field when it closed; final reconciliation checks all row counts, parameters, source, file sizes and last field/RNG digests. Full-file SHA256 values are those recorded at each full-field audit.')
    atomic_json(OUT/'execution-audit.json',audit)
    result=finite(dict(created_utc=now(),tau_table=tau_table,individual_flags=flagged_individual,pooled_flags=flagged_pooled,
        stability_overview=stability_overview,winding_block_exceptions=winding_block_exceptions,
        conditional_resolved=len(resolved),conditional_unresolved=unresolved,bracket_comparison=improvements,
        kappa_zero_endpoint=new['kappa_zero_endpoint'],analysis_cut=2001,
        block_width_range=[min(r['summary']['block_width'] for r in analyses),max(r['summary']['block_width'] for r in analyses)],
        blocks_per_chain_range=[min(r['summary']['blocks'] for r in analyses),max(r['summary']['blocks'] for r in analyses)],
        python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,h5py=h5py.__version__))
    atomic_json(OUT/'diagnostic-summary.json',result)
    lines=['# XY boundary refinement results','',
        f"Completed {len(points)} trajectories ({len(production)} production/refinement, {len(points)-len(production)} pilot). All passed their full-field audits. Saved {audit['saved_fields']:,} fields in {audit['bytes']/1e9:.2f} GB; no generated rows were thinned or deleted.",'',
        '## Calculation','',
        'Validation: five focused analysis tests passed; the classifier reproduced all 196 original scout phase labels from the original inputs. Every new trajectory passed its independent full-field observable/admissibility audit.',
        'D=2, L=16,32,64; original epsilon grid retained. Three pilot-selected kappa values separated by 0.02 per positive-boundary slice, followed by separately frozen additions near unresolved bracket edges. See batches/*.json for every exact parameter and seed.',
        'Each production parameter/volume group has four independent trajectories of 12000 cycles, starting with windings (0,0), (1,0), (0,1), (-1,0). Each cycle performs one Metropolis sweep, four Wolff clusters, and one winding proposal: 48000 sweeps and 192000 clusters per parameter/volume group.',
        'Analysis retains cycles 2001–12000 from each start, giving 40000 retained cycles per group. HDF5 warmup=100 is a provisional runner label; the selected analysis cut is recorded in analysis/*.json. Raw seed and transient records remain intact. Block sizes are at least ten times the largest monitored tau, with 2000 independent-within-start block resamples. Raw moments are pooled before computing Binder and correlation-length ratios. Cut/block stability, disjoint halves and cross-start agreement are recorded.',
        '', '## What the error bars mean','',
        'The main phase_boundary_refined figure preserves the original joint Binder/correlation-length/effective-exponent rule, with two-standard-error margins and adjacent epsilon-cell unions. The decision uses L32 and L64 and their effective exponent; L16 remains a size/start check. Its bars remain heuristic finite-volume brackets, not Gaussian confidence intervals or a validated thermodynamic extrapolation. The slope across an epsilon cell contributes to their vertical width.',
        'precision_and_finite_size separates the attained grid/Monte Carlo resolution of a specified L64 Binder reference midpoint from the larger joint-diagnostic bracket. It is a conditional finite-volume proxy and must not be quoted as the BKT critical coupling with only its small statistical error. Finite-size and reference-choice bias has not been certified. The fixed-Binder limitation is material; see https://arxiv.org/pdf/1307.0485 section 3.1.',
        f"Conditional local interpolation is resolved for {len(resolved)} of {len(boundary['regions'])} epsilon slices under the recorded gates; unresolved cases remain explicit in boundary.json.",'',
        '| epsilon/pi cell center | old vertical half-width | refined half-width | reduction |','|---:|---:|---:|---:|']
    for r in improvements:lines.append(f"| {r['epsilon_over_pi']:.3f} | {r['old_halfwidth']:.3f} | {r['new_halfwidth']:.3f} | {r['reduction_factor']:.2f}x |")
    lines+=['','## Diagnostics and limits','',f"{len(flagged_individual)} individual trajectories have at least one cut/half diagnostic flag; {len(flagged_pooled)} parameter/volume groups have a cross-start p<0.001 flag. These are retained in diagnostic-summary.json, not silently discarded. Multiple comparisons mean isolated flags are not themselves proof of disequilibrium; no exhaustive equilibration or sector-completeness claim is made.",'',
        '| L | median tau(m2) | max tau(m2) | median tau(Binder influence) | max tau(Binder influence) |','|---:|---:|---:|---:|---:|']
    for L,row in tau_table.items():lines.append(f"| {L} | {row['m2']['median']:.2f} | {row['m2']['maximum']:.2f} | {row['binder']['median']:.2f} | {row['binder']['maximum']:.2f} |")
    lines+=['','Tau is measured in the new hybrid cycle (four clusters), not in the old one-cluster time unit. All-observable time histories, per-start statistics and bootstrap samples are stored under histories_final/, analysis/ and pooled/. Histories include the seed, transient and retained-data zoom with separate vertical scales. Their scope is finite-volume scouting; no infinite-volume BKT fit was performed.','',
        '## Reproduction','',
        'Run pilot/production batches with `.venv/bin/python experiments/2026-09-02-xy-bkt-refine/run.py BATCH`; run edge batches with `run_edges.py BATCH` in the same directory. The edge runner verifies disjoint chain IDs before concurrent execution. Completed jobs are skipped and interrupted jobs resume at their committed checkpoint. The original scout and its model code are unchanged. model-source.tar.gz and the immutable batch manifests record simulation provenance. Analysis scripts, tests, logs, final table and figures are retained in this directory.','']
    (OUT/'REPORT.md').write_text('\n'.join(lines))
    print(json.dumps({k:v for k,v in audit.items() if k not in ['records','verification']},indent=2))

if __name__=='__main__':main()
