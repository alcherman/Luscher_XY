"""Reconcile every trajectory and analysis row, then write the scouting report."""
from datetime import datetime
import json
import sys
import h5py
import numpy as np
from run import OUT, points, STEPS, atomic_json, sha256, now, _provenance


def main():
    manifest=json.loads((OUT/'manifest.json').read_text())
    state=json.loads((OUT/'status.json').read_text())
    assert state['state']=='complete' and state['completed']==588 and not state['failures']
    assert manifest['source_sha256']==_provenance()['source_sha256']
    rows=json.loads((OUT/'results.json').read_text())
    by_id={r['point']['id']:r for r in rows}
    assert len(rows)==len(by_id)==588 and len({p['seed'] for p in points()})==588
    total_bytes=0;jobs=[]
    for p in points():
        job=json.loads((OUT/'jobs'/f"{p['id']}.json").read_text());r=by_id[p['id']]
        assert job['state']=='complete' and job['audit']['passed']
        assert job['point']==r['point']==p and r['summary']['cut']==2001
        assert r['analysis_source_sha256']==sha256(OUT/'analyze.py')
        assert [s['last'] for s in job['segments']]==list(range(1000,10001,1000))
        path=OUT/'chains'/f"{p['id']}.h5"
        assert path.stat().st_size==job['audit']['bytes']
        with h5py.File(path,'r') as f:
            assert f.attrs['committed_rows']==STEPS+1
            assert json.loads(f.attrs['provenance'])['source_sha256']==manifest['source_sha256']
            a=f['magnetization_squared'][2001:];b=f['magnetization_fourth'][2001:]
            assert len(a)==r['summary']['retained']==8000
            assert abs(1-b.mean()/(2*a.mean()**2)-r['summary']['values']['binder'])<1e-13
            assert abs(a.mean()-r['summary']['values']['m2'])<1e-13
            R2=f['structure_zero'][2001:].mean()/f['structure_min'][2001:].mean()-1
            if R2>=0:
                R=np.sqrt(R2)/(2*p['L']*np.sin(np.pi/p['L']))
                assert abs(R-r['summary']['values']['xi_over_L'])<1e-13
            else:assert r['summary']['values']['xi_over_L'] is None
        total_bytes+=path.stat().st_size;jobs.append(job)
    audit=dict(passed=True,completed_utc=now(),chains=588,generated_cycles=588*STEPS,
        audited_fields=588*(STEPS+1),retained_rows=588*8000,total_hdf5_bytes=total_bytes,
        source_sha256=manifest['source_sha256'],unique_seeds=588,complete_checkpoint_segments=5880,
        final_verification_runtime=_provenance(),python_executable=sys.executable,
        per_field_audits='Each jobs/*.json records its full-field audit and complete HDF5 SHA256.',
        final_checks='All committed row counts, source hashes, parameter/seed identities, checkpoint segments, sizes, Binder and xi/L centers reconciled.')
    atomic_json(OUT/'execution-audit.json',audit)
    phase=json.loads((OUT/'phase-brackets.json').read_text())
    flags=[dict(id=r['point']['id'],flags=r['flags']) for r in rows if r['flags']]
    atomic_json(OUT/'flagged-points.json',flags)
    elapsed=(datetime.fromisoformat(state['finished_utc'])-datetime.fromisoformat(state['started_utc'])).total_seconds()
    exponent=json.loads((OUT/'kappa_zero_exponents.json').read_text())
    point_table=[]
    for e in [500,510,520,530,540,565,580,600,620]:
        cells=[]
        for L in [16,32,64]:
            s=by_id[f'e{e:04d}_k0000_L{L:03d}']['summary']
            cells.append(f"{s['values']['binder']:.6f} ± {s['errors']['binder']:.6f}")
        point_table.append(f"| {e/1000:g} | "+' | '.join(cells)+' |')
    e53=[r for r in exponent if r['epsilon_over_pi']==.53]
    eta_text='; '.join(f"{r['sizes'][0]}→{r['sizes'][1]}: {r['eta_eff']:.4f} ± {r['error']:.4f}" for r in e53)
    v53=[]
    for L in [16,32,64]:
        s=by_id[f'e0530_k0000_L{L:03d}']['summary']
        v53.append(f"L={L}: {s['vortex_present_cycles']} retained configurations containing vortices")
    bracket_table=[]
    for r in phase['regions']:
        bracket_table.append(f"| {r['x_lower']:g}–{r['x_upper']:g} | {r['kappa_lower']:g}–{r['kappa_upper']:g} | {r['kind']} |")
    blocks=[r['summary']['blocks'] for r in rows]
    tau_rows=[]
    for L in [16,32,64]:
        group=[r for r in rows if r['point']['L']==L]
        cells=[]
        for key in ['m2','binder','xi_ratio','link_cosine']:
            t=[r['summary']['taus'][key] for r in group if r['summary']['taus'][key] is not None]
            cells.append(f"{np.median(t):.2f} / {max(t):.2f}")
        tau_rows.append(f"| {L} | "+' | '.join(cells)+' |')
    report=f'''# Two-dimensional XY scouting scan

All **588 trajectories completed**, covering the user's exact 28-value epsilon grid, seven kappa values and L=16,32,64. Each trajectory contains 10,000 hybrid cycles and its cold initial configuration. All **{audit['audited_fields']:,} saved fields** passed independent physical-observable recomputation. Execution took {elapsed/60:.1f} minutes with four concurrent single-thread Numba jobs, including per-chain audits. Full HDF5 retention occupies {total_bytes/1e9:.2f} GB; previous experiments were untouched.

## Boundary plot

![Scouting boundary regions](phase_boundary_scout.png)

The black points are representative centers of boundary regions; the bars span their horizontal and vertical brackets. They are **finite-volume scouting brackets, not one-sigma confidence intervals or a fitted infinite-volume BKT curve**. Gray points remain unresolved. The zero-coupling endpoint is represented by its horizontal bracket, without drawing a zero-kappa transition segment through the quasi-ordered region to its left.

For this initial map, `phase.py` uses the envelope of measured critical Binder and xi/L values in the ordinary and constraint-only XY models on L=32,64, with two-standard-error margins and a 32→64 effective-exponent check. L=16 supplies an additional displayed size comparison. These finite-volume observable references come from [Bietenholz et al., Tables 2 and 3](https://arxiv.org/pdf/1307.0485). **Assuming their envelope usefully brackets the mixed action is a heuristic**, not an established bound. No literature parameter endpoint is forced into the result. The raw finite-size data are the primary deliverable; precision reference matching and larger-volume checks remain for refinement.

At each epsilon column, the bracket lies between the largest kappa labelled massive-like and the smallest labelled quasi-ordered-like. An interior point represents an adjacent pair of epsilon columns: its vertical bar encloses both columns' brackets. This assumes a single continuous, locally monotone boundary within that cell. Sampling diagnostic flags prevent a point from establishing a phase label. Broad bars include unresolved grid rows; half the nominal spacing is not assigned automatically. These brackets do not quantify all finite-size bias.

| Epsilon/pi interval | Kappa interval | Region |
|---|---|---|
{chr(10).join(bracket_table)}

Machine-readable results: [phase-brackets.csv](phase-brackets.csv) and [phase-brackets.json](phase-brackets.json). Figure exports: [PDF](phase_boundary_scout.pdf), [SVG](phase_boundary_scout.svg).

## The region just above epsilon/pi = 1/2

![Magnification at zero coupling](near_threshold.png)

At kappa=0, the new 0.51–0.54 points show nearly size-independent nonzero Binder values and effective magnetic exponents below 1/4. This supports a quasi-ordered region extending beyond the vortex-admissibility threshold. It is finite-size evidence, and does not by itself rule out a much larger finite correlation length or measure a thermodynamic endpoint precisely.

The L=64 trajectory at epsilon/pi=0.54 has a flagged difference between the two retained-half means of squared magnetization, so it does not establish the lower endpoint bracket. The final conservative bracket is **0.53–0.60**. Its displayed midpoint, 0.565, is a grid representative and must not be read as a new precision determination agreeing with the literature endpoint.

| Epsilon/pi | U, L=16 | U, L=32 | U, L=64 |
|---|---:|---:|---:|
{chr(10).join(point_table)}

Errors in this table are conditional one-standard-deviation joint-block bootstrap errors. At epsilon/pi=0.53, the effective exponents are **{eta_text}**. Vortex-containing configurations were actually sampled there: {'; '.join(v53)}. These are counts of configurations, not independent vortex-creation events. Thus the quasi-order-like signal at this point is not solely an artifact of observing no vortices. At 0.51 no vortices were observed in these short trajectories; that does not imply a prohibition or supply a reliable upper bound. At exactly 0.50, zero plaquette vorticity is imposed by the strict constraint and was checked on every saved field.

The [full zero-kappa comparison](kappa_zero_focus.png) also shows xi/L, vortex density and the crossover towards rapidly decaying finite-size magnetization. The [Binder scan](binder_scan.png) and [correlation-length scan](xi_over_L_scan.png) display all kappa rows and all three volumes.

The [Binder heatmaps](binder_heatmap.png) and [xi/L heatmaps](xi_over_L_heatmap.png) give complementary views of the whole parameter plane.

## Sampling and analysis scope

The exact schedule is one complete Metropolis sweep, one Wolff cluster and one winding proposal per cycle, with fixed proposal half-width pi. All starts are cold; there is one seed per parameter/volume pair. No refinement runs or extra replicas were added. The model uses S = -kappa times the sum of positive-direction link cosines, strict absolute principal link angle less than epsilon, and periodic square lattices.

The primary retained interval is **cycles 2,001–10,000**, giving 8,000 measurements per trajectory and {audit['retained_rows']:,} total retained rows. Initial and full histories were generated; representative cold relaxations and full histories were visually inspected, with quantitative cut and half comparisons across the dataset. This does not certify every trajectory's equilibrium or initialization independence. The original chain warmup attribute is only a provisional label; `analysis-provenance.json` records the actual cut used here.

Tau_int uses the convention 1/2 plus the positive-lag autocorrelation sum and a nonperiodic FFT initial-positive monotone paired estimator. Block widths are at least 20 cycles and approximately ten times the largest estimated tau among squared magnetization, fourth moment, Binder influence, correlation-length-ratio influence, mean link cosine and vortex density. The cosine is monitored because action density is identically zero at kappa=0. Winding diagnostics are saved separately. There are {min(blocks)}–{max(blocks)} blocks per primary estimate.

The analysis uses 2,000 joint balanced contiguous-block resamples, retaining every post-cut row, and recalculates nonlinear ratios each time. Negative estimated structure-factor ratios are left unresolved, rather than clipped to a zero correlation length. The ratio is formed from ensemble-mean structure factors averaged over the two directions. Changing the cut to 1,000 or 4,000 discarded cycles, doubling block widths, and comparing disjoint retained halves provides sensitivity checks. Half differences beyond three combined standard errors are flagged; overlapping cut differences use a conservative two-combined-error diagnostic. These flags are diagnostics across many comparisons, not individually calibrated discoveries of nonstationarity.

**{len(flags)} trajectories have diagnostic flags**, listed in [flagged-points.json](flagged-points.json). Absence of a flag does not establish equilibration; independent starts and replicas remain necessary for a precision result. These tau estimates cannot exclude slower processes that the chains did not sample.

Median / maximum estimated tau_int, in hybrid cycles:

| L | Squared magnetization | Binder influence | xi ratio influence | Link cosine |
|---|---:|---:|---:|---:|
{chr(10).join(tau_rows)}

## Verification and reproducibility

The focused sampler/writer suite passed **155 tests**, including buffered versus single-row trajectory equality, exact resume, and replay after injected partial-batch write failures. **Three additional analysis tests passed**, checking the nonlinear bootstrap, tail inclusion and invalid correlation-length handling. These logs are [implementation-tests.txt](implementation-tests.txt) and [analysis-tests.txt](analysis-tests.txt). The full legacy suite was not run.

Each job independently recomputed magnetization moments, Fourier modes, cosine/action density, plaquette vortices, every axial winding and maximum link angle from every saved field. It also checked strict admissibility, canonical finite angles, RNG key sets, work counters and nonnegative timing fields. Full-file SHA256 hashes are in the per-job records. The final [execution audit](execution-audit.json) reconciles all 588 source hashes, seeds, file sizes and row counts, all 5,880 progress segments, and all retained Binder and correlation-length centers. These numerical consistency checks do not establish a thermodynamic phase transition.

All raw fields, observables, RNG states, work and timing data, including discarded equilibration, remain in `chains/`; unfiltered logs remain in `logs/`. [manifest.json](manifest.json), [status.json](status.json), [results.json](results.json), [observables.csv](observables.csv), [analysis-provenance.json](analysis-provenance.json), `analysis/*_bootstrap.npz`, `histories/`, and the simulation source archive preserve the evidence. See [README.md](README.md) for reproduction commands.
'''
    (OUT/'REPORT.md').write_text(report)
    print(json.dumps(audit,indent=2))


if __name__=='__main__':main()
