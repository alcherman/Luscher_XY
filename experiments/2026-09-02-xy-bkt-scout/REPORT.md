# Two-dimensional XY scouting scan

All **588 trajectories completed**, covering the user's exact 28-value epsilon grid, seven kappa values and L=16,32,64. Each trajectory contains 10,000 hybrid cycles and its cold initial configuration. All **5,880,588 saved fields** passed independent physical-observable recomputation. Execution took 27.4 minutes with four concurrent single-thread Numba jobs, including per-chain audits. Full HDF5 retention occupies 93.79 GB; previous experiments were untouched.

## Boundary plot

![Scouting boundary regions](phase_boundary_scout.png)

The black points are representative centers of boundary regions; the bars span their horizontal and vertical brackets. They are **finite-volume scouting brackets, not one-sigma confidence intervals or a fitted infinite-volume BKT curve**. Gray points remain unresolved. The zero-coupling endpoint is represented by its horizontal bracket, without drawing a zero-kappa transition segment through the quasi-ordered region to its left.

For this initial map, `phase.py` uses the envelope of measured critical Binder and xi/L values in the ordinary and constraint-only XY models on L=32,64, with two-standard-error margins and a 32→64 effective-exponent check. L=16 supplies an additional displayed size comparison. These finite-volume observable references come from [Bietenholz et al., Tables 2 and 3](https://arxiv.org/pdf/1307.0485). **Assuming their envelope usefully brackets the mixed action is a heuristic**, not an established bound. No literature parameter endpoint is forced into the result. The raw finite-size data are the primary deliverable; precision reference matching and larger-volume checks remain for refinement.

At each epsilon column, the bracket lies between the largest kappa labelled massive-like and the smallest labelled quasi-ordered-like. An interior point represents an adjacent pair of epsilon columns: its vertical bar encloses both columns' brackets. This assumes a single continuous, locally monotone boundary within that cell. Sampling diagnostic flags prevent a point from establishing a phase label. Broad bars include unresolved grid rows; half the nominal spacing is not assigned automatically. These brackets do not quantify all finite-size bias.

| Epsilon/pi interval | Kappa interval | Region |
|---|---|---|
| 0.53–0.6 | 0–0 | zero_kappa_endpoint |
| 0.6–0.62 | 0–0.6 | interior_grid_region |
| 0.62–0.64 | 0–0.6 | interior_grid_region |
| 0.64–0.66 | 0–1 | interior_grid_region |
| 0.66–0.68 | 0.4–1 | interior_grid_region |
| 0.68–0.7 | 0.4–0.8 | interior_grid_region |
| 0.7–0.72 | 0.6–0.8 | interior_grid_region |
| 0.72–0.74 | 0.6–1.2 | interior_grid_region |
| 0.74–0.76 | 0.6–1.2 | interior_grid_region |
| 0.76–0.78 | 0.6–1 | interior_grid_region |
| 0.78–0.8 | 0.8–1 | interior_grid_region |
| 0.8–0.82 | 0.8–1 | interior_grid_region |
| 0.82–0.84 | 0.8–1.2 | interior_grid_region |
| 0.84–0.86 | 0.8–1.2 | interior_grid_region |
| 0.86–0.88 | 0.8–1.2 | interior_grid_region |
| 0.88–0.9 | 0.8–1.2 | interior_grid_region |
| 0.9–0.92 | 0.8–1.2 | interior_grid_region |
| 0.92–0.94 | 1–1.2 | interior_grid_region |
| 0.94–0.96 | 0.8–1.2 | interior_grid_region |
| 0.96–0.98 | 0.8–1.2 | interior_grid_region |
| 0.98–1 | 1–1.2 | interior_grid_region |
| 1–1 | 1–1.2 | unconstrained_endpoint |

Machine-readable results: [phase-brackets.csv](phase-brackets.csv) and [phase-brackets.json](phase-brackets.json). Figure exports: [PDF](phase_boundary_scout.pdf), [SVG](phase_boundary_scout.svg).

## The region just above epsilon/pi = 1/2

![Magnification at zero coupling](near_threshold.png)

At kappa=0, the new 0.51–0.54 points show nearly size-independent nonzero Binder values and effective magnetic exponents below 1/4. This supports a quasi-ordered region extending beyond the vortex-admissibility threshold. It is finite-size evidence, and does not by itself rule out a much larger finite correlation length or measure a thermodynamic endpoint precisely.

The L=64 trajectory at epsilon/pi=0.54 has a flagged difference between the two retained-half means of squared magnetization, so it does not establish the lower endpoint bracket. The final conservative bracket is **0.53–0.60**. Its displayed midpoint, 0.565, is a grid representative and must not be read as a new precision determination agreeing with the literature endpoint.

| Epsilon/pi | U, L=16 | U, L=32 | U, L=64 |
|---|---:|---:|---:|
| 0.5 | 0.494463 ± 0.000191 | 0.494587 ± 0.000236 | 0.494501 ± 0.000213 |
| 0.51 | 0.493980 ± 0.000231 | 0.493985 ± 0.000219 | 0.493840 ± 0.000229 |
| 0.52 | 0.493600 ± 0.000229 | 0.493572 ± 0.000220 | 0.493074 ± 0.000272 |
| 0.53 | 0.493315 ± 0.000225 | 0.493029 ± 0.000272 | 0.492718 ± 0.000281 |
| 0.54 | 0.492574 ± 0.000250 | 0.492053 ± 0.000294 | 0.492445 ± 0.000292 |
| 0.565 | 0.491169 ± 0.000286 | 0.490784 ± 0.000319 | 0.490699 ± 0.000410 |
| 0.58 | 0.489306 ± 0.000331 | 0.489332 ± 0.000385 | 0.488967 ± 0.000469 |
| 0.6 | 0.485401 ± 0.000648 | 0.485331 ± 0.000570 | 0.481132 ± 0.000917 |
| 0.62 | 0.477631 ± 0.001011 | 0.465764 ± 0.001798 | 0.433016 ± 0.003294 |

Errors in this table are conditional one-standard-deviation joint-block bootstrap errors. At epsilon/pi=0.53, the effective exponents are **16→32: 0.2239 ± 0.0062; 32→64: 0.2229 ± 0.0070**. Vortex-containing configurations were actually sampled there: L=16: 2 retained configurations containing vortices; L=32: 10 retained configurations containing vortices; L=64: 37 retained configurations containing vortices. These are counts of configurations, not independent vortex-creation events. Thus the quasi-order-like signal at this point is not solely an artifact of observing no vortices. At 0.51 no vortices were observed in these short trajectories; that does not imply a prohibition or supply a reliable upper bound. At exactly 0.50, zero plaquette vorticity is imposed by the strict constraint and was checked on every saved field.

The [full zero-kappa comparison](kappa_zero_focus.png) also shows xi/L, vortex density and the crossover towards rapidly decaying finite-size magnetization. The [Binder scan](binder_scan.png) and [correlation-length scan](xi_over_L_scan.png) display all kappa rows and all three volumes.

The [Binder heatmaps](binder_heatmap.png) and [xi/L heatmaps](xi_over_L_heatmap.png) give complementary views of the whole parameter plane.

## Sampling and analysis scope

The exact schedule is one complete Metropolis sweep, one Wolff cluster and one winding proposal per cycle, with fixed proposal half-width pi. All starts are cold; there is one seed per parameter/volume pair. No refinement runs or extra replicas were added. The model uses S = -kappa times the sum of positive-direction link cosines, strict absolute principal link angle less than epsilon, and periodic square lattices.

The primary retained interval is **cycles 2,001–10,000**, giving 8,000 measurements per trajectory and 4,704,000 total retained rows. Initial and full histories were generated; representative cold relaxations and full histories were visually inspected, with quantitative cut and half comparisons across the dataset. This does not certify every trajectory's equilibrium or initialization independence. The original chain warmup attribute is only a provisional label; `analysis-provenance.json` records the actual cut used here.

Tau_int uses the convention 1/2 plus the positive-lag autocorrelation sum and a nonperiodic FFT initial-positive monotone paired estimator. Block widths are at least 20 cycles and approximately ten times the largest estimated tau among squared magnetization, fourth moment, Binder influence, correlation-length-ratio influence, mean link cosine and vortex density. The cosine is monitored because action density is identically zero at kappa=0. Winding diagnostics are saved separately. There are 29–400 blocks per primary estimate.

The analysis uses 2,000 joint balanced contiguous-block resamples, retaining every post-cut row, and recalculates nonlinear ratios each time. Negative estimated structure-factor ratios are left unresolved, rather than clipped to a zero correlation length. The ratio is formed from ensemble-mean structure factors averaged over the two directions. Changing the cut to 1,000 or 4,000 discarded cycles, doubling block widths, and comparing disjoint retained halves provides sensitivity checks. Half differences beyond three combined standard errors are flagged; overlapping cut differences use a conservative two-combined-error diagnostic. These flags are diagnostics across many comparisons, not individually calibrated discoveries of nonstationarity.

**14 trajectories have diagnostic flags**, listed in [flagged-points.json](flagged-points.json). Absence of a flag does not establish equilibration; independent starts and replicas remain necessary for a precision result. These tau estimates cannot exclude slower processes that the chains did not sample.

Median / maximum estimated tau_int, in hybrid cycles:

| L | Squared magnetization | Binder influence | xi ratio influence | Link cosine |
|---|---:|---:|---:|---:|
| 16 | 3.22 / 6.30 | 2.11 / 5.29 | 3.33 / 5.48 | 2.10 / 4.80 |
| 32 | 4.07 / 13.86 | 2.62 / 6.93 | 4.30 / 15.75 | 2.21 / 6.84 |
| 64 | 4.68 / 24.51 | 3.09 / 10.91 | 5.02 / 27.22 | 2.33 / 10.83 |

## Verification and reproducibility

The focused sampler/writer suite passed **155 tests**, including buffered versus single-row trajectory equality, exact resume, and replay after injected partial-batch write failures. **Three additional analysis tests passed**, checking the nonlinear bootstrap, tail inclusion and invalid correlation-length handling. These logs are [implementation-tests.txt](implementation-tests.txt) and [analysis-tests.txt](analysis-tests.txt). The full legacy suite was not run.

Each job independently recomputed magnetization moments, Fourier modes, cosine/action density, plaquette vortices, every axial winding and maximum link angle from every saved field. It also checked strict admissibility, canonical finite angles, RNG key sets, work counters and nonnegative timing fields. Full-file SHA256 hashes are in the per-job records. The final [execution audit](execution-audit.json) reconciles all 588 source hashes, seeds, file sizes and row counts, all 5,880 progress segments, and all retained Binder and correlation-length centers. These numerical consistency checks do not establish a thermodynamic phase transition.

All raw fields, observables, RNG states, work and timing data, including discarded equilibration, remain in `chains/`; unfiltered logs remain in `logs/`. [manifest.json](manifest.json), [status.json](status.json), [results.json](results.json), [observables.csv](observables.csv), [analysis-provenance.json](analysis-provenance.json), `analysis/*_bootstrap.npz`, `histories/`, and the simulation source archive preserve the evidence. See [README.md](README.md) for reproduction commands.
