# Proposed two-parameter BKT scan

Status: the authorized scouting run completed all 588 trajectories on 2026-09-02 (local time). Results and limitations are in `experiments/2026-09-02-xy-bkt-scout/REPORT.md`. The user expanded the epsilon grid below and requested all three volumes. Refinement and multiple-replica production remain future stages.

## Target and working update

Periodic square lattices, D = 2, L = 16, 32, 64, with

$$
S=-\kappa\sum_{x,\mu}\cos\theta_\mu(x),\qquad
|\theta_\mu(x)|<\epsilon,\qquad x_\epsilon=\epsilon/\pi.
$$

The requested plot has horizontal coordinate x_epsilon and vertical coordinate kappa. The domain is 0.5 <= x_epsilon <= 1 and 0 <= kappa <= 1.2.

Use hybrid updates as the campaign default: one complete Metropolis sweep, one Wolff cluster, and one winding proposal per fixed cycle. Specify this explicitly in each manifest rather than relying on the existing generic CLI default. Freeze the Metropolis proposal width during production; if tight constraints make the existing width inefficient, tune in a separate discarded pilot and record the selected value. Do not stop a cycle adaptively when accumulated cluster size first exceeds a volume. Parallelize independent parameter/volume/seed jobs, with a bounded worker count and no nested CPU oversubscription.

## Known endpoints and what they do not imply

For the ordinary square-lattice XY model, the reported critical inverse temperature in our normalization is kappa_c = 1.11996(6). [Komura and Okabe, arXiv:1210.6116](https://arxiv.org/abs/1210.6116).

For the pure angle-constraint action, the reported critical angle is epsilon_c = 1.77521(57), or x_epsilon,c = 0.56507(18). This is the zero-vortex-penalty row of Table 1, and hence the kappa = 0 model here; the other rows of that paper vary an explicit vortex penalty and are not points on our mixed cosine/constraint curve. [Bietenholz et al., arXiv:1212.0579, Table 1](https://arxiv.org/pdf/1212.0579).

Use these values to place calibration runs and assess bias, not to force our measured curve through prescribed endpoints. The expected quasi-ordered region extends beyond the vortex-forbidden threshold: admitting vortices at x_epsilon > 0.5 does not immediately disorder the system. Within the proposed nonnegative-kappa domain, the transition curve is expected to begin near x_epsilon = 0.565 and rise towards the ordinary-XY endpoint. Do not draw a fictitious kappa_c = 0 transition segment to the left of this endpoint.

## Stage 1: coarse reconnaissance

Use the Cartesian grid

```text
epsilon/pi = 0.50, 0.51, 0.52, 0.53, 0.54, 0.565, 0.58,
             0.60, 0.62, 0.64, ..., 0.98, 1.00
kappa      = 0.00, 0.20, 0.40, 0.60, 0.80, 1.00, 1.20
```

This is 28 epsilon values and 196 parameter pairs. Run L = 16, 32, 64 at every pair, one independent 10,000-cycle pilot chain per pair and volume: 588 trajectories and 5,880,000 cycles. Use exactly the user's grid, including its omission of 0.55. These short single-seed runs are reconnaissance, not final estimates. Prioritize the full kappa = 0 axis, then fill the remaining rows. Existing ordinary-XY data and literature endpoints provide context, but are not substituted for campaign data or used to force a measured boundary through preset locations.

Review the coarse magnetic finite-size trends and identify broad candidate transition intervals. A massive system with correlation length much greater than 64 can look quasi-ordered on these volumes, so unresolved intervals remain unresolved. The no-vortex line x_epsilon = 0.5 is a useful implementation control; the added points 0.51-0.54 directly probe the region above this threshold. The executable manifest, full trajectories and logs are in `experiments/2026-09-02-xy-bkt-scout/`. Four independent single-thread CPU jobs run concurrently; all physical observables are independently recomputed from every saved field before a job is marked complete.

## Stage 2: locate and refine the boundary

At fixed x_epsilon, refine each candidate kappa interval from spacing 0.20 to 0.05, then to 0.01-0.02 where justified by the statistical and finite-size resolution. Near kappa = 0, scan horizontally in x_epsilon and refine to spacings 0.005 and, if justified, 0.002-0.0025. Insert extra x_epsilon slices where the apparent curve bends; aim initially for about 10-12 transition locations rather than uniform precision throughout the plane.

Every refined neighborhood uses all three volumes. Generate four independent chains per volume and production parameter point, starting with 10,000 cycles each and extending checkpoints as needed. Pilot chains may contribute only if their production kernel was already fixed, their metadata are compatible, and they pass the same analysis checks. Keep separate labels for tuning, reconnaissance, and production.

Use at least a cold start and an independently prepared admissible start approached from the opposite side of the candidate transition. Unconstrained random spins are invalid at tight epsilon. All prepared starts require equilibration at the target parameters. Include admissible nonzero-winding starts in representative mixing checks; different RNG seeds with identical cold starts alone do not test sector dependence.

After inspecting complete histories including the initial field, select a retained plateau and estimate autocorrelations. Monitor squared magnetization, fourth moment, the Binder and correlation-length influence series, mean link cosine, vortex density, and winding histories. At kappa = 0, the action is identically zero and cannot diagnose mixing; use the fluctuating link cosine. At x_epsilon = 0.5, zero vorticity is a constraint rather than an equilibration test.

Use blocks at least roughly ten times the largest reliably estimated autocorrelation time among the nonconstant monitored observables. Check stability against longer blocks and later cuts; an unresolved long autocorrelation tail calls for more data rather than a nominal block width. As an initial statistical target, request roughly 2,000 pooled effective samples per volume/parameter point for each of the two primary influence series, with at least 50 sufficiently long blocks in each chain. These are minimum diagnostics rather than guarantees of a useful phase-location error. Extend until parameter-location errors meet the chosen resolution or finite-size uncertainty dominates. If independent starts disagree, diagnose equilibration before pooling them.

## Transition estimator

Measure jointly, from ensemble means,

$$
U_L=1-\frac{\langle |m|^4\rangle}{2\langle |m|^2\rangle^2},\qquad
R_L=\xi_{2\mathrm{nd}}/L,
$$

$$
R_L=\frac{1}{2L\sin(\pi/L)}
\sqrt{\frac{\langle F_0\rangle}{\langle (F_{p_x}+F_{p_y})/2\rangle}-1}.
$$

The current chain format already stores all inputs. Average the two directions at the level of structure factors and recompute both nonlinear observables within every block-bootstrap draw. Also inspect each direction separately as a symmetry check. Susceptibility scaling, including its BKT logarithmic correction, and vortex/winding diagnostics are supporting evidence.

Use joint matching of (U_L, R_L) to a well-determined critical finite-size reference, with a fitted nonuniversal length rescaling, as the initial quantitative method. Carry reference uncertainty and its valid size range through the fit. Do not extrapolate a published reference parameterization outside its tested range without validation. Do not fit several unconstrained logarithmic correction coefficients to only three volumes. Hasenbusch develops the joint-matching approach and explicitly shows the need for subleading Binder corrections. [arXiv:0804.1880, Sections 6-7](https://arxiv.org/pdf/0804.1880).

The constrained endpoint is a particularly important check. Bietenholz et al. found that matching their accessible-size Binder and correlation-length results to asymptotic spin-wave values suggested epsilon = 1.756(2), inconsistent with their preferred 1.775(1) determination. Thus even a small statistical error can accompany a significant finite-size bias. Their result is a warning against a fixed universal Binder threshold or uncorrected spin-wave matching, not a reason to abandon these observables. [arXiv:1307.0485, Section 3.1](https://arxiv.org/pdf/1307.0485).

Validate the estimator at both known endpoints, compare fits including and excluding L = 16, and test reasonable reference/correction choices. The remaining two volumes cannot support a free extrapolation; use the reduced fit as a sensitivity check. If the inferred location shifts appreciably or the reference matching fails, add L = 128 at a few nearby points. If that is insufficient, keep the location provisional rather than declaring a thermodynamic transition from L <= 64. A second route is a massive-side correlation-length fit, but only once convergence in L at fixed parameters is demonstrated; otherwise it is not an infinite-volume correlation length.

The usual cosine-only helicity estimator cannot be used with an active hard constraint: a twist also changes the admissible domain. At kappa = 0 that naive estimator would incorrectly give zero stiffness. A dynamical-twist/free-energy-curvature implementation would be a separate validation project if needed. It is not a prerequisite for the coarse scan, and is not currently supplied by the runner. The constraint-action paper above describes dynamical boundary conditions for this purpose.

## Horizontal and vertical error bars

The user clarified that the requested bars represent localization of the transition within the two-dimensional grid. Although each simulation's parameters are exact, finite grid resolution leaves an uncertain boundary region. Use bars spanning a supported rectangle in (x_epsilon, kappa); a representative midpoint is a plotting convention. Half a grid spacing is appropriate only when the neighboring sampled points actually bracket the boundary. Wider unresolved regions need wider bars. No diagonal-path construction is needed.

Label boundary-location brackets explicitly; they are not automatically one-sigma statistical errors. Propagate sampling uncertainty and include finite-size/model sensitivity, either by widening the bracket or showing a separately identified band. Show unresolved portions as provisional, and do not interpolate a smooth curve through a failed mixing or finite-size check. Display literature endpoints distinctly from new simulation estimates. Grid refinement can follow after this scouting run.

## Organization and retention

Create an immutable campaign manifest specifying exact float parameters, lattice geometry, update schedule, seeds, initialization provenance, proposal width, source hashes and diagnostic thresholds. Use separate directories for calibration, coarse runs and refinement rounds; a scheduler state file should atomically checkpoint each completed chunk and its timing. Parallelize complete independent chains. Save unfiltered run logs, every generated field, observables, RNG states, update counts, cluster sizes and timings, including the initial field and discarded equilibration history.

Keep analysis scripts separate from simulation. Repeated plotting, cut changes, blocking and fits should operate on existing files. Persist the phase-location table with both coordinates, statistical intervals, covariance, finite-size sensitivity, volumes, retained sample counts, effective sample sizes and status. Produce PDF/SVG/PNG figures and raw bootstrap locations. The L48 epsilon = pi timing comparison does not establish throughput or mixing across this constrained scan; measure those in the pilot before quoting a completion time.

## September 2 refinement run

The user requested smaller vertical resolution after the completed coarse scout. The refinement campaign is in `experiments/2026-09-02-xy-bkt-refine/`.

It keeps the requested epsilon grid and all three volumes. Pilot L64 runs choose three directly simulated kappa values at spacing 0.02 near each candidate boundary. Production uses four independent seeds and prepared windings (0,0), (1,0), (0,1), (-1,0), with 12000 cycles per chain and four fixed Wolff clusters per hybrid cycle. Each parameter/volume group therefore has 48000 Metropolis sweeps and 192000 clusters before the analysis cut. Missing nearby grid points may be added in a separately frozen batch after reviewing production diagnostics. Additional kappa=0 groups check epsilon/pi=0.54,0.565,0.58.

The narrow precision diagnostic locates a specified finite-volume Binder-reference midpoint. It must be labeled as a conditional proxy, not an infinite-volume BKT determination. Show its sampling/grid precision alongside the unchanged joint U/R/eta scouting bracket; do not relabel the latter as a statistical standard error or discard finite-size ambiguity to obtain shorter bars. All trajectories and all generated rows remain on disk, including pilots and initial transients.

Completed September 3 at 05:53 UTC: 1200 production/refinement trajectories (100 parameter pairs on three volumes, four starts each) and 22 design pilots. Every full-field audit and final record reconciliation passed. On the 20 comparable epsilon cells, the median reduction of the vertical half-width is a factor 3.54. Refined interior half-widths range from 0.04 to 0.12, with median 0.06; the steep low-kappa portion remains wider. The kappa=0 scouting endpoint is bracketed by epsilon/pi=0.54 and 0.58. Exact tables, diagnostics, the unchanged-rule comparison, and the main figure are in the refinement campaign's REPORT.md and phase-brackets.csv. This remains finite-volume scouting, with no thermodynamic extrapolation.

## September 3 low-kappa refinement

The subsequent request to refine both grid directions near small kappa is recorded in `experiments/2026-09-03-xy-bkt-low-kappa/`. Two completed batches add epsilon/pi=0.575 and 0.590 and tighten the 0.580 slice: seven parameter pairs, three volumes, four independent starts, 32000 cycles each. All 84 full-field audits and final reconciliation passed. Each parameter/volume group retains 120000 cycles after the supported 2000-cycle cut, three times the previous count. The lowest positive-curve half-widths are now 0.07–0.08, still above the requested 0.04–0.05 target; the target is not marked achieved. The measured fields occupy 42.69 GB. Further work is limited by disk space: about 25.5 GiB remains free against the pre-existing 25 GiB runner reserve. The concrete next four-point proposal needs about 23.7 GiB more free space; its manifest is a proposal, not a scheduled batch. Raw fields, observables, RNG states and all transients remain intact.

## September 3 kappa-zero endpoint refinement

The endpoint-specific request is completed in `experiments/2026-09-03-xy-bkt-kappa-zero/`. New points at epsilon/pi=0.550,0.570,0.575, all at kappa=0 and L=16,32,64, each use two independent 7000-cycle hybrid starts. The 2000-cycle cut leaves 10000 pooled cycles per parameter/volume group. The unchanged joint U/R/eta rule labels 0.550 quasi-ordered-like, 0.570 massive-like, and 0.575 unresolved. Thus the horizontal scouting bracket narrows from [0.54,0.58] to [0.55,0.57], halving the plotted half-width from 0.02 to 0.01. This is a finite-volume grid/diagnostic bracket, without thermodynamic extrapolation.

All 18 field audits and final reconciliation passed; 126,018 fields (2.006 GB) are preserved with measurements, initial/transient rows, RNG, update work and timings. No per-chain or pooled consistency flags were found. Five analysis tests passed and 292 prior classifications were reproduced. The updated full curve and endpoint diagnostics are in phase_boundary_refined.png and endpoint_diagnostics.png; REPORT.md and phase-brackets.csv contain the numerical results. This modest batch fits the remaining storage reserve; the separate low-kappa vertical-precision target remains unfinished.
