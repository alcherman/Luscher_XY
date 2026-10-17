# Infinite volume XY phase boundary campaign design

The target is the infinite-volume critical coupling beta_c(epsilon) in dimensions
2, 3, and 4, at fixed epsilon/pi values spaced by 0.02. Only the positive-beta
transition branch is included. Each accepted location should have a vertical
uncertainty half-width at most 0.02 as a target, with wider bars or unresolved
slices allowed within the deadline. The output is one plot per dimension, plus
a combined comparison, with vertical error bars only and a numerical table.

This replaces the earlier finite-volume scouting target. The user authorized
larger boxes where necessary, subject to a practical budget. Results are due
Friday October 16, 2026 at 3 p.m. Central. During weekdays 9 a.m. to 5 p.m.
Central, at most 16 single-thread simulation workers run; outside those hours,
up to 18 may run. Begin with modest dimension-specific sizes; no 64^4 run is
planned. A fit using small boxes is an extrapolation candidate; test its bias
before accepting the requested uncertainty.

Dimension 5 is an optional last-priority extension only if the dimensions
2, 3, and 4 deliverables are ready and enough time remains before the deadline.
It needs its own above-upper-critical-dimension finite-size treatment and
endpoint scouting. Do not assume the d=4 logarithmic ansatz applies. The runner
keeps d=5 disabled until `enable_d5` is switched on after checking the priority
deliverables, resource use, and time forecast.

## Existing measurements and candidate epsilon slices

| Dimension | Candidate positive branch on the requested grid | Existing input |
| --- | --- | --- |
| 2 | 0.58, 0.60, ..., 1.00 | All 22 slices have measurements at N=16,32,64 in the September 2 refinement; the September 3 low-coupling campaign adds longer measurements at 0.58. |
| 3 | 0.80, 0.82, ..., 1.00 | Zero-coupling measurements at N=8,12,16, with cutoff angles 2.40 to 2.60 radians. Positive-beta scans remain necessary. |
| 4 | 0.88, 0.90, ..., 1.00 | Zero-coupling measurements at N=8,12,16 and selected N=20 points. Positive-beta scans remain necessary. |

The candidate lists follow the existing scouting endpoints and are not final
infinite-volume branch classifications. Check neighboring requested slices
0.56 in d=2 and 0.86 in d=4, and the first candidate 0.80 in d=3, before fixing
the branch membership. A zero-coupling endpoint is useful for locating the
branch, but is not inserted as a positive-beta transition or plotted with an
artificially exact horizontal coordinate.

The existing d=2 finite-volume Binder-reference proxy at epsilon/pi=0.60 is
0.208678 at N=16, 0.249938 at N=32, and 0.272344 at N=64. The last size change
is 0.022405. These proxy values establish observable finite-size sensitivity;
they do not determine the infinite-volume critical coupling or its error.
See `experiments/2026-09-02-xy-bkt-refine/boundary.json` and `boundary.py`.

## Extrapolation methods

### October 14 user authorized adaptive configuration retention

The user explicitly authorized this policy for existing and future campaign
runs. It supersedes the earlier requirement to retain every lattice field and
the prohibition on new thinning or compaction. Keep every observable, winding,
RNG, timing and provenance record, along with full logs and earlier analyses.

Generate and independently audit the complete trajectory before conversion.
The detached `retain_fields.py` service checks the original file hash against
that full-field audit, writes a temporary archive, and verifies every retained
field and every other dataset byte for byte before replacing the source.
Durable prepared records, atomic replacement and recovery events preserve the
audit lineage. Temporary copies and outstanding simulations share the 25 GiB
disk reserve. Retention consumes one worker slot during conversion.

Retain the initial/final fields and every published 100-cycle checkpoint.
Estimate autocorrelation from complete post-warmup and late-half histories of
the relevant moments, influence series, energy, vortices and winding. Keep four
dense diagnostic stretches, each at least 64 cycles or 50 times the largest
measured autocorrelation time, and retain other fields at a stride of ten times
that estimate. This reduces measured dependence without certifying independence.

The field dataset retains its logical shape for native checkpoint replay.
Unstored fields are explicit NaN fill; exact saved indices are recorded in its
attributes. Consumers must consult these indices. Production archives are
immutable: resume a checkpoint from a working copy. Exact checkpoint replay,
observable preservation and interrupted-replacement recovery are tested.
Original closure hashes describe the original full files; current archive
hashes and physical byte counts are recorded separately in retention records.

### October 14 zero beta endpoint measurements

The user requested a beta=0 transition measurement for d=2,3,4. At fixed
beta=0 scan epsilon, estimate the infinite-volume epsilon_c/pi, and plot
that point at beta=0 exactly with horizontal Monte Carlo and tested
sensitivity bars only. The other fixed-epsilon transition measurements
retain vertical bars only. The prior beta half-width target does not specify
an epsilon precision; use practical local refinement and report attained
endpoint uncertainty honestly.

Use eta_eff and xi/L in 2D and the fixed Binder quarter definition in 3D/4D.
Increasing epsilon weakens ordering, reversing the response relative to
increasing beta. Dimension-specific finite-size coordinates can be applied
to epsilon/pi as the thermal parameter, conditional on a smooth local
mapping; preserve interpolation, minimum-size, correction, cut/block and
start sensitivity. Reuse the September zero-coupling archives at their
actual epsilon values, including the radian-grid 3D values without rounding.
New endpoint points have campaign_role=zero_beta_endpoint and an independent
seed namespace. They must not enter fixed-epsilon positive-beta groupings.
Preserve earlier plots and create updated plot/table/report versions.

### October 14 user directed 2D method update

The user resumed 2D work and selected the effective magnetic exponent
eta_eff and the second moment correlation length ratio R=xi_2nd/N.
This supersedes Binder based primary estimators for 2D only. The completed
3D/4D Binder quarter deliverables remain intact. Preserve the 316 unstarted
2D Binder quarter jobs in their deferred state; do not launch that old queue.

Start with compatible September data at the exact requested epsilon grid.
Compute eta_eff(N,2N)=-log2(<m^2>_(2N)/<m^2>_N), retaining covariance with
R and between overlapping size pairs through the joint moment bootstrap.
Use the limiting BKT targets eta=1/4 and R=0.7506912 (square periodic boxes;
Bietenholz et al., arXiv:1307.0485, equations 1.4 and 3.5) as local crossing
definitions, with BKT logarithmic size extrapolation and tested variants.
These are infinite-size targets, not asserted finite-box critical values.
Keep separate observable extrapolations and their disagreement visible;
any shared-limit fit must preserve covariance and test identifiability.
No straight line in 1/N and no Binder criterion replaces those checks.

First calibrate N=128/256 at epsilon/pi=.60,.80,1.00 using locally informed
beta neighborhoods and independent starts. Measure generation plus audit
cost, mixing diagnostics, storage and crossing coverage before committing
to full larger-box coverage. Maintain the Friday deadline, practical .02
uncertainty target, wider/unresolved allowance, full retention, detached
operation, 25GiB reserve, and daytime two-core reservation.

### October 13 user-directed simplification: fixed Binder thresholds

The user requested a simpler primary estimator: define a finite-box location
by a fixed Binder threshold, then extrapolate these locations with size.
This supersedes the requirements below to use joint reference matching or a
global observable-scaling fit as the primary estimator. Those methods remain
optional cross-checks; rejection of their broad beta windows does not invalidate
a threshold-crossing estimator.

Use the same convention in all dimensions, U=1-<m^4>/(2<m^2>^2), and start
with U=0.25 (equivalently Q=2/3). Bracket and interpolate this crossing locally
at each size rather than fitting the whole crossover to a straight line.
Preserve independent-start and block-bootstrap uncertainty, and propagate
shared-data covariance between any threshold variants. Compare U=0.20 and
0.30 where available as economical checks on finite-size drift.

For continuous 3D transitions, a fixed interior threshold has leading shift
beta_U(N)=beta_c+A*N^(-1/nu), with corrections to be checked through minimum
size and simple correction variants. In 4D XY the leading shift instead scales
as N^(-2)*(log N)^(-1/10), with possible logarithmic corrections. In 2D BKT,
thresholds on the disordered side of the limiting critical Binder value can
approach beta_BKT with leading inverse-square-logarithmic drift; do not use
3D power-law extrapolation. The square periodic XY critical value in this U
convention is about 0.4909, so U=0.25 is on the appropriate side. Not every
threshold in the 2D algebraically ordered phase identifies the BKT point.

The core check is whether a common limiting coupling is stable to the tested
size/threshold/shift variants at the requested practical accuracy. Report
wider bars or unresolved slices when it is not. No full crossover model or
free fit of universal critical Binder values is required. Existing data,
frozen manifests, rejected-fit diagnostics and raw retention remain intact.

The paragraphs below describe the earlier, more elaborate baseline.

For d=2, use joint BKT matching of the Binder ratio and second-moment
correlation-length ratio to a calibrated critical reference. Fit the
nonuniversal length rescaling, check corrections and the dependence on the
minimum included size, and restrict reference matching to its validated size
range. A polynomial in 1/N or a fixed Binder threshold is not a substitute.
The published ASOS reference in Hasenbusch equations 46 and 47 is calibrated
on reference sizes 96 through 1024; mapped sizes outside that range require a
different validated reference. Reuse N=16,32,64 data if the size policy allows.

For d=3, test continuous O(2) finite-size scaling using joint magnetic Binder
and correlation-length data. Compare leading correction fits with fits that
include a subleading correction or omit the smallest sizes. Constrain scaling
exponents using an independently sourced universality reference and propagate
their uncertainty. Keep covariance between observables and between any crossing
estimates that share a box size. Do not count such crossings as independent.

For d=4, include the logarithmic corrections at the upper critical dimension.
An initial Binder-ratio ansatz is equation 13 of Lv et al.: a thermal scaling
term, a leading inverse-square-root logarithmic correction, and a subleading
log(log N)/log N correction. Test identifiability and minimum-size stability.
Determine whether the data support the continuous-transition assumption before
using that ansatz for the constrained action.

Calibrate each estimator at epsilon/pi=1 against the ordinary XY endpoint,
without fixing the fit to the published coupling. Also test a constrained
slice near the zero-coupling endpoint, where cutoff-induced corrections may
differ. Successful calibration at the ordinary endpoint alone does not certify
the entire curve.

## Measurement and acceptance procedure

1. Reanalyse compatible existing measurements before choosing new couplings.
   At repeated coordinates, keep the longer independent run as primary input
   and preserve every older record and its provenance.
2. Measure beta neighborhoods adaptively at each exact epsilon slice rather
   than scanning a full rectangular grid. Begin with several box sizes and
   independent cold and admissible winding starts. Freeze each batch manifest
   and runner configuration before execution.
3. Inspect full histories, including the initial transient. Select a supported
   cut and use a stratified contiguous-block bootstrap within independent
   chains. Monitor autocorrelations, independent-start agreement, later cuts,
   doubled blocks, and retained halves for both primary influence series.
4. Fit the infinite-volume critical coupling. Retain fit diagnostics,
   bootstrap failures, covariance, and sensitivities to the beta window,
   interpolation, correction model, and minimum box size. Additional statistics
   address sampling error; additional box sizes address extrapolation stability.
5. Use a 95 percent block-bootstrap statistical interval as the initial
   reporting convention. Expand its endpoints by the tested model and
   interpolation sensitivity envelope. Record the statistical interval and
   sensitivity contributions separately. An envelope does not prove absence
   of all finite-size bias; acceptance also requires calibrated and stable fits.
6. Mark the precision target attained only when both distances from the central
   estimate to the reported interval endpoints are at most 0.02 and the
   diagnostics support the extrapolation. Report a stable extrapolation with
   wider bars when needed. If extrapolation remains unvalidated, report the
   slice as unresolved with its measured information; do not force it onto the
   completed curve.

## Data retention and outputs

Preserve simulation histories, observables, fields, RNG states, checkpoints,
provenance, batch manifests, and full unfiltered runner logs under the existing
retention rules. No new thinning or compaction is included in this design.
Before launching a batch, estimate its disk requirement from the actual record
schema and its work from a measured pilot, and respect the existing runner
space reserve. Historical timing is context, not a guaranteed current runtime.

Save PNG, PDF, and SVG plots and a CSV containing dimension, epsilon/pi,
beta_c, lower and upper vertical bounds, statistical interval, size set,
fit method, sensitivity measures, and diagnostic status. Epsilon is a fixed
input at each plotted point: use yerr only, without horizontal bars or the
earlier adjacent-epsilon-cell union. Do not fabricate completed plots from
the zero-coupling endpoint scouts.

## References and validation evidence

- Existing d=2 inference limits: `experiments/2026-09-02-xy-bkt-refine/boundary.py`
  and `experiments/2026-09-03-xy-bkt-low-kappa/REPORT.md`.
- Existing d=3 and d=4 endpoint evidence: the respective September 3 campaign
  README and REPORT files under `experiments/`.
- [Hasenbusch, The Binder Cumulant at the Kosterlitz Thouless Transition](https://arxiv.org/pdf/0804.1880),
  section 6 and equations 46 and 47.
- [Lv et al., Finite size Scaling of O(n) Systems at the Upper Critical Dimensionality](https://arxiv.org/html/1909.10347v2),
  section III.1 and equation 13.
- [Bietenholz et al., Topological Lattice Actions for the 2d XY Model](https://arxiv.org/pdf/1307.0485),
  section 3, for the limitations of asymptotic spin-wave matching.

The current preparation checked the saved d=2 pooled coordinates, proxy size
dependence, and the published d=2 and d=4 extrapolation prescriptions. No new
simulations or infinite-volume fits have been run for this campaign. The d=3
exponent reference, executable estimator, numerical calibration, and batch
resource measurements remain to be completed before production.


## Outside-bracket field reduction authorized October 14, 2026

The user explicitly authorized removing non-checkpoint configurations for runs clearly outside transition brackets, while preserving all observable/RNG/winding/timing histories, logs, provenance and exact checkpoints. This supersedes historical full-field and dense-diagnostic requirements for the eligible subset only. Numerical run.py and model kernels remain frozen.

The retention daemon refreshes a frozen, source-hashed eligibility manifest when analyzed pointers or completed audited jobs change. For every contributing size and observable (including BOTH sizes of eta pairs), it requires identifiable >=99% bootstrap brackets in all saved cut/block/threshold/interpolation comparisons. It protects the union of measured local windows, one measured neighbor on each side, and a full window-width margin. Unresolved curves, unknown coordinates, and material start/warmup/pooling flags are protected. Low-block-count flags alone do not remove eligibility; all histories and exact checkpoints remain complete. Endpoints are checked separately in epsilon/pi at exactly beta zero, with exact coordinates. The current endpoint comparisons protect their entire measured scans, so no endpoints qualify initially.

Eligible archives retain phi at every 100-cycle checkpoint plus first and final rows. Missing phi rows are explicit unallocated NaNs. All other datasets are compared byte-for-byte; existing full-audit provenance and any previous compaction record are preserved. A second reduction verifies the current compact SHA against the prior committed ledger. Durable prepared records distinguish crashes before replacement (restore prior ledger, preserve abandoned attempt) from after replacement (finish journal commit). The existing shared disk-admission and trajectory locks apply.

Initial stricter manifest: 348 eligible completed jobs (138 D3,210 D4),341.19 GiB current physical files at manifest creation. This is candidate storage, not projected or realized savings. D2 and endpoints remain protected at this stage. Eligibility refresh for future completed runs is automatic, but no run is classified solely by distance from an infinite-volume central estimate. Immutable hash-named manifests preserve every eligibility generation. No historical archive outside this campaign is modified.

Validation: 14 retention, checkpoint/replay, failure/recovery, shared admission and SSD routing checks passed. Test log: logs/outside-bracket-retention-final-tests-20261014.log.


## October 14: practical N512 stopping point and active 5D priority

The user explicitly capped all further 2D simulations at N=512 and requested attained results, then 5D. The two completed N1024 probes remain preserved and excluded from the primary curve. No additional 2D precision campaign is planned. `control.json` carries d2_max_N=512; the storage-aware coordinator blocks pending larger 2D launches without changing run.py or kernels. Seven routing/admission checks pass. All316oldBinderD2jobs remain deferred.

The independently reviewed attained delivery is `d2-d3-d4-attained-N512-with-endpoints-20261014`, including per-dimension and combined PNG/PDF/SVG, CSV, unresolved-locations.csv and attained-phase-boundary-report.md. N512endpoint roots retain99.95–100% draws in the tested cut/block variants. This is conditional coverage, not certified total accuracy. The q1 BKT benchmark discrepancy and wider descriptive unions remain explicit.

5D is enabled and active. All54initial N4/6/8,4000cycle calibration trajectories passed audits:0.9509worker-hours in4.912minutes,11.61effective workers,38.35GiB full audited files. Coarse interpolation is not accepted as final precision. The independently reviewed54-trajectory refinement (27physical coordinates) at N6/8/10,8000cycles is running in d5-local-bracket-and-N10-refinement-20261014. It narrows two positive-beta slices and the beta-zero endpoint and budgets271.13GiB full transient storage. No duplicate batches.

Use analyze_d5_calibration.py and dedicated d5_scaling.py. The leading periodic zero-mode hypothesis is N^-5/2 with N^-3 stress correction, not the4D logarithmic coordinate. Three-size correction fits havezeroDOF and are diagnostic. Four analytic/covariance/zero-mode tests pass. Sources: Chen–Dohm cond-mat/9809394 and Wittmann–Young1410.5296.

The measured pilot supports a planning allowance12–29h for qualified attained5D estimates versus about48h toFriday, not a runtime bound or ±.02 guarantee. The 168-trajectory baseline corresponds to84physical coordinates; after54trajectory reuse,114trajectories(57coordinates) remain. Conservative remaining fullstorage690.28GiB is marginally larger than post-refinement headroom before further actual recovery. Check actualN10 runtime, coverage, retained bytes and fresh reservation-aware space before freezing the full grid. Uniform8kremaininggrid needs1378.62GiB and does not fit current uncompactedcapacity. Current54run stage fits comfortably. See runtime-forecast-d5-calibration-20261014.md/json.

Preserve all histories/checkpoints/logs and numerical kernels. User-authorized adaptive and outside-bracket field retention continue under shared admission. Use atomic JSON updates for control and live pointers. Prioritize the remaining D5 slices on exact.02epsilon/pi increments and its endpoint, with honest uncertainty/flags, then updated combined unconnected MC points against the full tree criterion. Do not restore interpolation lines. AtFridaydeadline stopnew launches,checkpoint safely,deliverattained/unresolvedresults andpause; earlierpause onlyafterD5deliverablesready.
