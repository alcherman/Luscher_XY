# Attained 5D phase boundary and combined comparison

All282 frozen5D trajectories passed full-field closure audits. There are six positive-beta slices q=epsilon/pi=.90,.92,.94,.96,.98,1.00 and a beta=0 scan, with N=6,8,10,12 and independent starts. N4 pilot data are preserved but excluded from the uniform baseline. The primary contains31linear finite-size roots, all100% bootstrap-valid, with62linear/PCHIP roots total. Audits establish recorded-data correctness, not mixing or infinite-volume accuracy.

Five positive-beta slices and the endpoint have conditional baseline size fits; q=1 fails that diagnostic. Its nominal extrapolation is shown as an OPEN SQUARE only and is not an accepted critical coupling. No larger-box campaign is required for this practical attained delivery. Rejected q1 leading/correction fits and the exploratory minN8 alternative remain visible in the complete stress ledger.

## Locations and uncertainty

Positive rows give beta, the beta=0 row gives epsilon/pi. These are descriptive attained estimates under the method below.

| Scan | Baseline location | Conditional MC95 | Descriptive stress envelope | Status |
|---|---:|---|---|---|
| 0.90 | 0.015824 | [0.015313, 0.016245] | [0.011339, 0.017857] | conditional |
| 0.92 | 0.065534 | [0.065149, 0.065888] | [0.062180, 0.067005] | conditional |
| 0.94 | 0.109835 | [0.109372, 0.110166] | [0.107201, 0.114060] | conditional |
| 0.96 | 0.153482 | [0.152886, 0.153936] | [0.149725, 0.156966] | conditional |
| 0.98 | 0.192668 | [0.192052, 0.193062] | [0.189252, 0.198050] | conditional |
| 1.00 | 0.229017 | [0.228766, 0.229222] | [0.228520, 0.232217] | rejected-fit diagnostic |
| beta=0 | 0.893197 | [0.893033, 0.893324] | [0.891906, 0.894762] | conditional |

The uniform reference is locally bracketed linear U=.25 followed by beta_N=beta_infinity+a*N^-2.5 (or epsilon_N at beta=0), N>=6. This periodic above-upper-critical-dimension hypothesis is not the4D logarithmic law. O(2) quartic zero-mode Binder U*=1-pi/4=.2146: U=.25 is a shifted ordered-side fixed criterion, not the universal finite-box critical value. An additiveN^-3 correction is a stress hypothesis. Source references are Chen–Dohm cond-mat/9809394 and Wittmann–Young1410.5296, as documented in the campaign design.

Dark bars are conditional2000-draw block-bootstrap intervals, conditional on central-data local-window selection, autocorrelation/block estimates and the size ansatz. For q1 that ansatz is rejected and its interval is only diagnostic. Pale bars unite intervals from9source-pinned snapshots: baseline/fullwindow/cut2001/cut3001/doubledblocks/U.20/U.30/even/oddstarts; each tests linear/PCHIP, N>=6 or8 leading, and N>=6 correction. Only identifiable, positive-DOF fits with all required sizes enter. Rejected fits and flagged inputs remain INCLUDED in this stress envelope and explicitly flagged in d5-stress-ledger.json. Missing sizes, insufficient bootstrap retention, rank failures and zeroDOF fits are separately recorded. The envelope includes amplified Monte Carlo variance, threshold dependence and rejected-model diagnostics: it is neither pure systematic uncertainty nor a certified total95% interval. A small envelope does not demonstrate totalaccuracy. No model is selected by bestp. Even/odd draws share helper seeds and are not paired for difference significance.

## Sampling and interpolation limits

There are144pooled coordinate groups, with group flag counts{'fewer_than_50_blocks_per_start': 74, 'start_disagreement_Q': 2, 'start_disagreement_binder': 2}. The smallest retained block count per start is7. Selected-root flags and all-measured flags are preserved; absence of start/half flags does not certify equilibration. Long autocorrelations leave few blocks for some N12starts, limiting the reliability of small conditional intervals. Across the nine variants, all 28 N>=6 linear roots remain identifiable (minimum valid fraction .998), but variant root flags include Binder half-history shifts, insufficient half-history blocks and start disagreement. These warnings remain material qualifications; see d5-variant-sampling-flag-summary.json for per-variant counts. The complete first/last-half, per-start, pooled, autocorrelation and block diagnostics remain in the pinned analysis snapshots. Shared beta=0 q=.90 sources used on both scan axes are correlated evidence, not additional independent replication.

The final matched midpoints changed finite-size roots by as much as .003584 relative to the earlier two-knot windows; this comparison also changes measured sources and is not a pure interpolation-bias estimate. Most local positive-beta linear/PCHIP comparisons still use two knots and are identical by construction. The full measured-window PCHIP comparison adds a shape sensitivity check, but does not certify interpolation coverage or exclude curvature between samples. Local3SD windows are data-selected. All unresolved/tail losses and skipped groups in every variant remain explicit; successful common-root comparisons are not all-root validation.

## Comparison and prior dimensions

The five conditionally fitted positive-beta points lie above the tree prediction by beta differences .00619 to .01165. Their descriptive stress envelopes remain above that prediction. The beta-zero endpoint is epsilon/pi=.893197 versus tree .898458, a difference -.005261; its envelope [.891906,.894762] is below the tree endpoint. These are absolute descriptive comparisons within the tested families, not calibrated significance or guaranteed total accuracy. The q1 nominal difference .005409 is diagnostic only because its reference size fit rejects.

Dashed theory solves the full supplied (2D-1)v1=1 criterion, v1 the normalized exponential-cosine integral over[-epsilon,epsilon]. It is the regular-tree instability prediction, not an exact hypercubic critical curve. The supplied PDF definition is onpage2 and orangebox onpage5; itsSHA is recorded in analysis.json. MC points are unconnected. Fixed-epsilon points have vertical bars only; beta=0 endpoints horizontal bars only. Open D5square denotes rejected-fit diagnostic, not a validated extrapolation.

The42 attained D2/D3/D4 rows and their errors are copied exactly, with no new2Danalysis. See prior-d2-d3-d4-qualifications.md. In particular D2q1's dark interval excludes1.1199; its wider descriptive envelope includes it, so certified absolute BKT accuracy is not claimed. D2 uses N<=512 in its primary, with512only on representative slices and endpoint. Completed1024cost probes remain preserved/excluded. Prior unresolved neighbors remain in prior-unresolved-locations.csv. The fixed D5q=.88 positive slice was not measured; no absence inference is made from missing finite-box crossings.

## Provenance and preservation

analysis.json pins all9snapshot paths and hashes, the prior delivery hash and theory PDF hash. d5-finite-size-locations.csv contains root windows/status/flags; d5-stress-ledger.json identifies every included or excluded stress fit. All histories, logs, provenance and exact checkpoints remain intact. User-authorized field-only retention runs after full audits. No numerical simulation kernels or running environment were changed. This folder remains pending independent scientific and visual review until delivery-status.json says reviewed_attained_with_qualifications.
