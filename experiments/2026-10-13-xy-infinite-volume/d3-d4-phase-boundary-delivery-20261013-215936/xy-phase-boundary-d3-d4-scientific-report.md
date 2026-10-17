# XY phase boundary in three and four dimensions

The fixed Binder quarter estimator gives 10 conditional 3D and 7 conditional 4D transition locations on the requested epsilon grid. The analysis considers 788 audited 3D and 728 audited 4D trajectories, including endpoint and diagnostic runs. The largest tested sensitivity half-width is 0.01659 in beta. These are attained estimates with explicit sensitivity bars; total error coverage at the requested 0.02 has not been established.

## Transition estimates

| d | epsilon/pi | beta | Conditional MC 95 percent interval | Tested sensitivity union |
|---|---|---|---|---|
| 3 | 0.80 | unresolved neighbor | | |
| 3 | 0.82 | 0.06111 | [0.05990, 0.06235] | [0.04808, 0.06514] |
| 3 | 0.84 | 0.12138 | [0.12042, 0.12232] | [0.10479, 0.12469] |
| 3 | 0.86 | 0.17744 | [0.17534, 0.17910] | [0.16401, 0.18347] |
| 3 | 0.88 | 0.22667 | [0.22525, 0.22820] | [0.21462, 0.23153] |
| 3 | 0.90 | 0.27223 | [0.27037, 0.27361] | [0.26202, 0.28086] |
| 3 | 0.92 | 0.31544 | [0.31453, 0.31614] | [0.31071, 0.32173] |
| 3 | 0.94 | 0.35453 | [0.35306, 0.35610] | [0.34368, 0.36591] |
| 3 | 0.96 | 0.39080 | [0.38945, 0.39209] | [0.37988, 0.39761] |
| 3 | 0.98 | 0.42360 | [0.42230, 0.42482] | [0.41712, 0.42861] |
| 3 | 1.00 | 0.45497 | [0.45408, 0.45565] | [0.45014, 0.45849] |
| 4 | 0.86 | unresolved neighbor | | |
| 4 | 0.88 | 0.04758 | [0.04713, 0.04788] | [0.04526, 0.05179] |
| 4 | 0.90 | 0.09765 | [0.09715, 0.09810] | [0.09034, 0.09856] |
| 4 | 0.92 | 0.14569 | [0.14538, 0.14595] | [0.14357, 0.14713] |
| 4 | 0.94 | 0.18886 | [0.18805, 0.18946] | [0.18396, 0.19707] |
| 4 | 0.96 | 0.22996 | [0.22972, 0.23011] | [0.22807, 0.23451] |
| 4 | 0.98 | 0.26777 | [0.26732, 0.26810] | [0.26495, 0.27349] |
| 4 | 1.00 | 0.30176 | [0.30146, 0.30210] | [0.29908, 0.30560] |

## Estimator and finite size analysis

The convention is U = 1 - <m^4>/(2<m^2>^2). At each size, weighted monotone regression pools only nonmonotone measured Binder values, and local interpolation brackets U=0.25. The same procedure acts on each bootstrap draw. Monotonicity is a declared assumption, not an independent validation. The primary interpolation is linear; PCHIP is a sensitivity check. No rejected broad observable fit supplies a critical coupling.

The leading shift coordinate is N^(-1/0.67169) in 3D and N^(-2)(ln N)^(-1/10) in 4D. The primary size cuts are N>=6 in 3D and N>=10 in 4D. The uniform 4D cut was chosen after inspecting the pilots and small-box residuals; its goodness-of-fit p values are exploratory and do not have a confirmatory selection interpretation. All-size fits, further size cuts, and a subleading correction remain in the sensitivity calculations. This uses modest-box extrapolation and does not prove the asymptotic law.

Independent review of the earlier N24 holdout check found observed-minus-predicted quarter-root differences of -0.000034, -0.000224, and +0.000266 at 4D epsilon/pi=0.88,0.92,0.96, respectively, each with a conditional bootstrap interval including zero. The comparison supports small sensitivity within tested sizes, without certifying the thermodynamic limit.

## Monte Carlo uncertainty and tested sensitivity

The default analysis discards cycles through 1000 and uses 1000 independent-coordinate, contiguous-block bootstrap draws, stratified by independent starts and preserving joint observables. Block widths respond to measured autocorrelations. Dark plot bars show the primary conditional MC 95 percent interval. Light bars are the union of intervals across the saved threshold, interpolation, minimum-size, correction, transient-cut and block-length variants. The saved variant grid is not full factorial: transient-cut and block-length checks use U=0.25. The light union includes rejected exploratory fits and is a descriptive stress test, not a calibrated confidence interval. Sensitivity and MC uncertainty are kept as separate numerical columns.

Many groups retain fewer than 50 blocks per start. Initial start disagreement at 3D epsilon/pi=.92,N16,beta=.293 and 4D .98,N10,beta=.261 remains recorded. New independent 16000-cycle starts have descriptively agreeing late halves; window checks do not certify mixing. Removing those below-crossing coordinates leaves the central primary local-linear roots unchanged; a 4D PCHIP root shifts by 0.000190. Bootstrap influence from those deletions was not independently quantified.

Missing positive-beta quarter crossings at 3D .80 and 4D .86 are reported as unresolved branch neighbors. Their zero-beta finite-size behavior is consistent with an ordered side, but does not alone certify absence of a positive thermodynamic transition. No horizontal uncertainty bars or artificial endpoint coordinates are introduced.

## Validation and retained evidence

Every contributing trajectory passed independent field/observable recomputation and its closure audit. The frozen numerical sources remain unchanged. Exact checkpoint replay and focused XY tests were validated before production; nine focused threshold/calibration/storage tests passed during final analysis development. Correct generation is separate from demonstrated mixing or extrapolation accuracy.

The contributing raw HDF5 files occupy 1.358 TiB. Raw histories, fields, RNG states, checkpoints, all earlier snapshots and full unfiltered logs remain retained on the verified Office_SSD route. D2 remains deferred at the user request and d5 remains disabled.

Source snapshots:
- [threshold-curve-snapshot-20261013-215755](/Users/acherman/research/Luscher_XY/experiments/2026-10-13-xy-infinite-volume/threshold-curve-snapshot-20261013-215755/analysis.json)
- [threshold-curve-snapshot-20261013-215837-cut501-block1](/Users/acherman/research/Luscher_XY/experiments/2026-10-13-xy-infinite-volume/threshold-curve-snapshot-20261013-215837-cut501-block1/analysis.json)
- [threshold-curve-snapshot-20261013-215837-cut2001-block1](/Users/acherman/research/Luscher_XY/experiments/2026-10-13-xy-infinite-volume/threshold-curve-snapshot-20261013-215837-cut2001-block1/analysis.json)
- [threshold-curve-snapshot-20261013-215837-cut1001-block2](/Users/acherman/research/Luscher_XY/experiments/2026-10-13-xy-infinite-volume/threshold-curve-snapshot-20261013-215837-cut1001-block2/analysis.json)

Reproducible outputs:
- [Numerical table](/Users/acherman/research/Luscher_XY/experiments/2026-10-13-xy-infinite-volume/d3-d4-phase-boundary-delivery-20261013-215936/transition-locations.csv)
- [Analysis and all tested variants](/Users/acherman/research/Luscher_XY/experiments/2026-10-13-xy-infinite-volume/d3-d4-phase-boundary-delivery-20261013-215936/phase-boundary-analysis.json)
- [Combined plot PDF](/Users/acherman/research/Luscher_XY/experiments/2026-10-13-xy-infinite-volume/d3-d4-phase-boundary-delivery-20261013-215936/combined-transition-locations.pdf)
