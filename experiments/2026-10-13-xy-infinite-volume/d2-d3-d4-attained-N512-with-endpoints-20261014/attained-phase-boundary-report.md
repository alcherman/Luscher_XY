# Attained XY phase-boundary estimates through N=512 in 2D

The current 2D estimates support a substantial separation from the supplied tree-level curve within the tested extrapolation family. At epsilon/pi=1 the baseline is beta_c=1.114958, while the full tree-level criterion gives 0.707541; the difference is 0.407418. The descriptive sensitivity interval [1.108818,1.147607] remains well above that prediction. This is a comparison of attained estimates and tested models, not certified total confidence coverage.

## Attained data and observables

There are 22 positive-beta 2D slices at epsilon/pi=.58,.60,...,1.00, 17 reviewed positive-beta 3D/4D points, and three preliminary beta-zero endpoints. 2D uses eta_eff(N,2N)=-log2(m2(2N)/m2(N)) and the second-moment correlation-length ratio. Their limiting BKT targets are 0.25 and 0.7506912; they are not exact finite-box critical values. Matching N=512 neighborhoods exist on representative slices .60, .80,1.00 and at beta=0; the remaining new positive-beta slices reach N=256 with compatible smaller archived anchors. Nothing here asserts uniform N=512 coverage.

The user requested an N=512 cap and a practical stopping point on October 14. Completed N=1024 cost probes are preserved and excluded from the attained primary curve. No further 2D runs above N=512 will launch. No further 2D precision campaign is planned before shifting priority to 5D.

## Monte Carlo intervals versus model sensitivity

Dark bars are conditional bootstrap 95% intervals under the selected local-window and zero-offset logarithmic fit. Pale bars are unions of tested minimum-size, logarithmic-offset and separate/joint observable fit intervals. They include amplified Monte Carlo noise of weak extrapolations and are neither pure systematic errors nor certified total 95% intervals. The baseline q=1 Monte Carlo interval [1.112671,1.117301] still excludes the known 1.1199 benchmark, although the wider descriptive union includes it. Small current-chain cut/block shifts do not settle this finite-size model discrepancy. The scientific comparison above does not require claiming that the small dark intervals certify absolute BKT precision.

The updated positive-beta sensitivity study has 22 successful common-root comparisons with four correlated variants each; it does not validate every root. At q=.60 the latecut local-window comparison omits R256 and eta128:256. An unpromoted full measured-window diagnostic restores 99.9% and 99.75% root retention; this is a window coverage check, not added statistics. September bootstrap anchors retain their original cuts, so current-chain cut checks do not validate archived warmup. Individual flags and all earlier analyses remain available.

## Beta-zero endpoints

The 2D primary endpoint is epsilon_c/pi=0.558567, conditional Monte Carlo 95% [0.555898,0.561744], descriptive union [0.536293,0.568551]. Adding wider N=512 anchors at .546/.572 closes the largest-size bracket tails: R512 and eta256:512 retain 99.95–100% draws across default/current-only, early-cut, late-cut and doubled-block checks. This is coverage within the chosen bootstrap analysis, not total-accuracy certification. Size cuts, offsets and shared-source covariance remain part of the qualifications. Fixed-epsilon 2D endpoint neighbors .54 and .56 have no attained positive-beta estimates (see unresolved-locations.csv); a missing positive finite-box crossing does not establish absence of a thermodynamic positive-beta transition.

3D/4D positive-beta estimates are reused from the reviewed delivery. Their Binder-quarter criterion is followed by dimension-appropriate size extrapolation. Endpoint interpolation/threshold/cut/block checks are saved separately. The 4D endpoint leading fit has p about .036 and remains exploratory; correction and size-cut unions are descriptive. Do not portray any of these unions as calibrated total intervals. Unresolved 3D .80 and 4D .86 neighbors remain explicit in the source analysis.

## Figures and numerical table

Use combined-MC-vs-tree-level.png/pdf/svg or the per-dimension equivalents. MC points are unconnected, fixed-epsilon points have vertical bars only, and beta=0 diamonds have horizontal bars only. The dashed theory is (2D-1)v1(beta,epsilon)=1 from the supplied PDF, page 5; v1 is its normalized exponential-cosine integral. It is a tree instability prediction, not an asserted exact finite-dimensional hypercubic curve. The CSV preserves all 42 attained points and reported interval components.

## Preservation and next priority

All scalar, winding, RNG and timing histories, provenance, complete logs and exact checkpoints remain intact. User-authorized adaptive/configuration-only reductions operate after full closure audits. Original full-file hashes and current archive hashes are separately recorded. Completed N=1024 probes are retained. New 5D work begins with reviewed 54-trajectory N=4, 6 and 8 cost/coverage calibration, then a measured forecast; its above-upper-critical-dimension analysis will not substitute the 4D logarithmic formula.

Sources: experiments/2026-10-13-xy-infinite-volume/mc-vs-tree-level-20261014-141630; latest-d2-eta-correlation.json; latest-zero-beta-endpoints.json; latest-d2-current-chain-sensitivity.json; latest-endpoint-sensitivity.json; latest-d3-d4-delivery.json. Numerical model and run.py remain frozen. Twelve endpoint/GLS checks and seven size-cap/storage-admission checks pass; prior retention and replay checks remain preserved.
