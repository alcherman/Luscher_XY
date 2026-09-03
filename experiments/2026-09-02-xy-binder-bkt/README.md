# XY Binder experiment — L = 4, 8, 16, 32

Completed 28 cold-start chains on periodic square lattices, at κ = 0.80, 0.95,
1.05, 1.12, 1.20, 1.35, 1.50 and ε = π. Each chain has exactly 10,000 sweeps
plus its initial configuration at sweep zero. The action is S = −κ Σ cos θ,
counting positive-direction bonds once. The L32 extension adds seven chains;
the original 21 chains and their numerical estimates are unchanged.

![Binder cumulant](binder_vs_kappa.png)

[Main PDF](binder_vs_kappa.pdf) · [Main SVG](binder_vs_kappa.svg) ·
[Near-transition zoom](binder_zoom.png) · [CSV table](binder.csv) ·
[Original analysis before L32](before-L32/README.md)

The vertical line is the external infinite-volume reference κ ≈ 1.12 from
[Komura and Okabe](https://arxiv.org/abs/1210.6116), not a fit to these data.
Connecting lines guide the eye.

## Observable and run setup

With m the spatial mean of exp(iφ),

$$
U_L=1-\frac{\langle |m|^4\rangle}{2\langle |m|^2\rangle^2}.
$$

A sweep attempts one symmetric Metropolis update per site (proposal half-width π),
followed by one symmetric global winding proposal. Three chains run concurrently,
using serial Numba within each chain. The L32 seeds are 2026090422–2026090428.
Every configuration, observable and RNG state is saved without thinning.

All estimates discard the cold seed and sweeps 1–500, leaving 9,500 measurements.
For L32 this cut was selected after examining the cold-start relaxation and full
histories, available in [thermalization_L32.png](thermalization_L32.png),
[overview_L32.png](overview_L32.png) and the per-point plots in `histories/`.
Visible initial relaxation occurs within a few hundred sweeps. This observation
is not proof that every slow mode is equilibrated. The generation-time HDF5
`analysis_start_row=101` is provisional; the actual analysis starts at row **501**.

## Results

Errors are estimated one-standard-error block-bootstrap uncertainties.
For L32, six of seven points have fewer than 20 blocks, so the quoted errors
are provisional and do not establish precision equilibrium results.

| κ | L = 4 | L = 8 | L = 16 | L = 32 |
|---:|---:|---:|---:|---:|
| 0.80 | 0.43830 ± 0.00278 | 0.36757 ± 0.00640 | 0.18468 ± 0.02177 | 0.06295 ± 0.02721 |
| 0.95 | 0.47489 ± 0.00119 | 0.46248 ± 0.00282 | 0.42770 ± 0.01166 | 0.37823 ± 0.02676 |
| 1.05 | 0.48248 ± 0.00081 | 0.48238 ± 0.00109 | 0.47878 ± 0.00270 | 0.48296 ± 0.00381 |
| 1.12 | 0.48712 ± 0.00057 | 0.48823 ± 0.00081 | 0.49083 ± 0.00068 | 0.49080 ± 0.00132 |
| 1.20 | 0.49043 ± 0.00050 | 0.49221 ± 0.00049 | 0.49367 ± 0.00071 | 0.48850 ± 0.00509 |
| 1.35 | 0.49343 ± 0.00029 | 0.49487 ± 0.00035 | 0.49581 ± 0.00041 | 0.49682 ± 0.00032 |
| 1.50 | 0.49542 ± 0.00022 | 0.49674 ± 0.00016 | 0.49712 ± 0.00021 | 0.49689 ± 0.00055 |

The L32 values extend the strong size dependence at low κ and lie near the
smaller-lattice values around κ = 1.12. The apparent dip at κ = 1.20 comes with
a large uncertainty and a late magnetization excursion. It is not evidence for
a physical nonmonotonicity. No transition fit or thermodynamic extrapolation was made.

## L32 autocorrelations and uncertainty

Autocorrelation times are estimated on post-cut data using a nonperiodic FFT
and initial-positive monotone paired correlation sums, with
Var(mean) ≈ 2 τ_int Var(series)/n. The largest estimate among |m|², |m|⁴, the
linearized Binder estimator, action density and vortex density sets a target
block width of at least 10 τ_int. Global magnetization orientation and winding
are retained as separate diagnostics; orientation does not set the block size
for the rotation-invariant Binder ratio.

At L32, balanced contiguous blocks include **all 9,500 post-cut measurements**;
block sizes differ by at most one sweep and all meet the target width. This
avoids excluding a long remainder, which could otherwise discard up to 1,510
measurements. In each of 4,000 bootstrap resamples, the joint moment sums and
block lengths are resampled together before evaluating the ratio. The reported
point estimate uses all retained measurements. The original L4–L16 estimates
retain their previous complete-block method and values; their unused remainders
remain in the raw files. The balanced-block option does not rerun those estimates.

| κ (L32) | Estimated τ_max (sweeps) | Actual block width | Blocks |
|---:|---:|---:|---:|
| 0.80 | 16.8 | 169–170 | 56 |
| 0.95 | 159.8 | 1900 | 5 |
| 1.05 | 106.7 | 1187–1188 | 8 |
| 1.12 | 117.4 | 1187–1188 | 8 |
| 1.20 | 235.0 | 2375 | 4 |
| 1.35 | 122.3 | 1357–1358 | 7 |
| 1.50 | 62.6 | 633–634 | 15 |

The κ = 1.20 chain has a pronounced late excursion in |m|², accompanied by winding
updates near the end. It is fully included. Only four blocks remain, so this
trajectory cannot establish the frequency or equilibrium weight of such excursions.
The κ = 0.95 estimate also has only five blocks. These limitations are not cured
by the number of bootstrap resamples.

Cuts discarding 500, 1,000 and 2,000 sweeps shift L32 estimates by at most 0.45 times
the primary error. These are correlated comparisons. Doubling the target block
width changes L32 errors by factors about 0.64–1.17, while leaving only 2–28 blocks.
Some four-times-width checks and quarter-window estimates have fewer than two
blocks and are explicitly recorded as unavailable, not as zero error.
See `stability.json`, `cut_stability.png`, and `autocorrelation_blocks.png`.

There are no independent repeat chains or hot starts. At L32, no global winding
proposal was accepted in the retained window at κ = 1.12, 1.35 or 1.50. Individual
loop winding fluctuations alone do not establish sampling of rare global sectors.
The plotted errors exclude unresolved equilibration and rare-sector bias.
The earlier L4–L16 limitations remain documented in `before-L32/README.md`.

## Verification and evidence

`audit.py` independently recomputed moments, action density and maximum bond angle
from **all 280,028 configurations**, and checked the grid, row counts, cold seeds,
canonical angle range and strict admissibility. All checks passed. It also verified
that every original HDF5 checksum is unchanged and that sampler-source fingerprints
match across all 28 chains. Results and checksums are in `data-audit.json`.
The only exclusions at ε = π are exactly antipodal bonds, a measure-zero set.

A targeted analysis check verified inclusion of an artificial terminal excursion,
the minimum balanced-block width, all 28 unique estimates, unchanged original
results and use of all 9,500 L32 measurements (`analysis-checks-l32.log`). The sampler
code was unchanged for this extension; the earlier focused suite had 46 passes and
one thread-parity skip (`validation.log`) and was not rerun for this data addition.

The original generation driver and manifest are preserved. L32 generation is
recorded in `extend_l32.py`, `manifest-l32.json`, `status-l32.json` and `run-l32.log`.
The combined analysis uses `dataset.py` and records its source fingerprints and
cuts in `analysis-provenance.json`. Raw chains live in `chains/`, full histories in `histories/`, and bootstrap samples
in the experiment directory. Original analysis outputs are archived in `before-L32/`.

## Reproduction

From the repository root, set:

```sh
export MPLBACKEND=Agg
export MPLCONFIGDIR=/private/tmp/luscher-mpl
export NUMBA_CACHE_DIR=/private/tmp/luscher-numba
export NUMBA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
```

The L32 generation command was:

```sh
.venv/bin/python experiments/2026-09-02-xy-binder-bkt/extend_l32.py > experiments/2026-09-02-xy-binder-bkt/run-l32.log 2>&1
```

These commands only analyze stored data; the size selection preserves the original
L4–L16 estimates when updating L32:

```sh
.venv/bin/python experiments/2026-09-02-xy-binder-bkt/analyze.py --histories --sizes 32
.venv/bin/python experiments/2026-09-02-xy-binder-bkt/analyze.py --cut 501 --sizes 32 --include-tail
.venv/bin/python experiments/2026-09-02-xy-binder-bkt/audit.py
```

To reconstruct all estimates from raw files, analyze `--sizes 4 8 16 --cut 501`
without `--include-tail`, then use the L32 analysis command above. The original
`experiment.py` remains unchanged as a generation record; its provisional analysis
routine is superseded by `analyze.py`.

Autocorrelation follow-up: [size dependence, ACFs and full time table](AUTOCORRELATION.md).
