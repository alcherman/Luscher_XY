# XY Binder experiment — 2026-09-02

Completed 21 independent cold-start chains on periodic square lattices, with
L = 4, 8, 16 and κ = 0.80, 0.95, 1.05, 1.12, 1.20, 1.35, 1.50.
Each chain has exactly 10,000 sweeps plus its initial configuration at sweep zero.
The action is S = −κ Σ cos θ, counting each positive-direction bond once, with
ε = π. The exclusion of exactly antipodal bonds has zero measure in the usual XY ensemble.

![Binder cumulant](binder_vs_kappa.png)

[Main plot PDF](binder_vs_kappa.pdf) · [Main plot SVG](binder_vs_kappa.svg) ·
[Near-transition zoom](binder_zoom.png) · [CSV table](binder.csv)

The vertical line is an external infinite-volume reference, κ ≈ 1.12, from
[Komura and Okabe, arXiv:1210.6116](https://arxiv.org/abs/1210.6116), who report
1.11996(6). It was not fitted to these data. Lines between points only guide the eye.

## Estimator and uncertainty

With m the spatial mean of exp(iφ), we compute

$$
U_L = 1 - \frac{\langle |m|^4\rangle}{2\langle |m|^2\rangle^2}.
$$

One sweep visits every site once with a symmetric Metropolis proposal of half-width π,
then attempts one symmetric global winding update. Three chains ran concurrently,
using serial Numba within each chain; seeds and source provenance are recorded.

The plots in `histories/` retain every stored observable from sweep zero through
10,000. The three `thermalization_*.png` grids show the first 500 sweeps across
all 21 points. After inspecting the visible relaxation of |m|², energy and vortex
density, the common analysis cut was set to discard the cold seed and sweeps 1–500.
This leaves 9,500 measurements per chain; no measurements were thinned or deleted.
The HDF5 `analysis_start_row=101` attribute is the generation-time provisional
setting; **the final analysis uses row 501**, as recorded in `analysis-provenance.json`.

Autocorrelation times are estimated on the retained series by a nonperiodic FFT
and an initial-positive, monotone paired correlation sum, in the convention
Var(mean) ≈ 2 τ_int Var(series)/n. Blocks are at least ten times the largest
estimated τ_int among |m|², |m|⁴, the linearized Binder estimator, energy and vortex
density. This is an empirical estimate, not an upper bound on mixing time.
Global magnetization orientation is a separate diagnostic and does not determine
the block width for this rotation-invariant estimator.

Block widths range from 39 to 521 sweeps, giving 18–243 complete blocks. The point
estimates and errors use complete blocks (9,240–9,500 measurements); the remaining
0–260 measurements stay in the raw files. The two moments are bootstrapped jointly,
then combined into the ratio in each of 4,000 resamples. Errors below and in the
plots are bootstrap standard deviations, approximately one standard error.
Per-point τ_int, widths, block counts, used counts, acceptance rates and seeds
are in `results.json` and `binder.csv`; all bootstrap samples are saved.

## Results

| κ | L = 4 | L = 8 | L = 16 |
|---:|---:|---:|---:|
| 0.80 | 0.43830 ± 0.00278 | 0.36757 ± 0.00640 | 0.18468 ± 0.02177 |
| 0.95 | 0.47489 ± 0.00119 | 0.46248 ± 0.00282 | 0.42770 ± 0.01166 |
| 1.05 | 0.48248 ± 0.00081 | 0.48238 ± 0.00109 | 0.47878 ± 0.00270 |
| 1.12 | 0.48712 ± 0.00057 | 0.48823 ± 0.00081 | 0.49083 ± 0.00068 |
| 1.20 | 0.49043 ± 0.00050 | 0.49221 ± 0.00049 | 0.49367 ± 0.00071 |
| 1.35 | 0.49343 ± 0.00029 | 0.49487 ± 0.00035 | 0.49581 ± 0.00041 |
| 1.50 | 0.49542 ± 0.00022 | 0.49674 ± 0.00016 | 0.49712 ± 0.00021 |

At κ = 0.80 and 0.95, U decreases substantially with lattice size. Near and above
κ = 1.12 the values are close on the full vertical scale. The zoom resolves
remaining size dependence, including U increasing with L on the upper branch.
These finite lattices do not exhibit exact curve collapse and do not provide a
precision determination of the BKT point.

## Stability and limits

Cuts discarding 500, 1,000 and 2,000 sweeps and block widths of 1×, 2× and 4× the
primary width are recorded in `stability.json`; see `cut_stability.png` and
`autocorrelation_blocks.png`. The largest cut-induced shift is 1.17 times the
primary error, at L = 16, κ = 1.35, when discarding 2,000 sweeps. These are correlated
estimates from the same chain. Doubling the block width changes the estimated
standard errors by factors 0.82–1.10; larger widths also reduce the number of blocks.

The slowest point, L = 16, κ = 0.95, has estimated τ_max = 52.1 sweeps and only
18 primary blocks of width 521. Its uncertainty is particularly approximate.
Four consecutive quarter estimates are retained: the first and last differ by
2.90 combined standard errors at L = 4, κ = 1.35 and 2.41 at L = 8, κ = 1.12.
The full histories show fluctuations without an obvious persistent late trend,
but the quarter differences prevent treating the diagnostics as a blanket
stationarity certificate. There are no independent repeat chains or hot starts.

At L = 16 and κ ≥ 1.20 no global winding proposal was accepted in the retained
window. Individual loop winding values do fluctuate through local vortex events;
this alone does not establish sampling of rare global winding sectors. These
results are a short exploratory scan, with estimated statistical errors that
exclude any unresolved equilibration or rare-sector bias. No finite-size scaling
fit or thermodynamic extrapolation was performed.

## Verification and retained evidence

`audit.py` independently recomputed magnetization moments, action density and maximum
link angle directly from **all 210,021 saved configurations**, checked the canonical
angle interval, strict admissibility, parameter grid, cold seeds and row counts.
All checks passed; the largest absolute observable residual was 8.9e-16.
`data-audit.json` includes SHA-256 checksums of every HDF5 file.
The focused XY test suite passed 46 tests with one thread-parity test skipped
under NUMBA_NUM_THREADS=1 (`validation.log`). This validates the recorded checks,
not equilibrium sampling at every point.

`chains/` preserves every field, observable and RNG state, including warmup.
`manifest.json`, `status.json`, HDF5 provenance, `run.log` and
`analysis-provenance.json` record the run and analysis. The earlier low-κ scan was
terminated and retained separately in `../2026-09-02-xy-binder`; it contributes no
points to these plots.

## Reproduce from the repository root

The existing files are complete. The disk-only commands below regenerate plots,
estimates and validation without adding sweeps:

```sh
export MPLBACKEND=Agg
export MPLCONFIGDIR=/private/tmp/luscher-mpl
export NUMBA_CACHE_DIR=/private/tmp/luscher-numba
export NUMBA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
.venv/bin/python experiments/2026-09-02-xy-binder-bkt/analyze.py --histories
.venv/bin/python experiments/2026-09-02-xy-binder-bkt/analyze.py --cut 501
.venv/bin/python experiments/2026-09-02-xy-binder-bkt/audit.py
```

The original generation command, with the same environment, was:

```sh
.venv/bin/python experiments/2026-09-02-xy-binder-bkt/experiment.py run > experiments/2026-09-02-xy-binder-bkt/run.log 2>&1
```

Use `analyze.py` for the final analysis. The generation script is preserved byte
for byte as run, including its earlier provisional analysis routine.
