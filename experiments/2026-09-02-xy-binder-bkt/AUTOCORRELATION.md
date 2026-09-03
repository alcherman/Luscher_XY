# Autocorrelation times versus lattice size

This is a disk-only analysis of the existing 28 chains. No new sweeps were run.
All estimates use sweeps 501–10,000, after the previously selected thermalization cut.

![Size dependence](autocorrelation_vs_L.png)

Times are integrated autocorrelation estimates in full-lattice sweeps:

$$
\tau_{\mathrm{int},A}=\frac12+\sum_{t\geq1}\rho_A(t),\qquad
\operatorname{Var}(\overline A)\simeq\frac{2\tau_{\mathrm{int},A}}{n}\operatorname{Var}(A).
$$

Independent samples therefore have τ_int = 1/2 in this convention. The estimates
use the initial-positive, monotone paired correlation sum implemented in
`experiment.py`. This is a finite-run estimator, not a bound on the slowest
relaxation time. No uncertainty bars on τ or exponential-time fits were computed.

The Binder column uses the linearized combination of the two moments that controls
the uncertainty of the ratio. With a_t = |m_t|² and b_t = |m_t|⁴, its series is

$$
h_t=\frac{\overline b}{\overline a^3}(a_t-\overline a)
-\frac{b_t-\overline b}{2\overline a^2}.
$$

It does not assign a Binder cumulant to each individual configuration. The maximum
in the right panel is over h, |m|², |m|⁴, action density and vortex density. It is
the blocking diagnostic, not the maximum over every observable of the Markov chain.

## Binder-relevant integrated times, in sweeps

| κ | L = 4 | L = 8 | L = 16 | L = 32 |
|---:|---:|---:|---:|---:|
| 0.80 | 3.1 | 6.1 | 9.7 | 7.5 |
| 0.95 | 2.8 | 7.6 | 52.1 | 84.7 |
| 1.05 | 2.6 | 6.5 | 22.6 | 65.3 |
| 1.12 | 2.8 | 6.9 | 10.0 | 54.9 |
| 1.20 | 3.4 | 5.9 | 13.6 | 168.7 |
| 1.35 | 3.1 | 6.9 | 22.8 | 31.8 |
| 1.50 | 2.6 | 4.7 | 11.2 | 45.0 |

## Interpretation

- At κ = 0.80, the Binder time stays below about 10 sweeps and the L32 estimate
  is smaller than L16. The measured size dependence flattens; there is no evidence
  here for continuing L² growth at that coupling.
- Near the transition, there is strong slowing as L grows. For example, at κ = 1.05,
  the |m|² times for L = 8, 16, 32 are 8.7, 30.7, 106.7 sweeps. At κ = 1.12 they
  are 8.0, 17.8, 117.4. The Binder ratio often decorrelates faster because some
  fluctuations cancel between numerator and denominator.
- A roughly quadratic increase is a useful comparison near criticality. Local
  Metropolis dynamics for the square-lattice XY model have been found compatible
  with z ≈ 2 [Michel, Mayer and Krauth, Fig. 5](https://arxiv.org/pdf/1508.06541).
  Our dashed lines are arbitrary-normalization L² guides, not fits. These small
  lattices and short single chains do not establish z or deviations from it.
- At L32, κ = 1.20 is exceptional: τ_B ≈ 169 and τ_max ≈ 235 sweeps, with a large
  late magnetization excursion. Four primary blocks do not establish the rate of
  such excursions. Its long times should not be interpreted as a new critical point.
- Global magnetization orientation is much slower: Re(m) and Im(m) estimates reach
  roughly 1,000–2,000 sweeps on larger lattices. They are poorly determined by
  9,500 measurements. A common rotation leaves |m|² and the Binder ratio unchanged,
  so these orientation times do not directly set the Binder error.

![Near-critical autocorrelation functions](autocorrelation_functions_k1.12.png)

At κ = 1.12 the L32 |m|² ACF has a long positive tail extending several hundred
sweeps. The Binder influence ACF also shows long-lag structure, and is not well
summarized by assuming a single exponential. Finite-run noise limits interpretation
of the detailed long-lag shapes. The plotted ACFs are normalized covariance
estimates; no physical-time dynamics is inferred from Monte Carlo sweep time.

## Reliability and reproducibility

Six of the seven L32 points have only 4–15 primary analysis blocks. Even estimates
with more blocks have no calculated τ uncertainty. Recomputing after discarding
500, 1,000 or 2,000 sweeps gives, for example, τ_B = 54.9, 50.5, 46.5 and
τ_max = 117.4, 105.6, 87.2 at L32, κ = 1.12. These same-chain cut checks are useful
sensitivity diagnostics, not independent error bars or an equilibration proof.

`autocorrelation_times.csv` contains all monitored times, seeds, cuts and block counts.
`autocorrelation_cut_checks.json` contains the additional cut checks, and
`autocorrelation_functions_k1.12.csv` contains the plotted ACF values. The reporting
script reproduces all stored primary times and checks its FFT covariance at four
lags against direct products. `autocorrelation-report-provenance.json` records the
source and input-results hashes. Existing Binder estimates and raw chains are unchanged.

From the repository root:

```sh
MPLBACKEND=Agg MPLCONFIGDIR=/private/tmp/luscher-mpl NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python experiments/2026-09-02-xy-binder-bkt/autocorrelation_report.py
```
