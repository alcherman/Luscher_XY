# D = 2, L = 48 cluster experiment

Completed the same seven κ values at ε = π with S = −κ Σ cos(dφ). Each of seven independent cold-start chains has 10,000 saved update cycles: one Wolff cluster followed by one symmetric global winding proposal. These cycles are not Metropolis volume sweeps.

Generation took 124.5 seconds of wall time with three worker processes. All 70,007 configurations, including the cold seeds, full fields, observables, RNG states and cluster work, are retained. Every saved field passed an independent action/magnetization/vortex/winding audit. The previous Binder results were read only and their checksum is unchanged.

![Binder comparison](binder_vs_kappa.png)

## Analysis

After inspecting cold-start histories, the primary cut discards cycles 1–2,000 plus the seed, retaining cycles 2,001–10,000 (8,000 measurements). The κ = 0.80 energy and vortex relaxation motivated the common conservative cut; its remaining long correlation tail prevents a trustworthy error estimate. A visible initial plateau is not an equilibration proof.

The cumulant is U = 1 − mean(|m|⁴)/(2 mean(|m|²)²). Errors use 4,000 joint contiguous-block bootstrap draws, recomputing the ratio in each draw. Balanced blocks retain all 8,000 post-cut samples. The minimum width is max(20, ceil(10 × τ_max)), where τ_max covers |m|², |m|⁴, Binder influence, action density and vortex density. The τ estimator is the same nonperiodic FFT initial-positive monotone paired estimator used previously; independent samples have τ_int = 1/2.

| κ | Binder U | Bootstrap SE | τ_int(|m|²), cycles | τ_max, cycles | Blocks | Mean cluster / volume |
|---:|---:|---:|---:|---:|---:|---:|
| 0.80 | 0.058716 | not estimated | 17.8 | 654.7 | 1 | 0.0104 |
| 0.95 | 0.303013 | 0.009094 | 9.3 | 31.1 | 25 | 0.0978 |
| 1.05 | 0.478367 | 0.001095 | 5.5 | 10.9 | 72 | 0.2623 |
| 1.12 | 0.490977 | 0.000408 | 4.8 | 9.9 | 80 | 0.3271 |
| 1.20 | 0.494279 | 0.000243 | 5.1 | 9.6 | 83 | 0.3720 |
| 1.35 | 0.496460 | 0.000195 | 5.0 | 8.8 | 90 | 0.4201 |
| 1.50 | 0.497453 | 0.000112 | 5.4 | 12.1 | 65 | 0.4681 |

At κ = 0.80, only one block meets the conservative width rule, so the plotted hollow marker is a provisional central estimate without an error bar. Its τ estimates are poorly constrained by this run; particularly τ_energy ≈ 655 cycles is not a precise determination. A longer run is needed for this point under the chosen uncertainty protocol. The other points have 25–90 primary blocks.

For κ ≥ 0.95, changing the discarded prefix to 500, 1,000 or 4,000 cycles changes U by at most about 1.2 primary standard errors in these same-chain comparisons. Blocks doubled and quadrupled give broadly comparable errors, but the largest blocks leave as few as six blocks at κ = 0.95. These checks are not independent replicates and do not prove equilibration.

Near κ = 1.12, τ_int(|m|²) ≈ 4.84 cycles, with mean cluster membership about 32.7% of the lattice. Their product is about 1.58 lattice volumes of cluster membership. This work measure omits whole-field validation, random tapes, winding proposals, measurements and I/O. There is no same-L Metropolis comparison here, so this is not a measured CPU speedup.

## Artifacts

- [Magnified Binder comparison](binder_zoom.png)
- [Cold-start relaxation](thermalization.png) and [full histories](overview.png)
- [Cluster/autocorrelation diagnostics](cluster_diagnostics.png)
- [Cut and block sensitivity](stability.png)
- [L48 numerical results](results.json), [all five sizes](combined-results.json), and [cut/block checks](stability.json)
- [Run manifest](manifest.json), [timing and runtime provenance](status.json), [analysis provenance](analysis-provenance.json), and [field audit/checksums](data-audit.json)

From the repository root:

```sh
MPLBACKEND=Agg MPLCONFIGDIR=/private/tmp/luscher-mpl NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python experiments/2026-09-02-xy-binder-bkt-l48-cluster/run.py
MPLBACKEND=Agg MPLCONFIGDIR=/private/tmp/luscher-mpl NUMBA_NUM_THREADS=1 .venv/bin/python experiments/2026-09-02-xy-binder-bkt-l48-cluster/analysis.py --histories --cut 2001
MPLBACKEND=Agg MPLCONFIGDIR=/private/tmp/luscher-mpl NUMBA_NUM_THREADS=1 .venv/bin/python experiments/2026-09-02-xy-binder-bkt-l48-cluster/audit.py
```
