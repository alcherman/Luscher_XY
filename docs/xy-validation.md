# XY implementation validation — 2026-09-02

This is infrastructure validation, not a phase-diagram measurement. The local
runtime is Python 3.14.6, NumPy 2.5.2, Numba 0.67.0, and pytest 9.1.1.

## Correctness

The new test files are `test/test_xy.py` and `test/test_xy_scan.py`:
**47 passed** after correcting the convention to `S = -kappa*sum(cos(theta))`.
They check:

- `S_cold = -kappa*D*N**D`, the small-gradient quadratic coefficient `kappa/2`,
  the measured action density, independent explicit-loop action comparisons,
  compact-angle invariance, strict bounds including the pi/2 vortex threshold,
  and tiny link differences that an unstable wrapping formula could erase.
- Multi-dimensional local updates, odd periodic coloring, the two distinct
  bonds at N=2, float conversion, and noncontiguous input arrays.
- Identical final fields and acceptance counts under shared proposal tapes in
  the Python oracle, serial Numba and parallel Numba; parallel replay with
  different CPU thread counts.
- A two-site Metropolis acceptance threshold calculated from the full action.
- The two-site equilibrium moment of cos(theta), whose density is proportional
  to `exp(2*kappa*cos(theta))` on `(-epsilon,epsilon)`. At kappa=0 this moment is
  `sin(epsilon)/epsilon`; at kappa=1.7 it is checked against quadrature. Each
  fixed-seed test uses 2,000 warmup and 10,000 retained sweeps, with 200-sweep
  blocks and a six-block-standard-error gate. This is a bounded statistical
  regression, not proof of mixing on large constrained lattices.
- Configuration copying, HDF5 ensemble continuation, exact checkpoint/restart,
  recovery of uncommitted trailing rows, mismatch rejection, retained warmup,
  winding/vorticity measurements, and a command-line smoke test.

Run these with:

```sh
MPLCONFIGDIR=/private/tmp/luscher-mpl NUMBA_CACHE_DIR=/private/tmp/luscher-numba \
  NUMBA_NUM_THREADS=4 .venv/bin/python -m pytest test/test_xy.py test/test_xy_scan.py -q
```

Before the normalization correction, the selected existing-component regression
command passed **1,938 tests** (not rerun for the correction):

```sh
MPLCONFIGDIR=/private/tmp/luscher-mpl NUMBA_CACHE_DIR=/private/tmp/luscher-numba \
  NUMBA_NUM_THREADS=4 .venv/bin/python -m pytest \
  test/test_batch.py test/test_batch_h5.py test/test_batch_dtype.py \
  test/test_lattice_tables.py test/test_lattice_kernels.py test/test_lattice.py \
  test/test_autocorrelation_fallback.py test/test_autocorrelation_time.py \
  test/test_origin.py test/test_streaming.py -q
```

The full legacy suite was started and interrupted after **274 passes** in about
three minutes, during its repeated Villain/Worldline simulation tests. It was
not completed. The selected regression above exercises the shared lattice,
storage, and analysis components. These pass counts overlap and should not be
added together as a count of distinct tests.

A separate four-point CLI smoke scan before the correction exercised kappa=0,1 and epsilon=1.4,1.8
on a 3x3 lattice, followed by extension from 6 to 8 sweeps. All four output files
resumed successfully. The existing Ensemble API also measured ActionDensity,
Links and Spin_Spin on a generated XY ensemble.

## Performance

`benchmarks/xy.py` writes every completed repeat to `docs/xy-benchmark.json` by
atomic replacement, retaining timing, load and parity records. The benchmark
warms the implementations before timing and includes the public replay wrapper's
validation/copy overhead; random tape generation, HDF5 output, global winding
proposals and measurements are outside the timed region.

The retained benchmark JSON records the original `-kappa/2` convention and its
source digests; it has not been regenerated after the normalization correction.
That local four-thread run measured about **5.59 ms per serial sweep**
and **2.66 ms per parallel sweep** on a 32^3 lattice, a factor of about **2.1**.
All repeated final-state digests agreed. These are fixed-tape microbenchmarks on
a shared machine, not throughput or statistical-efficiency guarantees for a
full production scan. Compiled serial sweeps were faster than parallel sweeps
on the two smaller benchmark cases. Exact timings and provenance are in the JSON.

## Remaining scientific scope

No large-volume equilibration, sector-completeness, critical-slowing-down study,
finite-size scaling, or phase boundary has been established. Global twist moves
supplement local updates; their existence does not certify ergodicity. GPU
acceleration is outside this implementation.
