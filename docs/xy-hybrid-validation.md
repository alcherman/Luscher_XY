# XY hybrid cycles and timing — 2026-09-02

`run_chain(..., update="hybrid")` now performs a full Metropolis sweep,
then a fixed `clusters_per_step` Wolff batch, then the optional winding proposal.
Each component has an independent reproducible random stream. The same action
and strict admissibility constraint are used throughout; the cluster component
requires kappa >= 0. Fixed composition of invariant kernels preserves the target
distribution and needs no final combined acceptance test.

The HDF5 record retains local acceptance, every cluster size and bond-test count,
cumulative cluster work, all physical fields/observables, and all RNG states.
Optional `record_timing=True` saves per-cycle wall times for Metropolis, cluster,
winding, their combined update, measurement, and the cycle with record writes.
Cycle timings exclude periodic checkpoint flushes; end-to-end experiment times
must be recorded separately. Timings are performance observations, so replay
reproduces the physical trajectory but not the elapsed seconds.

```
NUMBA_NUM_THREADS=2 MPLCONFIGDIR=/private/tmp/luscher-mpl .venv/bin/python -m pytest \
  test/test_xy.py test/test_xy_cluster.py test/test_xy_scan.py -q
```

Result: **141 passed in 5.34 seconds**. The new tests compare every hybrid field,
acceptance count, cluster size, cumulative work count, and RNG state to an explicit
Metropolis/cluster/winding composition for both backends and a nontrivial cutoff.
They verify uninterrupted-versus-resumed trajectories, timing independence for
all three methods, CLI creation, and negative-coupling rejection before file
creation. Existing cluster balance, strict-boundary, dtype, and interrupted-write
tests remain in the focused suite.

The experiment in `experiments/2026-09-02-xy-hybrid-k080` compares Metropolis,
one-sweep/one-cluster hybrid cycles, and batches of 100 clusters at D=2, L=48,
kappa=.80, epsilon=pi. These are fixed proposed schedules, not a tuning search.
Its raw timing records, full data, and analysis provide the performance evidence;
unit tests establish no performance or equilibration claim.
