# XY cluster validation — 2026-09-02

Implementation: `ClusterUpdate` in `supervillain/generator/xy.py`, compiled
growth in `supervillain/generator/_xy_cluster.py`, and independent full-graph
oracle in `supervillain/generator/reference_implementation/xy_cluster.py`.
The public API and scan examples are in [xy.md](xy.md).

## Target and balance argument

The target has action `S = -kappa * sum_positive_bonds cos(dphi)` and strict
principal-angle constraint `abs(dphi) < epsilon`. The cluster generator accepts
only XY actions with nonnegative finite kappa. Each positive-direction link is
a distinct edge, including both edges of a periodic two-site direction.

For a fixed spin-space reflection, let `w` be the current allowed pair weight
and `w'` its weight after reflecting one endpoint. The activation probability is
`p = 1 - min(1, w'/w)`. Thus an inadmissible single-end reflection forces a bond;
otherwise its exponent is `max(0, 2*kappa*(r.s_x)*(r.s_y))`. Simultaneously
reflecting both endpoints preserves their pair weight. The probabilities for
an unactivated cluster-boundary edge have forward/reverse ratio `w'/w`.
Multiplying those factors gives the Boltzmann weight ratio for the completed
cluster flip. The same seed sites generate the reverse cluster, so the seed
factor cancels. This is a dimension-independent balance argument in exact
arithmetic, with a uniformly selected seed and reflection normal.

Neither that argument nor the finite test suite proves exact detailed balance
on the finite set of floating-point angle values. Both implementations test
the full candidate with the strict cutoff before mutation and raise explicitly
if rounding violates it. Such exceptions stop a run; they are not silent
Metropolis rejections. Near activation thresholds the independent algebraic
forms may produce different rounding.

## Executed checks

```
NUMBA_NUM_THREADS=2 .venv/bin/python -m pytest \
  test/test_xy.py test/test_xy_cluster.py test/test_xy_scan.py -q
```

After both review fixes and the injected-write regression: **135 focused tests
pass**. This includes:

- Complete transition-matrix enumeration on five small fixed-reflection
  orbits, verifying normalization, detailed balance, and stationarity against
  independently computed full-action weights; no sampling-error tolerance.
- Independent full-graph oracle versus compiled growth in dimensions 1–4,
  odd extents, two-site parallel edges, unconstrained and constrained cases.
- Zero-coupling mandatory bonds, strict equality at both pi/2 and pi, very
  small epsilon, large finite coupling, invalid inputs, nonmutation, and
  explicit numerical failure before field/counter mutation.
- A prescribed nonlocal cluster that changes torus winding from one to zero
  at epsilon=pi/2, without introducing forbidden vortices.
- Exact seeded backend replay; generator HDF5 round-trip and continuation;
  scan interruption/replay with both backends and optional winding proposals;
  all per-cluster work records and cumulative counters; changed-schedule
  rejection; CLI and pre-file parameter-validation checks.

A repository-wide `pytest -q` run was deliberately interrupted in the lengthy
legacy observable tests after **144 passed in 145.78 seconds**. This is not a
full-suite pass. Retained evidence:
[focused test log](validation-logs/xy-cluster-2026-09-02/focused-tests.txt) and
[interrupted full-suite log](validation-logs/xy-cluster-2026-09-02/full-tests-interrupted.txt).
[Runtime versions and final simulation-source fingerprint](validation-logs/xy-cluster-2026-09-02/runtime-provenance.json)
are also retained.
No production phase scan was run.

## Independent adversarial reviews

The first fresh GPT5.6 Sol review completed with no remaining material findings
after the fix below. Its independent script constructs neighbors from coordinate
tuples and reflections from complex spins. Five nonzero-angle, branch-crossing
orbits enumerate 800 graphs and 3,144 seeded moves; the maximum weighted-flux
asymmetry was 4.44e-16. The primary agent reran the retained
[script](validation-logs/xy-cluster-2026-09-02/independent-orbits.py);
[results](validation-logs/xy-cluster-2026-09-02/independent-orbits.txt) agree.
The second fresh GPT5.6 Sol review found no remaining cluster-specific blocker
after the dtype fix below. It independently rechecked all eight dtype/backend
combinations and enumerated three nonzero-angle, branch-wrapping three-spin
orbits; the maximum weighted-flux error was 1.11e-16.
It independently injected a genuine mid-row HDF5 write failure, leaving unequal
dataset lengths, and verified exact replay of every dataset after `run_chain`
resume. Complete generator/ensemble HDF5 snapshots also restored correctly.
The injected-write case is now a permanent scan test. The review's integration
[script](validation-logs/xy-cluster-2026-09-02/integration-review.py) and
[output](validation-logs/xy-cluster-2026-09-02/integration-review.txt) are retained.

The first reviewer found a Python-oracle cancellation error at kappa=1e18,
reflection normal angle `pi/2 + 1e-9`, and a cold three-site chain: subtracting
two nearly equal cosines erased an order-one activation exponent. The oracle
now uses the stable two-sine identity for their difference. A regression checks
both backends with the supplied tape; the focused suite passes after the fix.

The second reviewer identified an existing XY action input-type bug that also
affected cluster validation: unsigned neighbor subtraction underflowed before
conversion to float. For `uint64([0, 1, 0])`, kappa=1 and epsilon=1.2, the action
was -2.4753024434 instead of -2.0806046117. `XY.links` now converts to float64
before subtracting. Eight regression cases cover nonuniform signed/unsigned
integer fields, correct links and action, unchanged input, and agreement with
float64 cluster inputs in both backends. Villain code is unchanged.

Known inherited limitation found in that review: generic `Ensemble.extend_h5`
appends configurations without refreshing the saved generator RNG/counters
(`supervillain/h5/extendable.py`). Reloading an appended generic ensemble is
therefore not an exact checkpoint. That unrelated generic API was not changed;
the guide now explicitly directs streaming runs to `run_chain`, and bounds
generic HDF5 continuation claims to complete `Ensemble.to_h5` snapshots.

## Scope and interpretation

One generator step is one cluster. The runner records after a fixed number of
clusters, optionally followed by a winding proposal. It does not stop recording
cycles when accumulated membership reaches the volume. All cluster sizes and
growth bond-test counts are retained, along with cumulative membership/volume.
These are work diagnostics, not replacements for the fixed observation clock.

Cluster growth is compiled but serial. There is no new parallel cluster kernel
or GPU implementation. The oracle builds the entire graph; even the compiled
backend uses volume-sized random tapes and validates the full field. Work
counts exclude these overheads. No speedup or improved autocorrelation time is
claimed by this validation. Tight constraints can produce global reflections
that leave energy and magnetization magnitude unchanged, so mixing must be
measured separately for the physics regime of interest.

The supplied Hasenbusch--Meyer paper (Phys. Lett. B 241, 238–242, 1990;
DOI [10.1016/0370-2693(90)91286-K](https://doi.org/10.1016/0370-2693(90)91286-K))
was inspected in text and rendered pages 238–239. It describes an embedded-Ising
multicluster algorithm. The new code uses the related single-cluster Wolff
construction with mandatory constraint bonds; its performance is not inferred
from that paper's dynamic exponent.
