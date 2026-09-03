# Admissible XY simulations

The XY implementation uses the existing periodic hypercubic `Lattice(D, N)`:
any positive integer dimension, one common extent `N >= 2` in every direction.
The action and measure are

$$
S=-\kappa\sum_{x,\mu}\cos\theta_\mu(x),
$$

$$
\theta_\mu(x)=\operatorname{Arg}\left(e^{i[\varphi(x+\hat\mu)-\varphi(x)]}\right),
$$

$$
\prod_x d\varphi_x\;e^{-S}\prod_{x,\mu}\mathbf{1}_{|\theta_\mu(x)|<\epsilon}.
$$

Each positive-direction bond appears once. In particular, the two bonds of a
periodic two-site direction are distinct. Relative to the convention
`S = -beta * sum(cos(theta))`, **beta = kappa**. The small-gradient expansion is

$$
S=-\kappa D N^D+\frac{\kappa}{2}\sum_{x,\mu}\theta_\mu(x)^2
-\frac{\kappa}{24}\sum_{x,\mu}\theta_\mu(x)^4+\cdots.
$$

The code accepts finite real
`kappa`, including zero and negative values, and `0 < epsilon <= pi` in radians.
Zero epsilon has no admissible configurations under a strict inequality.
Spins are stored as float64 angles in `[0, 2*pi)`. `action.canonical(phi)` converts
other finite angle representatives to that interval. Action evaluation is
periodic in those representatives; `valid(configuration)` additionally checks
the stored angle interval. Bounds are strict without a tolerance; exactly
antipodal bonds are rejected even at `epsilon = pi`.

## Python API

```python
import numpy as np
from supervillain import Lattice, XY, Ensemble
from supervillain.generator.xy import SiteUpdate, WindingUpdate
from supervillain.generator.combining import Sequentially
from supervillain.xy import measure

action = XY(Lattice(D=3, N=16), kappa=1.0, epsilon=1.6)
local = SiteUpdate(action, seed=123, backend="numba", parallel=True)
generator = Sequentially((local, WindingUpdate(action, seed=456)))
ensemble = Ensemble(action).generate(100, generator, start="cold")
ensemble.measure(["ActionDensity", "Links", "Spin_Spin"])
diagnostics = measure(action, ensemble.configuration[-1]["phi"])
```

`backend="python"` selects the intentionally slow oracle: it independently
recomputes the entire action for every trial spin, in explicit Python loops.
`backend="numba", parallel=False` selects compiled serial sweeps.
The parallel backend uses a serial sequence of independent site colors and
`prange` within each color. Odd extents use more than two colors. Incoming and
outgoing bonds are both checked before accepting any trial. No `fastmath` is used.
Proposals are symmetric uniform offsets with default width `epsilon`, followed
by Metropolis acceptance. The width stays fixed throughout a run.

All random proposals and acceptance uniforms are drawn outside parallel kernels.
`local.sweep(phi, proposals, uniforms)` replays arrays indexed by flat site;
it returns a fresh field. The same random tape is used to test the oracle,
compiled serial, and compiled parallel algorithms. Floating-point action
subtraction in the oracle can differ from local differences at a Metropolis
threshold; agreement tests avoid treating roundoff-scale differences as physics.
Numba's rules about independent `prange` iterations are documented in its
[parallel execution guide](https://numba.readthedocs.io/en/stable/user/parallel.html).

## Admissibility-aware Wolff clusters

```python
from supervillain.generator.xy import ClusterUpdate

action = XY(Lattice(D=3, N=16), kappa=0.5, epsilon=1.6)
cluster = ClusterUpdate(action, seed=123, backend="numba")
ensemble = Ensemble(action).generate(100, cluster, start="cold")
print(cluster.report())
```

`ClusterUpdate` accepts **XY actions with kappa >= 0 only**. Villain actions are
rejected. Both backends perform exactly one Wolff cluster per `step`; spatial
dimension does not change the two-component spin reflection. Choose a uniformly
random seed site and a random normal vector
$r=(\cos\alpha,\sin\alpha)$, and reflect spins by
$R s=s-2(r\cdot s)r$, or $\varphi'=2\alpha+\pi-\varphi$ modulo $2\pi$.
For a bond joining a cluster site to an unvisited neighbor, set

$$
p_{xy}=\begin{cases}
1, & |\operatorname{Arg}(e^{i(\varphi_y-\varphi'_x)})|\ge\epsilon,\\
1-\exp\!\left[-\max\!\left(0,2\kappa(r\cdot s_x)(r\cdot s_y)\right)\right],
& \text{otherwise}.
\end{cases}
$$

All bond decisions use the original configuration. An activated bond adds the
neighbor to the cluster. After growth finishes, reflect every member. Internal
bonds keep their cosine and absolute angle difference; forbidden boundary bonds
have been activated with certainty. The move therefore needs no final Metropolis
test. The bond rule is $p=1-\min(1,w'/w)$ for the combined cosine and hard-cutoff
pair weight. Its usual boundary-bond balance argument is independent of spatial
dimension. At zero coupling only mandatory constraint bonds remain.

The supplied Hasenbusch--Meyer paper,
[Critical exponents of the 3D XY model from cluster update Monte Carlo](https://doi.org/10.1016/0370-2693(90)91286-K),
uses an embedded-Ising **multicluster** update. This implementation uses the
related **single-cluster Wolff** construction and adds the hard-cutoff bond rule;
the paper's measured dynamic exponent is not a performance prediction for it.

`backend="python"` constructs the entire random bond graph independently,
evaluating cosine energy differences directly. `backend="numba"` uses compiled
queue-based growth and spin projections. Cluster growth is serial in both;
independent chains can be distributed across workers. No GPU support or parallel
cluster kernel is supplied. Validation, random-tape generation, and the final
strict constraint check do work proportional to the whole lattice even for a
small cluster.

`cluster.cluster(phi, seed_site, reflection_angle, uniforms)` replays a supplied
random tape and returns a fresh field. The normal angle is `reflection_angle`;
`uniforms` has shape `(D, N**D)`, with one draw per positive bond indexed by its
direction and flat origin site. Both distinct bonds in a two-site periodic
direction are retained. No `fastmath` is used. Near the cutoff, both single-end
reflections are checked. If rounding makes the completed reflected field
inadmissible, `FloatingPointError` stops the update before field or counter
mutation. No tolerance relaxes the constraint and no silent rejection is added.
These guards do not establish exact detailed balance on floating-point states.

The generator records `last_cluster_size`, `last_bond_tests`, and cumulative
`clusters`, `flipped_sites`, and `bond_tests`. `equivalent_sweeps` is accumulated
cluster membership divided by volume, a work diagnostic. Use a **fixed number
of clusters** between measurements: stopping whenever accumulated cluster size
reaches the volume introduces a state-dependent observation schedule. Compare
autocorrelations per fixed cycle alongside mean work and wall time.

Complete `Ensemble.to_h5` snapshots preserve the generator's RNG and counters
for continuation. The existing generic `Ensemble.extend_h5` path appends rows
without refreshing the saved generator state; reloading that appended file is
not an exact checkpoint. Use `run_chain` below for streaming, resumable cluster
runs, or write a complete new ensemble snapshot.

Clusters can change global winding even below the vortex threshold. Very tight
constraints can instead produce nearly whole-lattice reflections, which leave
energy and magnetization magnitude unchanged. Check histories, sector mixing,
and autocorrelations; cluster acceptance alone is not an efficiency diagnostic.

## Parameter scans and restart

A local `.venv` contains the runtime dependencies. From the repository root:

```sh
.venv/bin/python -m supervillain.xy --D 2 --N 16 32 \
  --kappa 0 1 2 --epsilon 1.4 1.5707963267948966 1.8 3.141592653589793 \
  --replicas 4 --sweeps 10000 --warmup 2000 --seed 12345 \
  --threads 4 --checkpoint-every 10 --output runs/xy-pilot
```

This is an example launch command, not a validated choice of volume, run length,
or thermalization cut. Grid points and replicas run sequentially; each local
sweep uses the requested CPU threads. Use `--serial` or `--backend python` to
select the alternatives. `--interval-phi` changes the fixed proposal width.
The first Numba call compiles; benchmark warmed kernels separately from startup.
Parallel thread launch overhead can dominate small lattices.

For cluster runs, add `--update cluster --clusters-per-step 1` and use a fresh
output directory. `--sweeps` then counts saved cycles, each containing the fixed
number of cluster flips followed by an optional global winding proposal. The
same choice is available as `run_chain(..., update="cluster", clusters_per_step=1)`.
`--interval-phi` is rejected for clusters. CPU thread settings affect Metropolis
only; cluster metadata records `parallel=false` and `time_unit="cluster_cycle"`.
Cluster files additionally retain every cluster size and growth bond-test count
within each cycle, their cumulative totals, and accumulated membership/volume.
The seed row has zero work. These counts do not include oracle whole-graph work,
validation, RNG generation, or winding proposals. The checkpoint format is now
`supervillain.xy.v2`; historical files remain untouched and are not resumable
under the new source fingerprint.

`--update hybrid --clusters-per-step 1` performs one full Metropolis sweep,
then the fixed number of Wolff clusters, then the optional winding proposal.
It uses independent RNG streams for all three generators, restores them on
resume, and retains both local acceptance and cluster work records. The
Metropolis part accepts `--interval-phi` and `--serial` as usual. The hybrid
mode requires nonnegative kappa, as does its cluster component.

`--record-timing` additionally saves wall-clock seconds per cycle for the
Metropolis, cluster and winding parts, their combined update, measurement,
and the cycle including record writes. `cycle_seconds` excludes periodic
checkpoint flushes; total run wall time should also be retained for performance
comparisons. Timings do not affect random draws and are not bitwise replayable
on resume. Model fields, measurements, RNG states and work counts still are.

`--buffer-rows 100` batches HDF5 writes while retaining every generated row,
including every field and RNG state. Each checkpoint flushes the pending rows
before advancing its committed marker. The default is one row; buffered writes
and checkpoint flushes are excluded from per-cycle timing, so retain total wall
time for throughput comparisons. Buffering changes neither random draws nor
the saved scientific data, and interrupted batches replay from the last commit.

Every chain has its own HDF5 file and child seed. Files retain the cold or twisted
initial state at row 0, all warmup configurations, every subsequent configuration,
all measured values, acceptance counts, and both RNG states. Parameters, runtime
versions, thread count, Git revision/dirty status, and a source digest are stored.
`sweeps` is the total target number of updates including warmup; there are
`sweeps + 1` saved rows. The suggested analysis interval starts at `warmup + 1`.

Repeat the command with `--resume` to restore each existing chain's last committed
RNG checkpoint, finish interrupted points, and create points not yet started.
Increasing `--sweeps` extends a chain. Keep the grid order, root seed, replica
count, and all other parameters identical. Mismatched parameters or source/runtime
versions fail closed. The HDF5 writer flushes data before advancing the committed
row count; incomplete trailing rows are recomputed from the saved RNG state.
This is an application checkpoint, not a guarantee against filesystem corruption.
No existing file is overwritten without an explicit resume request.

For a single chain use `supervillain.xy.run_chain(path, D=..., N=..., kappa=...,
epsilon=..., ...)` directly. `--initial-winding 1 0` (in two dimensions) starts
from a smooth unit twist if admissible; use a separate output directory for each
choice of initial sector. `--no-winding` disables the extra global twist proposal;
it does not force fixed winding, since cluster moves can also change the sector.

## Vortices, winding, and observables

Plaquette vorticity is the oriented circulation of principal differences divided
by `2*pi`. Four strictly bounded edges imply `|circulation| < 4*epsilon`.
Consequently **elementary vortices are forbidden for epsilon <= pi/2**, including
equality, and are kinematically allowed above that threshold. This statement is
about allowed configurations, not their equilibrium density or where a phase
transition lies. In dimensions above two, the dual vortex defects extend in
`D-2` dimensions. In one dimension there are no plaquettes.

Smooth fields can have nonzero winding around the torus without any local
vorticity. A local chain can fail to move between winding sectors at tight
admissibility. The additional global proposal chooses an axis and sign uniformly,
adds `+/- 2*pi*x_axis/N`, and applies the full constraint and Metropolis test.
It obeys the same target measure and can cross sectors in one move. This proposal
does not prove ergodicity or ensure adequate sector mixing; monitor its acceptance
and compare independently initialized sectors. At very tight bounds, some sectors
cannot exist even as smooth twists.

`measure(action, phi)` and the scan files provide:

- `mean_link_cosine` and `action_density = -kappa*D*mean_link_cosine`.
- Complex magnetization components, magnitude, `m2`, and `m4` (dataset names
  `magnetization_squared` and `magnetization_fourth`).
- `structure_zero = V*m2` and `structure_min` for each axial minimal momentum.
- `vortex_density`: mean absolute integer plaquette charge, averaged over all
  orientations and sites (zero by convention in D=1).
- `winding`: an integer for every straight periodic loop, shape
  `(D, N**(D-1))`. With vortices, different parallel loops need not agree.
- Maximum absolute link angle and per-sweep local/global acceptance counts.

Useful dimensionless finite-size observables formed from **ensemble means** are
the O(2) Binder cumulant `U = 1 - <m4>/(2*<m2>**2)` and `xi/N`, where
`xi = sqrt(<structure_zero>/<structure_min[a]> - 1)/(2*sin(pi/N))` for axis `a`.
The spin susceptibility in the unconnected convention is `V*<m2>`.
These estimators are intended for the ferromagnetic, small-momentum regime;
negative kappa may require observables at the ordering wavevector instead.
Crossings, drifts and scaling across N are needed to locate a transition; the
vortex admissibility threshold alone does not locate it.

The usual cosine-action helicity estimator is deliberately not supplied: under
admissibility the domain also changes with an imposed twist, so differentiating
only the smooth cosine weight omits contributions from the constraint boundary.

Before quoting errors or a transition curve, inspect histories including the
initial field, choose a stable thermalized interval, estimate retained-data
autocorrelation times, and bootstrap common block means (start around
`10*tau_int`). Recompute ratios within each bootstrap draw. Check cut and block
size stability, independent seeds, and sector occupancy. A constant vorticity
series below the threshold is an exact constraint and gives no mixing diagnostic.
The saved raw data support these analyses; this initial infrastructure does not
automatically certify equilibration or produce a phase boundary.

## Validation

See [the cluster validation record](xy-cluster-validation.md),
[the original validation record](xy-validation.md), and [benchmark repeats](xy-benchmark.json)
for measured results and the bounds of the checks actually performed.

```sh
NUMBA_NUM_THREADS=4 .venv/bin/python -m pytest test/test_xy.py test/test_xy_cluster.py test/test_xy_scan.py -q
NUMBA_NUM_THREADS=4 .venv/bin/python -m pytest -q
```

Tests exercise normalization, compact angles and strict boundaries, multiple
dimensions, two-site bonds, odd periodic coloring, shared-tape oracle parity,
thread-count reproducibility, HDF5 continuation, and analytic small-system checks.
The scan tests also exercise interrupted checkpoints and vortex/winding diagnostics.
Only bounded verification simulations have been run during implementation.
