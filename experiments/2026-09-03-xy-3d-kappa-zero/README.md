# Three-dimensional admissible XY scout at zero coupling

This campaign measures the O(2) Binder cumulant at `kappa = 0` on periodic
`L^3` lattices.  The epsilon grid is 2.40, 2.45, 2.50, 2.55, and 2.60 radians;
the repeated 2.50 in the request is interpreted as 2.55.  Sizes are
`L = 8, 12, 16`.

Each parameter point has two independent 10,000-cycle trajectories, initialized
in winding sectors `(0,0,0)` and `(1,0,0)`.  A cycle consists of one Metropolis
sweep, four epsilon-aware Wolff clusters, and one global winding proposal.  The
runner checkpoints every 100 cycles and independently recomputes the saved
observables from the fields before marking a trajectory complete.

The Binder estimator is

$$
U_L = 1 - \frac{\langle |m|^4\rangle}
                 {2\langle |m|^2\rangle^2}.
$$

The analysis uses a visually checked equilibration cut, integrated
autocorrelation times, and a stratified contiguous-block bootstrap across the
two starts.  Blocks are at least ten times the largest monitored autocorrelation
time.

```bash
MPLBACKEND=Agg MPLCONFIGDIR=/private/tmp/luscher-mpl \
NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
.venv/bin/python experiments/2026-09-03-xy-3d-kappa-zero/design.py

MPLBACKEND=Agg MPLCONFIGDIR=/private/tmp/luscher-mpl \
NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
.venv/bin/python experiments/2026-09-03-xy-3d-kappa-zero/run.py scout --workers 4

MPLBACKEND=Agg MPLCONFIGDIR=/private/tmp/luscher-mpl \
.venv/bin/python experiments/2026-09-03-xy-3d-kappa-zero/analyze.py --histories --cut 2001
```
