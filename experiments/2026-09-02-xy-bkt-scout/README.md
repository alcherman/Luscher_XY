# Admissible XY scouting scan

Authorized grid: epsilon/pi = 0.50, 0.51, 0.52, 0.53, 0.54, 0.565, 0.58, and 0.60 through 1.00 in steps of 0.02; kappa = 0 through 1.2 in steps of 0.2; periodic D=2, L=16,32,64. There are 588 independently seeded cold-start chains, each with 10,000 hybrid cycles plus its initial field. One cycle is one Metropolis sweep, one Wolff cluster, and one winding proposal. The Metropolis proposal half-width stays at pi.

The immutable [manifest](manifest.json) contains exact parameters, seeds and source hashes. [status.json](status.json) records overall progress; `jobs/` holds atomic per-chain progress, segment timings, final checksums and full-field audit results. All stdout/stderr is retained in `run.log` and `logs/`. No row or physical column is thinned or discarded from `chains/`.

`run.py` executes four independent single-thread Numba jobs concurrently, checkpoints each chain every 100 cycles and its progress every 1,000 cycles. HDF5 writes are buffered in batches of 100, including every field, observable, RNG state and update count. A full independent recomputation of all physical observables from all 10,001 saved fields must pass before a job is called complete. Resume restores the last committed random state and replays any uncommitted suffix. The source snapshot preserves all Python files in the simulation package, the focused tests, runner and manifest.

`analyze.py` reads completed files only. Histories precede selection of the analysis cut. The scouting analysis retains cycles 2,001-10,000; the chain's original warmup attribute is only a provisional label, and does not override this cut. It bootstraps joint contiguous blocks, retaining every post-cut row, and recalculates Binder and correlation-length ratios on each resample. All sampling errors are conditional single-chain errors, without finite-size or initialization uncertainty. Cut, block, and disjoint-half diagnostics remain separate from any equilibrium claim.

`phase.py` produces conservative **finite-volume scouting brackets**. Its explicit heuristic compares Binder and xi/L against an envelope of published critical finite-volume observables, with an effective-exponent cross-check. This is not the precision reference-matching analysis proposed for a later refinement stage. In particular, it does not force either endpoint's parameter value, fit a thermodynamic transition, or treat grid brackets as one-sigma confidence intervals. The final report describes the assumptions and unresolved points.

From the repository root:

```sh
# Execute or resume the exact immutable scouting grid.
.venv/bin/python experiments/2026-09-02-xy-bkt-scout/run.py --workers 4

# Disk-only diagnostic and analysis commands.
.venv/bin/python experiments/2026-09-02-xy-bkt-scout/analyze.py --histories
.venv/bin/python experiments/2026-09-02-xy-bkt-scout/analyze.py --cut 2001 --plots
.venv/bin/python experiments/2026-09-02-xy-bkt-scout/plot_threshold.py
.venv/bin/python experiments/2026-09-02-xy-bkt-scout/phase.py
```

The complete implementation test log is [implementation-tests.txt](implementation-tests.txt), and nonlinear estimator checks are in [analysis-tests.txt](analysis-tests.txt). `simulation-source.tar.gz` and `simulation-source-files.json` preserve the executed source independently of the uncommitted checkout.
