# Boundary refinement

Requested September 2, 2026: refine the kappa grid and increase simulation work until vertical grid/sampling errors are comparable to horizontal resolution.

The original complete scout is preserved in `../2026-09-02-xy-bkt-scout/`.

Completed September 3, 2026, 05:53 UTC: 1200 production/refinement trajectories at 100 parameter pairs and three volumes, plus 22 design pilots. All 1222 full-field audits and final record reconciliation passed. See [REPORT.md](REPORT.md), [phase_boundary_refined.png](phase_boundary_refined.png), [conservative_comparison.png](conservative_comparison.png), and [phase-brackets.csv](phase-brackets.csv).

- Pilot: 22 independent L64 trajectories, 10000 hybrid cycles, used only to choose neighborhoods.
- Production: directly simulate three kappa values at spacing 0.02 near each candidate boundary, with L=16,32,64 and four independent seeds per point.
- Each production trajectory has 12000 cycles. Each cycle has one Metropolis sweep, four fixed Wolff clusters, and one winding proposal. Thus each parameter/volume point has 48000 sweeps and 192000 clusters before the analysis cut.
- Starts have winding (0,0), (1,0), (0,1), (-1,0). The number of clusters is fixed, not a stopping rule based on the current field.
- Additional kappa=0 checks at epsilon/pi=0.54,0.565,0.58.
- The initial field and the full field after every hybrid cycle are retained, together with every measurement, RNG checkpoint, cluster size, work count and timing. Measurements use a fixed cycle cadence. Each completed chain gets an independent full-field physical-observable audit and SHA256 digest. Progress is checkpointed every 100 cycles, and job metadata every 1000 cycles.
- Four independent serial Numba workers per batch. Disjoint batches ran concurrently, using up to eight numerical workers; no GPU.

Run designs were frozen in `batches/*.json` before generation. Pilot choices are data-dependent, and pilots are excluded from production inference. The selected 2000-cycle cut is supported by the inspected cold/prepared-start histories and recorded cut/block checks; subsequent analysis retains all remaining rows with blocks at least ten times the largest measured integrated autocorrelation time. Independent starts and cut/block sensitivity are checked separately. Diagnostic flags remain in `diagnostic-summary.json`; the two cross-start flags occur at L16, which is a size/start check rather than a phase-classification input.

The old phase bars were conservative grid brackets using a heuristic finite-volume envelope of critical observables. Fine coupling spacing and longer runs can reduce grid/sampling uncertainty but cannot certify the infinite-volume transition. Report any finite-size/estimator spread separately and retain the original conservative diagnostic as a comparison. Do not call a fixed-Binder crossing a thermodynamic BKT determination.

Edge batches can run alongside production only through `run_edges.py`, which rejects overlapping chain IDs across the frozen manifests and budgets storage for all pending batches. Worker counts are recorded in each batch status file. Final history plots distinguish initial, full, and retained ranges in `histories_final/`.

Rebuild the complete analysis from the existing chains with `.venv/bin/python experiments/2026-09-02-xy-bkt-refine/finish_analysis.py`. It analyzes only closed files and does not start simulations. The main plot preserves the original scouting-bracket rule; the median vertical half-width is now 0.06, with a median reduction factor of 3.54 over the 20 comparable original cells. The remaining widths include finite-volume diagnostic ambiguity and the slope across each epsilon cell. Five conditional Binder-proxy interpolations remain unresolved under their stricter gates; they are not substituted for the main brackets.
