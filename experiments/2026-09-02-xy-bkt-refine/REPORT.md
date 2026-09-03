# XY boundary refinement results

Completed 1222 trajectories (1200 production/refinement, 22 pilot). All passed their full-field audits. Saved 14,621,222 fields in 236.22 GB; no generated rows were thinned or deleted.

## Calculation

Validation: five focused analysis tests passed; the classifier reproduced all 196 original scout phase labels from the original inputs. Every new trajectory passed its independent full-field observable/admissibility audit.
D=2, L=16,32,64; original epsilon grid retained. Three pilot-selected kappa values separated by 0.02 per positive-boundary slice, followed by separately frozen additions near unresolved bracket edges. See batches/*.json for every exact parameter and seed.
Each production parameter/volume group has four independent trajectories of 12000 cycles, starting with windings (0,0), (1,0), (0,1), (-1,0). Each cycle performs one Metropolis sweep, four Wolff clusters, and one winding proposal: 48000 sweeps and 192000 clusters per parameter/volume group.
Analysis retains cycles 2001–12000 from each start, giving 40000 retained cycles per group. HDF5 warmup=100 is a provisional runner label; the selected analysis cut is recorded in analysis/*.json. Raw seed and transient records remain intact. Block sizes are at least ten times the largest monitored tau, with 2000 independent-within-start block resamples. Raw moments are pooled before computing Binder and correlation-length ratios. Cut/block stability, disjoint halves and cross-start agreement are recorded.

## What the error bars mean

The main phase_boundary_refined figure preserves the original joint Binder/correlation-length/effective-exponent rule, with two-standard-error margins and adjacent epsilon-cell unions. The decision uses L32 and L64 and their effective exponent; L16 remains a size/start check. Its bars remain heuristic finite-volume brackets, not Gaussian confidence intervals or a validated thermodynamic extrapolation. The slope across an epsilon cell contributes to their vertical width.
precision_and_finite_size separates the attained grid/Monte Carlo resolution of a specified L64 Binder reference midpoint from the larger joint-diagnostic bracket. It is a conditional finite-volume proxy and must not be quoted as the BKT critical coupling with only its small statistical error. Finite-size and reference-choice bias has not been certified. The fixed-Binder limitation is material; see https://arxiv.org/pdf/1307.0485 section 3.1.
Conditional local interpolation is resolved for 17 of 22 epsilon slices under the recorded gates; unresolved cases remain explicit in boundary.json.

| epsilon/pi cell center | old vertical half-width | refined half-width | reduction |
|---:|---:|---:|---:|
| 0.610 | 0.300 | 0.110 | 2.73x |
| 0.630 | 0.300 | 0.090 | 3.33x |
| 0.650 | 0.500 | 0.090 | 5.56x |
| 0.670 | 0.300 | 0.080 | 3.75x |
| 0.690 | 0.200 | 0.080 | 2.50x |
| 0.710 | 0.100 | 0.080 | 1.25x |
| 0.730 | 0.300 | 0.060 | 5.00x |
| 0.750 | 0.300 | 0.070 | 4.29x |
| 0.770 | 0.200 | 0.060 | 3.33x |
| 0.790 | 0.100 | 0.050 | 2.00x |
| 0.810 | 0.100 | 0.050 | 2.00x |
| 0.830 | 0.200 | 0.060 | 3.33x |
| 0.850 | 0.200 | 0.050 | 4.00x |
| 0.870 | 0.200 | 0.050 | 4.00x |
| 0.890 | 0.200 | 0.040 | 5.00x |
| 0.910 | 0.200 | 0.050 | 4.00x |
| 0.930 | 0.100 | 0.050 | 2.00x |
| 0.950 | 0.200 | 0.040 | 5.00x |
| 0.970 | 0.200 | 0.040 | 5.00x |
| 0.990 | 0.100 | 0.040 | 2.50x |

## Diagnostics and limits

21 individual trajectories have at least one cut/half diagnostic flag; 2 parameter/volume groups have a cross-start p<0.001 flag. These are retained in diagnostic-summary.json, not silently discarded. Multiple comparisons mean isolated flags are not themselves proof of disequilibrium; no exhaustive equilibration or sector-completeness claim is made.

| L | median tau(m2) | max tau(m2) | median tau(Binder influence) | max tau(Binder influence) |
|---:|---:|---:|---:|---:|
| 16 | 0.95 | 1.28 | 0.71 | 0.90 |
| 32 | 1.07 | 1.37 | 0.79 | 1.09 |
| 64 | 1.19 | 1.54 | 0.86 | 1.08 |

Tau is measured in the new hybrid cycle (four clusters), not in the old one-cluster time unit. All-observable time histories, per-start statistics and bootstrap samples are stored under histories_final/, analysis/ and pooled/. Histories include the seed, transient and retained-data zoom with separate vertical scales. Their scope is finite-volume scouting; no infinite-volume BKT fit was performed.

## Reproduction

Run pilot/production batches with `.venv/bin/python experiments/2026-09-02-xy-bkt-refine/run.py BATCH`; run edge batches with `run_edges.py BATCH` in the same directory. The edge runner verifies disjoint chain IDs before concurrent execution. Completed jobs are skipped and interrupted jobs resume at their committed checkpoint. The original scout and its model code are unchanged. model-source.tar.gz and the immutable batch manifests record simulation provenance. Analysis scripts, tests, logs, final table and figures are retained in this directory.
