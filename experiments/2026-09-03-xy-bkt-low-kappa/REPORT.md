# Low-kappa refinement results

Completed 84 trajectories; every full-field audit and the final reconciliation passed. Saved 2,688,084 fields (42.69 GB), with all measurements, RNG records, update work, timing and initial/transient rows retained.

The first batch adds epsilon/pi=0.575 and 0.590 and direct kappa values 0.10 and 0.18 at epsilon/pi=0.580. A follow-up at (epsilon/pi,kappa)=(0.590,0.16) established the missing massive-side edge. Each parameter/volume group has four independent 32000-cycle trajectories, or 128000 Metropolis sweeps and 512000 clusters. The selected 2000-cycle cut leaves 120000 pooled cycles, three times the previous refined count. Inspected initial/retained histories, cut/block checks and independent-start comparisons support this cut; no exhaustive equilibration or thermodynamic BKT claim is made.

Five analysis tests passed, and the classifier reproduced all 286 previous refined phase labels from unchanged previous inputs. The model package and previous campaign inputs are unchanged.

The full curve combines the new results with the unchanged earlier campaigns. At repeated parameter/volume coordinates, the new longer estimate replaces the earlier estimate; earlier raw records remain intact. The same joint U/R/eta two-SE classification and adjacent-epsilon-cell union define the bars. They are finite-volume scouting brackets, not Gaussian confidence intervals.

| epsilon/pi cell | kappa bracket | vertical half-width | target <= 0.05 |
|---|---|---:|---|
| 0.575–0.580 | 0.040–0.200 | 0.080 | not yet |
| 0.580–0.590 | 0.100–0.260 | 0.080 | not yet |
| 0.590–0.600 | 0.160–0.300 | 0.070 | not yet |
| 0.600–0.620 | 0.220–0.440 | 0.110 | not yet |
| 0.620–0.640 | 0.340–0.520 | 0.090 | not yet |
| 0.640–0.660 | 0.420–0.600 | 0.090 | not yet |

Diagnostic flags: 1 individual trajectories and 0 pooled groups. See diagnostic-summary.json and the individual analysis records. Flags are preserved, not silently removed.

The current disk budget restricts further refinement. Where the table misses the target, additional targeted grid points/statistics or a better controlled finite-volume estimator remain necessary. No smaller total error is claimed from an auxiliary fixed-Binder crossing.

Artifacts: phase_boundary_refined.png/pdf/svg, low_kappa_comparison.png/pdf/svg, phase-brackets.csv, pooled.json, histories_final/, and the frozen batches, raw chains and audit records.
