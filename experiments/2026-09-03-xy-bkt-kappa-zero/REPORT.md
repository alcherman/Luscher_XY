# Kappa-zero endpoint refinement

Completed 18 trajectories on L=16,32,64 at epsilon/pi=0.550,0.570,0.575, all at kappa=0. Each coordinate/volume uses two independent 7000-cycle starts (cold and winding (1,0)). Every cycle contains one Metropolis sweep, four fixed Wolff clusters and one winding proposal. The 2000-cycle cut leaves 10000 pooled cycles per coordinate/volume.

The previous horizontal bracket was [0.540, 0.580]. The updated bracket is [0.550, 0.570], with plotting midpoint 0.5600 and half-width 0.0100. Its width is smaller by a factor 2.00.

| New epsilon/pi | Joint U/R/eta label | eta from L32,64 |
|---:|---|---:|
| 0.550 | quasi_ordered_like | 0.23779 +/- 0.00388 |
| 0.570 | massive_like | 0.26186 +/- 0.00431 |
| 0.575 | unresolved | 0.25857 +/- 0.00431 |

The same joint U/R/eta two-standard-error rule is used, including unresolved coordinates. These are finite-volume scouting brackets; they are not Gaussian confidence intervals or a thermodynamic extrapolation. Existing higher-statistics estimates at other coordinates remain in use. This small added batch has lower statistics per coordinate than the longer low-kappa campaign.

All 18 full-field audits and the final reconciliation passed: 126,018 saved fields, 2.006 GB, with initial/transient fields, physical observables, RNG, work and timing records retained. Five analysis tests passed. All 292 prior phase classifications were reproduced from unchanged input files. The model package is unchanged.

Diagnostic flags: 0 individual trajectories and 0 pooled groups. Blocks per trajectory range from 250 to 250. See diagnostic-summary.json and per-chain analysis for autocorrelations, alternate cuts, doubled blocks, retained-half checks and independent-start comparisons. Initial/full/retained histories are saved in histories_final/. Visual inspection of the epsilon/pi=0.570 initial and retained histories supports the cut; doubled-block standard errors range from 0.862 to 1.119 times the primary estimates. See thermalization-review.json.

Artifacts: phase_boundary_refined.png/pdf/svg, endpoint_diagnostics.png/pdf/svg, phase-brackets.csv, endpoint-observables.csv, pooled.json, execution-audit.json, and complete frozen manifests and raw chains. The separate target of shrinking all low-kappa vertical bars remains unfinished; this batch targets only the kappa=0 horizontal bracket.
