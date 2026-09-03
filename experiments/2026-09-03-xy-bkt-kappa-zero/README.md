# Kappa-zero endpoint refinement

Requested September 3: add a few points at kappa=0 to reduce the horizontal uncertainty of the phase-curve endpoint.

Frozen design: epsilon/pi=0.550,0.570,0.575; L=16,32,64; two independent starts (cold and winding (1,0)), each with 7000 hybrid cycles. Each cycle contains one Metropolis sweep, four fixed Wolff clusters and one winding proposal. The provisional 2000-cycle cut leaves 10000 pooled retained cycles per parameter/volume group. This is a modest endpoint scouting batch, with less statistics per coordinate than the longer low-kappa refinement. It fits the remaining space while preserving the existing 25 GiB reserve.

Every initial and subsequent cycle field, physical observable, RNG state, cluster size/work record and timing is retained. The same independent full-field audit checks every generated field. The model/update code is unchanged. New seed streams are distinct from earlier campaigns. Earlier datasets and their uncertainty interpretation are preserved; the same joint U/R/eta two-standard-error rule updates only directly measured coordinates. These are finite-volume scouting brackets, not an infinite-volume BKT fit.

Run `.venv/bin/python experiments/2026-09-03-xy-bkt-kappa-zero/run.py endpoint --workers 4`. The exact manifest and unfiltered per-chain and batch logs are retained. No simulations are started by analysis scripts.
