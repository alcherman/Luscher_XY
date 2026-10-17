# Runtime forecast for the 2D eta and correlation length campaign

All 54 pilot trajectories passed their field audits. Generation and auditing completed in 12.86 minutes from manifest creation to last closure, totaling 2.477 worker-hours. They stored 140.86 GiB with full retention. The computational baseline remains on track for Friday October 16 at 3 p.m. Central; 54.9 hours remain at this forecast. This does not establish infinite-volume precision.

## Measured cost

| N | Trajectories | Median generation seconds | Median audit seconds |
|---|---|---|---|
| 64 | 18 | 26.66 | 3.44 |
| 128 | 18 | 87.11 | 10.53 |
| 256 | 18 | 323.48 | 38.57 |

Equivalent concurrent worker processes averaged 11.55 under the 18-worker nighttime cap. This sums process wall times, including simulation I/O and auditing; it is not a CPU utilization measurement. Applying the observed efficiency to the 16-worker weekday cap gives 10.27 equivalent workers. The pilot had no observed removable-volume permission interruption.

## Baseline remaining plan

There are 22 candidate requested slices epsilon/pi=0.58,0.60,...,1.00. Three have larger-box pilot data, leaving19. A baseline assumes three beta samples per slice, independent starts 2 and 3, and 12,000 cycles at N=128/256: 228 additional trajectories. If every slice also needs fresh matching-beta N=64 anchors, the count becomes 342. Prefer compatible September N=16/32/64 data and matching-beta N=64 anchors. The 54 completed trajectories are additional to those reused archives.

The currently frozen 24 bracket extensions add shallow-response coverage at 0.60 and a lower anchor at 0.80. Including those extensions, the measured per-size means imply 15.79–16.79 worker-hours and about 1.54–1.63 hours at the observed efficiency under the weekday cap. Budget 2–4 hours for baseline generation, allowing different cutoffs, I/O and batch tails. This is not a promise that three beta samples will suffice.

Reserve a further 6–18 hours as a working allowance for analysis, cut/block and start diagnostics, local refinement, and selected extra statistics or larger-size checks. The resulting 8–24-hour working allowance is not a measured completion bound. The principal risk is BKT inference and the shallow response near the low-epsilon end, rather than the measured simulation throughput.

## Coverage and storage limits

At epsilon/pi=0.60, xi/L changes only slightly across the pilot beta range, and eta_eff(128,256) does not yet bracket 0.25. Several bootstrap root tails therefore leave their measured windows. At 0.80 the N=128 correlation-ratio lower endpoint lies too close to its target. The dedicated 24-run extensions address these coverage failures before full larger-box commitment. The limiting targets are not finite-box critical values, and no pilot infinite-volume coupling is accepted solely from the small-box crossings.

Current free space is 2184.6 GiB, with 25 GiB reserved. Based on observed stored bytes, the remaining baseline plus 24 extensions needs about 910.4–954.7 GiB; frozen conservative estimates give 1115.3 GiB. A uniform fourfold statistics plan across all 22 slices would require approximately 3991.0 GiB of new files beyond the completed pilot and does not fit. Use targeted refinement and honest wider bars, retaining everything.

A possible sensitivity example of three slices at N=512, three beta values, two starts and 12,000 cycles would add 18 trajectories and roughly 427.2 GiB by volume scaling. Its runtime and mixing need a separate measured pilot; no broad N=512 campaign is committed.

The completed 3D/4D outputs remain intact. Optional 5D remains disabled until the 2D/3D/4D deliverables are ready and the remaining budget supports its separate analysis. No data deletion, thinning, compaction or environment change is part of this plan.

[Full forecast data](/Users/acherman/research/Luscher_XY/experiments/2026-10-13-xy-infinite-volume/runtime-forecast-d2-eta-correlation-pilot-2026-10-14.json)
