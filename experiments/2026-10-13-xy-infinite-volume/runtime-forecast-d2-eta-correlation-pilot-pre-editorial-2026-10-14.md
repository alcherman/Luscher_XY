# Runtime forecast for the 2D eta and correlation length campaign

All54 pilot trajectories passed their field audits. Generation and auditing completed in 12.86 minutes from manifest creation to last closure, totaling 2.477 worker-hours. They stored 140.86 GiB with full retention. The computational baseline remains on track for Friday October16 at3p.m. Central; 55.0 hours remain at this forecast. This does not establish infinite-volume precision.

## Measured cost

| N | Trajectories | Median generation seconds | Median audit seconds |
|---|---|---|---|
| 64 | 18 | 26.66 | 3.44 |
| 128 | 18 | 87.11 | 10.53 |
| 256 | 18 | 323.48 | 38.57 |

Equivalent concurrent worker processes averaged 11.55 under the18-worker nighttime cap. This sums process wall times, including simulation I/O and auditing; it is not a CPU utilization measurement. Applying the observed efficiency to the16-worker weekday cap gives 10.27 equivalent workers. The pilot had no observed removable-volume permission interruption.

## Baseline remaining plan

There are22 candidate requested slices epsilon/pi=.58,.60,...,1.00. Three have larger-box pilot data, leaving19. A baseline assumes three beta samples per slice, independent starts2/3, and12000 cycles atN128/256:228 additional trajectories. If every slice also needs fresh matching-betaN64 anchors, the count becomes342. Prefer compatible SeptemberN16/32/64 data and matching-betaN64 anchors. The54 completed trajectories are additional to those reused archives.

The currently frozen24 bracket extensions add shallow-response coverage at.60 and a lower anchor at.80. Including those extensions, the measured per-size means imply 15.79–16.79 worker-hours and about 1.54–1.63 hours at the observed efficiency under the weekday cap. Budget2–4 hours for baseline generation, allowing different cutoffs, I/O and batch tails. This is not a promise that three beta samples will suffice.

Reserve a further6–18 hours as a working allowance for analysis, cut/block and start diagnostics, local refinement, and selected extra statistics or larger-size checks. The resulting8–24-hour working allowance is not a measured completion bound. The principal risk is BKT inference and the shallow response near the low-epsilon end, rather than the measured simulation throughput.

## Coverage and storage limits

At epsilon/pi=.60, xi/L changes only slightly across the pilot beta range, and eta128:256 does not yet bracket.25. Several bootstrap root tails therefore leave their measured windows. At.80 theN128 correlation-ratio lower endpoint lies too close to its target. The dedicated24-run extensions address these coverage failures before full larger-box commitment. The limiting targets are not finite-box critical values, and no pilot infinite-volume coupling is accepted solely from the small-box crossings.

Current free space is 2208.0 GiB, with25GiB reserved. Based on observed stored bytes, the remaining baseline plus24 extensions needs about 910.4–954.7 GiB; frozen conservative estimates give 1115.3 GiB. A uniform fourfold statistics plan acrossall22 slices would require approximately 3991.0 GiB of new files beyond the completed pilot and does not fit. Use targeted refinement and honest wider bars, retaining everything.

A possible sensitivity example of three slices atN512, three beta values, two starts and12000 cycles would add18 trajectories and roughly 427.2 GiB by volume scaling. Its runtime and mixing need a separate measured pilot; no broadN512 campaign is committed.

The completed3D/4D outputs remain intact. Optional5D remains disabled until the2D/3D/4D deliverables are ready and the remaining budget supports its separate analysis. No data deletion, thinning, compaction or environment change is part of this plan.

[Full forecast data](/Users/acherman/research/Luscher_XY/experiments/2026-10-13-xy-infinite-volume/runtime-forecast-d2-eta-correlation-pilot-2026-10-14.json)
