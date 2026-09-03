# XY experiment storage review

Reviewed and executed 2026-09-03. The user subsequently authorized recommendations
2 and 3 below. The completed migration is recorded in
`storage-compaction-summary.json`, `storage-compaction.jsonl`, and the refinement
campaign's `pilot-deletion-manifest.json`.

## Completed migration

All 1,938 surviving trajectories were validated and rewritten. Every scalar and
topological observable remains present at every original cycle. Exact `phi` and
RNG states are retained at cycles divisible by 100, the standard cut-boundary
candidates, and both endpoints. Five explicitly listed representative trajectories
retain every field and RNG state. All retained arrays are losslessly compressed.

The 22 design-pilot HDF5 trajectories were deleted after their closure hashes,
manifests, job records, logs, production-selection summaries, and available
per-chain analyses were verified and recorded. They had already been excluded
from production inference.

The final reconciliation checked every compact file hash, the sparse selection
rule, matching field/RNG cycles, complete observable row counts, the five dense
representatives, and the absence and archival records of all pilots. It passed.
Surviving trajectories fell from 344.65 GiB to 8.01 GiB; pilot removal saved a
further 7.19 GiB. Total savings by file accounting are **343.83 GiB**. The
experiment directory now occupies about 8.8 GiB, and the filesystem reports
about 370 GiB available.

## Inventory

The eight experiment directories occupy 352.9 GiB. The filesystem has about
24 GiB free. There are 1,960 HDF5 trajectory files occupying 351.8 GiB.

| Campaign | Total GiB | Role |
|---|---:|---|
| `2026-09-02-xy-binder` | 0.05 | Early small-volume Metropolis test |
| `2026-09-02-xy-binder-bkt` | 0.94 | Earlier unconstrained Binder study |
| `2026-09-02-xy-binder-bkt-l48-cluster` | 1.31 | L=48 cluster-algorithm study |
| `2026-09-02-xy-bkt-refine` | 220.60 | Main three-volume refinement |
| `2026-09-02-xy-bkt-scout` | 87.52 | Broad parameter scout |
| `2026-09-02-xy-hybrid-k080` | 0.58 | Update-algorithm comparison |
| `2026-09-03-xy-bkt-kappa-zero` | 1.88 | Final kappa=0 endpoint refinement |
| `2026-09-03-xy-bkt-low-kappa` | 39.98 | Long low-kappa refinement |

Within the HDF5 datasets, the full spin fields `phi` occupy 319.1 GiB and the
per-loop winding arrays occupy 13.3 GiB. All remaining datasets together occupy
about 6.0 GiB. HDF5 allocation and metadata account for most of the difference
between those dataset sizes and the files' 351.8 GiB total.

The existing analysis reads the stored observable histories and winding arrays;
it excludes `phi` and `rng_state`. Thus the reported Binder cumulants,
correlation-length ratios, autocorrelations, bootstrap errors, and phase brackets
can be recomputed without the dense field history. Full fields remain useful for
inventing new observables, studying spatial vortex structure, reproducing the
original independent field audits, or restarting/replaying trajectories.

## Recommended order

1. **Losslessly rewrite completed HDF5 files with compression.** This sacrifices
   no information. Tests on 1,000-row samples from four representative L=64
   chains reduced the combined field and winding storage to 67.6--78.2% of its
   original size. Extrapolation suggests roughly 70--105 GiB of savings. The
   winding arrays alone compressed to 0.8--21.8% of their original size. Process
   one file at a time, verify all dataset values and attributes, record old and
   new hashes, and atomically replace only after verification. The largest source
   file is 1.04 GiB, so the current free space is sufficient for this procedure.

2. **Compact the dense field histories if substantially more space is wanted.**
   Preserve every observable at every cycle, but retain exact `phi` and RNG state
   only at the initial state, analysis-cut boundaries, the final state, and a
   fixed sparse cadence such as every 100 cycles. Preserve a small, explicitly
   selected set of complete trajectories for future field-level validation. A
   1-in-100 field cadence would remove about 316 GiB of the present `phi` payload
   before accounting for compression. It would preserve the present statistical
   analysis but would prevent arbitrary new field-level measurements on the
   discarded configurations.

3. **Discard raw design-pilot trajectories after retaining their manifests,
   summaries, bootstrap results, and selection records.** The 22 refinement
   pilots occupy 7.19 GiB. They were excluded from production inference and are
   unlikely to be needed as complete field histories.

4. **Discard raw trajectories superseded at identical parameter/volume
   coordinates, after retaining compact observables and analysis products.**
   Thirty scout chains (4.46 GiB) and twelve refinement chains (2.13 GiB) have
   later replacements at the same coordinates. They can still serve as
   independent replication checks, so their summaries should remain.

5. **Retire raw fields from the exploratory algorithm studies.** The early
   Binder, L=48 cluster, and hybrid-comparison campaigns use about 3 GiB in total.
   Their reports, manifests, source snapshots, derived tables, and representative
   histories are more likely to be useful than every saved field. These campaigns
   do not supply the final phase boundary.

6. **Drop per-cycle performance timings once benchmarking is finished.** The
   `*_seconds` datasets occupy about 1.1 GiB. They do not enter the physics
   analysis. Cluster-work counters and RNG states occupy less than 5 GiB
   combined, so deleting them alone has limited benefit.

Rendered histories, duplicate PNG/SVG/PDF figures, Python caches, and bootstrap
archives together are well below 1 GiB. Cleaning those first would add complexity
without materially relieving the storage pressure.

## Proposed preservation set

Always retain the frozen manifests and seeds, source archives and hashes, run and
audit logs, complete scalar observable series, production labels, analysis JSON,
bootstrap samples, phase tables, final figures, and reports. Any rewrite should
produce a conversion manifest mapping every old file hash to its replacement,
describe the exact field-selection rule, and validate the compact files before
the originals are removed.

The lowest-risk action is therefore lossless compression. If the aim is to free
hundreds of GiB, sparse retention of `phi` is the only change large enough to do
that; removing plots, caches, or scalar observables will not.
