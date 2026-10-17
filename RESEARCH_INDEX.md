# Research handoff: 17 October 2026

Start with [AGENTS.md](AGENTS.md) and [PYTHON_ENVIRONMENT.md](PYTHON_ENVIRONMENT.md).
Use the existing `.venv/bin/python`; do not synchronize or upgrade the environment
as a side effect of analysis. Historical manifests and source snapshots are
preserved evidence, not instructions to relaunch a campaign.

## Current results

- **3D negative-kappa continuation:**
  [reviewed delivery pointer](experiments/2026-10-17-xy-3d-negative-kappa/latest-exploratory-negative-kappa-delivery.json),
  [report](experiments/2026-10-17-xy-3d-negative-kappa/exploratory-negative-kappa-transition-report.txt),
  [full-range plot](experiments/2026-10-17-xy-3d-negative-kappa/exploratory-negative-kappa-transition-curve-full-range.png),
  and [independent review](experiments/2026-10-17-xy-3d-negative-kappa/exploratory-negative-kappa-scientific-review.json).
  The measured finite-box brackets at epsilon/pi = .52, .51, .505 are
  [-2.1,-2.0], [-4.4,-4.0], [-8.8,-8.0]. They support an exploratory
  approach toward one-half, without proving a thermodynamic asymptote or
  supplying calibrated confidence intervals. The q=.505 primary cohort uses
  an independently tested exact restricted-angle heatbath; the older trapped
  Metropolis cohort remains preserved separately. The exploratory goal is
  complete; no further simulations are planned under it.
- **Positive-coupling phase boundaries, D=2,3,4,5:**
  [attained combined delivery pointer](experiments/2026-10-13-xy-infinite-volume/latest-attained-d2-d3-d4-d5-delivery.json).
  Read its report and reviews before quoting numbers. Descriptive sensitivity
  unions are not certified total confidence intervals; rejected fits remain
  diagnostic. The full tree-level criterion supersedes the earlier straight-line
  approximation. D2's practical primary data are capped at N=512; the two larger
  probes are preserved separately.
- **3D scaling dimension and spin correlations:**
  [common-xi delivery](experiments/2026-10-15-xy-3d-spin-dimension/latest-common-xi-delivery.json),
  [four-size crossing plot](experiments/2026-10-15-xy-3d-spin-dimension/latest-four-size-xi-plot.json),
  [dimension readout plots](experiments/2026-10-15-xy-3d-spin-dimension/latest-dimension-readout-plots.json),
  [existing-data joint fits](experiments/2026-10-15-xy-3d-spin-dimension/latest-joint-fss-existing.json),
  and [spin-correlation analysis](experiments/2026-10-15-xy-3d-spin-dimension/latest-spin-correlations.json).
  These have exploratory qualifications; do not turn descriptive sensitivity
  ranges into a calibrated systematic standard deviation.
- **Collaborator report and talk figures:**
  [package guide](output/pdf/xy-collaborator-brief-2026-10-16/README.txt).
  The report and atlas summarize the earlier campaign; the October 17 negative-
  kappa continuation above is newer.

## Paused work and resources

The separate [spin follow-up control](experiments/2026-10-16-xy-3d-spin-followup/control.json)
remains user-paused with new launches disabled. Do not resume it from historical
manifests. The most recent allocation for the negative-kappa work was **four
CPU-heavy processes total**, implemented as at most three single-thread
simulations plus one analysis. Its simulations have finished. Respect newer
user resource instructions before launching anything.

## What Git contains, and where raw data live

Git contains scientific source, tests, frozen designs, job/audit records,
analysis snapshots and derived samples, provenance, reviews, numerical tables,
and presentation artifacts. Rebuildable cache directories, raw HDF5 files,
full runner logs, temporary files, lock files and machine-local `chains`
symlinks remain on disk outside Git. Nothing was deleted for this commit.

[RAW_DATA_LOCATIONS.json](RAW_DATA_LOCATIONS.json) records the current external
links and local raw-file directory counts. The bulk archive is documented in
[the verified archive record](experiments/2026-10-16-mc-data-hdd-archive/README.md).
Its manifests/checksums are versioned. Recent negative-kappa raw files remain on
Office_SSD. A Git clone alone is **not** the raw Monte Carlo archive.

Some historical JSON/plist provenance uses absolute paths for this Mac.
Preserve those records; use the relative links in this index for navigation.
On another machine, configure and verify raw-data locations before reanalysis;
never bootstrap historical launchd plists merely to open results.
