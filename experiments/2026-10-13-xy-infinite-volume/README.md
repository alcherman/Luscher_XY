# Infinite volume XY phase boundary campaign

Target: beta_c(epsilon) for d=2,3,4 at fixed epsilon/pi values spaced by 0.02
on the positive-beta branch, extrapolated to infinite volume. Aim for vertical
uncertainty half-width 0.02, allowing wider reported bars or unresolved slices
within the time budget. Deliver per-dimension and combined plots with vertical
error bars only, the numerical table, and a scientific uncertainty report by
Friday October 16, 2026 at 3 p.m. America/Chicago.

Dimension 5 is optional and last in priority. Enable it only after the d=2,3,4
deliverables are ready and a resource/time forecast supports an additional
campaign before the deadline. `enable_d5` is initially false; the runner refuses
5D boxes above N=12 without a separately reviewed budget.

The detailed estimator design is in
`../../docs/xy-infinite-volume-phase-boundary-design-2026-10-13.md`. Current
stage: additional endpoint size checks and constrained-slice reconnaissance.
The initial ordinary XY endpoint calibration completed. This initial
batch has 36 trajectories: dimensions 3 and 4, N=6,8,12, three couplings near
the respective unconstrained critical points, and two independent cold/winding
starts. Each trajectory has 4000 hybrid cycles. The 1000-cycle warmup label is
provisional; it does not certify equilibration. Every initial and transient row
is saved. These pilots do not establish final infinite-volume transition data.

## Detached operation and restart

Raw-data storage moved to `/Volumes/Office_SSD/Luscher_XY/2026-10-13-xy-infinite-volume/chains`.
The local `chains/` path is a symlink to that directory. All 36 completed files
matched their recorded SHA-256 hashes after copying; the originals remain in
`chains-internal-copy-2026-10-13/`. `storage-routing.json` records every verified
file and the coordinator wrapper hash. Keep Office_SSD connected while jobs
run. The SSD wrapper refuses new launches if the volume is unmounted or the
raw-data link differs from this route.

The detached service now runs `coordinator_ssd.py`, which imports the unchanged
frozen `run.py` and routes its disk-space guard to the raw-data filesystem.
The original launch configuration is retained as `runner.internal-storage.plist`.
Source hashes and scientific trajectory parameters remain unchanged. Two
storage-routing tests passed, covering external-filesystem selection, refusal
when unmounted, and refusal of an unexpected raw-data path.

`runner.plist` defines a user macOS launchd job. Once bootstrapped into
`gui/501`, the coordinator and simulation workers run independently of Codex.
Updating or restarting the Codex app does not terminate them. The plist stays
in this campaign directory; no system service or privileged daemon is installed.
The hourly Codex heartbeat inspects progress and advances validated adaptive
batches after the app is available. Numerical jobs do not depend on that
heartbeat while a batch is running.

The coordinator waits for additional frozen manifests under `batches/` until
the deadline, so follow-up batches do not require another service. It uses an
exclusive coordinator lock and each trajectory has its own exclusive lock.
After a coordinator interruption it identifies and adopts live workers by
their exact command and trajectory ID. An interrupted worker resumes the same
file from its last committed field/RNG checkpoint; it does not reseed.

Use `launchctl print gui/501/org.luscherxy.phase-boundary-20261013` to inspect
the service. Do not bootstrap a duplicate. `status.json` gives counts and the
current limits. `jobs/` records each trajectory and field-audit closure hash.
`logs/coordinator.log`, `logs/coordinator.stderr.log`, and per-trajectory logs
retain full unfiltered output. Numerical histories are in `chains/`.

The app-restart guarantee is process independence, not protection from a
machine shutdown or power loss. Graceful worker termination closes HDF5 and
allows exact checkpoint continuation. Abrupt OS termination can require HDF5
recovery; preserve files and diagnose rather than deleting or replacing them.

## Resources and retention

The launchd classification changed from Background to Standard on October 13
after runnable workers were observed receiving about one third of a CPU each.
The original background classification imposed macOS resource restrictions
more conservatively than intended for this user-requested compute campaign.
The coordinator was restarted while independently session-separated workers
continued; it adopted the live workers. New workers inherit Standard
scheduling. The old launch configuration is retained as
`runner.background-scheduling.plist`; `scheduling-change-2026-10-13.json`
records the before/after settings and adopted PIDs. `Nice=5` is retained, and
office-hour worker caps still leave two cores available. This changes resource
scheduling only, not the frozen numerical kernel or simulation parameters.
Report measured worker CPU utilization and throughput rather than assuming
the classification change guarantees a performance gain. The CPU reports
Apple M5 Max with six Super and twelve Performance cores.

The machine has 18 physical cores and 64 GiB memory. Workers use one Numba and
BLAS thread each. During weekdays 9 a.m. to 5 p.m. Central, at most 16 workers
run. Outside office hours the ceiling is 18, subject to `control.json` and
observed memory/IO limits. If the ceiling decreases while jobs run, surplus
process groups are paused without being killed, then resumed when slots return.
Pause/resume events are logged; elapsed timings can include a pause.

The first calibration workers measured about 160 to 220 MiB resident memory
each. `max_workers` is now 18, automatically capped at 16 during office hours.
Continue checking memory and IO behavior as sizes increase. `launch_new=false` stops new
launches while letting current trajectories finish. No new trajectory is
launched after the Friday 3 p.m. deadline. A trajectory already running may
finish after that cutoff; the heartbeat should checkpoint and stop work if
finishing it would interfere with delivery.

Each trajectory retains all scalar and topological observables, every spin
field and RNG state, initial/transient rows, timings and update counters. There
is no new thinning or compaction. A conservative dense-record estimate gates
launches against currently live reservations and a 25 GiB free-space reserve.
Disk usage and process memory must be monitored before adding larger boxes.
The runner refuses newly frozen 4D boxes above N=24; do not launch 64^4 boxes.

Frozen manifests include source and environment-file hashes. Never edit the
runner or hashed model sources after freezing a batch. To fix a runner bug,
stop new launches, let current writes close safely, preserve the old runner
and manifest, and introduce a separately versioned replacement. Estimator code
may evolve independently, with its versions and input hashes recorded.

## Validation

Operating notes: a two-second current-utilization sample after the scheduling
change reported 86.71 percent user CPU, 2.72 percent system CPU and 10.55
percent idle for the whole machine. This is a descriptive snapshot, not an
independent throughput-improvement measurement. Its full output is retained
in `logs/cpu-utilization-sample-2026-10-13.log`. Worker process-table averages
alone should not be treated as the current whole-machine utilization.

The Mac's configured idle-sleep timer was one minute. A separate detached
`keep-awake.plist` job now runs caffeinate with an idle-system-sleep assertion
through Friday October 16 at 3:15 p.m. Central, so numerical work can proceed
unattended. It does not alter display-sleep or screen-lock settings. The
assertion was verified with pmset. Remove the keep-awake launchd job early if
the campaign completes; otherwise its timeout releases the assertion. A
machine shutdown, power loss or explicitly requested system sleep can still
interrupt computation and must be handled by checkpoint recovery.

`test/test_xy_infinite_volume_campaign.py` checks the office-hour boundaries,
independent d=2,3,4 field recomputation, rejection of deliberately corrupted
observables, and exact split-run checkpoint continuation of all physical and
RNG records. The focused XY suite is run before freezing and launch. Future
analysis must test its statistical and finite-size fitting methods separately;
correct trajectory generation alone does not validate an extrapolation.

Before launch, the focused XY suite passed 170 tests with one skip. The final
campaign-specific run passed 9 tests, including termination of a separate
worker process and exact recovery of all physical and RNG records compared
with an uninterrupted trajectory. Full validation logs are under `logs/`.
The service was bootstrapped on October 13, and launchd reported it running;
the coordinator's parent process is PID 1. No Codex app restart was performed
as part of verification.

## Next analysis work

On October 13 the first SSD-routed detached workers reached HDF5 open but did
not create files. A sampled worker stack was waiting in the OS open call;
macOS privacy logs attributed external-volume requests to the background
Python executable. Computer Use reported the Mac locked. The user was asked
to unlock it and approve the Python removable-volume prompt if presented.
The user approved access, all 16 active workers created their SSD data files,
and `launch_new` was reenabled. System Settings was then inspected directly:
Privacy & Security > Files & Folders > python3.14 > Removable Volumes is ON.
That standing permission applies to the current background Python runtime.
Keep the existing runtime fixed during the campaign. Before a future campaign
changes runtime or storage routing, test one small write through the actual
background launcher while the Mac is unlocked, then verify progress before
launching the full batch. Do not count an alive but blocked process as
productive throughput. The first SSD-routed larger-box generation wall times
include the permission wait; label or exclude that wait from cost forecasts
without altering the raw timing records. The verified SSD copies and preserved
internal originals remain readable.

`analyze_calibration.py` has analysed the 18 original parameter/volume groups,
saved joint U/R/Q bootstrap samples, covariance, histories, cut/block tests,
and start comparisons in `calibration-analysis/`. The first broad linear
finite-size fits are rejected by chi-square: 175.4/5 in 3D and 646.1/5 in 4D.
The nominal ratios of fit coefficients must not be quoted as critical
couplings. Subleading correction and minimum-size tests with only three
sizes are explicitly unidentifiable. These observations call for narrower
coupling windows or validated histogram reweighting and the queued size
checks. No infinite-volume critical coupling has been accepted. Two synthetic
fit tests passed, checking recovery of known roots and rejection of
underidentified correction models; that is not an empirical calibration PASS.

The initial runtime forecast has been reported to the user and independently
reviewed. Its arithmetic passed; cautious first-pass feasibility before Friday
is supported, but a guarantee of complete validated +/-0.02 curves is not.
The 12–36-hour estimate is a planning allowance, not an observed runtime bound.
`runtime-forecast-2026-10-13.md` now uses the SSD capacity. The original
internal-disk forecast remains in `runtime-forecast-internal-storage-2026-10-13.md`.
Consult `milestones.json` before repeating the first-forecast notification.

Two new frozen batches have been added: 18 endpoint size-check trajectories
(N=16 in 3D and N=10,16 in 4D), and 192 reconnaissance trajectories on the
requested constrained epsilon slices. The scout uses N=6,8, two starts and
4000 cycles at beta=0,0.25,0.50 in 3D or beta=0,0.20,0.40 in 4D. These broad
anchors locate coupling neighborhoods and measure constrained-slice costs;
they are not final critical-coupling fits. Analyse closed groups while batches
continue, then choose refined beta neighborhoods and additional sizes from
actual diagnostics. All initial/transient fields and RNG states are retained.

As soon as the 36 calibration trajectories finish and pass audits, prepare
the requested full d=2,3,4 runtime forecast and notify the user whether the
campaign is on track for Friday. Use measured generation plus audit costs,
actual throughput under the 16 daytime and 18 otherwise schedule, stored
bytes, reused d=2 data, and remaining deadline time. Give explicit trajectory
counts, beta samples, size sets, independent starts and cycle counts for the
baseline, plus scenarios for longer statistics and larger boxes. Include
analysis/validation time and the disk reserve. Distinguish sampling completion
from validated infinite-volume precision. Save a descriptive forecast report;
revise it when measured cost or required scope meaningfully changes.

The hourly heartbeat should first verify the 36 pilot closures and produce
history, cut, block-size and independent-start diagnostics. Fit the ordinary
3D and 4D endpoints and compare with literature without anchoring the fits to
those couplings. Three pilot sizes are insufficient to freely fit a 4D
critical ratio and two separate correction amplitudes plus the critical
coupling: explicitly check rank/identifiability rather than reporting a small
optimizer covariance. Add selected N=16 measurements for size-cut and
correction-model tests after using the measured pilot cost for the forecast.

In parallel, reanalyse the existing d=2 pooled measurements at the 22 requested
slices using joint BKT matching, reference-range guards, covariance and a
minimum-size sensitivity check. A fit that maps the reference outside its
validated range is unresolved. Prioritize estimator validation, then the full
grid, then uncertainty reduction, and optional d=5 last. Every added batch
must use this frozen runner and its unchanged model source hashes; keep the
analysis implementation separate and versioned. Do not confuse pilot field
audit success with demonstrated infinite-volume uncertainty.

On October 13 at about 15:56 Central, all 246 initial trajectories had closed
and passed their independent field audits. `scout.py` summarized 96 paired-start
groups; ten groups have diagnostic flags (mostly low-beta correlation-ratio
bootstrap samples where the square-root argument becomes negative). Undefined
ratios are written as JSON null, and the bootstrap NPZ retains the original
values; no data are removed. The first summary serialization attempt failed
on a NaN; its complete log and temporary summary remain preserved.

The frozen `beta-neighborhood-scout` batch adds 244 trajectories: four interior
beta values at N=6,8 with two independent starts on 15 finite-box design
brackets, plus d3 epsilon/pi=.80 beta0 branch checks at N=12,16. The ordinary
endpoint N8 R merely selects simulation neighborhoods, and is not used as a
fixed critical ratio in an infinite-volume fit. The .80 slice was not assigned
a positive critical coupling. Estimated full-retention storage is 24.78 GiB.

`endpoint-narrow-window-validation` adds 30 trajectories: d3 beta=.453,.457
at N=8,12,16; d4 beta=.300,.304 at N=8,10,12,16; and d4 beta=.302 at N=20,
all with two starts and 4000 cycles. Existing midpoint measurements are reused.
These measurements address beta-window curvature and empirical histogram
overlap, while N20 provides another size for correction-model/minimum-size
checks. Their conservative storage estimate is 33.01 GiB. The N20 planning
cost is about55 minutes per chain from volume-scaled measured N16 wall time,
including any interruptions; this is not a claim of productive throughput
or adequate mixing. Details are in `endpoint-narrow-window-budget.json`.
No infinite-volume result is accepted by these design steps. During active
analysis the coordinator temporarily reserves an extra core, then restores
the daytime16/outside18 worker schedule.

All 54 initial endpoint/size-check trajectories have now been analysed as
27 paired groups. Adding sizes did not repair the broad linear beta-window
fits: leading all-size chi-square/dof is330.02/8 (3D),1893.01/11 (4D).
Near-critical 4D groups include fewer-than-50-block flags; narrower-window
measurements do not remove the need to review mixing and extend statistics.
The nominal roots remain rejected and must not enter final tables or plots.

The user simplified the primary estimator to a fixed U=1/4 crossing in each
box, followed by size extrapolation (design addendum takes precedence over
the earlier joint-matching requirements). `threshold_calibration.py` now
applies this to all 84 completed epsilon/pi=1 calibration trajectories,
without reweighting or beta extrapolation. See
`binder-quarter-calibration-2026-10-13/` for CSV, JSON, bootstrap records and
PNG/PDF/SVG plots. Direct brackets exist at d3 N12,16 and d4 N6,8,10,12,16;
d3 N6,8 and d4 N20 are unbracketed and omitted, not invented. Leading
intercepts are0.4556049 (3D, two-size provisional) and0.3018244 (4D,
conditional MC95[0.3014724,0.3021165], chi2/dof3.402/3). Interpolation,
size-cut/correction and cut/block variants are recorded separately. Black
rings flag sparse independent blocks. Four synthetic tests passed and an
independent review reproduced the intercepts and d4 chi-square. Next use
threshold scouting rather than requiring acceptance of the broad linear
observable fits. Bracket missing sizes economically when needed.

At the October13 16:50 heartbeat, the detached coordinator was reused and
three new frozen batches were advanced:120 constrained N12/16 threshold
brackets (102.80GiB conservative estimate),180 N6/8 local bracket refinements
(18.12GiB), and54 d2 pilot trajectories at epsilon/pi=.60,.80,1, N16/32/64
(5.25GiB). All retain complete fields/RNG/transients/logs. d3 .80 has no
positive-beta small-box quarter bracket; it remains a branch-check case.
Small-box root predictions only design future couplings; in particular the
coarse d3 .90 bracket[0,.25] is being refined before accepting its crossing.

All542 earlier trajectories, including endpoint bracket extensions, passed
closure audits. The refined endpoint snapshot is in
`binder-quarter-calibration-refined-2026-10-13/`: all4 d3 sizes6/8/12/16
and all6 d4 sizes6/8/10/12/16/20 now bracket quarter. Leading conditional
intercepts are .4555147 (MC95 [.4538244,.4568697],chi2/dof3.204/2) and
.3018292 (MC95 [.3015945,.3020288],chi2/dof3.403/4). An independent
review reproduced these numbers. Model sensitivities and low-block flags
remain separate from MC error; valid-draw fractions are retained. The
doubled-block d3 joint crossing vector retains98.4% of draws, so its
interval is explicitly conditional on complete draws. No total-error
precision certification is asserted. Zero-dof variants cannot test fit
quality. Six focused threshold/storage tests passed after metadata fixes.
The original user-shown calibration snapshot remains in its original folder.

Follow the updated fixed-threshold automation prompt and design addendum.
Raw files remain on the verified SSD route; source hashes of all9 frozen
manifest batches agree, all896 planned IDs are unique, and worker caps
restore18 (16 during office hours) after brief analysis reservations.

October13 17:48 heartbeat: all896 planned trajectories had passed field
audits. The measured d2 pilot cost was modest, but q=.60 remained above
quarter even at beta0,N64. A threshold-choice question was asked. The user
then explicitly deferred2D;316 pending d2 jobs were marked `deferred`, and
none were running. No raw files or manifests were removed. `control.json`
now declares deferred_dimensions=[2], the adaptive designer skips these
dimensions, and the heartbeat prompt preserves this user steering. Resume
2D only when the user returns to it; do not let earlier scope override this.

D3/d4 continue. Frozen adaptive followup406 includes36 d2 local-refinement
trajectories already completed before the deferral, plus370 d3/d4 followups.
Unstarted expanded d2 scout306 and size probe10 remain deferred. An extra
72 d4 size-sensitivity trajectories (N10,N20,three beta anchors,two starts,
4000cycles on six constrained slices) test simple finite-size approximations.
Their conservative storage estimate is256.62GiB; roughly5.8 productive
worker-hours from actual endpoint timings, with up to4x planning overhead
allowed. This fits the SSD reserve and multi-day budget; no64^4 boxes.

The first conditional full-curve snapshot is
`threshold-curve-snapshot-20261013-175611/`:10 d3 preliminary locations and
2 d4 conditional locations,4 rejected leading d4 size fits,2 insufficient
size cases. It is a diagnostic preview, not a final accepted curve. Four
rejected d4 leading laws must not be reported as accepted critical values.
Omitted-small-size/log-correction checks and finer beta brackets are ongoing.
Weighted monotonic regression regularizes noisy Binder responses (declared
shape assumption) consistently in central data and bootstrap draws; raw
central roots/pooling shifts retained. Maximum central pooling shift was
2.01 SD. An independent review found coherent conditional inference, and
requested explicit conditional coverage, full contributing flags, independent
epsilon seeds and size-cut masks. These changes are implemented. Bootstrap
helper before seed revision is preserved, caches are versioned by its hash,
and later curve snapshots copy their analysis sources. Nine focused
threshold/calibration/storage tests passed. No model or generator source
hash changed, and all13 frozen manifests retain1690 unique IDs.

October13 18:38 user status request: all1374 nondeferred trajectories
closed with passed audits and no failures. D2 remains deferred316jobs.
Updated reviewed snapshot `threshold-curve-snapshot-20261013-183856/`
contains10 conditional d3 locations (.82..1) and6 conditional d4 locations
(.90..1). Finer beta sampling removes the earlier d4 .92-.98 leading-fit
rejections; this does not establish final total-error accuracy. D4 .88
remains strongly rejected (chi2/dof80.3/1), and five of6 constrained N20
slices are not yet quarter-bracketed. Nominal rejected points are diagnostic
crosses only, not accepted critical couplings. Shifted beta windows and
local refinements134trajectories are now queued (279.17GiB conservative
retention estimate), plus12 sampling-boundary coverage runs (6.78GiB) and
4 targeted independent-start diagnostics (5.52GiB). The latter use new
replicas2/3 with16000cycles at d3 .92 N16 beta=.293 and d4 .98 N10
beta=.261, where start disagreement occurs below the primary brackets.
No root bias is established; review histories/cuts/blocks before final claims.
These targeted followups remain within the multi-day time and SSD budget.
The adaptive designer now skips deferred dimensions before opening raw files
and accepts groups with additional independent replicas. All original
snapshots, raw data, source versions and bootstrap caches remain retained.

October13 18:52 heartbeat: shifted-window runs continue without failures;
2D remains deferred. The curve CLI now supports cut/block/threshold options
and an explicit frozen source snapshot. Three variants (cuts501/2001 and
doubled blocks atcut1001) used exactly1284 original d3/d4 trajectories.
`compare_threshold_stability.py` refits the intersection of bracketed sizes
for17/18 slice cases (d3 .80 unresolved). Maximum leading-intercept shift
is.001970 at d3 .94; conditional d4 .90..1 shifts are at most.000172.
Doubled-block MC SD ratios range.897-1.097; retained complete-draw fractions
are at least.998. An independent review reproduced these numbers. These
are overlapping-data descriptive checks, not mixing or total-error PASS.
The combined report is `threshold-cut-block-stability-2026-10-13.json`.

`check_independent_starts.py` retains full-history displays and running
Binder estimates at the two flagged below-crossing coordinates, with
cut/block bootstrap diagnostics. The initial r0/r1 Binder discrepancies
persist aftercut2001 (gaps .11859 and.11202); do not attribute them solely
to transient behavior. New independent16000cycle r2/r3 runs are still
pending. Pending counts now consult frozen manifests, not only materialized
job files; an earlier diagnostic snapshot undercounted unstarted replicas.
All snapshots remain retained.

Two modest d3 batches add8 branch-check trajectories at epsilon/pi=.80,
N24/32,beta0/.010 (7.31GiB estimate), and18 representative N24 threshold
measurements at epsilon/pi=.84,.92,1 (10.19GiB). These test endpoint
membership and size-law sensitivity rather than requiring full Binder
crossover fits. No2D job is resumed or analysed. After the analysis core
reservation the coordinator cap restored18 for nighttime; weekday16
remains automatic. Nine focused threshold/calibration/storage tests pass.

October13 19:52 heartbeat: all1550 nondeferred trajectories closed with passed
field audits. Snapshot `threshold-curve-snapshot-20261013-195243/` includes
all constrained d4N20 quarter brackets and d3N24 at .84,.92,1. An independent
review reproduced primary-fit checks. Strict leading size fits at d4 .88,
.92,.96 remain rejected; the latter two maximum residuals are about .0011
and .0014 in beta, small compared with the requested .02. Rejection is
retained rather than converting nominal roots into accepted couplings.
The full tested interval unions span less than .02 in each plotted slice,
but these descriptive unions have no total-error coverage guarantee.

Four new16000-cycle independent-start diagnostics finished. Last-half U
estimates are .09476/.08909 in d3 and .10336/.10063 in d4 (per-start MC SD
about .013-.019). Window checks are in the latest `long-start-window-check-*`
folder; apparent late agreement does not prove equilibration. The initial
discrepant starts and flags remain retained. Independent leave-coordinate-out
review finds the central local-linear quarter roots unchanged, with d4
PCHIP sensitivity .000190. This does not certify all bootstrap effects zero.

Eight d4 beta0 neighbor checks at epsilon/pi=.86,N8/12/16/20 and36 local
N20/N24 sensitivity trajectories at .88,.92,.96 are now queued (21.36 and
358.52GiB conservative estimates). The latter narrow coarse N20 chords and
add N24 size checks. Approximately2.8TiB SSD free space remains, far above
the25GiB reserve. There is no observed deadline risk; no new forecast
precision is asserted. All20 manifests retain identical numerical source
hashes and1910 unique IDs. D2 remains deferred316jobs; d5 remains disabled.

October13 20:37 user status request: all1594 nondeferred trajectories closed
with passed field audits, no failures. Snapshot `threshold-curve-snapshot-20261013-203715/`
incorporates44 followups. Narrowed N20 beta chords remove strict leading-fit
rejections at d4 .92 (chi2/dof4.085/5, includingN24) and .96(1.259/4).
This supports the local interpolation update; it is not proof of total
accuracy or equilibration. D4 .88 remains rejected(17.45/4), with N24
quarter still unbracketed. N24 .96 needs a lower-beta bracket. Ten further
N24 runs extend .88 to .048/.051, refine .92 at .146, and extend .96
to .228/.229 (133.02GiB conservative estimate). All raw records retained.
D4 .86 beta0 Binder rises from .3871 atN8 to .4911 atN20; the missing
positive finite-box quarter bracket does not by itself settle branch
membership. Approximately2.47TiB free SSD space remains, above reserve;
no current workload evidence puts Friday at risk. D2 remains deferred.

October13 20:49 heartbeat: ten N24 jobs remain active with advancing
checkpoints and no failures; no duplicate runner launched. Updated cut/block
checks use1504 identical frozen IDs from snapshot203715. Report
`threshold-cut-block-stability-203715-2026-10-13.json` compares17 common-size cases: max
beta shift .001390(d3 .82), maxd4 .000148(.94); doubled-block SD ratios
.864-1.060. Independent review reproduced these values. Strict leading-fit
acceptance is sensitive at d4 .90(double-block p=.00704) and .94(cut501
p=.00432), so do not claim all checks pass. D4 .92 N24 crossing validity
is .973/.985 atcuts501/2001; common-size comparisons exclude that size.
Pending N24 .146 may help interpolation; explicitly check lower-beta
coverage after closure. These checks do not certify mixing or total error.

New `threshold-curve-snapshot-20261013-203715/finite-size-panels/` contains
per-dimension PNG/PDF/SVG displays of all measured quarter locations,
vertical MC bars, and conditional/rejected leading size fits. Figure source
and provenance are retained; labels and layout were visually checked.
Nine focused threshold/calibration/SSD tests pass. D2 remains deferred;
d5 stays disabled. No material workload or deadline-risk change observed.

October13 22:00 active d3/d4 deliverables completed: `d3-d4-phase-boundary-delivery-20261013-220038/`
contains per-dimension and combined PNG/PDF/SVG curves with vertical bars
only, finite-size panels, numerical CSV, full variant/bootstrap provenance,
and `xy-phase-boundary-d3-d4-scientific-report.md`. Independent review
reproduced17conditional estimates and all19CSVrows. Largest tested
sensitivity halfwidth .01659175; total accuracy coverage is not certified.
Uniform postpilot d4 minimumN10 is disclosed, with all-size and other
variants retained in descriptive bars. D3 .80 andd4 .86 remain unresolved
branch neighbors. The final lower anchor restores100percent N24.92 crossing
validity across all four cut/block snapshots. Cosmetic rerender was checked
numerically identical to the independently reviewed delivery draft.
All1606closures passed audits; source counts are788d3 and728d4, plus90
previously completedd2 records. All316pendingd2 jobs stay deferred, d5disabled.
No active orpending jobs remain. New launches stopped and heartbeat paused
while awaiting the user return to2D. No data, histories or logs removed.

October14 user resumed2D usingeta_eff and second-moment xi/L;
`d2-user-resumption-eta-correlation-20261014.json` supersedes deferral
for neweta/R work only. Old316Binder-quarter jobs remain state=deferred
and blocked_binder_dimensions=[2] guards that adaptive designer.
September archives have matching physicalmodel SHA256
ef507b3cc18bc126e297c6ecf529e27d6b90b64c528d2d02820f12234fc232c8.
Current reuse snapshot in latest-d2-eta-correlation.json uses735selected
coordinates:291refined,6longerlow-kappa,438scoutfallbacks. There are105
conditional finite-size crossings outof110 possibleobservable/size cases.
Targets eta=.25 andxi/L=.7506912 are infinite-size critical limits,
not finite-box criticalvalues. Jointcachedm2/R draws retain covariance.
Localbracket plusadjacentknots avoids censoring roots on remotelyundefined
massive-phase ratios; earlier snapshots/logs are retained. Source flags
and local-versus-full shifts are propagated; max centralshift .000455.
No infinite-volume2D result has been accepted from these crossings.

Frozen pilot54 trajectories spans .60,.80,1.00 at N64/128/256, three
matched beta values and independentfreshreplicas2/3,12000cycles. Matching
N64 anchors make exactbeta64:128 and128:256 eta pairs possible; they
supplement36larger-box trajectories. Full-retention estimate164.55GiB.
Detached coordinator resumed without duplication, SSDguard active,
weekday16/outside18worker limits retained. Keepawake expiry remains
October16 at15:15Central. Hourly heartbeat resumed with updated2Dscope.
Once54closures pass, measure generation plusauditcost, storage, mixing
and crossingcoverage, then report2Dfull-campaign forecast againstFriday
before broad commitment. Completed3D/4Ddeliverables remain intact;
d5disabled. No cleanup, thinning, or environment change performed.

October14 2D eta/R pilot milestone: all54 trajectories closed and passed
full-field audits in12.863minutes, totaling2.47682process worker-hours and
140.859GiB of full retained data. Runtime forecast is
`runtime-forecast-d2-eta-correlation-pilot-2026-10-14.md` with numerical
JSON and independent arithmetic review. Baseline remaining19slices assumes
3beta samples,2starts,12000cycles,N128/256 with optional matchingN64
anchors:228-342newtrajectories, plus24coverageextensions. Planning
allowance2-4h generation and6-18h analysis/refinement fits roughly55h
remaining; these allowances are not a measured completion bound.
Observed storage baseline plusextensions910-955GiB, conservative1115GiB;
uniform4xstatistics would require3991GiB and exceedcurrentfree space.
All existing raw evidence and unfiltered logs remain retained.

The current two-observable reconstruction has direct matched-beta eta
pairs at64:128 and128:256. Lowepsilon.60 responses are shallow and
eta128:256 remains unbracketed; .80N128 R loweranchor is near target.
Dedicated24runs extend .60 tobeta.221/.341/.381 and .80 to.881 at
N64/128/256. No oldBinder-quarter job is thawed. The new joint-moment
bootstrap preserves sharedm2/R samples and sourceflags; exploratoryGLS
uses the full rootcovariance (endpointmaxoffdiagonalcorrelation.533).
Its zero-offsetendpointestimate1.11306 has conditionalMC95
[1.10996,1.11607], missing known1.1199 despite chi2/dof3.235/4.
Offset sensitivity givescentral1.10939-1.12776, so no validated
infinite-volume2D precision is claimed. Fixedoffsetfamilies are exploratory
stress tests; do not force the benchmarkcoupling. Eight focused analytical
and crossing tests pass. Completebrackets and sensitivity checks before
a broad grid commitment; optional5D remains disabled.

October14 user addedbeta0 endpointmeasurements forD2/D3/D4, tobe
shown with horizontalerrorbars inepsilon/pi ratherthanverticalbeta bars.
Frozen `zero-beta-endpoint-pilot-20261014.json` contains64trajectories
(18D2,16D3,30D4), conservative217.20GiB fullretention. D2useseta/R,
D3/D4Binderquarter. Dedicated endpointrole and independentseednamespace
2026101403 preventoffgrid epsilon aliasing intoexistingpositivebeta
curves; analysis/adaptivefilters explicitlyskipthisrole. Existingzero
couplingdata andactualradian-grid epsilons remainimmutableinputs.
Beta0vertical coordinateisexact; endpointMC/modelsensitivity is horizontal
in epsilon/pi units. Practicalprecisionisreportedratherthanconverting
thebeta .02target intoan unspecifiedepsilon target. Independentreview
verifiedrole,points,seeds,storage andconditionalthermalcoordinate mapping.
Earlier3D/4Ddeliveriesstaypreserved; updatedversionsfollowafterendpoint
analysis. Optional5Dremainsdisabled. Old316D2Binderquarterjobsremain
deferred, whileauthorizedeta/R andendpointwork continues.

October 14 endpoint bracket refinement: first snapshot showed incomplete D2 bootstrap tails, a missing D4 N10 lower bracket, and coarse D3 epsilon interpolation. Frozen additional 40 trajectories (12 D2, 8 D3, 20 D4), estimated 141.99 GiB, use fresh starts and exact epsilon coordinates. Dedicated analyzer and source-preserving snapshots are available via latest-zero-beta-endpoints.json. Independent review confirmed orientation/covariance/scaling and identified two fixes now applied: decimal-safe cache names and common masks after size selection. Eleven relevant tests pass. Endpoint infinite-volume estimates remain preliminary pending completed brackets, historical reuse and sensitivity checks.

October 14 08:47 Central: independently reviewed remaining 2D eta/R coverage batch frozen: 290 new runs over 19 non-pilot epsilon slices, 3 beta samples per slice, 12000 cycles, fresh starts2/3. Reuses26 clean multistart N64 coordinates; fresh counts N64=62,N128=114,N256=114. Conservative1016.02GiB; measured pilot projection15.23worker-hours and871.89GiB actual storage. Current storage budget leaves approximately844GiB beyond25GiB reserve after remaining endpoint runs and this conservative batch. Finite pilot root coverage passes >=.99; infinite-volume accuracy remains conditional on BKT size/offset/observable checks. This extends local coverage, not a precision claim. No model/runner/environment change, old316Binder jobs remain deferred.

October14 09:40: endpoint pilot64/64 audits passed; six refinement jobs still pending. Legacy zero-beta adapter integrates64 archived coordinate groups at exactepsilon, preserving joint D2 bootstrap draws and historicalD3 radians/pi. Duplicatecoordinate policy prefers current pooled starts; archived Binder first2000of4000draws used without modifying archives. Legacy included only defaultcut2001/block1; alternative cut/block snapshots are current-only. All endpoint extrapolations remain preliminary. Completed low-epsilon grid slices require24 bracket extensions (q.58 beta.020/.280,q.62 .300,q.66 .640,matchingN64/128/256,freshreplicas4/5,12000cycles,73.13GiBconservative). Incompleteq.92 tails are not yet refined because its neighboring grid runs were still pending in the analysis snapshot. Tests11pass.

October14 10:08 user-requested preliminary phase-boundary plots now includeD2/D3/D4 beta-zero endpoints. Pointer latest-preliminary-all-dimensions.json;42points,12graphicsPNG/PDF/SVG,andnumericalCSV/JSON. All290baseline plus24targetedD2 bracket runs and104endpoint runs closed with passedaudits. D2minN32 fixedoffset0 sharedeta/R GLS is exploratory; sensitivityunion usesminN32/64,offsets0/1/2/4,joint/separateobservables andincludespoorexploratoryfits. EndpointprimaryD3min12,D4min10,D2min32C0joint. Diamondendpoint hashorizontalepsilon/pi bars only; everypositive-beta point hasverticalbars only. DarkconditionalMC95 versuslightdescriptivesensitivityunion lacks totalcoveragecertification. Final scientificacceptance remainsoutstanding, especiallyD2 log/size checks andendpointlocal-window/start/cut/block sensitivity; no finaldelivery orD5enablement implied. Graphicsvisuallychecked and42interval/axisassignments verified.

October14 user plot revision: show unconnected MC extrapolated points plus leading largeD prediction, (1-epsilon/pi)+beta/2=1/(2D), hence beta=2(epsilon/pi-1)+1/D onpositivebeta branch. Newpointerlatest-MC-vs-large-D.json; original42MCvalues/errorbarsunchanged,12PNG/PDF/SVGfiles. Preserve this preference for updatedfinalfigures. Endpointhorizontalbarsremain. WideD2endpointpaleunion isdrivenbyextrapolationchoiceandamplifiedMCvariance ineta-onlyfits; itisnotpurecentralshift orcalibratedtotalCI. PrimaryMCwidth.007785versusunionwidth.042587; centralvariantspan.008685. NoMCconnectinglines remain; theoreticalendpointsq=.75,.833333,.875 andbeta(q=1)=1/D arechecked.

October14 10:50 validation phase: seven endpoint sensitivity snapshots test local3SD versusfullscan, current-only cuts501/2001/3001, doubleblock, and Binderthresholds.20/.25/.30. Endpoint sensitivity summary uses common current roots at identicalminimumsizes andmodels: maxcut/block epsilon/pi shift D2.001950,D3.0000608,D4.0001111. Earlier cut4001 attempt encounteredemptyD4arrays; failure/source/log preserved, retentionguardadded,cut3001 replaces it. No fieldaudit failed. ThresholdshiftsD3/4below.001 are descriptive, not totalprecision certification. D2positivecurrent-chain variations hold archivedanchorsfixed; commonrootC0GLS maxbetashift.011379(q.58),.010186(q.60). All11focusedtests pass. Further scientificreview requested. Sixrun N512D2cost/observablepilot (162.77GiBconservative) independentlyreviewedandfrozen atq.60beta.281,q.80beta.921,q1beta1.111;matchesexistingN256, freshstarts6/7,12000cycles. Measureactualcost/mixingbeforemoreN512; a singlebeta perslice cannotlocatecrossings. FinalD2andendpointdeliverablesremainincomplete,D5disabled,largeDpredictionandunconnectedMCplotpreferencepreserved.

Independent sensitivity review reproduced headlines and requested explicit incompleteness labeling: q.66 latecut anddoubleblock retain.9855/.9895 completejointrootdraws. Summarynowreports successfulvariantcounts andincompletecomparisons, storesendpointinputflags/jointvalidfractions, andguardsemptycommonrootsets. New6run matchedlowerbeta.500 anchor atq.66 Ns64/128/256,freshstarts6/7,12000cycles,18.28GiB frozen toextendmeasuredtails. Reanalysebaselineandthreecurrentchainvariants aftertheseclose; do not silentlyconditionaway invaliddraws. Reviewed reportpointers latest-endpoint-sensitivity.json andlatest-d2-current-chain-sensitivity.json supersedethefirstsummarycopies, whichremainpreserved.

October14 USER AUTHORIZED retention change for existing/future campaign runs, superseding previous no-thinning/fullfield retention rules. retain_fields.py and detachedretention.plist service preserveall observable/winding/RNG/provenance histories andexact100cyclecheckpoints. Fullfieldclosureaudit+originalSHA precede byteverifiedsparsechunkconversion, NaNfillforunsavedfields, retainindicesexplicit. Fourdense diagnosticwindows(max64,50tau), remainingstrideceil10maxmeasuredtau acrosspostwarmup/latehalf relevantseries; noindependenceguarantee. Durablepreparedjournal+directoryfsync+atomicreplace+recovery event, actualtrajectorylocks,25GiBtemporaryguardincludingactivejobreservation. Nativecheckpointcontinuation/histories/recovery tests11pass. Archiveimmutable; continuationfromworkingcopy, changedarchivefailsreconciliation. FirstrealN512 file23.61->1.78GiB(~92%reduction). Servicecountsonecore, productioncap15day/17otherwisewhileconversionbusy; futurefullgeneration+fullaudit stilloccurbeforeconversion. LatestN512 sixpilotcomplete, R512 andeta256:512 nearlimitingtargets atonebeta; nofullN512crossingsyet. TwelveextraN512bracket runs arefrozenmatchingoldN256betaanchors,325.55GiBtransientfullupperbudget. N1024notyetjustified/frozen; decideafterN512brackets andsize/offset/benchmarkvalidation. See retention-records/events/logs, latestD2snapshotd2-eta-correlation-existing-20261014-113051.

October14 11:46: reanalysed all current-chain cut/block variants after the q.66 loweranchor. All22 slices now have allfour identifiable/common-valid variantfits; q.66 baselinebeta∞.574883 andmaxshift.00738347. Globalmaxshift stays.0113788atq.58. These remain correlated sensitivity checks, not bias/precision certification. Sourcepointer latest-d2-current-chain-sensitivity.json nowreferences114555summary; previousreports preserved. Adaptive retention tests now12pass, including exactnativecontinuation at interiorcheckpoint100 and final600 fromworkingcopies. SixN512pilotarchives converted,130.99GiBfreed; old4Darchiveconversion and12newN512bracketruns continue. Useraskedwhetherthinningalreadyactive; confirmeddetachedservice, notsubagent.

October14 12:12 shared storage admission fix: user correctly identified independent generation/retention space-check race. Currentbatch ample, butfutureadmissionnowserializesunderstorage-admission.lock. Retention reserves unfinishedtemporarycopy, coordinator snapshotsactivejob growth beforefreshphysicalfree and avoidsdouble-counting thosejobs in frozen schedulingloop. Newpendingjobs reservefull footprints; savingscountonlyafterreplacement. Sharedtemporary/globalcopy locks, staleleasecleanupunderlocks, durableleaseupdates.18tests passed includinginterprocessserialization andsimulatedactive/temp growth; independentreviewno blocker. Two wrappers reloaded whilelaunchgated andactive0; newPID69566 publishesphysical_disk_free_bytes,retention_reserved_bytes,active_reserved_bytes separately. Launchgate reopened, no numericalkernel/run.py changes. All12N512bracket runs nowcomplete; latestD2snapshot includesN512R andeta256:512 crossings. FurtherFSSreviewpending, noN1024batch.

October14 user corrected comparison theory using constrained_xy_large_d_note.pdf: page5orangebox (z-1)v1=1,page1z=2D,page2v1 normalizedintegral exp(beta*cosθ)cosθ over[-epsilon,epsilon]. The earlierlinearformula was anaccidentalapproximation; finalcomparisonfiguresmustusethefull treecriterion. New tree_level_prediction.py usesGL128quadrature+brent,independentquad/sinc/Bessel/residual tests6pass. latest-MC-vs-tree-level.json points to122036plot set; latestMCdatarefreshedincludescompletedN512crossings. Exacttreeendpoint epsilon/pi=.7253845,.8262494,.8720063; beta(epsilon=pi)=.7075406,.4082772,.2886803. These aretree-levelpredictions, notexactfiniteDhypercubictransitionclaims. TwelvePNG/PDF/SVGgraphics,CSV/JSONandoriginalPDFSHA retained. PDFcombinedrender visuallychecked, MCpointsremainunconnectedandbeta0bars horizontalonly. Priorplots preserved.

October14 12:49 size-validation phase: validate_d2_heldout_size.py trains onlysmallerboxes andpredictsR512/eta256 roots, preserving sharedsource bootstrap covariance. C0/min32 atq.80 underpredicts by.008849/.008634,jointnominalp.02038; higheroffset/minimumsize variantsremainplausible. Descriptiveheldoutcheck, notconfirmatorysignificance/modelselection orproofN1024required. Independentreview reproducedandapproved limitedcostprobes. New8run batch d2-eta-correlation-targeted-N1024-cost-and-N512-endpoint-probes-20261014: twoN1024 q.80beta.921, plus sixN512beta0epsilon/pi.553/.558/.565,matchingcompletedN512/N256respectively,fresh8/9,12000cycles. Conservativefulltransient378.97GiB (216.20+162.77), about509GiBheadroombeyond25GiBreserveandcurrentretentionbudget atdesign. N1024forecastallowance1.5–3.5h,notactualcost; measurebeforefurthercommitments. NoN1024crossingfromonebeta; endpointcoveragecheckedafteraudit. Fiveprediction/GLStests pass. Adaptivefieldretentionandsharedadmissionactive;29convertedarchiveshadfreed444.65GiB. D5remainsdisabled; finalprioritydeliverablesstillpendingBKT/endpointreview.


## Outside-bracket field reduction authorized October 14, 2026

The user explicitly authorized removing non-checkpoint configurations for runs clearly outside transition brackets, while preserving all observable/RNG/winding/timing histories, logs, provenance and exact checkpoints. This supersedes historical full-field and dense-diagnostic requirements for the eligible subset only. Numerical run.py and model kernels remain frozen.

The retention daemon refreshes a frozen, source-hashed eligibility manifest when analyzed pointers or completed audited jobs change. For every contributing size and observable (including BOTH sizes of eta pairs), it requires identifiable >=99% bootstrap brackets in all saved cut/block/threshold/interpolation comparisons. It protects the union of measured local windows, one measured neighbor on each side, and a full window-width margin. Unresolved curves, unknown coordinates, and material start/warmup/pooling flags are protected. Low-block-count flags alone do not remove eligibility; all histories and exact checkpoints remain complete. Endpoints are checked separately in epsilon/pi at exactly beta zero, with exact coordinates. The current endpoint comparisons protect their entire measured scans, so no endpoints qualify initially.

Eligible archives retain phi at every 100-cycle checkpoint plus first and final rows. Missing phi rows are explicit unallocated NaNs. All other datasets are compared byte-for-byte; existing full-audit provenance and any previous compaction record are preserved. A second reduction verifies the current compact SHA against the prior committed ledger. Durable prepared records distinguish crashes before replacement (restore prior ledger, preserve abandoned attempt) from after replacement (finish journal commit). The existing shared disk-admission and trajectory locks apply.

Initial stricter manifest: 348 eligible completed jobs (138 D3,210 D4),341.19 GiB current physical files at manifest creation. This is candidate storage, not projected or realized savings. D2 and endpoints remain protected at this stage. Eligibility refresh for future completed runs is automatic, but no run is classified solely by distance from an infinite-volume central estimate. Immutable hash-named manifests preserve every eligibility generation. No historical archive outside this campaign is modified.

Validation: 14 retention, checkpoint/replay, failure/recovery, shared admission and SSD routing checks passed. Test log: logs/outside-bracket-retention-final-tests-20261014.log.


### October 14 13:48: N512 endpoint coverage and tested-analysis update

Six N512 beta-zero probes closed with full-field audits passed. Generation took1073–1083s and audit145–147s each. The local-window fallback now retains an additional neighboring knot when a3SDanchor is unavailable;12focused analysis tests pass and independent review verified the sign convention and coverage correction. Earlier snapshots remain intact. The largest roots remain unresolved: R512epsilon/pi .557984 with98.25% valid draws, eta256:512 .557141 with98.15%. Their current-chain cut/block variants retain93.7–98.7%; they are excluded from accepted common-root fits. This is not evidence of attained infinite-volume endpoint accuracy.

Four reviewed wider anchors at epsilon/pi .546/.572,N512,replicas10/11,12000cycles are frozen in d2-eta-correlation-N512-endpoint-bracket-extensions-20261014. MatchingN256 data already exist. Conservative transientfull storage108.52GiB, with~880GiBheadroom beyond reservations and25GiBreserve at freeze. Do not duplicate; reanalyse after audited closure before further size decisions.

Updated positiveD2 current-chain sensitivity includes completedN512 brackets:22slice comparisons, four successful common-root variants each, maxbeta shift.011379 atq.58. These are not all-root coverage checks: q.60 latecut excludes R256 andeta128:256. Diagnostic use of the full existing measured beta window restores99.9% and99.75% draws respectively (q060-latecut-full-window-diagnostic-20261014.json). This full-window check is not promoted and does not increase precision; inspect window selection before authorizing newq.60 sampling. Theq1 baseline1.114958,MC95[1.112671,1.117301] still misses1.1199 despite small cut/block shifts; preserve model-dependence qualification.

All348 initial outside-bracket reductions completed, saving323.58GiB additionally; total cumulative retention saving889.12GiB at13:43. This milestone was reported. Existing twoN1024probes continue; no furtherN1024 orD5 commitment.


## October 14: practical N512 stopping point and active 5D priority

The user explicitly capped all further 2D simulations at N=512 and requested attained results, then 5D. The two completed N1024 probes remain preserved and excluded from the primary curve. No additional 2D precision campaign is planned. `control.json` carries d2_max_N=512; the storage-aware coordinator blocks pending larger 2D launches without changing run.py or kernels. Seven routing/admission checks pass. All316oldBinderD2jobs remain deferred.

The independently reviewed attained delivery is `d2-d3-d4-attained-N512-with-endpoints-20261014`, including per-dimension and combined PNG/PDF/SVG, CSV, unresolved-locations.csv and attained-phase-boundary-report.md. N512endpoint roots retain99.95–100% draws in the tested cut/block variants. This is conditional coverage, not certified total accuracy. The q1 BKT benchmark discrepancy and wider descriptive unions remain explicit.

5D is enabled and active. All54initial N4/6/8,4000cycle calibration trajectories passed audits:0.9509worker-hours in4.912minutes,11.61effective workers,38.35GiB full audited files. Coarse interpolation is not accepted as final precision. The independently reviewed54-trajectory refinement (27physical coordinates) at N6/8/10,8000cycles is running in d5-local-bracket-and-N10-refinement-20261014. It narrows two positive-beta slices and the beta-zero endpoint and budgets271.13GiB full transient storage. No duplicate batches.

Use analyze_d5_calibration.py and dedicated d5_scaling.py. The leading periodic zero-mode hypothesis is N^-5/2 with N^-3 stress correction, not the4D logarithmic coordinate. Three-size correction fits havezeroDOF and are diagnostic. Four analytic/covariance/zero-mode tests pass. Sources: Chen–Dohm cond-mat/9809394 and Wittmann–Young1410.5296.

The measured pilot supports a planning allowance12–29h for qualified attained5D estimates versus about48h toFriday, not a runtime bound or ±.02 guarantee. The 168-trajectory baseline corresponds to84physical coordinates; after54trajectory reuse,114trajectories(57coordinates) remain. Conservative remaining fullstorage690.28GiB is marginally larger than post-refinement headroom before further actual recovery. Check actualN10 runtime, coverage, retained bytes and fresh reservation-aware space before freezing the full grid. Uniform8kremaininggrid needs1378.62GiB and does not fit current uncompactedcapacity. Current54run stage fits comfortably. See runtime-forecast-d5-calibration-20261014.md/json.

Preserve all histories/checkpoints/logs and numerical kernels. User-authorized adaptive and outside-bracket field retention continue under shared admission. Use atomic JSON updates for control and live pointers. Prioritize the remaining D5 slices on exact.02epsilon/pi increments and its endpoint, with honest uncertainty/flags, then updated combined unconnected MC points against the full tree criterion. Do not restore interpolation lines. AtFridaydeadline stopnew launches,checkpoint safely,deliverattained/unresolvedresults andpause; earlierpause onlyafterD5deliverablesready.

### October 14 15:03: D5 measured-window correction

The first narrower q=.92 data invalidate treating the coarse pilot's linear root as a precision result. At N6/N8 the quarter crossing is between beta .060 and .070 (linear locations .064510/.064330); all six audited N10 samples at .070/.090/.110 are above U=.25, with U(.070)=.493305. The independently reviewed frozen d5-q092-lower-bracket-refinement-20261014 adds 18 trajectories (nine coordinates) at beta .060/.064/.068, N6/8/10, 8000 cycles, fresh replicas6/7. Its90.38GiB full-transient budget left534GiB beyond recorded active/copy reservations and25GiB reserve at freeze. This is a bracket correction, not mixing certification; the N8 beta .070 start disagreement remains recorded. No N12 q=.92 commitment until the lower bracket is measured.

analyze_d5_calibration.py now accepts cut/factor/target/window/start-parity variants, --no-promote and --source-analysis to pin an immutable analyzed source-job cohort. It records local/full windows and shifts, selected source IDs, unresolved omitted sizes, and first/last retained-half Binder values with measured-width influence-function block errors. Four-SE half flags are descriptive only. Even/odd-start draws share helper random seeds and must not be paired for start-difference inference. Compare cuts1001/2001/3001 on common sources; cut4001 on old4000cycle chains changes composition and must be separately labeled. Eight focused tests pass in logs/d5-local-analysis-tests-20261014.log; independent review found no implementation blocker. Local windows are selected on central data and bootstrap coverage remains conditional on that selection. All original snapshots, histories, checkpoints and full logs remain intact.

### October 14 15:16: physical HDF5 chunk accounting correction

All54 D5refinements passed. Actual352.43GiB exceeded the old271.13GiB estimate because N10 autochunk padding expands four phi axes10to12. The new storage_admission.estimated_full_bytes counts padded phi/winding allocation plus margin and is shared by coordinator/retention growth reservations. N10,8k corrected17.824GiB ceiling exceeds observed15.472GiB. Seventeen retention/admission/routing checks pass; independent review verified geometry and no estimator recursion. Retention queue now rescans after each successful copy to see newly closed large archives, preserving all eligibility/exactness rules. See d5-storage-forecast-chunk-correction-20261014.md/json. The114-trajectory structural remainder is now793.68GiB, not690.28GiB; do not freeze it wholesale or assume future thinning. Numerical sources still match frozen hashes. Coordinator reloaded toPID95610 with all six detached lower-bracket workers preserved; launchgate temporarily closed while retention reload waits for its current verified copy. Reopen only after the single refreshed service is verified.

### October 14 15:21: lower bracket complete and N12 stage running

All18 lower-bracket trajectories passed audits. The126-source primary snapshot152036 has100% valid quarter-root draws at q=.92,N6/8/10: beta .064611/.065290/.065066. This establishes coverage of that measured window, not total extrapolation accuracy. The q1 minN6 leading-size fit has p~1e-17 and is now explicitly rejected_size_fit, retaining its diagnostic values. The cause is not isolated: the new midpoint tests interpolation while N12 tests additional size leverage. Nine D5 analysis/oracle tests pass, including explicit rejection of an inconsistent size ansatz.

The reviewed14-run d5-N12-q1-endpoint-cost-and-coverage-20261014 is now running: twelveN12,4000cycle q1(.224/.228/.232) and beta0epsilon/pi(.894/.896/.898) probes, plus twoN10,8000cycle q1beta.228 anchors. Corrected full-transient ceiling183.823GiB is written into the new frozen manifest, and its refreshed design left333GiB beyond previous live/pending/copy/25GiBreserve at freeze. No N12q92 or remainingfourpositive-epsilon batches are frozen. Inspect actualN12cost/coverage/mixing/storage before subsequent commitments. NewcoordinatorPID95610 andretentionPID95868 were verified as single instances, currentlaunchgate reopened; allworkerssurvivedreload and no failures were observed.

The original forecast_d5_calibration.py retains historical hardcoded old storage prose; do not rerun it unchanged. Use d5-storage-forecast-chunk-correction-20261014.md/json for the current793.68GiB structural remainder and staged-budget limitations. Retention is processing the completedN1024probes first by physical size; those remain preserved but excluded from the primary2Dcurve. Advance remainingD5slices(.90/.94/.96/.98) after the current stage; finalize local/full/interpolation/threshold/cut/block/start/size comparisons on pinned sources and update combined MC/tree figures. Neither favorable fit p nor bracket completion is total-coverage certification.

### October 14 16:08: measured N12 cost, remaining grid, and targeted local refinements

All14 N12/midpoint trajectories passed full-field audits; frozen numerical hashes still match. N12,4000cycle generation took834–1134s plus119–146s audit, producing128.61GiB across12 files. The twoN10,8000cycle midpoint files generated in876–882s plus130–131s audit. See d5-N12-cost-and-coverage-update-20261014.md and d5-N12-cost-and-next-stage-forecast-20261014.json. N12q1 has a conditionalquarterroot.229562 betweenbeta.228/.232; its lowblock flags remain. N12beta0 atq.894/.896/.898 is unresolved because allpooledU lie below.25. Do not treat omission of that size as endpoint acceptance.

The independently reviewed d5-remaining-low-boxes-and-N12-bracket-coverage-20261014 is frozen with78 new trajectories/39coords, physical ceiling400.322GiB. It comprises66positive runs onq.90/.94/.96/.98 atN6/8/10 plus6N12q92 and6N12lower-epsilon endpoint anchors. It reuses10 existingauditedbeta0 starts on3exactq.90coordinates rather than launching6 redundant runs. Analyzer exact_zero_beta_fixed_scan_anchor accepts only beta0, declaredq=.90 andmilli900; allother endpointcoordinates stayexcluded and actualradians/roles preserved. Sharedsources acrossaxes must not be treated as independent evidence. Ten focused tests pass in logs/d5-exact-anchor-reuse-analysis-tests-20261014.log.

Two concrete stages were also independently reviewed and frozen: d5-q1-local-curvature-check-20261014 adds4trajectories/2coords atq1,beta.230,N10/12,4000cycles,fresh12/13 (42.777GiB); d5-zero-beta-local-window-refinement-20261014 adds12trajectories/6coords atbeta0,epsilon/pi.893/.894,N6/8/10,4000cycles,fresh14/15 (45.845GiB). These target specific existing interpolation limitations; they do not establish why the rejected q1sizefit failed or certify mixing. Including the78run stage,94newtrajectories are committed; do notduplicate them.

Fresh physical-budget snapshots include all prior pending jobs, livegrowth, temporarycopy and25GiBreserve. The12endpointdesign left326.8GiB after allknownwork/reserve. The optional/uncommitted24N12trajectories for four remainingq slices must wait for observed coverage and another freshbudgetcheck. The forecast JSON's58GiB extra-N12scenario margin predates the12endpointaddition; do not use it as a current launchauthorization. At the recorded recovery check425committedarchives had freed1410.54GiB; all425current sizes matchedledgers, with145.75GiB recoveredfrom11D5archives. This is actualrecovery, notassumedfuturethinning. Allhistories/checkpoints/logs remain complete. Measured N12cost remainscompatiblewith the Friday planning scenario, without guaranteedtotalprecision. Continue source-pinned local/full/PCHIP/threshold/cut/block/start/size review and retain rejected/zeroDOF/unresolved qualifications in finalD5andcombinedMC/tree deliverables. D2/D3/D4attaineddelivery remainsfrozen; no newD2precisionwork.

### October 14 17:08: final matched N12 stage and detached sensitivity analysis

All234 preceding D5 trajectories passed full-field audits. Primary snapshot165416 has27/27 measured linear quarter crossings with100% valid draws; this establishes bracket coverage only. The q1 minN6 leading fit is rejected, and minN8 remains exploratory. Low-block flags remain. Two-knot local windows make linear and PCHIP identical and do not test interpolation bias.

The independently reviewed d5-final-N12-and-matched-midpoints-20261014 is frozen/running:48 trajectories at24 coordinates,4000cycles,fresh16/17. At q=.90/.94/.96/.98 the N12 windows are respectively .013/.019/.026, .101/.107/.113, .144/.150/.156, .183/.189/.196; matching middle values are measured atN6/8/10. Corrected full-field ceiling388.037GiB left178.27GiB beyond all prior commitments, copy reservations and25GiB reserve at freeze. Do not duplicate. No further uniform production stage is planned if these attained results suffice.

Exactly one finite detached analysis service org.luscherxy.d5-final-sensitivity-20261014 was verified running, PID11601. run_d5_final_sensitivity.py --wait --parallel3 waits for all282 frozenD5 trajectories and passed audits, then runs a primary plus eight source-pinned variants (full window,cut2001,cut3001,block2,U20,U30,even,odd). State is in d5-final-sensitivity-manager-status.json; full log is logs/d5-final-sensitivity-daemon-20261014.log. It survives Codex app restart. Do not duplicate manual analyses while active. Snapshots_complete_pending_scientific_review does not mean final deliverables are complete. Cohort-specific atomic bootstrap caches prevent interrupted cache publication;11analysis/oracle and3manager guard tests pass. Independent review found no blocker. Missing sizes, rejected and zeroDOF fits remain explicit incomplete/diagnostic comparisons. Final D5/combined plots,tables,report and independent scientific review still follow.

### October 14 18:00: all 282 D5 simulations audited; final analysis and draft assembly

All 282 frozen D5 trajectories passed full-field audits, with no failed jobs. Primary snapshot173811 has 31/31 linear finite-size roots retaining 100% bootstrap draws. Five positive-beta slices and the beta-zero endpoint have conditional minN6 leading fits; q1 minN6 leading and correction fits reject and remain diagnostic. The uniform baseline will not be switched to minN8 merely for a better p-value. Draft plots use an open square for the rejected q1 extrapolation, never an accepted coupling. Descriptive pale stress envelopes retain identifiable positive-DOF rejected/flagged variants, with explicit exclusions for missing sizes or unidentifiable cases. Local two-knot linear/PCHIP agreement does not validate interpolation; the matched midpoint update moved finite roots by up to .003584 and also changes sources.

The original single sensitivity manager is actively computing the cut/block variants, then thresholds and start parities. A separate finite, lock-protected launchd service org.luscherxy.d5-delivery-builder-20261014 waits for all nine successful snapshots and creates a pending-review folder through build_attained_d5_delivery.py. It does not mark scientific completion. Check logs/d5-delivery-builder-daemon-20261014.log, latest-draft-d5-delivery.json and manager status. Six combined manager/reporting guard tests pass. The builder copies all 42 prior rows unchanged, preserves missing-family/coverage records and flags, and requires intervals to contain plotted centers. No new production batch is planned. Final visual/scientific review, delivery promotion, keep-awake stop and heartbeat pause remain pending.

### October 14 18:28: reviewed attained D5 and combined delivery complete

All 282 D5 simulation trajectories and nine source-pinned sensitivity snapshots completed. The reviewed delivery is d2-d3-d4-d5-attained-with-tree-20261014-182254, referenced by latest-attained-d2-d3-d4-d5-delivery.json and latest-MC-vs-tree-level.json. It preserves all 42 prior D2/D3/D4 rows exactly, adding six conditional D5 extrapolations and one explicitly rejected q1 diagnostic (49 plotted points). Independent review reproduced rows/envelopes, checked hashes/CSV, all28 larger-size linear roots across9variants (minimumvalid .998), error orientations and q1 open-square labeling. Largest descriptive halfwidth .005382 is not certified totalaccuracy. Variant half-history/start/low-block flags remain material and preclude a mixing PASS. See report, full stress ledger, complete coverage and independent review records in the delivery folder. No further production is necessary for this attained scope.

Numerical checks finished18:17. Initial export stalled reopening protected Downloads PDF for a redundant sourcehash; the blocked figure process was terminated, partial folder181718/full log preserved. Corrected export182254 reuses the already verified immutable sourcehash. No simulation or numerical-analysis failure. Combined/D5/size-diagnostic PDFs rendered and viewed; size labels corrected. control.launch_new=false; deferred jobs preserved. Stop temporary keep-awake and pause heartbeat now that priority deliverables are reviewed. Authorized field retention continues independently with all histories/checkpoints/logs intact.

Completion actions verified: control.launch_new=false, finite analysis/export processes exited, campaign keep-awake service absent, heartbeat PAUSED. Retention continues independently under the authorized evidence-preserving policy. Six final reporting/manager guard tests pass.

### October 14: user-authorized parallel retention, four-core headroom

The user requested all but four cores for thinning. Hardware reports6Super+12Performance,18total, with no efficiency cluster; use up to14 single-thread worker processes, leaving four cores worth of capacity to the OS. Exact core placement remains managed by macOS. No simulation launches are enabled and6D remains on hold. Field selection, all histories/checkpoints/logs, source/compact hashes, byte verification and crash recovery are unchanged.

The retention manager now owns a spawned14-process pool with unique in-flight archive IDs. Workers hold shared retention-copy locks and exclusive archive locks; exclusive maintenance waits for all active copies. Temporary growth has per-archive leases, summed with all other copies and live simulation growth under the shared admission lock before fresh free-space is read.25GiB reserve remains protected. Temporary archives are created only after admission with unique attempt IDs. Ordinary/recovery event appends are serialized. Parent-only worker accounting prevents child completion from restoring simulation capacity prematurely.24 concurrency/preservation/routing tests pass, including simultaneous real conversions, same-file exclusion, post-replacement recovery, exact checkpoint continuation and near-limit admission with writer growth. See parallel-retention-deployment-report-20261014.json and logs/parallel-retention-final-tests-20261014.log.

Reload waits for the current serial copy to finish under exclusive maintenance lock, verifies the old daemon has exited, then starts exactly one replacement retention service with --workers14. The unused idle simulation coordinator is stopped so its old single-lease accounting cannot survive. For future retention reloads, write retention-stop-request.json, let the parent drain all workers, verify exit, then reload; exclusive maintenance remains available for transaction-safe inspection. Source snapshots are retained under retention-implementation-snapshots. Do not claim14x speedup; inspect retention-status.json and process/SSD measurements after deployment.

Parallel rollout verified at19:13Central: old serial converter finished normally; one replacement daemonPID28116 and14 worker processes active. Snapshot aggregateCPU about13.3cores andRSS about5GiB on64GiBhost; no speedup claim. Fourteencopy leases reserved73.0GiB remaining, physicalfree484.5GiB, leaving386.4GiB beyondleases+25GiBreserve. No active/pending simulations; idle simulationcoordinator stopped. Scientific heartbeat remainsPAUSED and6Donhold. See deploymentreport/log and retention-status.json for live progress.

### October 14 21:18: temporary-copy accounting race fixed and retention restarted

The14-worker converter completed139archives, then drained and stopped19:37 after a FileNotFoundError while admitting d2_e074_b0809_N256_r03_s012000. Status checked at21:15 revealed the halt. The named archive remains byte-identical to its original full-auditSHA; all626committed current file sizes match their ledgers. Old logs preserved exceptionrepr only, so the exact historical line is unavailable. A reproducible race was found in peers’ temporary-file exists()+stat(): another worker can complete an atomic rename outside the admission lock between those calls. The repair uses one stat and treats FileNotFound as zero bytes written, conservatively reserving the full lease until its owner releases it. No storage safety check or data-retention policy is weakened. Future failures save full remote tracebacks.

Twenty-five preservation, concurrentcopy/recovery, admission and routing tests pass, including a regression reproducing the stale-positive-exists observation. The existing inactive14-worker service was kickstarted after verifying no old retention workers remained. No new simulations;6D remains on hold and scientific heartbeat remains paused. See parallel-retention-temp-rename-race-recovery-20261014.json and logs/parallel-retention-temp-rename-race-tests-20261014.log. Earlier2–4h completion estimate must account for the roughly100minute stoppage; recompute from resumed throughput rather than counting halted walltime as productive.

### October 15: authorized finer 3D beta-zero spin-dimension scan

The user requested a finer epsilon scan to resolve interpolation in the 3D beta=0 spin-dimension measurement. The proposed two-CFT motivation was explicitly corrected by the user: this O(2) study compares with the standard Wilson-Fisher O(2) CFT, without fixing Delta or nu in the fit. The completed D2/D3/D4/D5 phase-boundary delivery stays fixed. No D2, D5, D6 or additional N1024 work is authorized by this new scan. All316 old BinderD2 jobs remain deferred. SSD originals must remain; the archive-deletion authorization is still revoked.

Reviewed frozen batch d3-zero-beta-spin-fine-epsilon-20261015-stage1 contains82 new trajectories at41 exact coordinates, sizes8/12/16/24/32,12000cycles,two fresh RNG starts with explicit zero/unit-x winding. Coordinate and seed encoding uses integer1e-5 units; role=zero_beta_endpoint excludes these from rounded positive-beta grouping. Full physical storage ceiling96.44GB leaves2.56TB beyond25GiB reserve and known commitments, with no hypothetical thinning credit. Design, reviews, logs and live state are in ../2026-10-15-xy-3d-spin-dimension/.

One finite detached service org.luscherxy.d3-spin-fine-20261015 runs run_fine_scan.py, which uses the existing single global coordinator and immutable simulation kernels. It waits for full closure audits, closes the launch gate, stops the idle coordinator and runs fine-only, source-pinned cut/block/threshold/parity/window variants plus a separate pooled diagnostic. fine-scan-status.json distinguishes simulations, analyses, failures and snapshots_complete_pending_scientific_review. Do not duplicate this service or launch its82 jobs again. It survives Codex app restarts. Friday15Central remains the safety launch deadline for this initial finite stage.

The existing retention daemon was drained through retention-stop-request.json during this small finite scan so its automatic capacity writer cannot override four-core headroom. Its PID34826 exited cleanly at16:58 with no active copies/failures. max_workers14 was restored after exit; the finite scan manager performs the already-authorized adaptive field retention after all new trajectories audit. All histories, RNG, checkpoints, provenance and full logs stay intact. Current storage/admission guards remain active. Production and analysis outputs do not establish mixing, interpolation accuracy or CFT precision until scientific review.

### October 15 17:15: finer beta-zero 3D spin scan attained and reviewed

All82 new trajectories passed full closure audits (1.55934worker-hours generation+audit,79.615GB originalfullfiles; about9minutes elapsed). New analysis independently reviewed at ../2026-10-15-xy-3d-spin-dimension/zero-beta/latest-fine-delivery.json, folderfine-delivery-20261015-171155. It pins15 source-preserving snapshots and a correlated same-source grid-decimation diagnostic. Conditional baseline moment-consistent,U=.25,L8–32 gives Delta_eff=.505400±.007578stat (one bootstrapSD), MC95[.489177,.518262]. SYSTEMATIC ERROR IS UNCONTROLLED: sizecutsL12/L16 shiftcentralto.530758/.550138, weakcorrectionmodelspermitlargerchanges. Do notinvent±.04sys,claimtwo reliableCFTdigits,orclaimWFexclusionfromthebaselineMCinterval. Interpolationcentersnowspan.5054–.5081 descriptively; decimationparity1same-size shift−.00903 hasshared-drawMC95[−.02256,.00773],notestablishedbias/bound. U20/U30stresscoverageisstillincomplete.

User asked what controlled systematic/two-digit precision would take. Recommendation only: reviewed L48/L64 pilot withfineepsilonbrackets,independentstarts,cost/mixing/sizecut/thresholdchecks (~1–3h initialconditionalallowance),thenhigherstatistics/largersizesonlyifwarranted. A1–3dayprecisionATTEMPTisnotapromiseof±.005combinedcoverage. No L48/L64/L96 production is frozen or launched; a question about cost is not authorization. New manager has closed launch gate and stopped idle coordinator; its remaining work, if any, is finite user-authorized retention ofauditednewarchives, not furtherMC. SSD originals remainheld.

### October 15 17:21: user-authorized N48/N64 spin-dimension pilot running

Direct user authorization: go ahead with48^3/64^3pilot, assess systematics before committing to much longer runs. Frozen d3-zero-beta-spin-N48-N64-pilot-20261015 has28trajectories=14exactcoordinates,20000cycles,two fresh zero/unit-x winding starts(r30/31). N48q=.80075:.00025:.80225; N64q=.80050:.00025:.80200. Windows cover prior-root stress predictions, not a presumednu/Delta. Correctedfulltransient1.01435TB left1.633TB beyondallknownpending/live/copygrowthand25GiBreserve atdesign. Initial1–3hcalendar/17.2workerhforecast isconditional volume scaling; actualcost/mixing/I/Omustbemeasured.

Detached service org.luscherxy.d3-spin-larger-pilot-20261015 runs ../2026-10-15-xy-3d-spin-dimension/run_larger_pilot.py. Verifiedone managerPID3267andonecoordinatorPID3271at17:20,14N48jobsactive/14N64pending. Neverduplicate. The previous82-runfinite manager finisheditsretentionwithnofailures beforethestudylockwasacquired. Resourcecap14single-threadsimulations leaves4coreheadroom. Coordinatorwrappernowhonorscontrol.launch_allowed_ids; onlythese28pendingIDsareallowed. Frozen numericalrun.py/modelkernelsunchanged.18focusedanalysis/topology/routing/admissiontests pass.

Manager waits all28fullaudits,closeslaunchgate/stopsidlecoordinator, pins combined110sources(fine82+pilot28), executes11cut/block/threshold/fullwindow/startparity/large-onlysnapshots,streamscompletewindinghistoriesforactualaccepted-eventdenominators/topologytau/constant-axis/start/cutdiagnostics, thenproducesaPRELIMINARYreport/PNG/PDF/SVGviareport_larger_pilot.py. Itsatomicstateis larger-pilot-status.json; fullstdout/stderr,resourceRSS/CPUsamplesandworkerlogsarepreserved. Topologytauexceedingmomenttaurequiresuncertaintyreview;flagsdonotautomaticallyenlargebootstrapblocks. No good-lookingfit/audit/acceptedmoveprovesmixingorcontrolledsystematics. Additionalmin24/32leadingandcorrectionpower.5/.8/1/1.5/2stressfitskeepomitted/zeroDOF/rejectedqualifications.

Finalindependentscientificreviewmustexamineheadline/noise/absolutescale/scope/coveragebeforeacceptance. .00025epsilonstepsmaystillbetoocoarseinlargerboxes; noautomaticrefinement,longerrunsorN96areauthorized. Afteralldiagnosticextraction,finite2-workeradaptivefieldretentionkeepsallhistories/RNG/checkpoints/logs/provenanceandperformsatomicverifiedconversion. SSDoriginalsremainheld,nocopyrelocation/deletion. ThepilotandanalysisareindependentofCodexapprestarts.

### October 15: common-coupling xi/L quotient analysis, existing data only

User questioned the spin-dimension systematic estimate and requested common-coupling xi/L. The earlier arbitrary size-cut/correction-fit spread is not a calibrated sigma_sys and must not be quoted as one. New ../2026-10-15-xy-3d-spin-dimension/common_xi_quotients.py uses common-epsilon xi/L crossings and m2 quotients, joint within-start bootstrap of m2/m4/mabs/S0/F, shared coordinate draws across overlapping pairs, exact piecewise polynomial crossing roots, and explicit invalid/multiple-root masks. Only <=8-ulp coordinate differences are pooled; original source epsilons remain recorded. Numerical kernels and production launch gate are unchanged; no new MC.

Source-pinned baseline/fourfold-block snapshots contain150 completed audited trajectories including older coarse epsilon anchors. latest-common-xi-delivery.json points to common-xi-delivery-20261015-195403 (study directory), with report, PNG/PDF/SVG, table, native per-size interpolation brackets and provenance. Linear Delta_eff for12:24,16:32,32:64 is .52323/.52247/.52209 with conditional SD .00633/.00845/.00706. These are finite-pair estimates, not a final asymptotic Delta or total uncertainty. N32 still interpolates .8000:.8015; finer N64 knots do not refine that native interval. 24:48 remains unresolved: no accepted linear crossing, PCHIP draw retention .665 baseline/.682 block4. Do not equate the initial direct-R central crossing diagnostic with moment-consistent crossing coverage. Diagnostic fixed-omega.789 correction fits retain omitted-pair/zeroDOF/rejection/common-mask qualifications. Sixteen focused oracle/analysis tests pass; independent review is recorded separately in delivery review.json. No further production follows automatically.

### October 15 20:03: authorized common-coupling grid/bracket refinement running

Direct user authorization: improve coarse epsilon grid and resolve24:48 common-xi crossing. Frozen d3-zero-beta-common-xi-crossing-refinement-20261015 has120trajectories=60newcoordinates, plus reuse of two clean coordinate groups atN24/32q.8015. N6/8/12/16/24/32 receive common .000125-spaced inner knots aroundq.8005-.80125 with .80025/.8015 guards,12000cycles,fresh40/41zero/unit-xstarts. N48 addsq.80025/.8005/.800625/.800875/.801125; N64 adds.800625/.800875/.801125,20000cycles. No largerboxes. LongN64/N48trajectoriesfirst. Reviewed fulltransient636.475GBleft1.843TBbeyond25GiBreserve; no assumedfuturethinning. CorrectedmatchedD3beta0medianforecast11.178workerh,0.798ideal14workerh,1–3hcalendarplanningallowance incltail/audits/analysis, notcoverageorprecisionguarantee.

Finite launchd service org.luscherxy.d3-common-xi-refinement-20261015, plistandrunner under study directory. Verified manager20946/coordinator20950,14running106pending at20:03. Initial preproduction startup hit a directory-versus-analysis-file provenance path error; no jobs launched in that attempt. Fixed path normalization, original traceback preserved, regression test added;30focused tests nowpass. Exactly one active manager/coordinator. Frozennumericalrun.py/modelkernelsunchanged. Previouspilotmanager exited/lockfree/retentioncomplete. Newmanager uses singleton+studylock,120-IDallowlist,14single-threadcap,Friday15deadline andfullclosureaudits; keepsallhistories/checkpoints/logs andSSDoriginalhold.

Afterclosure it pins EXACTold150+new120sourcecohort, computes common-xi primary thenblock4/latecut4001/even/odd variants, exports plot/table/report, and separately checks24:48linear/PCHIP conditionalcrossingcoverage inprimary/block4/latecut. Incompletecoverage publishes crossing_still_needs_measured_refinement, notcompletedscience. Independent scientific review of attainedfits/starts/cuts/interpolation remainsrequired; no artificial sigma_sys orautomaticprecisionclaim. Only afteranalysis do2retentionworkers apply already-authorized field policy. Status common-xi-refinement-status.json, unfiltereddaemon/coordinator/workerlogs preserved; workers surviveCodexrestart. No unrelatedproductionorarchivecleanup.

### October15 21:38: refined common-coupling scan audited; local-response follow-up

All120newtrajectories passed full audits by20:52:52, about50minutes/10.4027worker-hours. Manager thenfailed beforeanalysis because progress-state variable previous shadowed its prior-analysis Path. Repaired to distinct prior_analysis/previous_states; existingidle service restarted21:33, noMC rerun. Originaltracebacks/logs retained. Integrationtestnow exercises running->complete->analysis and exactold150+new120cohort.31combinedprioranalysis/admission/manager tests, and25focusednewlocalresponse/analysis tests passed (see currentlogs for actualchecks).

Five270/135source-pinned pointwise snapshotsclosed. Fine nativeN32/N64knots .80075:.800875 removeoldN32coarsegap;32:64rawlinearDelta=.525933±.005717stat/PCHIP.526122±.006109. Startparity.536232±.006881 versus.513416±.008229 isdescriptive2.1combinedSD, notmixingproof/failure. Finepointwiseinterpolation exposes noise/multiplecrossings:24:48linearvalid.9705,PCHIP.9285; managercorrectlyretains crossing_still_needs_measured_refinement for THATestimator. Allunresolvedtablelabels/lateparitydetailsareexplicit inupdatedreviewedpointwise213724delivery.

Analysis-only local_xi_response.py tests PREDECLARED physicalwindows inner.800625:.801000 andbroad.800500:.801500,degree1/2. Marginal WLSensembleR/m2fits use matchingknots,scaledcoordinates,fixedcentralweights, originaljointbootstrapdraws andsharedcoordinatecovariance; explicit rootfailurecauses retained. This changesestimatorfromrawmomentknotinterpolation andisnotanunannouncedrootchoice. Latestlocal213702,independentlyaccepted_with_qualifications: broadquadratic24:48Delta_eff=.515981±.002921stat,q=.800756472;32:64=.519027±.003260stat,q=.800742754. Bothretain2000/2000inall5pooled/block/late/parityvariants. This is conditional local-modelcoverage, notall-model resolution: innerquadratic24retains.97 andremainsunresolved;N6broadquadratic.6885unresolved; broadlinear missesN48/N64curvatureandisexplicitlyrejected-response. Goodpdoesnotcertifybias.

refined-common-coupling-summary.txtandfinite-size-diagnostics.jsoninthe local213702folderretainallconditions. Correlatedfixedomega.789sizefitsillustrate .523209±.003462stat(minpair8) versus.518340±.006158(min12),innerlinear .520570±.006399; no preferredasymptoticCFTresult, invented numericalsigma_sys orcertifiedtotalprecision. This follow-up required no new MC; no further simulationbatchwaslaunched. Preserveallhistories/checkpoints/logs andSSD-originalhold; finite2workeradaptivefieldretentioncontinuesseparately.

### October15: four-size correlation-length overlay

Userrequested allL24/32/48/64 onONEplot, superseding the pairwise presentation. Study latest-four-size-xi-plot.json points to four-size-xi-crossing-20261015-223504: PNG/PDF/SVG, plotted-data/provenance and independentreview. Exact same8measuredphysicalknots and independentquadraticresponsefits; no commoncrossingimposed. Allcoefficients/sourcehashesmatchreviewedlocalanalysis; rawpoints/fitbandsconditionalMC95, nottotalerrors. PDF rendered/checked. This is presentationonly, no new MC or fit-method change.

### October15: scaling-dimension read-off plots

UserrequestedanalogousplotreadDelta. Study latest-dimension-readout-plots.json points to dimension-readout-plots-20261015-224446: log-logRMSmagnetizationvsL24/32/48/64atEXACTmeasuredq=.80075, slope=-Delta_eff=.519421±.007947stat,MC95[.503978,.535052],chi2/dof=.345895/2; this isfixedrepresentativeqnearxi-crossings, notepsilon_c ornewthermodynamicDelta. NoDeltaimposed; shared-bootstrapGLS matchesindependentpower1oracle. Companionplotshowscommon-qDeltaeffcurvesfor12:24,16:32,24:48,32:64frommeasuredratios/reviewedmagnetizationresponses. Greystripisdescriptiveunionoffinite-xi-paircrossingintervals, notinfinite-volumeepsilon_cCI. Allplottedbars/bandsconditionalMC95,finite-sizecorrectionsremain. Two-pagePDF/PNG/SVG,data/sourcehashes,independentreviewacceptedqualified; no newMC.
