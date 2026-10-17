# XY campaign runtime and storage forecast

As of October 13 at 02:51 PM Central, all 36 calibration
trajectories completed and passed their full-field audits. The Friday October
16 at 3 p.m. Central deadline is 72.1 hours away.

The current resource evidence supports a first pass and targeted size checks
before Friday. A conditional 12 to 36 hour working allowance covers sampling,
adaptive refinement, fitting and diagnostic review; it is a planning range,
not a guaranteed duration or evidence that every slice will have a validated
infinite-volume uncertainty of 0.02. Storage and extrapolation validation are
more restrictive than the initial CPU forecast. Wider stable error bars and
unresolved slices remain allowed deliverables.

## Measured calibration costs

Each trajectory has 4000 cycles. Times below include generation and the
independent field audit. All six trajectories per row are actual events,
including different couplings and starts; these are descriptive ranges, not
independent performance-replication confidence intervals.

| Dimension | N | Trajectories | Median seconds | Range seconds | GiB per trajectory |
| --- | --- | --- | --- | --- | --- |
| 3 | 6 | 6 | 53.3 | 52.9–53.8 | 0.013 |
| 3 | 8 | 6 | 60.5 | 60.4–60.8 | 0.025 |
| 3 | 12 | 6 | 113.1 | 109.4–130.7 | 0.069 |
| 4 | 6 | 6 | 133.3 | 107.8–134.0 | 0.069 |
| 4 | 8 | 6 | 203.4 | 109.5–203.7 | 0.190 |
| 4 | 12 | 6 | 358.1 | 304.5–402.0 | 0.845 |

The batch ran from 01:50:55 PM to
01:59:25 PM Central, 8.5 minutes elapsed. Its summed
generation and audit times are 1.47 worker-hours, for
10.41 effective concurrent workers over the entire batch.
The configured maximum was 16; the lower achieved number includes launch
overhead, shared resource contention and the tail of six larger trajectories.
The closed data total is 7.27 GiB.

## Explicit sampling scenario

Reuse the existing d=2 data on all 22 requested slices from epsilon/pi=0.58
to 1.00 at N=16,32,64. Reanalysis and validated BKT matching are needed;
having these measurements does not imply that the infinite-volume estimates
are already determined.

For d=3, use 11 candidate slices 0.80,0.82,...,1.00. For d=4, use 7 candidate
slices 0.88,0.90,...,1.00. Branch membership at the first slices and neighboring
zero-coupling endpoint cells still needs validation. Three beta points per
slice, N=6,8,12, two independent starts and 4000 cycles per start give
324 trajectories, of which the 36 calibration trajectories are
already available. Couplings must be chosen adaptively to cover each transition;
three samples are an initial neighborhood, not a commitment to stop there.

Add N=16 on three representative slices per dimension. In d=4, also add N=10
on those three slices to check identifiability and the minimum-size sensitivity
of the two-correction fit. These checks add 54 trajectories.
Together the scenario has 342 new d=3,4
trajectories. Scaling measured N=12 costs by site count for unmeasured sizes,
and linearly with cycles, gives 18.7 additional worker-hours,
or 1.8 wall hours at the observed pilot
throughput. Both scaling assumptions need replacement with actual measurements;
constrained slices and larger boxes can cost more than the ordinary endpoint.
Allow roughly 2 to 6 elapsed hours for this first sampling pass rather than
treating the central throughput calculation as a guarantee.

As a provisional d=2 allowance, four selected slices at N=128, three beta
points, two starts and 10000 cycles require 24 new trajectories. Their current
runtime has not been measured. These are targeted follow-ups only after
existing-data diagnostics identify where larger boxes are necessary.

## Storage ceiling and refinement scenarios

Current usable disk space above the 25 GiB reserve is 174.8 GiB.
The projected new d=3,4 records occupy about 104.6 GiB; the
runner's conservative launch estimates total 124.5 GiB. The
optional d=2 allowance has a conservative estimate of 35.6 GiB.
Their combined conservative estimate is 160.1 GiB,
leaving about 14.7 GiB for adaptive additions,
other outputs and estimation error above the reserve.

Doubling or quadrupling statistics everywhere would approximately double or
quadruple the dense-record storage as well as work. A 2x or 4x blanket extension
of this scenario does not fit the current disk budget. Repeated new independent
trajectories consume additional storage, too. Refine only slices diagnosed as
limiting the result, use measured resource forecasts before every new batch,
and retain all raw histories; no thinning, deletion or compaction is proposed.

The 12 to 36 hour delivery allowance reserves substantial time beyond the
2 to 6 hour first sampling pass for estimator implementation, fit calibration,
thermalization and autocorrelation review, bootstrap fits, additional starts,
coupling refinement and figure verification. It is conditional on targeted
additions fitting the storage ceiling. If extrapolation bias remains material,
report wider supported intervals or unresolved slices by the deadline.
Dimension 5 stays last in priority and disabled until the d=2,3,4 deliverables
are ready and a separate resource/time forecast supports it.

## Evidence and limits

Inputs are the frozen `batches/calibration.json`, all 36 matching `jobs/*.json`,
their per-trajectory logs, and the current `status.json`. Current chain file
sizes were checked against all 36 audit records. The companion JSON preserves
the numbers and scenario parameters. The audits establish recorded-field and
observable correctness; they do not establish equilibration, independent-start
agreement, fit identifiability or infinite-volume precision. Those checks are
the next scientific stage. No elapsed runtime for the unmeasured full grid or
larger boxes is being reported as an observed result.
