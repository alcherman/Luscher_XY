"""Measured pilot costs and explicit, conditional campaign resource scenarios."""
from collections import defaultdict
from datetime import datetime
import importlib.util
import json
from pathlib import Path
import statistics

OUT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("frozen_xy_runner", OUT / "run.py")
run = importlib.util.module_from_spec(spec)
spec.loader.exec_module(run)


def main():
    jobs = [json.loads(p.read_text()) for p in (OUT / "jobs").glob("*.json")]
    manifest = json.loads((OUT / "batches/calibration.json").read_text())
    ids = {p["id"] for p in manifest["points"]}
    jobs = [j for j in jobs if j["point"]["id"] in ids]
    assert len(jobs) == 36 and all(j["state"] == "complete" and j["audit"]["passed"] for j in jobs)
    grouped = defaultdict(list)
    for j in jobs:
        p = j["point"]
        path = OUT / "chains" / f"{p['id']}.h5"
        assert path.stat().st_size == j["audit"]["bytes"]
        grouped[(p["D"], p["N"])].append(j)
    costs = {}
    for (D, N), js in grouped.items():
        times = [j["attempt_generation_seconds"] + j["audit_seconds"] for j in js]
        costs[(D, N)] = dict(D=D, N=N, trajectories=len(js), cycles=4000,
                            median_seconds=statistics.median(times),
                            min_seconds=min(times), max_seconds=max(times),
                            mean_bytes=statistics.mean(j["audit"]["bytes"] for j in js))
    first = min(datetime.fromisoformat((OUT / "logs" / f"{j['point']['id']}.log").read_text().splitlines()[0].split()[0]) for j in jobs)
    last = max(datetime.fromisoformat(j["updated"]) for j in jobs)
    wall = (last - first).total_seconds()
    worker_seconds = sum(j["attempt_generation_seconds"] + j["audit_seconds"] for j in jobs)
    effective_workers = worker_seconds / wall
    def projection(D, N, cycles):
        anchor = costs[(D, N)] if (D, N) in costs else costs[(D, 12)]
        scale = (N / anchor["N"])**D * cycles / 4000
        return anchor["median_seconds"] * scale, anchor["mean_bytes"] * scale
    # Candidate slices are not yet validated infinite-volume branch membership.
    slices = {3: 11, 4: 7}
    baseline_time = baseline_bytes = gate_bytes = 0.
    baseline_count = 0
    for D, count in slices.items():
        for N in [6, 8, 12]:
            n = count * 3 * 2
            t, b = projection(D, N, 4000)
            baseline_time += n * t
            baseline_bytes += n * b
            gate_bytes += n * run.estimated_bytes(run.point(D, 100, 300, N, 0, 4000))
            baseline_count += n
    # Three representative slices per dimension get a fourth size. A fifth
    # modest size in 4D enables an identifiable two-correction size-cut test.
    sensitivity_time = sensitivity_bytes = sensitivity_gate = 0.
    sensitivity_count = 0
    for D, Ns in [(3, [16]), (4, [10, 16])]:
        for N in Ns:
            n = 3 * 3 * 2
            t, b = projection(D, N, 4000)
            sensitivity_time += n * t
            sensitivity_bytes += n * b
            sensitivity_gate += n * run.estimated_bytes(run.point(D, 100, 300, N, 0, 4000))
            sensitivity_count += n
    pilot_bytes = sum(j["audit"]["bytes"] for j in jobs)
    new_seconds = baseline_time + sensitivity_time - worker_seconds
    new_bytes = baseline_bytes + sensitivity_bytes - pilot_bytes
    new_gate = gate_bytes + sensitivity_gate - manifest["estimated_bytes"]
    # Optional allowance, not a frozen batch or a claimed measured current cost.
    d2_allowance = [run.point(2, 60, 250, 128, r, 10000) for _ in range(4) for _beta in range(3) for r in range(2)]
    d2_gate = sum(run.estimated_bytes(p) for p in d2_allowance)
    current = datetime.now(run.CENTRAL)
    status = json.loads((OUT / "status.json").read_text())
    usable = status["disk_free_bytes"] - run.RESERVE
    result = dict(created=current.isoformat(), deadline=run.DEADLINE.isoformat(),
                  remaining_hours=(run.DEADLINE - current).total_seconds() / 3600,
                  calibration=dict(trajectories=36, audits_passed=36, start=first.isoformat(),
                                   end=last.isoformat(), wall_seconds=wall,
                                   worker_seconds=worker_seconds, effective_workers=effective_workers,
                                   bytes=pilot_bytes, costs=list(costs.values())),
                  baseline=dict(candidate_slices=slices, beta_points=3, sizes=[6,8,12],
                                independent_starts=2, cycles_per_start=4000,
                                trajectories_including_calibration=baseline_count),
                  selected_size_checks=dict(slices_per_dimension=3, d3_sizes=[16], d4_sizes=[10,16],
                                            beta_points=3, starts=2, cycles=4000, new_trajectories=sensitivity_count),
                  projected_new_d3_d4=dict(trajectories=baseline_count + sensitivity_count - 36,
                                          worker_hours=new_seconds / 3600,
                                          wall_hours_at_observed_throughput=new_seconds / effective_workers / 3600,
                                          bytes=new_bytes, conservative_gate_bytes=new_gate),
                  optional_d2_allowance=dict(slices=4, beta_points=3, sizes=[128], starts=2,
                                            cycles=10000, trajectories=24, conservative_gate_bytes=d2_gate,
                                            current_runtime_unmeasured=True),
                  usable_disk_bytes=usable,
                  conditional_delivery_range_hours=[12,36],
                  conditional_delivery_range_note="Planning allowance for adaptive refinement, fitting, diagnostic review and selected d2 runs; not a statistical prediction or a precision guarantee.")
    run.atomic_json(OUT / "runtime-forecast-2026-10-13.json", result)
    rows = "\n".join(f"| {c['D']} | {c['N']} | 6 | {c['median_seconds']:.1f} | {c['min_seconds']:.1f}–{c['max_seconds']:.1f} | {c['mean_bytes']/2**30:.3f} |" for _, c in sorted(costs.items()))
    text = f"""# XY campaign runtime and storage forecast

As of {current.strftime('%B %d at %I:%M %p Central')}, all 36 calibration
trajectories completed and passed their full-field audits. The Friday October
16 at 3 p.m. Central deadline is {result['remaining_hours']:.1f} hours away.

All timing measurements below are at epsilon/pi=1, the ordinary XY endpoint;
constrained-slice runtimes remain unmeasured. The current resource evidence
supports a first pass and targeted size checks before Friday. A conditional
12 to 36 hour working allowance covers sampling,
adaptive refinement, fitting and diagnostic review; it is a planning range,
not a guaranteed duration or evidence that every slice will have a validated
infinite-volume uncertainty of 0.02. Extrapolation validation is the principal
scientific uncertainty. Storage figures below apply to the filesystem currently
checked by the coordinator. Wider stable error bars and
unresolved slices remain allowed deliverables.

## Measured calibration costs

Each trajectory has 4000 cycles. Times below include generation and the
independent field audit. All six trajectories per row are actual events,
including different couplings and starts; these are descriptive ranges, not
independent performance-replication confidence intervals.

| Dimension | N | Trajectories | Median seconds | Range seconds | GiB per trajectory |
| --- | --- | --- | --- | --- | --- |
{rows}

The batch ran from {first.strftime('%I:%M:%S %p')} to
{last.strftime('%I:%M:%S %p Central')}, {wall/60:.1f} minutes elapsed. Its summed
generation and audit times are {worker_seconds/3600:.2f} worker-hours, for
{effective_workers:.2f} effective concurrent workers over the entire batch.
The configured maximum was 16; the lower achieved number includes launch
overhead, shared resource contention and the tail of six larger trajectories.
The closed data total is {pilot_bytes/2**30:.2f} GiB.

## Explicit sampling scenario

Reuse the existing d=2 data on 22 baseline candidate slices from epsilon/pi=0.58
to 1.00 at N=16,32,64. Reanalysis and validated BKT matching are needed;
having these measurements does not imply that the infinite-volume estimates
are already determined.

For d=3, use 11 candidate slices 0.80,0.82,...,1.00. For d=4, use 7 candidate
slices 0.88,0.90,...,1.00. Branch membership at the first slices and neighboring
zero-coupling endpoint cells still needs validation. Three beta points per
slice, N=6,8,12, two independent starts and 4000 cycles per start give
{baseline_count} trajectories, of which the 36 calibration trajectories are
already available. Couplings must be chosen adaptively to cover each transition;
three samples are an initial neighborhood, not a commitment to stop there.

Add N=16 on three representative slices per dimension. In d=4, also add N=10
on those three slices to check identifiability and the minimum-size sensitivity
of the two-correction fit. These checks add {sensitivity_count} trajectories.
Together the scenario has {baseline_count + sensitivity_count - 36} new d=3,4
trajectories. Scaling measured N=12 costs by site count for unmeasured sizes,
and linearly with cycles, gives {new_seconds/3600:.1f} additional worker-hours,
or {new_seconds/effective_workers/3600:.1f} wall hours at the observed pilot
daytime throughput. This takes no credit for the permitted increase to 18
workers outside office hours; overnight throughput has not yet been measured.
Both scaling assumptions need replacement with actual measurements;
constrained slices and larger boxes can cost more than the ordinary endpoint.
Allow roughly 2 to 6 elapsed hours for this first sampling pass rather than
treating the central throughput calculation as a guarantee.

As a provisional d=2 allowance, four selected slices at N=128, three beta
points, two starts and 10000 cycles require 24 new trajectories. Their current
runtime has not been measured. These are targeted follow-ups only after
existing-data diagnostics identify where larger boxes are necessary.

## Storage ceiling and refinement scenarios

Current usable disk space above the 25 GiB reserve is {usable/2**30:.1f} GiB.
The projected new d=3,4 records occupy about {new_bytes/2**30:.1f} GiB; the
runner's conservative launch estimates total {new_gate/2**30:.1f} GiB. The
optional d=2 allowance has a conservative estimate of {d2_gate/2**30:.1f} GiB.
Their combined conservative estimate is {(new_gate+d2_gate)/2**30:.1f} GiB,
leaving about {(usable-new_gate-d2_gate)/2**30:.1f} GiB for adaptive additions,
other outputs and estimation error above the reserve.

Doubling or quadrupling statistics everywhere would approximately double or
quadruple the dense-record storage as well as work. A 2x or 4x blanket extension
of this scenario would not fit the original internal-disk budget. The external
SSD route, when active, removes that particular ceiling; remeasure writing and
audit throughput there before treating internal-disk timing as an SSD timing
measurement. Repeated new independent
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
"""
    (OUT / "runtime-forecast-2026-10-13.md").write_text(text)
    print(json.dumps({k: result[k] for k in ["remaining_hours", "calibration", "projected_new_d3_d4", "usable_disk_bytes"]}, indent=2))


if __name__ == "__main__":
    main()
