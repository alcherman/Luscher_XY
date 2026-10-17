"""History plots, blocked Binder analysis, and final 4D scout plot."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np

from run import (OUT, EPSILON_RATIOS, LARGE_VOLUME_RATIOS, SIZES,
                 atomic_json, now, sha256)

sys.path.insert(0, str(OUT.parent / "2026-09-02-xy-binder-bkt"))
from experiment import tau_int


PLOT_SIZES = SIZES + [20]
COLORS = {8: "#0072B2", 12: "#D55E00", 16: "#009E73", 20: "#CC79A7"}
MARKERS = {8: "o", 12: "s", 16: "^", 20: "D"}
RATIO_LABELS = {
    0.5: r"$1/2$",
    2 / 3: r"$2/3$",
    5 / 6: r"$5/6$",
    0.84: r"$0.84$",
    0.85: r"$0.85$",
    0.86: r"$0.86$",
    0.87: r"$0.87$",
    7 / 8: r"$7/8$",
    11 / 12: r"$11/12$",
    23 / 24: r"$23/24$",
    1.0: r"$1$",
}


def completed():
    result = []
    for path in (OUT / "jobs").glob("*.json"):
        job = json.loads(path.read_text())
        if job.get("state") == "complete":
            result.append(job["point"])
    return sorted(result, key=lambda p: (p["epsilon_over_pi"], p["L"], p["replica"]))


def read(spec):
    with h5py.File(OUT / "chains" / f"{spec['id']}.h5", "r") as f:
        assert f.attrs["committed_rows"] == spec["steps"] + 1
        return {key: value[:] for key, value in f.items() if key not in {"phi", "rng_state"}}


def finite(value):
    if isinstance(value, dict):
        return {key: finite(item) for key, item in value.items()}
    if isinstance(value, list):
        return [finite(item) for item in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def binder(m2, m4):
    return float(1 - np.mean(m4) / (2 * np.mean(m2) ** 2))


def chain_diagnostics(data, cut):
    m2 = data["magnetization_squared"][cut:]
    m4 = data["magnetization_fourth"][cut:]
    influence = (m4.mean() / m2.mean() ** 3 * (m2 - m2.mean())
                 - (m4 - m4.mean()) / (2 * m2.mean() ** 2))
    winding2 = (data["winding"][cut:].astype(float) ** 2).mean(axis=(1, 2))
    monitored = dict(
        m2=m2, m4=m4, binder=influence,
        link_cosine=data["mean_link_cosine"][cut:],
        vortex=data["vortex_density"][cut:], winding_squared=winding2,
    )
    return {key: tau_int(values) for key, values in monitored.items()}


def block_sums(values, width):
    n_blocks = len(values) // width
    if n_blocks < 2:
        return None, None
    edges = np.linspace(0, len(values), n_blocks + 1, dtype=int)
    return np.add.reduceat(values, edges[:-1], axis=0), np.diff(edges)


def pooled_summary(group, cut, draws=4000, width_factor=1, seed_offset=0):
    chains = [(spec, read(spec)) for spec in group]
    diagnostics = [chain_diagnostics(data, cut) for _, data in chains]
    widths = []
    arrays = []
    for (spec, data), taus in zip(chains, diagnostics):
        finite_taus = [value for value in taus.values() if value is not None and np.isfinite(value)]
        width = width_factor * max(20, int(np.ceil(10 * max(finite_taus))))
        widths.append(width)
        arrays.append(np.column_stack((data["magnetization_squared"][cut:],
                                       data["magnetization_fourth"][cut:])))
    pooled = np.concatenate(arrays)
    central = binder(pooled[:, 0], pooled[:, 1])
    rng = np.random.default_rng(np.random.SeedSequence(
        [2026090404, int(round(1000 * group[0]["epsilon_over_pi"])), group[0]["L"], cut,
         width_factor, seed_offset]
    ))
    blocks = [block_sums(values, width) for values, width in zip(arrays, widths)]
    samples = np.full(draws, np.nan)
    if all(sums is not None for sums, _ in blocks):
        for draw in range(draws):
            total = np.zeros(2)
            count = 0
            for sums, lengths in blocks:
                indices = rng.integers(len(sums), size=len(sums))
                total += sums[indices].sum(axis=0)
                count += int(lengths[indices].sum())
            means = total / count
            samples[draw] = 1 - means[1] / (2 * means[0] ** 2)
    valid = samples[np.isfinite(samples)]
    error = float(valid.std(ddof=1)) if len(valid) > 1 else None
    per_start = [binder(values[:, 0], values[:, 1]) for values in arrays]
    return finite(dict(
        binder=central, error=error, cut=cut,
        retained_per_start=[len(values) for values in arrays],
        block_widths=widths,
        blocks_per_start=[len(sums) if sums is not None else 0 for sums, _ in blocks],
        enough_blocks=all(sums is not None and len(sums) >= 20 for sums, _ in blocks),
        taus=diagnostics, per_start_binder=per_start,
        start_difference=abs(per_start[0] - per_start[1]),
        local_acceptance=[float(data["local_accepted"][cut:].mean() / spec["L"] ** 3)
                          for spec, data in chains],
        mean_cluster_fraction=[float(data["cluster_sizes"][cut:].mean() / spec["L"] ** 3)
                               for spec, data in chains],
        winding_acceptance=[float(data["winding_accepted"][cut:].mean())
                            for _, data in chains],
    )), samples


def make_histories(specs):
    folder = OUT / "histories"
    folder.mkdir(exist_ok=True)
    datasets = {spec["id"]: read(spec) for spec in specs}
    for stop, label in ((2001, "thermalization"), (10001, "full")):
        fig, axes = plt.subplots(len(EPSILON_RATIOS), len(SIZES), figsize=(13, 11), squeeze=False)
        for row, epsilon_ratio in enumerate(EPSILON_RATIOS):
            for column, L in enumerate(SIZES):
                ax = axes[row, column]
                for replica, linestyle in ((0, "-"), (1, "--")):
                    candidates = [spec for spec in specs
                                  if spec["epsilon_over_pi"] == epsilon_ratio
                                  and spec["L"] == L and spec["replica"] == replica]
                    if not candidates:
                        continue
                    values = datasets[candidates[0]["id"]]["magnetization_squared"][:stop]
                    x = np.arange(len(values))
                    width = 50
                    n_blocks = (len(values) - 1) // width
                    ax.plot(x[1:1 + n_blocks * width].reshape(-1, width).mean(axis=1),
                            values[1:1 + n_blocks * width].reshape(-1, width).mean(axis=1),
                            linestyle, color=COLORS[L], lw=0.9,
                            label=f"start {replica}")
                    ax.plot(0, values[0], ".", color=COLORS[L], ms=4)
                ax.set_title(rf"$\epsilon/\pi={epsilon_ratio:.3f}$, $L={L}$", fontsize=9)
                ax.grid(alpha=0.15)
                if row == len(EPSILON_RATIOS) - 1:
                    ax.set_xlabel("update cycle")
                if column == 0:
                    ax.set_ylabel(r"$|m|^2$")
        axes[0, 0].legend(fontsize=8)
        fig.suptitle(r"4D, $\kappa=0$: 50-cycle magnetization means; seed shown at cycle 0")
        fig.tight_layout(rect=(0, 0, 1, 0.975))
        fig.savefig(folder / f"magnetization_{label}.png", dpi=140)
        plt.close(fig)
    large = [spec for spec in specs if spec["L"] == 20]
    if large:
        fig, axes = plt.subplots(len(LARGE_VOLUME_RATIOS), 1, figsize=(8, 5.5), squeeze=False)
        for row, ratio in enumerate(LARGE_VOLUME_RATIOS):
            ax = axes[row, 0]
            for replica, linestyle in ((0, "-"), (1, "--")):
                spec = next(item for item in large
                            if item["epsilon_over_pi"] == ratio and item["replica"] == replica)
                values = datasets[spec["id"]]["magnetization_squared"]
                width = 50
                n_blocks = (len(values) - 1) // width
                x = np.arange(len(values))
                ax.plot(x[1:1 + n_blocks * width].reshape(-1, width).mean(axis=1),
                        values[1:1 + n_blocks * width].reshape(-1, width).mean(axis=1),
                        linestyle, color=COLORS[20], lw=0.9, label=f"start {replica}")
                ax.plot(0, values[0], ".", color=COLORS[20], ms=4)
            ax.set_title(rf"$\epsilon/\pi={ratio:.3f}$, $L=20$")
            ax.set_ylabel(r"$|m|^2$")
            ax.grid(alpha=0.15)
        axes[-1, 0].set_xlabel("update cycle")
        axes[0, 0].legend(fontsize=8)
        fig.suptitle(r"4D, $\kappa=0$: $L=20$ 50-cycle magnetization means")
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        fig.savefig(folder / "magnetization_L20.png", dpi=140)
        plt.close(fig)


def analyze(specs, cut):
    rows = []
    stability = []
    for epsilon_ratio in EPSILON_RATIOS:
        sizes = SIZES + ([20] if epsilon_ratio in LARGE_VOLUME_RATIOS else [])
        for L in sizes:
            group = [spec for spec in specs if spec["epsilon_over_pi"] == epsilon_ratio and spec["L"] == L]
            if len(group) != 2:
                raise RuntimeError(f"Expected two starts at epsilon/pi={epsilon_ratio}, L={L}")
            summary, samples = pooled_summary(group, cut)
            epsilon_id = int(round(1000 * epsilon_ratio))
            np.savez_compressed(OUT / f"bootstrap_q{epsilon_id:04d}_L{L:03d}.npz", samples=samples)
            checks = []
            for other_cut in sorted({1001, cut, 4001}):
                factors = (1, 2) if other_cut == cut else (1,)
                for factor in factors:
                    check, _ = pooled_summary(group, other_cut, draws=1000,
                                              width_factor=factor, seed_offset=1)
                    checks.append(check)
            flags = []
            if not summary["enough_blocks"]:
                flags.append("fewer_than_20_blocks_per_start")
            if summary["error"] is not None and summary["start_difference"] > 3 * np.sqrt(2) * summary["error"]:
                flags.append("independent_start_difference")
            for check in checks:
                if check["cut"] == cut:
                    continue
                if None not in (summary["error"], check["error"]):
                    if abs(summary["binder"] - check["binder"]) > 2 * np.hypot(summary["error"], check["error"]):
                        flags.append("cut_sensitive_binder")
            row = dict(epsilon=float(np.pi * epsilon_ratio), epsilon_over_pi=epsilon_ratio,
                       L=L, binder=summary["binder"], error=summary["error"],
                       flags=sorted(set(flags)), diagnostics=summary)
            rows.append(row)
            stability.append(dict(epsilon_over_pi=epsilon_ratio, L=L, checks=checks))
            print(f"epsilon/pi={epsilon_ratio:.3f}, L={L}: U={summary['binder']:.6f} "
                  f"+/- {summary['error']:.6f}; widths={summary['block_widths']}; "
                  f"blocks={summary['blocks_per_start']}; flags={flags}", flush=True)

    atomic_json(OUT / "results.json", rows)
    atomic_json(OUT / "stability.json", stability)
    with (OUT / "binder.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=("epsilon_over_pi", "epsilon", "L", "binder", "error", "flags"))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: ";".join(value) if key == "flags" else value
                             for key, value in row.items() if key in writer.fieldnames})
    atomic_json(OUT / "analysis-provenance.json", dict(
        updated_utc=now(), cut=cut, bootstrap_draws=4000,
        analysis_sha256=sha256(Path(__file__)),
        estimator="ratio of pooled raw moments from two equal-length independent starts",
        bootstrap="stratified contiguous-block bootstrap within each start",
        block_rule="per start: max(20, ceil(10*largest monitored tau_int)); all retained rows used",
        monitored="m2, m4, Binder influence, mean link cosine, vortex density, winding squared",
    ))
    make_plot(rows)
    return rows


def make_plot(rows):
    plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.2})
    fig, (ax, zoom) = plt.subplots(
        1, 2, figsize=(11.5, 5.3), sharey=True,
        gridspec_kw={"width_ratios": [1.65, 1]},
    )
    for axis in (ax, zoom):
        for L in PLOT_SIZES:
            group = [row for row in rows if row["L"] == L]
            axis.errorbar([row["epsilon_over_pi"] for row in group],
                          [row["binder"] for row in group],
                          yerr=[row["error"] for row in group],
                          color=COLORS[L], marker=MARKERS[L], ms=6, lw=1.5,
                          capsize=3, label=f"L={L}")
        axis.set_xlabel(r"admissibility cutoff $\epsilon/\pi$")
    ax.set_ylabel(r"$U_L=1-\langle |m|^4\rangle/(2\langle |m|^2\rangle^2)$")
    ax.legend(frameon=False)
    broad_ticks = [0.5, 2 / 3, 5 / 6, 7 / 8, 11 / 12, 23 / 24, 1.0]
    ax.set_xticks(broad_ticks, [RATIO_LABELS[value] for value in broad_ticks])
    ax.tick_params(axis="x", labelrotation=30)
    ax.set_xlim(0.48, 1.02)
    zoom_ticks = [5 / 6, 0.84, 0.85, 0.86, 0.87, 7 / 8]
    zoom.set_xticks(zoom_ticks, [RATIO_LABELS[value] for value in zoom_ticks])
    zoom.tick_params(axis="x", labelrotation=45)
    zoom.set_xlim(0.828, 0.88)
    zoom.set_title("transition-region detail", fontsize=11)
    fig.suptitle(r"4D admissible XY model at $\kappa=0$")
    fig.tight_layout()
    for extension in ("png", "pdf"):
        fig.savefig(OUT / f"binder_vs_epsilon.{extension}", dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cut", type=int)
    parser.add_argument("--histories", action="store_true")
    args = parser.parse_args()
    specs = completed()
    expected = 2 * (len(EPSILON_RATIOS) * len(SIZES) + len(LARGE_VOLUME_RATIOS))
    if len(specs) != expected:
        raise RuntimeError(f"Expected {expected} completed trajectories, found {len(specs)}")
    if args.histories:
        make_histories(specs)
    if args.cut is not None:
        analyze(specs, args.cut)


if __name__ == "__main__":
    main()
