"""Reproduce the requested epsilon=pi, 21-point Binder experiment.

Run from the repository root with .venv/bin/python. No data are thinned or
deleted. Generation logs and one atomic status record per completed point
supplement the scan writer's per-100-sweep RNG checkpoints.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import hashlib
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
OUT = Path(__file__).resolve().parent
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

import h5py
import matplotlib.pyplot as plt
import numpy as np

from supervillain.xy import run_chain

SIZES = (4, 8, 16)
KAPPAS = (.80, .95, 1.05, 1.12, 1.20, 1.35, 1.50)
COLORS = {4: "#0072B2", 8: "#D55E00", 16: "#009E73"}
MARKERS = {4: "o", 8: "s", 16: "^"}


def points():
    for i, (L, kappa) in enumerate((L, k) for L in SIZES for k in KAPPAS):
        yield dict(L=L, kappa=kappa, seed=2026090401 + i,
                   file=f"chains/L{L:02d}_k{kappa:.2f}.h5")


def atomic_json(path, value):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temp.replace(path)


def run_point(point):
    t = perf_counter()
    path = OUT / point["file"]
    run_chain(path, D=2, N=point["L"], kappa=point["kappa"], epsilon=np.pi,
              sweeps=10_000, warmup=100, seed=point["seed"], backend="numba",
              parallel=False, winding=True, checkpoint_every=100,
              resume=path.exists())
    with h5py.File(path, "r") as f:
        assert int(f.attrs["committed_rows"]) == 10_001
        assert np.array_equal(f["sweep"][:], np.arange(10_001))
    return point | dict(seconds=perf_counter() - t, completed_utc=datetime.now(timezone.utc).isoformat())


def run():
    (OUT / "chains").mkdir(exist_ok=True)
    manifest = dict(action="S=-kappa*sum_positive_links(cos(theta))", D=2,
                    epsilon=float(np.pi), sweeps_per_point=10_000,
                    initialization="cold, phi=0", local_update="serial Numba, proposal width pi",
                    global_update="one symmetric winding proposal per sweep",
                    provisional_warmup=100, analysis_cut="selected after inspecting histories",
                    retained_data="initial field and every sweep, all observables and RNG states",
                    workers=3, points=list(points()),
                    experiment_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    atomic_json(OUT / "manifest.json", manifest)
    record = dict(started_utc=datetime.now(timezone.utc).isoformat(), completed=[])
    atomic_json(OUT / "status.json", record)
    with ProcessPoolExecutor(max_workers=3) as pool:
        jobs = {pool.submit(run_point, p): p for p in points()}
        for future in as_completed(jobs):
            result = future.result()
            record["completed"].append(result)
            atomic_json(OUT / "status.json", record)
            print(f"{len(record['completed']):2d}/21 complete: L={result['L']:2d} "
                  f"kappa={result['kappa']:.2f} {result['seconds']:.1f}s", flush=True)
    record["finished_utc"] = datetime.now(timezone.utc).isoformat()
    atomic_json(OUT / "status.json", record)


def read(point):
    with h5py.File(OUT / point["file"], "r") as f:
        n = int(f.attrs["committed_rows"])
        return {k: f[k][:n] for k in f if k not in ("phi", "rng_state")}


def style():
    plt.rcParams.update({"font.size": 11, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.dpi": 170})


def histories():
    style()
    directory = OUT / "histories"
    directory.mkdir(exist_ok=True)
    # All stored physical observables and acceptance counts, including row zero.
    for p in points():
        data = read(p)
        keys = [k for k in data if k != "sweep"]
        fig, axes = plt.subplots((len(keys) + 2) // 3, 3, figsize=(15, 11), squeeze=False)
        for ax, key in zip(axes.flat, keys):
            y = data[key].reshape(len(data["sweep"]), -1)
            ax.plot(data["sweep"], y, lw=.4, alpha=.6)
            ax.set_title(key, fontsize=9)
            ax.set_xlabel("sweep (seed at 0)", fontsize=8)
        for ax in list(axes.flat)[len(keys):]:
            ax.axis("off")
        fig.suptitle(f"XY histories: L={p['L']}, κ={p['kappa']:.2f}, ε=π")
        fig.tight_layout()
        fig.savefig(directory / f"L{p['L']:02d}_k{p['kappa']:.2f}.png")
        plt.close(fig)
    for key, label in (("magnetization_squared", r"$|m|^2$"),
                       ("action_density", r"$S/L^2$"), ("vortex_density", "vortex density")):
        fig, axes = plt.subplots(7, 3, figsize=(12, 15), sharex=True, squeeze=False)
        for p in points():
            row, col = KAPPAS.index(p["kappa"]), SIZES.index(p["L"])
            ax = axes[row, col]
            data = read(p)
            t, y = data["sweep"][:501], data[key][:501]
            ax.plot(t, y, lw=.65, color=COLORS[p["L"]], alpha=.8)
            ax.plot(0, y[0], "ko", ms=3)
            ax.set_title(f"L={p['L']}, κ={p['kappa']:.2f}", fontsize=10)
            ax.grid(alpha=.15)
        for ax in axes[-1]:
            ax.set_xlabel("sweep (cold seed at 0)")
        fig.suptitle(f"Initial relaxation: {label}", fontsize=16)
        fig.tight_layout(rect=(0, 0, 1, .97))
        fig.savefig(OUT / f"thermalization_{key}.png")
        plt.close(fig)


def tau_int(x):
    """Nonperiodic FFT ACF, Geyer initial-positive paired window (tau >= 1/2)."""
    x = np.asarray(x, dtype=float)
    x = x - x.mean()
    n = len(x)
    size = 1 << (2 * n - 1).bit_length()
    z = np.fft.rfft(x, n=size)
    covariance = np.fft.irfft(z * z.conj(), n=size)[:n] / np.arange(n, 0, -1)
    if covariance[0] <= 0:
        return None
    rho = covariance / covariance[0]
    pair_sum = 0.0
    previous = np.inf
    for lag in range(0, n - 1, 2):
        pair = rho[lag] + rho[lag + 1]
        if pair <= 0:
            break
        pair = min(pair, previous)
        pair_sum += pair
        previous = pair
    return float(max(.5, pair_sum - .5))


def binder_and_error(m2, m4, width, *, seed, replicas=4000):
    blocks = len(m2) // width
    if blocks < 20:
        raise ValueError(f"Only {blocks} complete blocks for width {width}")
    used = blocks * width
    a = m2[:used].reshape(blocks, width).mean(axis=1)
    b = m4[:used].reshape(blocks, width).mean(axis=1)
    estimate = 1 - b.mean() / (2 * a.mean() ** 2)
    rng = np.random.default_rng(seed)
    samples = np.empty(replicas)
    for begin in range(0, replicas, 200):
        stop = min(replicas, begin + 200)
        index = rng.integers(blocks, size=(stop - begin, blocks))
        samples[begin:stop] = 1 - b[index].mean(axis=1) / (2 * a[index].mean(axis=1) ** 2)
    return float(estimate), float(samples.std(ddof=1)), blocks, used, samples


def analyze(cut):
    style()
    rows, stability = [], []
    for i, p in enumerate(points()):
        data = read(p)
        if len(data["sweep"]) != 10_001:
            raise ValueError("Every chain must have exactly 10,000 sweeps plus its seed")
        m2, m4 = data["magnetization_squared"][cut:], data["magnetization_fourth"][cut:]
        # Delta-method influence series for the Binder ratio, plus physical diagnostics.
        influence = m4.mean() / m2.mean() ** 3 * (m2 - m2.mean()) - .5 / m2.mean() ** 2 * (m4 - m4.mean())
        diagnostics = dict(m2=m2, m4=m4, binder_influence=influence)
        for name in ("action_density", "vortex_density", "magnetization_real", "magnetization_imag"):
            diagnostics[name] = data[name][cut:]
        diagnostics["winding_squared"] = (data["winding"][cut:].astype(float) ** 2).mean(axis=(1, 2))
        taus = {name: tau_int(values) for name, values in diagnostics.items()}
        largest = max(t for t in taus.values() if t is not None)
        width = max(20, int(np.ceil(10 * largest)))
        value, error, blocks, used, samples = binder_and_error(m2, m4, width, seed=2026090300 + i)
        np.save(OUT / f"bootstrap_L{p['L']:02d}_k{p['kappa']:.2f}.npy", samples)
        q16, q84 = np.quantile(samples, [.16, .84])
        row = p | dict(binder=value, standard_error=error, bootstrap_q16=float(q16),
                       bootstrap_q84=float(q84), first_retained_sweep=cut,
                       retained_sweeps=len(m2), blocked_sweeps=used,
                       remainder_sweeps=len(m2) - used, block_width=width,
                       blocks=blocks, tau_max=largest, tau=taus,
                       local_acceptance=float(data["local_accepted"][cut:].mean() / p["L"] ** 2),
                       winding_acceptance=float(data["winding_accepted"][cut:].mean()))
        rows.append(row)
        for extra_cut in (cut, 501, 1001, 2001):
            for factor in (1, 2, 4):
                u, e, nb, ns, _ = binder_and_error(
                    data["magnetization_squared"][extra_cut:],
                    data["magnetization_fourth"][extra_cut:], width * factor,
                    seed=100000 + i * 100 + extra_cut + factor, replicas=2000)
                stability.append(dict(L=p["L"], kappa=p["kappa"], cut=extra_cut,
                                      block_width=width * factor, binder=u, error=e,
                                      blocks=nb, used=ns))
        # Consecutive time quarters use the same block estimator to expose drift.
        for quarter, indices in enumerate(np.array_split(np.arange(len(m2)), 4)):
            u, e, nb, ns, _ = binder_and_error(m2[indices], m4[indices], width,
                                              seed=200000 + i * 10 + quarter, replicas=2000)
            stability.append(dict(L=p["L"], kappa=p["kappa"], quarter=quarter + 1,
                                  block_width=width, binder=u, error=e, blocks=nb, used=ns))
    atomic_json(OUT / "results.json", rows)
    atomic_json(OUT / "stability.json", stability)
    fields = [k for k in rows[0] if k != "tau"]
    with (OUT / "binder.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    fig, ax = plt.subplots(figsize=(8.4, 5.8))
    for L in SIZES:
        series = [r for r in rows if r["L"] == L]
        ax.errorbar([r["kappa"] for r in series], [r["binder"] for r in series],
                    yerr=[r["standard_error"] for r in series], color=COLORS[L],
                    marker=MARKERS[L], ms=6, lw=1.6, capsize=3, label=f"L = {L}")
    ax.axhline(0, color=".4", ls="--", lw=1)
    ax.set(xlabel=r"Coupling $\kappa$", ylabel=r"Binder cumulant $U_L=1-\langle |m|^4\rangle/(2\langle |m|^2\rangle^2)$",
           title=r"2D XY model: $S=-\kappa\sum_\ell\cos\theta_\ell$, $\epsilon=\pi$")
    ax.set_xticks(KAPPAS)
    ax.grid(alpha=.2)
    ax.legend(frameon=False)
    fig.text(.5, .025, f"10,000 sweeps per point • sweeps 1–{cut-1} excluded • error bars: 1σ block bootstrap",
             ha="center", fontsize=10, color=".35")
    fig.tight_layout(rect=(0, .055, 1, 1))
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(OUT / f"binder_vs_kappa.{suffix}")
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    for ax, L in zip(axes, SIZES):
        for cut_value, marker in zip((cut, 501, 1001, 2001), ("o", "s", "^", "D")):
            group = [r for r in stability if r["L"] == L and r.get("cut") == cut_value
                     and r["block_width"] == next(a["block_width"] for a in rows if a["L"] == L and a["kappa"] == r["kappa"])]
            ax.errorbar([r["kappa"] for r in group], [r["binder"] for r in group],
                        yerr=[r["error"] for r in group], marker=marker, ms=3, lw=.8,
                        alpha=.8, capsize=2, label=f"discard {cut_value-1}")
        ax.set(title=f"L = {L}", xlabel=r"$\kappa$")
        ax.axhline(0, c=".5", ls="--", lw=.7)
        ax.grid(alpha=.15)
    axes[0].set_ylabel("Binder cumulant")
    axes[-1].legend(fontsize=8, frameon=False)
    fig.suptitle("Sensitivity to the thermalization cut")
    fig.tight_layout()
    fig.savefig(OUT / "cut_stability.png")
    plt.close(fig)
    print(f"Analyzed {len(rows)} points; tau_max range {min(r['tau_max'] for r in rows):.2f}–"
          f"{max(r['tau_max'] for r in rows):.2f}; SE range {min(r['standard_error'] for r in rows):.4f}–"
          f"{max(r['standard_error'] for r in rows):.4f}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("run", "histories", "analyze"))
    parser.add_argument("--cut", type=int)
    args = parser.parse_args()
    if args.mode == "run":
        run()
    elif args.mode == "histories":
        histories()
    else:
        if args.cut is None or not 1 <= args.cut < 10_001:
            parser.error("analyze requires --cut specifying the first retained sweep")
        analyze(args.cut)
