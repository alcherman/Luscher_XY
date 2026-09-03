"""Disk-only Binder analysis. Inspect histories before specifying --cut.

The generation script is preserved as run. This analysis blocks the joint
rotation-invariant moments; slow rotation of the global magnetization is
recorded but does not set the block width for a rotation-invariant estimator.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import dataset as exp


def bootstrap(m2, m4, width, seed, draws=4000, include_tail=False):
    nblocks = len(m2) // width
    if nblocks < 2:
        return dict(binder=None, error=None, blocks=nblocks, used=0, enough_blocks=False)
    used = len(m2) if include_tail else nblocks * width
    pairs = np.column_stack((m2[:used], m4[:used]))
    # Balanced contiguous blocks retain the whole tail without creating a short block.
    edges = np.linspace(0, used, nblocks + 1, dtype=int)
    lengths = np.diff(edges)
    sums = np.add.reduceat(pairs, edges[:-1], axis=0)
    def ratio(x):
        return 1 - x[..., 1] / (2 * x[..., 0] ** 2)
    central = ratio(sums.sum(axis=0) / lengths.sum())
    rng = np.random.default_rng(seed)
    values = np.empty(draws)
    for begin in range(0, draws, 200):
        end = min(begin + 200, draws)
        indices = rng.integers(nblocks, size=(end - begin, nblocks))
        values[begin:end] = ratio(sums[indices].sum(axis=1) / lengths[indices].sum(axis=1)[:, None])
    return dict(binder=float(central), error=float(values.std(ddof=1)),
                q16=float(np.quantile(values, .16)), q84=float(np.quantile(values, .84)),
                blocks=nblocks, used=used, enough_blocks=nblocks >= 20,
                block_min_width=int(lengths.min()), block_max_width=int(lengths.max()),
                include_tail=include_tail,
                bootstrap_bias=float(values.mean() - central), samples=values)


def analyze(cut, sizes=None, include_tail=False):
    exp.style()
    selected = exp.SIZES if sizes is None else sizes
    rows = ([r for r in json.loads((exp.OUT / 'results.json').read_text()) if r['L'] not in selected]
            if sizes is not None and (exp.OUT / 'results.json').exists() else [])
    stability = ([r for r in json.loads((exp.OUT / 'stability.json').read_text()) if r['L'] not in selected]
                 if sizes is not None and (exp.OUT / 'stability.json').exists() else [])
    cuts = list(dict.fromkeys((cut, 501, 1001, 2001)))
    for i, point in enumerate(exp.points()):
        if point['L'] not in selected:
            continue
        data = exp.read(point)
        if not np.array_equal(data['sweep'], np.arange(10_001)):
            raise ValueError('Expected exactly 10,000 sweeps plus the cold seed')
        a, b = data['magnetization_squared'][cut:], data['magnetization_fourth'][cut:]
        # dU/d<a> = <b>/<a>^3; dU/d<b> = -1/(2<a>^2).
        influence = b.mean() / a.mean() ** 3 * (a - a.mean()) - (b - b.mean()) / (2 * a.mean() ** 2)
        series = dict(m2=a, m4=b, binder_influence=influence,
                      energy=data['action_density'][cut:], vortex=data['vortex_density'][cut:])
        taus = {name: exp.tau_int(x) for name, x in series.items()}
        largest = max(t for t in taus.values() if t is not None)
        width = max(20, int(np.ceil(10 * largest)))
        diagnostic_tau = {name: exp.tau_int(data[name][cut:])
                          for name in ('magnetization_real', 'magnetization_imag')}
        winding_squared = (data['winding'][cut:].astype(float) ** 2).mean(axis=(1, 2))
        diagnostic_tau['winding_squared'] = exp.tau_int(winding_squared)
        winding_changes = int(np.count_nonzero(np.any(np.diff(data['winding'][cut:], axis=0), axis=(1, 2))))
        result = bootstrap(a, b, width, seed=2026090601 + i, include_tail=include_tail)
        if result['binder'] is None:
            raise ValueError(f'{point}: insufficient blocks for any Binder error estimate')
        samples = result.pop('samples')
        np.save(exp.OUT / f"bootstrap_L{point['L']:02d}_k{point['kappa']:.2f}.npy", samples)
        row = point | result | dict(
            first_retained_sweep=cut, retained=len(a), remainder=len(a) - result['used'],
            block_width=width, tau_max=largest, tau=taus, diagnostic_tau=diagnostic_tau,
            winding_changes=winding_changes,
            local_acceptance=float(data['local_accepted'][cut:].mean() / point['L'] ** 2),
            winding_acceptance=float(data['winding_accepted'][cut:].mean()))
        rows.append(row)
        for other_cut in cuts:
            for factor in (1, 2, 4):
                check = bootstrap(data['magnetization_squared'][other_cut:],
                                  data['magnetization_fourth'][other_cut:], width * factor,
                                  seed=300000 + i * 10000 + other_cut + factor, draws=2000,
                                  include_tail=include_tail)
                check.pop('samples', None)
                stability.append(dict(L=point['L'], kappa=point['kappa'], cut=other_cut,
                                      block_factor=factor, block_width=width * factor) | check)
        for quarter, indices in enumerate(np.array_split(np.arange(len(a)), 4)):
            check = bootstrap(a[indices], b[indices], width, seed=400000 + i * 10 + quarter, draws=2000,
                              include_tail=include_tail)
            check.pop('samples', None)
            stability.append(dict(L=point['L'], kappa=point['kappa'], quarter=quarter + 1,
                                  block_width=width) | check)
        print(f"L={point['L']:2d} κ={point['kappa']:.2f}: U={row['binder']:.6f} ± {row['error']:.6f}; "
              f"estimated tau_max={largest:.1f}, width={width}, blocks={row['blocks']}", flush=True)
    rows.sort(key=lambda row: (row['L'], row['kappa']))
    cuts = sorted(set(r['cut'] for r in stability if 'cut' in r))
    exp.atomic_json(exp.OUT / 'results.json', rows)
    exp.atomic_json(exp.OUT / 'stability.json', stability)
    exp.atomic_json(exp.OUT / 'analysis-provenance.json', dict(
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        dataset_sha256=hashlib.sha256((exp.OUT / 'dataset.py').read_bytes()).hexdigest(),
        updated_sizes=list(selected), updated_cut=cut,
        updated_include_tail=include_tail,
        first_retained_sweep_by_L={str(r['L']): r['first_retained_sweep'] for r in rows},
        previous_analysis='before-L32/analysis-provenance.json',
        bootstrap_draws=4000, tau_method='nonperiodic FFT, initial-positive monotone paired ACF',
        block_rule='max(20, ceil(10*max(tau_m2,tau_m4,tau_Binder_influence,tau_energy,tau_vortex)))',
        errors='standard deviation of joint-block bootstrap ratio; no thinning',
        convention='U=1-mean(m4)/(2*mean(m2)**2)',
        note='global rotation and winding histories retained as separate diagnostics'))
    fields = list(dict.fromkeys(key for row in rows for key in row if key not in ('tau', 'diagnostic_tau')))
    with (exp.OUT / 'binder.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    plot(rows, stability, cuts, cut)


def draw_series(ax, rows, *, zoom=False):
    for L in exp.SIZES:
        group = [row for row in rows if row['L'] == L]
        if zoom:
            group = [row for row in group if row['kappa'] >= 1.05]
        ax.errorbar([r['kappa'] for r in group], [r['binder'] for r in group],
                    yerr=[r['error'] for r in group], color=exp.COLORS[L],
                    marker=exp.MARKERS[L], lw=1.5, ms=5.5, capsize=3, label=f'L = {L}')
    ax.axvline(1.12, color='.55', lw=1, ls=':', zorder=0)
    ax.grid(alpha=.15)
    ax.set_xlabel(r'Coupling $\kappa$')


def plot(rows, stability, cuts, cut):
    # The requested all-sizes plot, with a separate magnified view of the upper branch.
    fig, ax = plt.subplots(figsize=(8.6, 5.9))
    draw_series(ax, rows)
    ax.set_ylabel(r'$U_L = 1-\langle |m|^4\rangle/(2\langle |m|^2\rangle^2)$')
    ax.set_title(r'XY Binder cumulant on periodic $L\times L$ lattices', pad=15)
    ax.set_xticks(exp.KAPPAS)
    ax.set_ylim(min(0, min(r['binder']-r['error'] for r in rows)-.015), .515)
    ax.legend(frameon=False, loc='lower right')
    ax.text(.20, .06, r'$\epsilon=\pi,\quad S=-\kappa\sum_\ell\cos\theta_\ell$',
            transform=ax.transAxes, fontsize=11, color='.25')
    ax.text(1.13, .14, r'$\kappa_{\rm BKT}\simeq1.12$'+'\n(infinite-volume reference)', fontsize=9, color='.4')
    cut_by_L = {r['L']: r['first_retained_sweep'] - 1 for r in rows}
    if len(set(cut_by_L.values())) == 1:
        cut_label = f'discard first {next(iter(cut_by_L.values()))}'
    else:
        cut_label = '; '.join('discard ' + str(c) + ' for L=' + ','.join(str(L) for L in cut_by_L if cut_by_L[L] == c)
                              for c in sorted(set(cut_by_L.values())))
    fig.text(.5, .025, f'10,000 sweeps per point • {cut_label} • 1σ block bootstrap',
             ha='center', fontsize=9, color='.3')
    fig.tight_layout(rect=(0, .055, 1, 1))
    for suffix in ('png', 'pdf', 'svg'):
        fig.savefig(exp.OUT / f'binder_vs_kappa.{suffix}')
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    draw_series(ax, rows, zoom=True)
    subset = [r for r in rows if r['kappa'] >= 1.05]
    ax.set_ylim(min(r['binder']-r['error'] for r in subset)-.002, .501)
    ax.set_ylabel(r'Binder cumulant $U_L$')
    ax.set_title('Magnified view near and above the expected BKT transition')
    ax.legend(frameon=False, loc='lower right')
    fig.tight_layout()
    for suffix in ('png', 'pdf'):
        fig.savefig(exp.OUT / f'binder_zoom.{suffix}')
    plt.close(fig)
    fig, axes = plt.subplots(1, len(exp.SIZES), figsize=(4 * len(exp.SIZES), 4.5))
    for ax, L in zip(axes, exp.SIZES):
        for other_cut, marker in zip(cuts, ('o', 's', '^', 'D')):
            group = [r for r in stability if r['L'] == L and r.get('cut') == other_cut
                     and r.get('block_factor') == 1 and r['binder'] is not None]
            ax.errorbar([r['kappa'] for r in group], [r['binder'] for r in group],
                        yerr=[r['error'] for r in group], lw=.8, ms=3, marker=marker,
                        capsize=2, alpha=.85, label=f'discard {other_cut-1}')
        ax.set(title=f'L = {L}', xlabel=r'$\kappa$')
        ax.grid(alpha=.15)
    axes[0].set_ylabel('Binder cumulant')
    axes[-1].legend(fontsize=8, frameon=False)
    fig.suptitle('Thermalization-cut sensitivity (same-chain estimates)')
    fig.tight_layout()
    fig.savefig(exp.OUT / 'cut_stability.png')
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for L in exp.SIZES:
        group = [r for r in rows if r['L'] == L]
        axes[0].plot([r['kappa'] for r in group], [r['tau']['binder_influence'] for r in group],
                     marker=exp.MARKERS[L], c=exp.COLORS[L], label=f'L = {L}')
        axes[1].plot([r['kappa'] for r in group], [r['blocks'] for r in group],
                     marker=exp.MARKERS[L], c=exp.COLORS[L])
    axes[0].set(xlabel=r'$\kappa$', ylabel=r'$\tau_{\rm int}$ (sweeps)', title='Binder influence-series autocorrelation')
    axes[1].set(xlabel=r'$\kappa$', ylabel='Number of complete blocks', title='Blocks used for uncertainties')
    axes[0].legend(frameon=False)
    for ax in axes:
        ax.grid(alpha=.15)
    fig.tight_layout()
    fig.savefig(exp.OUT / 'autocorrelation_blocks.png')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--histories', action='store_true')
    parser.add_argument('--cut', type=int)
    parser.add_argument('--sizes', type=int, nargs='+', choices=exp.SIZES,
                        help='Analyze these sizes and preserve existing results for other sizes')
    parser.add_argument('--include-tail', action='store_true',
                        help='Use balanced blocks to include all measurements after the cut')
    args = parser.parse_args()
    if args.histories:
        exp.histories(args.sizes)
    elif args.cut is not None and 1 <= args.cut < 10001:
        analyze(args.cut, args.sizes, args.include_tail)
    else:
        parser.error('Use --histories or specify the first retained sweep with --cut')
