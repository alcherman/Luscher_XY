"""Disk-only size dependence and representative autocorrelation functions."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np

import dataset as ds

OUT = Path(__file__).resolve().parent


def correlation(x):
    x = np.asarray(x, dtype=float)
    x = x - x.mean()
    n = len(x)
    size = 1 << (2 * n - 1).bit_length()
    ft = np.fft.rfft(x, n=size)
    cov = np.fft.irfft(ft * ft.conj(), n=size)[:n] / np.arange(n, 0, -1)
    return cov / cov[0]


def series(data, cut):
    a = data['magnetization_squared'][cut:]
    b = data['magnetization_fourth'][cut:]
    influence = b.mean() / a.mean()**3 * (a-a.mean()) - (b-b.mean()) / (2*a.mean()**2)
    return dict(m2=a, m4=b, binder_influence=influence,
                energy=data['action_density'][cut:], vortex=data['vortex_density'][cut:])


def main():
    ds.style()
    rows = json.loads((OUT / 'results.json').read_text())
    records, stability = [], []
    for row in rows:
        record = {k: row[k] for k in ('L', 'kappa', 'seed', 'first_retained_sweep', 'retained', 'blocks', 'tau_max')}
        record.update({f'tau_{k}': v for k, v in row['tau'].items()})
        record.update({f'tau_{k}': v for k, v in row['diagnostic_tau'].items()})
        record['binder_effective_samples_estimate'] = row['retained'] / (2 * row['tau']['binder_influence'])
        records.append(record)
        data = ds.read(row)
        for cut in (501, 1001, 2001):
            times = {name: ds.tau_int(x) for name, x in series(data, cut).items()}
            stability.append(dict(L=row['L'], kappa=row['kappa'], first_retained_sweep=cut, tau=times))
            if cut == row['first_retained_sweep']:
                for name, tau in times.items():
                    np.testing.assert_allclose(tau, row['tau'][name], rtol=1e-12)
    with (OUT / 'autocorrelation_times.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    ds.atomic_json(OUT / 'autocorrelation_cut_checks.json', stability)

    colors = dict(zip(ds.KAPPAS, plt.get_cmap('tab10').colors))
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.7))
    handles = []
    for ax, key, title in zip(axes, ('binder_influence', 'max'),
            ('Binder estimator', 'Maximum used to choose block width')):
        for kappa in ds.KAPPAS:
            group = sorted((r for r in rows if r['kappa'] == kappa), key=lambda r: r['L'])
            ys = [r['tau_max'] if key == 'max' else r['tau'][key] for r in group]
            line, = ax.plot([r['L'] for r in group], ys, color=colors[kappa], lw=1.5,
                            label=fr'$\kappa={kappa:.2f}$')
            if ax is axes[0]:
                handles.append(line)
            for row, y in zip(group, ys):
                ax.plot(row['L'], y, 'o', ms=6, color=colors[kappa],
                        markerfacecolor=colors[kappa] if row['enough_blocks'] else 'white')
        grid = np.linspace(4, 32, 100)
        ax.plot(grid, grid**2 / 8, '--', c='.5', lw=1, zorder=0)
        ax.text(22, 22**2 / 8 * 1.2, r'$\propto L^2$', fontsize=10, color='.4')
        ax.set_xscale('log', base=2)
        ax.set_yscale('log')
        ax.set_xticks(ds.SIZES, labels=ds.SIZES)
        ax.set_yticks((2, 5, 10, 20, 50, 100, 200))
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.set(xlabel='Lattice size L', ylabel=r'$\tau_{\rm int}$ (sweeps)', title=title, ylim=(1.7, 300))
        ax.grid(which='major', alpha=.18)
    fig.suptitle('Integrated autocorrelation times versus lattice size', y=.98, fontsize=16)
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, .91), ncol=7,
               frameon=False, fontsize=10, columnspacing=1)
    fig.text(.5, .055, '9,500 post-cut sweeps per chain • open markers: fewer than 20 analysis blocks',
             ha='center', fontsize=10, color='.3')
    fig.text(.5, .017, 'Times are estimates without error bars; dashed lines are slope guides, not fits.',
             ha='center', fontsize=9, color='.4')
    fig.subplots_adjust(top=.76, bottom=.19, left=.075, right=.98, wspace=.24)
    for suffix in ('png', 'pdf', 'svg'):
        fig.savefig(OUT / f'autocorrelation_vs_L.{suffix}')
    plt.close(fig)

    # Raw normalized ACFs of two observables near the expected transition.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    acf_records = []
    for row in rows:
        if row['kappa'] != 1.12:
            continue
        s = series(ds.read(row), row['first_retained_sweep'])
        for ax, key, title in zip(axes, ('m2', 'binder_influence'),
                                  (r'Magnetization magnitude squared $|m|^2$', 'Binder influence series')):
            rho = correlation(s[key])
            # Independent direct-lag checks of the FFT normalization.
            centered = s[key] - s[key].mean()
            for lag in (0, 1, 10, 100):
                direct = np.mean(centered[:len(centered)-lag] * centered[lag:]) / np.mean(centered**2)
                np.testing.assert_allclose(rho[lag], direct, rtol=1e-11, atol=1e-14)
            ax.plot(np.arange(501), rho[:501], color=ds.COLORS[row['L']],
                    label=f"L = {row['L']}", lw=1.3)
            ax.set(title=title, xlabel='Lag (sweeps)', ylabel=r'Normalized autocorrelation $\rho(t)$',
                   xlim=(0, 500), ylim=(-.25, 1.03))
            for lag in range(501):
                acf_records.append(dict(L=row['L'], kappa=1.12, observable=key, lag=lag, rho=float(rho[lag])))
    for ax in axes:
        ax.axhline(0, color='.5', lw=.7)
        ax.grid(alpha=.15)
    axes[0].legend(frameon=False)
    fig.suptitle(r'Autocorrelation functions near the transition: $\kappa=1.12$', fontsize=15)
    fig.text(.5, .015, 'Post-cut data; long-lag wiggles are noisy finite-run estimates. No exponential fit.',
             ha='center', fontsize=9, color='.4')
    fig.tight_layout(rect=(0, .05, 1, .94))
    for suffix in ('png', 'pdf'):
        fig.savefig(OUT / f'autocorrelation_functions_k1.12.{suffix}')
    plt.close(fig)
    with (OUT / 'autocorrelation_functions_k1.12.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(acf_records[0]))
        writer.writeheader()
        writer.writerows(acf_records)
    ds.atomic_json(OUT / 'autocorrelation-report-provenance.json', dict(
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        results_sha256=hashlib.sha256((OUT / 'results.json').read_bytes()).hexdigest(),
        checks='All stored primary times reproduced; FFT ACF checked against direct lag products',
        new_monte_carlo_sweeps=0, tau_error_bars_estimated=False, dynamical_exponent_fitted=False))
    print('Saved size-dependence plots, near-critical ACFs, complete time table and cut checks. No new sweeps.')


if __name__ == '__main__':
    main()
