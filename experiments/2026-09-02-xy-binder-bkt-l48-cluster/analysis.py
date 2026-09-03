"""Disk-only L48 cluster histories, blocked Binder estimates, and comparison plots."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

from run import OUT, PREVIOUS, KAPPAS, STEPS, atomic_json, points
sys.path.insert(0, str(PREVIOUS))
import experiment as previous
from analyze import bootstrap
import h5py
import matplotlib.pyplot as plt
import numpy as np


def read(point):
    with h5py.File(OUT / point['file'], 'r') as f:
        assert f.attrs['committed_rows'] == STEPS+1
        assert all(len(ds) == STEPS+1 for ds in f.values())
        return {key: ds[:] for key, ds in f.items() if key not in ('phi', 'rng_state')}


def histories():
    previous.style()
    directory = OUT / 'histories'
    directory.mkdir(exist_ok=True)
    all_data = [read(p) for p in points()]
    for data, point in zip(all_data, points()):
        keys = [key for key in data if key != 'sweep']
        fig, axes = plt.subplots((len(keys)+2)//3, 3, figsize=(15, 17), squeeze=False)
        for ax, key in zip(axes.flat, keys):
            ax.plot(data['sweep'], data[key].reshape(STEPS+1, -1), lw=.35, alpha=.55)
            ax.set_title(key, fontsize=9)
            ax.set_xlabel('cluster cycle (cold seed at 0)', fontsize=8)
        for ax in list(axes.flat)[len(keys):]:
            ax.axis('off')
        fig.suptitle(f"L=48, κ={point['kappa']:.2f}, ε=π: all saved observables")
        fig.tight_layout()
        fig.savefig(directory / f"k{point['kappa']:.2f}.png")
        plt.close(fig)
    for stop, name in ((3001, 'thermalization'), (STEPS+1, 'overview')):
        fig, axes = plt.subplots(7, 3, figsize=(13, 16), squeeze=False)
        for row, (data, point) in enumerate(zip(all_data, points())):
            for col, (key, label) in enumerate((('magnetization_squared', r'$|m|^2$'),
                        ('action_density', r'$S/L^2$'), ('vortex_density', 'vortex density'))):
                ax = axes[row, col]
                t, y = data['sweep'][:stop], data[key][:stop]
                ax.plot(t, y, lw=.45, alpha=.55, color='#9467bd')
                bins = (len(t)-1)//100
                ax.plot(t[1:1+100*bins].reshape(-1,100).mean(axis=1),
                        y[1:1+100*bins].reshape(-1,100).mean(axis=1), lw=1, color='black')
                ax.plot(0, y[0], 'ko', ms=3)
                if col == 0:
                    ax.set_yscale('log')
                ax.set_title(f"κ={point['kappa']:.2f}: {label}", fontsize=10)
                ax.grid(alpha=.15)
        for ax in axes[-1]:
            ax.set_xlabel('cluster cycle: one cluster + one winding proposal')
        fig.suptitle('L=48 cold-start histories; black line = 100-cycle mean')
        fig.tight_layout(rect=(0,0,1,.975))
        fig.savefig(OUT / f'{name}.png')
        plt.close(fig)


def diagnostics(data, cut):
    a, b = data['magnetization_squared'][cut:], data['magnetization_fourth'][cut:]
    h = b.mean()/a.mean()**3*(a-a.mean()) - (b-b.mean())/(2*a.mean()**2)
    series = dict(m2=a, m4=b, binder_influence=h,
                  energy=data['action_density'][cut:], vortex=data['vortex_density'][cut:])
    return {key: previous.tau_int(value) for key, value in series.items()}


def analyze(cut):
    previous.style()
    rows, checks = [], []
    for i, point in enumerate(points()):
        data = read(point)
        taus = diagnostics(data, cut)
        tau = max(t for t in taus.values() if t is not None)
        width = max(20, int(np.ceil(10*tau)))
        a, b = data['magnetization_squared'][cut:], data['magnetization_fourth'][cut:]
        result = bootstrap(a, b, width, seed=2026090701+i, include_tail=True)
        samples = result.pop('samples', None)
        if samples is not None:
            np.save(OUT / f"bootstrap_L48_k{point['kappa']:.2f}.npy", samples)
        # A central ratio remains defined even when the run has too few blocks
        # to quote an uncertainty. Never substitute a zero error in that case.
        if result['binder'] is None:
            result['binder'] = float(1-b.mean()/(2*a.mean()**2))
        fraction = float(data['cluster_sizes'][cut:].mean()/48**2)
        row = point | result | dict(update='cluster', first_retained_cycle=cut,
                    retained=len(a), binder_samples=len(a), tau=taus, tau_max=tau, block_width=width,
                    mean_cluster_fraction=fraction,
                    tau_m2_cluster_work=taus['m2']*fraction,
                    total_cluster_work=float(data['cluster_equivalent_sweeps'][-1]),
                    winding_acceptance=float(data['winding_accepted'][cut:].mean()))
        rows.append(row)
        for other in sorted(set([cut, 501, 1001, 2001, 4001])):
            other_taus = diagnostics(data, other)
            base = max(20, int(np.ceil(10*max(t for t in other_taus.values() if t is not None))))
            for factor in (1,2,4):
                check = bootstrap(data['magnetization_squared'][other:], data['magnetization_fourth'][other:],
                                  base*factor, seed=2026090800+i*100+other+factor,
                                  draws=2000, include_tail=True)
                check.pop('samples', None)
                if check['binder'] is None:
                    x,y = data['magnetization_squared'][other:],data['magnetization_fourth'][other:]
                    check['binder'] = float(1-y.mean()/(2*x.mean()**2))
                checks.append(dict(kappa=point['kappa'], cut=other, factor=factor,
                                   block_width=base*factor, tau=other_taus) | check)
        for quarter, indices in enumerate(np.array_split(np.arange(len(a)), 4), 1):
            check = bootstrap(a[indices], b[indices], width,
                              seed=2026091000+10*i+quarter, draws=2000, include_tail=True)
            check.pop('samples', None)
            if check['binder'] is None:
                check['binder'] = float(1-b[indices].mean()/(2*a[indices].mean()**2))
            checks.append(dict(kappa=point['kappa'], quarter=quarter, cut=cut) | check)
        print(f"kappa={point['kappa']:.2f}: U={row['binder']:.6f} error={row['error']}; "
              f"tau_m2={taus['m2']:.1f} cycles, max_tau={tau:.1f}, blocks={row['blocks']}, "
              f"cluster_fraction={fraction:.4f}", flush=True)
    atomic_json(OUT / 'results.json', rows)
    atomic_json(OUT / 'stability.json', checks)
    old = json.loads((PREVIOUS / 'results.json').read_text())
    combined = [dict(r, update='metropolis') for r in old] + rows
    atomic_json(OUT / 'combined-results.json', combined)
    fields = [key for key in rows[0] if key != 'tau']
    with (OUT / 'binder.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)
    plot(combined, checks, cut)
    atomic_json(OUT / 'analysis-provenance.json', dict(
        cut=cut, bootstrap_draws=4000, bootstrap='joint contiguous balanced block means; all retained rows used',
        block_rule='max(20, ceil(10*max(tau_m2,tau_m4,tau_Binder_influence,tau_energy,tau_vortex)))',
        clock='fixed cycle: one cluster + one winding proposal',
        work_clock='cluster membership/volume; excludes full-field checks, RNG, I/O, winding and measurement work',
        analysis_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        previous_results_sha256=hashlib.sha256((PREVIOUS / 'results.json').read_bytes()).hexdigest(),
        previous_bootstrap_sha256=hashlib.sha256((PREVIOUS / 'analyze.py').read_bytes()).hexdigest(),
        previous_tau_sha256=hashlib.sha256((PREVIOUS / 'experiment.py').read_bytes()).hexdigest()))


def plot(rows, checks, cut):
    colors = {4:'#0072B2',8:'#D55E00',16:'#009E73',32:'#CC79A7',48:'#6A3D9A'}
    markers = {4:'o',8:'s',16:'^',32:'D',48:'P'}
    for zoom in (False, True):
        fig, ax = plt.subplots(figsize=(9.5,6.3))
        for L in colors:
            group = [r for r in rows if r['L']==L and (not zoom or r['kappa']>=1.05)]
            good = [r for r in group if r['error'] is not None]
            label = f'L = {L}' + (' (cluster)' if L==48 else ' (Metropolis)')
            ax.errorbar([r['kappa'] for r in good], [r['binder'] for r in good],
                        yerr=[r['error'] for r in good], color=colors[L], marker=markers[L],
                        lw=1.6, ms=6, capsize=3, label=label)
            bad = [r for r in group if r['error'] is None]
            ax.scatter([r['kappa'] for r in bad], [r['binder'] for r in bad],
                       facecolors='none', edgecolors=colors[L], marker=markers[L], s=60)
            if bad and not zoom:
                point = bad[0]
                ax.annotate('Provisional: too few blocks\nfor an error estimate',
                            xy=(point['kappa'],point['binder']),
                            xytext=(point['kappa']+.07,point['binder']+.07),
                            fontsize=9,color=colors[L],
                            arrowprops=dict(arrowstyle='-',color=colors[L],lw=.8))
        ax.set(xlabel=r'Coupling $\kappa$', ylabel=r'$U_L=1-\langle |m|^4\rangle/(2\langle |m|^2\rangle^2)$',
               title=r'2D XY Binder cumulant, $\epsilon=\pi$')
        ax.set_xticks([k for k in KAPPAS if not zoom or k>=1.05])
        ax.axvline(1.12, color='.6', ls=':', lw=1)
        ax.grid(alpha=.15)
        ax.legend(frameon=False, fontsize=9, loc='lower right')
        if not zoom:
            ax.set_ylim(min(-.02,min(r['binder']-(r['error'] or 0) for r in rows)-.015), .515)
        fig.text(.5,.018,f'L48: 10,000 cluster cycles, discard {cut-1}; L4–32: existing Metropolis data\n'
                 'Error bars: 1σ joint-block bootstrap; cluster cycles and Metropolis sweeps have different work costs',
                 ha='center',fontsize=9,color='.3')
        fig.tight_layout(rect=(0,.075,1,1))
        for extension in ('png','pdf','svg'):
            fig.savefig(OUT / (('binder_zoom' if zoom else 'binder_vs_kappa')+'.'+extension))
        plt.close(fig)
    fig, axes = plt.subplots(1,3,figsize=(14,4.5))
    group = [r for r in rows if r['L']==48]
    axes[0].plot(KAPPAS,[r['tau']['m2'] for r in group],'o-',label=r'$|m|^2$')
    axes[0].plot(KAPPAS,[r['tau']['energy'] for r in group],'s-',label='energy')
    axes[0].set(ylabel=r'$\tau_{int}$ (cluster cycles)',yscale='log')
    axes[0].legend(frameon=False)
    axes[1].plot(KAPPAS,[r['mean_cluster_fraction'] for r in group],'o-',color=colors[48])
    axes[1].set(ylabel='Mean cluster size / volume')
    axes[2].plot(KAPPAS,[r['blocks'] for r in group],'o-',color=colors[48])
    axes[2].axhline(20,color='.5',ls=':')
    axes[2].set(ylabel='Bootstrap blocks')
    for ax in axes:
        ax.set_xlabel(r'$\kappa$'); ax.grid(alpha=.15)
    fig.suptitle('L48 cluster diagnostics; no equal-L Metropolis benchmark in this run')
    fig.tight_layout(); fig.savefig(OUT/'cluster_diagnostics.png'); plt.close(fig)
    fig, axes = plt.subplots(1,2,figsize=(12,4.7))
    for other in sorted({r['cut'] for r in checks}):
        group = [r for r in checks if r['cut']==other and r.get('factor')==1 and r['error'] is not None]
        axes[0].errorbar([r['kappa'] for r in group],[r['binder'] for r in group],
                         yerr=[r['error'] for r in group],marker='o',lw=.8,capsize=2,label=f'discard {other-1}')
    group = [r for r in checks if r['cut']==cut and 'factor' in r and r['error'] is not None]
    for factor in (1,2,4):
        g = [r for r in group if r['factor']==factor]
        axes[1].plot([r['kappa'] for r in g],[r['error'] for r in g],'o-',label=f'block width ×{factor}')
    axes[0].set_ylabel('Binder cumulant'); axes[1].set_ylabel('Bootstrap standard error')
    for ax in axes:
        ax.set_xlabel(r'$\kappa$'); ax.grid(alpha=.15); ax.legend(frameon=False,fontsize=8)
    fig.suptitle('L48 cut and block sensitivity (correlated same-chain estimates)')
    fig.tight_layout(); fig.savefig(OUT/'stability.png'); plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--histories', action='store_true')
    parser.add_argument('--cut', type=int)
    args = parser.parse_args()
    if args.histories:
        histories()
    if args.cut is not None:
        if not 1 <= args.cut < STEPS:
            raise ValueError('cut must retain samples after the seed')
        analyze(args.cut)
