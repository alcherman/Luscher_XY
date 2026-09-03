"""Combined original and L32 manifests; disk-only data and plot helpers."""
import json

import matplotlib.pyplot as plt
import numpy as np

from experiment import OUT, KAPPAS, atomic_json, read, style, tau_int

SIZES = (4, 8, 16, 32)
COLORS = {4: '#0072B2', 8: '#D55E00', 16: '#009E73', 32: '#CC79A7'}
MARKERS = {4: 'o', 8: 's', 16: '^', 32: 'D'}


def points():
    rows = json.loads((OUT / 'manifest.json').read_text())['points']
    extension = OUT / 'manifest-l32.json'
    if extension.exists():
        rows += json.loads(extension.read_text())['points']
    keys = [(p['L'], p['kappa']) for p in rows]
    assert len(set(keys)) == len(keys)
    assert len(set(p['seed'] for p in rows)) == len(rows)
    return sorted(rows, key=lambda p: (p['L'], p['kappa']))


def histories(sizes=None):
    style()
    selected = SIZES if sizes is None else sizes
    directory = OUT / 'histories'
    directory.mkdir(exist_ok=True)
    for point in points():
        if point['L'] not in selected:
            continue
        data = read(point)
        keys = [k for k in data if k != 'sweep']
        fig, axes = plt.subplots((len(keys) + 2) // 3, 3, figsize=(15, 11), squeeze=False)
        for ax, key in zip(axes.flat, keys):
            ax.plot(data['sweep'], data[key].reshape(len(data['sweep']), -1), lw=.4, alpha=.6)
            ax.set_title(key, fontsize=9)
            ax.set_xlabel('sweep (cold seed at 0)', fontsize=8)
        for ax in list(axes.flat)[len(keys):]:
            ax.axis('off')
        fig.suptitle(f"XY histories: L={point['L']}, κ={point['kappa']:.2f}, ε=π")
        fig.tight_layout()
        fig.savefig(directory / f"L{point['L']:02d}_k{point['kappa']:.2f}.png")
        plt.close(fig)
    for L in selected:
        group = [p for p in points() if p['L'] == L]
        if not group:
            continue
        for end, name in ((2001, 'thermalization'), (10001, 'overview')):
            fig, axes = plt.subplots(len(group), 3, figsize=(12, 15), squeeze=False)
            for row, point in enumerate(group):
                data = read(point)
                for col, (key, label) in enumerate((('magnetization_squared', r'$|m|^2$'),
                        ('action_density', r'$S/L^2$'), ('vortex_density', 'vortex density'))):
                    ax = axes[row, col]
                    ax.plot(data['sweep'][:end], data[key][:end], lw=.5, c=COLORS[L], alpha=.7)
                    ax.plot(0, data[key][0], 'ko', ms=3)
                    ax.set_title(f"κ={point['kappa']:.2f}: {label}", fontsize=10)
                    ax.grid(alpha=.15)
            for ax in axes[-1]:
                ax.set_xlabel('sweep (cold seed at 0)')
            fig.suptitle(f'L = {L}: cold-start relaxation' if name == 'thermalization' else f'L = {L}: full histories')
            fig.tight_layout(rect=(0, 0, 1, .97))
            fig.savefig(OUT / f'{name}_L{L:02d}.png')
            plt.close(fig)
