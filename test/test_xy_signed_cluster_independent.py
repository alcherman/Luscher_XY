"""Independent signed-reflection balance and provenance regression checks."""
from pathlib import Path
import hashlib
import itertools
import json
import sys

import h5py
import numpy as np
import pytest

P = Path(__file__).resolve().parents[1] / 'experiments/2026-10-17-xy-3d-negative-kappa'
sys.path.insert(0, str(P))
import signed_cluster as signed
from supervillain import XY, Lattice
from supervillain.generator.xy import ClusterUpdate
from supervillain import xy


def test_frustrated_ring_hard_and_antiferromagnetic_bonds_exact_balance():
    # One hard ferromagnetic edge and three finite antiferromagnetic edges:
    # their preferred Ising signs cannot all be satisfied around this loop.
    kappa, epsilon = -.7, 1.8
    base = np.array([.2, .4, 1.3, 1.4])
    reflected = np.pi - base
    bits = np.array(list(itertools.product([0, 1], repeat=4)))
    fields = np.where(bits, reflected, base)
    wrap = lambda x: (x + np.pi) % (2*np.pi) - np.pi
    links = wrap(np.roll(fields, -1, axis=1) - fields)
    allowed = np.flatnonzero(np.all(abs(links) < epsilon, axis=1))
    energies = -kappa * np.cos(links[allowed]).sum(1)
    assert len(allowed) == 8 and np.ptp(energies) > .1
    action = XY(Lattice(1, 4), kappa, epsilon)
    update = signed.SignedCluster(action)
    transition = np.zeros((8, 8))
    positions = {int(state): i for i, state in enumerate(allowed)}
    saw_hard = saw_variable = False
    for state in allowed:
        field = fields[state]
        probs = []
        for i in range(4):
            j = (i+1) % 4
            single = reflected[i] if not bits[state, i] else base[i]
            if abs(wrap(single-field[j])) >= epsilon:
                prob = 1.
                saw_hard = True
            else:
                delta = kappa*(np.cos(field[i]-field[j])-np.cos(single-field[j]))
                prob = -np.expm1(-max(0., delta))
                saw_variable |= 0 < prob < 1
            probs.append(prob)
        probs = np.array(probs)
        variable = np.flatnonzero((probs > 0) & (probs < 1))
        for choices in itertools.product([False, True], repeat=len(variable)):
            active = probs == 1
            active[variable] = choices
            probability = np.where(active, probs, 1-probs).prod()
            tape = np.where(active, probs/2, (1+probs)/2)[None, :]
            for seed in range(4):
                out = np.asarray(update.cluster(action.canonical(field), seed, 0., tape)).ravel()
                target = int(np.argmin(np.max(abs(fields-out), axis=1)))
                np.testing.assert_allclose(out, fields[target], atol=2e-15, rtol=0)
                assert target in positions
                transition[positions[state], positions[target]] += probability/4
    assert saw_hard and saw_variable
    weights = np.exp(-energies); weights /= weights.sum()
    np.testing.assert_allclose(transition.sum(1), 1, atol=3e-13, rtol=0)
    flux = weights[:, None]*transition
    np.testing.assert_allclose(flux, flux.T, atol=3e-13, rtol=0)


@pytest.mark.parametrize('kappa', [0., .3, 2.])
def test_positive_coupling_native_supplied_randomness_identity(kappa):
    action = XY(Lattice(2, 3), kappa, 2.1)
    native, extension = ClusterUpdate(action), signed.SignedCluster(action)
    rng = np.random.default_rng(9145)
    phi = action.canonical(rng.uniform(0, .4, (3, 3)))
    for _ in range(20):
        seed = int(rng.integers(9)); alpha = rng.uniform(0, 2*np.pi)
        tape = rng.random((2, 9))
        a = native.cluster(phi, seed, alpha, tape)
        b = extension.cluster(phi, seed, alpha, tape)
        np.testing.assert_array_equal(np.asarray(a), np.asarray(b))
        assert native.last_cluster_size == extension.last_cluster_size
        assert native.last_bond_tests == extension.last_bond_tests
        phi = a


def test_positive_native_resume_and_tampered_extension_refuse_without_writes(tmp_path):
    path = tmp_path/'signed.h5'
    kw = dict(D=2, N=4, kappa=.2, epsilon=2.1, seed=419, warmup=2,
              update='hybrid', clusters_per_step=4, parallel=False,
              checkpoint_every=5, buffer_rows=5)
    signed.run_chain(path, sweeps=10, **kw)
    sha = lambda: hashlib.sha256(path.read_bytes()).hexdigest()
    before = sha()
    with pytest.raises((ValueError, RuntimeError)):
        xy.run_chain(path, sweeps=15, resume=True, **kw)
    assert sha() == before
    with h5py.File(path, 'r+') as h:
        provenance = json.loads(h.attrs['provenance'])
        provenance['signed_cluster_extension']['kernel_sha256'] = 'different'
        h.attrs['provenance'] = json.dumps(provenance)
    before = sha()
    with pytest.raises(ValueError, match='extension mismatch'):
        signed.run_chain(path, sweeps=15, resume=True, **kw)
    assert sha() == before
