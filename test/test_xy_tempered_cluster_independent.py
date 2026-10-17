"""Independent path-level balance oracle for the auxiliary tempered move."""
from pathlib import Path
import importlib.util
import itertools
import sys

import h5py
import numpy as np
import pytest

P = Path(__file__).resolve().parents[1]/'experiments/2026-10-17-xy-3d-negative-kappa'
sys.path.insert(0, str(P))
import tempered_cluster as t


def test_exact_accepted_rejected_path_matrix_nonvacuous():
    # Three finite states, fixed support and deliberately non-equilibrium legs.
    # Each component is reversible, but their bare composite is not at target.
    C = np.array([-1.7, .2, 1.4])
    target = -.65
    levels = np.array([-.35, 0., -.35])
    proposal = np.array([[0., .2, .15], [.2, 0., .25], [.15, .25, 0.]])
    kernels = []
    for k in levels:
        K = proposal*np.minimum(1., np.exp(k*(C[None, :]-C[:, None])))
        K[np.diag_indices(3)] = 1-K.sum(1)
        w = np.exp(k*C); w /= w.sum()
        np.testing.assert_allclose(w[:, None]*K, (w[:, None]*K).T, atol=2e-15)
        kernels.append(K)
    transition = np.zeros((3, 3))
    saw_rejection = saw_full_acceptance = False
    for start in range(3):
        for tail in itertools.product(range(3), repeat=len(levels)):
            path = [start, *tail]
            forward = np.prod([K[a, b] for K, a, b in zip(kernels, path[:-1], path[1:])])
            reverse = np.prod([K[b, a] for K, a, b in zip(kernels, path[:-1], path[1:])])
            # Palindrome makes these the physical reverse-path kernels too.
            oracle_log_ratio = target*(C[path[-1]]-C[start])+np.log(reverse/forward)
            actual = t.path_log_ratio(target, levels, C[path])
            assert actual == pytest.approx(oracle_log_ratio, abs=2e-15)
            accept = min(1., np.exp(actual))
            saw_rejection |= accept < .9
            saw_full_acceptance |= accept == 1.
            transition[start, path[-1]] += forward*accept
            transition[start, start] += forward*(1-accept)
    assert saw_rejection and saw_full_acceptance
    weights = np.exp(target*C); weights /= weights.sum()
    np.testing.assert_allclose(transition.sum(1), 1., atol=2e-15)
    np.testing.assert_allclose(weights[:, None]*transition,
                               (weights[:, None]*transition).T, atol=3e-15)
    bare = kernels[0]@kernels[1]@kernels[2]
    assert np.max(abs(weights[:, None]*bare-(weights[:, None]*bare).T)) > .001


def test_tampered_tempered_definition_refuses_before_file_writes(tmp_path):
    import hashlib
    import json
    file = tmp_path/'tempered.h5'
    kw = dict(D=2, N=4, kappa=-.25, epsilon=2.1, seed=102,
              warmup=2, update='hybrid', clusters_per_step=4,
              winding=True, parallel=False, checkpoint_every=5, buffer_rows=5)
    t.run_chain(file, sweeps=10, **kw)
    with h5py.File(file, 'r+') as h:
        meta = json.loads(h.attrs['provenance'])
        meta['tempered_cluster_extension']['attempt_probability'] = .9
        h.attrs['provenance'] = json.dumps(meta)
    before = hashlib.sha256(file.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match='Tempered extension mismatch'):
        t.run_chain(file, sweeps=15, resume=True, **kw)
    assert hashlib.sha256(file.read_bytes()).hexdigest() == before
