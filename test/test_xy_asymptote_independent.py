"""Independent isolated-producer routing, joint replay and source checks."""
from pathlib import Path
import hashlib
import importlib.util
import inspect
import json

import h5py
import numpy as np
import pytest

P = Path(__file__).resolve().parents[1]/'experiments/2026-10-17-xy-3d-negative-kappa/asymptote'
def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, P/filename)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def test_isolated_tiny_checkerboard_worker_and_fft_normalization(tmp_path, monkeypatch):
    runner = load('independent_asymptote_runner', 'runner.py')
    assert Path(inspect.getsourcefile(runner.audit)).resolve() == P/'audit.py'
    assert Path(runner.pt.__file__).resolve() == P/'pt.py'
    sources = runner.sources()
    point = runner.point(-8., .52, 6, 73, steps=6)
    point.update(kappas=[-8.], initial_winding=[0, 0, 0], initial_kind='checkerboard')
    for directory in ['jobs', 'chains']:
        (tmp_path/directory).mkdir()
    manifest = tmp_path/'batch.json'
    manifest.write_text(json.dumps(dict(points=[point], source_hashes=sources)))
    monkeypatch.setattr(runner, 'P', tmp_path)
    monkeypatch.setattr(runner, 'route', lambda: None)
    monkeypatch.setattr(runner, 'sources', lambda: sources)
    runner.worker(manifest, point['id'])
    job = json.loads((tmp_path/'jobs'/f"{point['id']}.json").read_text())
    assert job['state'] == 'complete' and job['audit']['passed']
    with h5py.File(tmp_path/'chains'/f"{point['id']}.h5") as h:
        spin = np.exp(1j*h['phi'][:, 0, 0])
        spectrum = abs(np.fft.fftn(spin, axes=(1, 2, 3)))**2/6**6
        np.testing.assert_allclose(h['staggered_squared'][:, 0], spectrum[:, 3, 3, 3], atol=2e-15)
        np.testing.assert_allclose(h['magnetization_squared'][:, 0], spectrum[:, 0, 0, 0], atol=2e-15)


def test_checkerboard_two_slot_odd_epoch_joint_replay(tmp_path):
    pt = load('independent_asymptote_pt', 'pt.py')
    kw = dict(N=6, epsilon=.52*np.pi, kappas=[-8., -4.], warmup=1,
              seed=997, initial_winding=[0, 0, 0], initial_kind='checkerboard', checkpoint_every=3)
    whole, interrupted = tmp_path/'a.h5', tmp_path/'b.h5'
    pt.run(whole, steps=9, **kw)
    with pytest.raises(RuntimeError, match='Injected crash'):
        pt.run(interrupted, steps=9, inject_crash_after_flush=5, **kw)
    with h5py.File(interrupted) as h:
        assert h.attrs['committed_rows'] == 4
    pt.run(interrupted, steps=9, resume=True, **kw)
    with h5py.File(whole) as a, h5py.File(interrupted) as b:
        for key in a:
            if not key.endswith('_seconds'):
                np.testing.assert_array_equal(a[key][:], b[key][:], err_msg=key)


def test_source_refusal_does_not_truncate_uncommitted_tail(tmp_path):
    pt = load('independent_asymptote_source_pt', 'pt.py')
    kw = dict(N=6, epsilon=.5*np.pi, kappas=[-4.], warmup=1,
              seed=1003, initial_winding=[0, 0, 0], initial_kind='checkerboard', checkpoint_every=2)
    file = tmp_path/'source.h5'; pt.run(file, steps=4, **kw)
    with h5py.File(file, 'r+') as h:
        source = json.loads(h.attrs['source_hashes']); source['replica_exchange_sha256'] = 'unreviewed'
        h.attrs['source_hashes'] = json.dumps(source)
        h['nematic_squared'].resize(7, axis=0)
    before = hashlib.sha256(file.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match='parameters/source mismatch'):
        pt.run(file, steps=6, resume=True, **kw)
    assert hashlib.sha256(file.read_bytes()).hexdigest() == before
