"""Check every stored field against direct NumPy moments and bond sums.

Disk-only: does not run a Markov chain or change the raw HDF5 files.
"""
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

OUT = Path(__file__).resolve().parent
manifest = json.loads((OUT / 'manifest.json').read_text())
assert hashlib.sha256((OUT / 'experiment.py').read_bytes()).hexdigest() == manifest['experiment_sha256']
records = []
for point in manifest['points']:
    path = OUT / point['file']
    residuals = dict(m2=0., m4=0., action_density=0., max_link_angle=0.)
    with h5py.File(path, 'r') as f:
        pars = json.loads(f.attrs['parameters'])
        assert pars['D'] == 2 and pars['N'] == point['L']
        assert pars['kappa'] == point['kappa'] and pars['epsilon'] == np.pi
        assert pars['seed'] == point['seed'] and pars['winding']
        assert int(f.attrs['committed_rows']) == 10001
        assert all(ds.shape[0] == 10001 for ds in f.values())
        np.testing.assert_array_equal(f['sweep'][:], np.arange(10001))
        assert np.all(f['phi'][0] == 0)
        for begin in range(0, 10001, 512):
            stop = min(10001, begin + 512)
            phi = f['phi'][begin:stop, 0]
            assert np.all(np.isfinite(phi)) and np.all((0 <= phi) & (phi < 2 * np.pi))
            spin = np.cos(phi) + 1j * np.sin(phi)
            m2 = np.abs(spin.mean(axis=(1, 2))) ** 2
            bonds = np.stack([np.roll(phi, -1, axis=axis) - phi for axis in (1, 2)], axis=1)
            angle = (bonds + np.pi) % (2 * np.pi) - np.pi
            assert np.all(np.abs(angle) < np.pi)
            expected = dict(magnetization_squared=m2, magnetization_fourth=m2 ** 2,
                            action_density=-point['kappa'] * np.cos(bonds).sum(axis=(1, 2, 3)) / point['L'] ** 2,
                            max_link_angle=np.abs(angle).max(axis=(1, 2, 3)))
            for (name, values), short in zip(expected.items(), residuals):
                actual = f[name][begin:stop]
                np.testing.assert_allclose(actual, values, atol=3e-14, rtol=3e-14)
                residuals[short] = max(residuals[short], float(np.max(np.abs(actual - values))))
        source_sha256 = json.loads(f.attrs['provenance'])['source_sha256']
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    records.append(point | dict(rows_checked=10001, max_absolute_residual=residuals,
                                source_sha256=source_sha256, hdf5_sha256=digest))
assert len(records) == 21
result = dict(status='passed', total_stored_configurations_checked=21 * 10001,
              checks='parameters, row counts, cold seed, canonical range, strict admissibility, moments, action, max link angle',
              records=records)
(OUT / 'data-audit.json').write_text(json.dumps(result, indent=2) + '\n')
print(f"Passed: {len(records)} chains, {result['total_stored_configurations_checked']} configurations.")
