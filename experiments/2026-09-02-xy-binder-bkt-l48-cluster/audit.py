"""Audit every saved field and observable in the completed L48 cluster data."""
import hashlib
import json
import h5py
import numpy as np
from run import OUT, PREVIOUS, STEPS, atomic_json, points


def checksum(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(4*1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def wrap(value):
    reduced = value % (2*np.pi)
    reduced = np.where(reduced >= np.pi, reduced-2*np.pi, reduced)
    return np.where((value >= -np.pi) & (value < np.pi), value, reduced)


def main():
    records = []
    for point in points():
        path = OUT / point['file']
        with h5py.File(path, 'r') as f:
            assert f.attrs['schema'] == 'supervillain.xy.v2'
            assert f.attrs['committed_rows'] == STEPS+1
            assert all(len(ds) == STEPS+1 for ds in f.values())
            params = json.loads(f.attrs['parameters'])
            assert params['D'] == 2 and params['N'] == 48 and params['epsilon'] == np.pi
            assert params['kappa'] == point['kappa'] and params['seed'] == point['seed']
            assert params['update'] == 'cluster' and params['clusters_per_step'] == 1
            assert params['winding'] and params['time_unit'] == 'cluster_cycle'
            assert np.array_equal(f['sweep'][:], np.arange(STEPS+1))
            assert np.array_equal(f['clusters_total'][:], np.arange(STEPS+1))
            assert np.all(f['phi'][0] == 0)
            sizes = f['cluster_sizes'][:,0]
            assert sizes[0] == 0 and np.all((sizes[1:] >= 1) & (sizes[1:] <= 48**2))
            assert np.array_equal(f['cluster_sites_total'][:], np.cumsum(sizes))
            assert np.array_equal(f['cluster_equivalent_sweeps'][:], np.cumsum(sizes)/48**2)
            assert np.array_equal(f['cluster_bond_tests_total'][:], np.cumsum(f['cluster_bond_tests'][:,0]))
            for start in range(0,STEPS+1,128):
                end = min(start+128,STEPS+1)
                a = f['phi'][start:end,0]
                assert np.isfinite(a).all() and np.all((a >= 0) & (a < 2*np.pi))
                dx,dy = wrap(np.roll(a,-1,axis=1)-a),wrap(np.roll(a,-1,axis=2)-a)
                assert np.all(np.abs(dx)<np.pi) and np.all(np.abs(dy)<np.pi)
                m = np.exp(1j*a).mean(axis=(1,2)); m2 = np.abs(m)**2
                cosine = (np.cos(dx).mean(axis=(1,2))+np.cos(dy).mean(axis=(1,2)))/2
                vortex = np.abs(np.rint((dx+np.roll(dy,-1,axis=1)-np.roll(dx,-1,axis=2)-dy)/(2*np.pi))).mean(axis=(1,2))
                wind = np.stack((np.rint(dx.sum(axis=1)/(2*np.pi)),np.rint(dy.sum(axis=2)/(2*np.pi))),axis=1)
                expected = dict(magnetization_real=m.real,magnetization_imag=m.imag,
                    magnetization_squared=m2,magnetization_fourth=m2*m2,
                    mean_link_cosine=cosine,action_density=-2*point['kappa']*cosine,
                    vortex_density=vortex,winding=wind)
                for key,value in expected.items():
                    np.testing.assert_allclose(f[key][start:end],value,rtol=2e-12,atol=2e-14,err_msg=key)
            for value in f['rng_state'][:]:
                state = json.loads(value)
                assert set(state) == {'local','winding'}
        records.append(point | dict(rows=STEPS+1,bytes=path.stat().st_size,sha256=checksum(path)))
        print(f"Audited all {STEPS+1} rows: kappa={point['kappa']:.2f}",flush=True)
    manifest = json.loads((OUT/'manifest.json').read_text())
    assert checksum(PREVIOUS/'results.json') == manifest['prior_results_sha256']
    atomic_json(OUT/'data-audit.json',dict(passed=True,points=records,total_rows=sum(r['rows'] for r in records),
        prior_results_unchanged=True,audit_source_sha256=checksum(__file_path := OUT/'audit.py')))


if __name__ == '__main__':
    main()
