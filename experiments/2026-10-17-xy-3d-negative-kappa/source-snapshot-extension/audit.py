"""Independent full-field audit adapted from validated legacy audit; Metropolis metadata and RNG keys."""
from pathlib import Path
import json,hashlib
import h5py,numpy as np
def sha256(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def audit(path, spec):
    """Independent full-field recomputation, generalized from prior 3D/4D audits."""
    D, N, steps = spec["D"], spec["N"], spec["steps"]
    volume = N**D
    spatial = tuple(range(1, D + 1))
    with h5py.File(path, "r") as f:
        pars = json.loads(f.attrs["parameters"])
        expected = dict(D=D, N=N, kappa=spec["beta"], epsilon=spec["epsilon"],
                        seed=spec["seed"], update="metropolis", parallel=False,
                        clusters_per_step=1, initial_winding=spec["initial_winding"],
                        winding=True, buffer_rows=100)
        for key, value in expected.items():
            assert pars[key] == value, (key, pars[key], value)
        assert f.attrs["committed_rows"] == steps + 1
        assert all(len(ds) == steps + 1 for ds in f.values())
        np.testing.assert_array_equal(f["sweep"][:], np.arange(steps + 1))
        coordinates = np.indices((N,) * D)
        initial = np.mod(sum(2 * np.pi * w * coordinates[a] / N for a, w in enumerate(spec["initial_winding"])), 2 * np.pi)
        np.testing.assert_array_equal(f["phi"][0, 0], initial)
        assert np.all((f["local_accepted"][:] >= 0) & (f["local_accepted"][:] <= volume))
        for key in f:
            if key.endswith("_seconds"):
                v = f[key][:]
                assert np.isfinite(v).all() and np.all(v >= 0)
        for state in f["rng_state"][:]:
            assert set(json.loads(state)) == {"local", "winding"}
        for first in range(0, steps + 1, 16):
            stop = min(first + 16, steps + 1)
            field = f["phi"][first:stop, 0]
            assert np.isfinite(field).all() and np.all((field >= 0) & (field < 2 * np.pi))
            links = np.stack([np.remainder(np.roll(field, -1, axis=a + 1) - field + np.pi, 2 * np.pi) - np.pi for a in range(D)], axis=1)
            assert np.all(np.abs(links) < spec["epsilon"])
            spin = np.exp(1j * field)
            magnetization = spin.mean(axis=spatial)
            m2 = np.abs(magnetization)**2
            cosine = np.cos(links).mean(axis=tuple(range(1, D + 2)))
            vortex_sum = np.zeros(stop - first)
            for a in range(D):
                for b in range(a + 1, D):
                    curl = links[:, a] + np.roll(links[:, b], -1, axis=a + 1) - np.roll(links[:, a], -1, axis=b + 1) - links[:, b]
                    vortex_sum += np.abs(np.rint(curl / (2 * np.pi))).mean(axis=spatial)
            winding = np.stack([np.rint(links[:, a].sum(axis=a + 1) / (2 * np.pi)).reshape(stop - first, -1) for a in range(D)], axis=1)
            minimum = []
            for a in range(D):
                shape = [1] * (D + 1)
                shape[a + 1] = N
                phase = np.exp(-2j * np.pi * np.arange(N) / N).reshape(shape)
                minimum.append(np.abs((spin * phase).sum(axis=spatial))**2 / volume)
            values = dict(magnetization_real=magnetization.real, magnetization_imag=magnetization.imag,
                          magnetization_abs=np.abs(magnetization), magnetization_squared=m2,
                          magnetization_fourth=m2 * m2, mean_link_cosine=cosine,
                          action_density=-D * spec["beta"] * cosine,
                          vortex_density=vortex_sum / (D * (D - 1) / 2), winding=winding,
                          structure_zero=volume * m2, structure_min=np.stack(minimum, axis=1),
                          max_link_angle=np.abs(links).max(axis=tuple(range(1, D + 2))))
            for key, value in values.items():
                np.testing.assert_allclose(f[key][first:stop], value, rtol=2e-12, atol=3e-13, err_msg=key)
    return dict(passed=True, rows=steps + 1, bytes=Path(path).stat().st_size, sha256=sha256(path))

