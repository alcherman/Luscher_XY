"""Independent joint-replica field, swap, history and checkpoint audits."""
from pathlib import Path
import json,hashlib
import h5py,numpy as np

def audit(path,spec,source_hashes):
    N=spec['N'];D=3;steps=spec['steps'];ks=np.asarray(spec['kappas']);R=len(ks);volume=N**3;spatial=(1,2,3)
    with h5py.File(path) as f:
        pars=json.loads(f.attrs['parameters']);assert pars['N']==N and pars['seed']==spec['seed'] and pars['epsilon']==spec['epsilon'] and pars['kappas']==spec['kappas'] and pars['initial_winding']==spec['initial_winding']
        assert json.loads(f.attrs['source_hashes'])==source_hashes
        assert f.attrs['committed_rows']==steps+1 and all(len(x)==steps+1 for x in f.values())
        np.testing.assert_array_equal(f['sweep'][:],np.arange(steps+1))
        sites=f['cluster_sizes'][:];tests=f['cluster_bond_tests'][:];assert np.all(sites[0]==0) and np.all((sites[1:]>=1)&(sites[1:]<=volume))
        np.testing.assert_array_equal(f['clusters_total'][:],np.repeat((4*np.arange(steps+1))[:,None],R,axis=1))
        np.testing.assert_array_equal(f['cluster_sites_total'][:],np.cumsum(sites.sum(2),axis=0));np.testing.assert_array_equal(f['cluster_bond_tests_total'][:],np.cumsum(tests.sum(2),axis=0))
        assert np.all((f['local_accepted'][:]>=0)&(f['local_accepted'][:]<=volume));assert np.all(np.isfinite(f['cycle_seconds'][:])) and np.all(f['cycle_seconds'][:]>=0)
        coords=np.indices((N,)*3);initial=np.mod(sum(2*np.pi*w*coords[a]/N for a,w in enumerate(spec['initial_winding'])),2*np.pi)
        np.testing.assert_array_equal(f['phi'][0],np.repeat(initial[None,None],R,axis=0))
        attempted=f['exchange_attempted'][:];accepted=f['exchange_accepted'][:];assert np.all((attempted==0)|(attempted==1)) and np.all((accepted==0)|(accepted==1)) and np.all(accepted<=attempted)
        for cycle in range(steps+1):
            state=json.loads(f['rng_state'][cycle]);assert len(state['slots'])==R and all(set(v)=={'local','cluster','winding'} for v in state['slots']);assert 'exchange' in state and state['next_exchange_parity']==cycle%2
            walker=f['walker_id'][cycle];np.testing.assert_array_equal(np.sort(walker),np.arange(R))
            if cycle==0:np.testing.assert_array_equal(walker,np.arange(R));continue
            before=f['pre_exchange_C'][cycle];after=before.copy();previous=f['walker_id'][cycle-1].copy();expected=np.zeros(R-1,int);expected[(cycle-1)%2::2]=1;np.testing.assert_array_equal(attempted[cycle],expected)
            for i in np.flatnonzero(expected):
                logratio=(ks[i]-ks[i+1])*(before[i+1]-before[i]);np.testing.assert_allclose(f['exchange_log_ratio'][cycle,i],logratio,rtol=1e-13,atol=1e-10)
                u=f['exchange_uniform'][cycle,i];assert 0<=u<1;assert bool(accepted[cycle,i])==(u==0 or np.log(u)<min(0.,logratio))
                if accepted[cycle,i]:after[i],after[i+1]=after[i+1],after[i];previous[i],previous[i+1]=previous[i+1],previous[i]
            np.testing.assert_allclose(f['post_exchange_C'][cycle],after,rtol=2e-12,atol=1e-9);np.testing.assert_array_equal(walker,previous)
        np.testing.assert_allclose(f['post_exchange_C'][:],3*volume*f['mean_link_cosine'][:],rtol=2e-12,atol=1e-9)
        for first in range(0, steps + 1, 2):
            stop = min(first + 2, steps + 1)
            field = f["phi"][first:stop, :, 0].reshape((-1,)+(N,)*D)
            assert np.isfinite(field).all() and np.all((field >= 0) & (field < 2 * np.pi))
            links = np.stack([np.remainder(np.roll(field, -1, axis=a + 1) - field + np.pi, 2 * np.pi) - np.pi for a in range(D)], axis=1)
            assert np.all(np.abs(links) < spec["epsilon"])
            spin = np.exp(1j * field)
            magnetization = spin.mean(axis=spatial)
            m2 = np.abs(magnetization)**2
            cosine = np.cos(links).mean(axis=tuple(range(1, D + 2)))
            vortex_sum = np.zeros((stop-first)*R)
            for a in range(D):
                for b in range(a + 1, D):
                    curl = links[:, a] + np.roll(links[:, b], -1, axis=a + 1) - np.roll(links[:, a], -1, axis=b + 1) - links[:, b]
                    vortex_sum += np.abs(np.rint(curl / (2 * np.pi))).mean(axis=spatial)
            winding = np.stack([np.rint(links[:, a].sum(axis=a + 1) / (2 * np.pi)).reshape((stop-first)*R, -1) for a in range(D)], axis=1)
            minimum = []
            for a in range(D):
                shape = [1] * (D + 1)
                shape[a + 1] = N
                phase = np.exp(-2j * np.pi * np.arange(N) / N).reshape(shape)
                minimum.append(np.abs((spin * phase).sum(axis=spatial))**2 / volume)
            values = dict(magnetization_real=magnetization.real, magnetization_imag=magnetization.imag,
                          magnetization_abs=np.abs(magnetization), magnetization_squared=m2,
                          magnetization_fourth=m2 * m2, mean_link_cosine=cosine,
                          action_density=-D * np.tile(ks,stop-first) * cosine,
                          vortex_density=vortex_sum / (D * (D - 1) / 2), winding=winding,
                          structure_zero=volume * m2, structure_min=np.stack(minimum, axis=1),
                          max_link_angle=np.abs(links).max(axis=tuple(range(1, D + 2))))
            for key, value in values.items():
                np.testing.assert_allclose(f[key][first:stop].reshape(value.shape), value, rtol=2e-12, atol=3e-13, err_msg=key)
    h=hashlib.sha256()
    with open(path,'rb') as stream:
        for b in iter(lambda:stream.read(1024*1024),b''):h.update(b)
    return dict(passed=True,rows=steps+1,slots=R,field_configurations=(steps+1)*R,bytes=Path(path).stat().st_size,sha256=h.hexdigest())
