"""Verify every saved configuration and its physical observables in the comparison."""
import hashlib
import json
from pathlib import Path
import h5py
import numpy as np
from run import OUT, ARMS, STEPS, atomic_json


def checksum(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(4*1024**2),b''):h.update(block)
    return h.hexdigest()


def wrap(a):
    reduced=a%(2*np.pi)
    reduced=np.where(reduced>=np.pi,reduced-2*np.pi,reduced)
    return np.where((a>=-np.pi)&(a<np.pi),a,reduced)


def main():
    rows=[]
    for arm in ARMS:
        path=OUT/'chains'/f"{arm['name']}.h5"
        with h5py.File(path,'r') as f:
            p=json.loads(f.attrs['parameters'])
            assert p['D']==2 and p['N']==48 and p['kappa']==.8 and p['epsilon']==np.pi
            assert p['update']==arm['update'] and p['seed']==arm['seed']
            assert p['clusters_per_step']==arm['clusters_per_step'] and p['record_timing']
            assert not p['parallel'] and p['winding']
            assert f.attrs['committed_rows']==STEPS+1
            assert all(len(ds)==STEPS+1 for ds in f.values())
            assert np.array_equal(f['sweep'][:],np.arange(STEPS+1))
            assert np.all(f['phi'][0]==0)
            if arm['update']!='metropolis':
                sizes=f['cluster_sizes'][:]
                assert np.all((sizes[1:]>=1)&(sizes[1:]<=2304))
                assert np.array_equal(f['clusters_total'][:],arm['clusters_per_step']*np.arange(STEPS+1))
                assert np.array_equal(f['cluster_sites_total'][:],np.cumsum(sizes.sum(1)))
                assert np.array_equal(f['cluster_equivalent_sweeps'][:],np.cumsum(sizes.sum(1))/2304)
                assert np.array_equal(f['cluster_bond_tests_total'][:],np.cumsum(f['cluster_bond_tests'][:].sum(1)))
            for key in f:
                if key.endswith('_seconds'):
                    assert np.isfinite(f[key][:]).all() and np.all(f[key][:]>=0)
            for start in range(0,STEPS+1,128):
                end=min(start+128,STEPS+1);a=f['phi'][start:end,0]
                assert np.isfinite(a).all() and np.all((a>=0)&(a<2*np.pi))
                dx,dy=wrap(np.roll(a,-1,axis=1)-a),wrap(np.roll(a,-1,axis=2)-a)
                assert np.all(np.abs(dx)<np.pi) and np.all(np.abs(dy)<np.pi)
                spin=np.exp(1j*a);m=spin.mean((1,2));m2=np.abs(m)**2
                cosine=(np.cos(dx).mean((1,2))+np.cos(dy).mean((1,2)))/2
                vortex=np.abs(np.rint((dx+np.roll(dy,-1,axis=1)-np.roll(dx,-1,axis=2)-dy)/(2*np.pi))).mean((1,2))
                winding=np.stack((np.rint(dx.sum(1)/(2*np.pi)),np.rint(dy.sum(2)/(2*np.pi))),axis=1)
                phase=np.exp(-2j*np.pi*np.arange(48)/48)
                minimum=np.stack((np.abs((spin*phase[None,:,None]).sum((1,2)))**2/2304,
                                  np.abs((spin*phase[None,None,:]).sum((1,2)))**2/2304),axis=1)
                expected=dict(magnetization_real=m.real,magnetization_imag=m.imag,
                    magnetization_abs=np.abs(m),magnetization_squared=m2,magnetization_fourth=m2*m2,
                    mean_link_cosine=cosine,action_density=-1.6*cosine,vortex_density=vortex,winding=winding,
                    structure_zero=2304*m2,structure_min=minimum,
                    max_link_angle=np.maximum(np.abs(dx).max((1,2)),np.abs(dy).max((1,2))))
                for key,value in expected.items():
                    np.testing.assert_allclose(f[key][start:end],value,rtol=2e-12,atol=3e-13,err_msg=key)
            expected_states={'local','winding','cluster'} if arm['update']=='hybrid' else {'local','winding'}
            for state in f['rng_state'][:]:assert set(json.loads(state))==expected_states
        rows.append(dict(arm=arm['name'],rows=STEPS+1,bytes=path.stat().st_size,sha256=checksum(path)))
        print(f"{arm['name']}: audited {STEPS+1} rows",flush=True)
    atomic_json(OUT/'data-audit.json',dict(passed=True,arms=rows,total_rows=sum(r['rows'] for r in rows),
                                         audit_source_sha256=checksum(Path(__file__))))


if __name__=='__main__':main()
