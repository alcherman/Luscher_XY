"""Independent energy-graph and exhaustive balance tests at negative coupling."""
from pathlib import Path
import sys,importlib.util
import numpy as np,pytest,h5py
P=Path(__file__).resolve().parents[1]/'experiments/2026-10-17-xy-3d-negative-kappa';sys.path.insert(0,str(P))
import signed_cluster as signed
from supervillain import XY,Lattice
from supervillain.generator.reference_implementation.xy_cluster import reference_cluster
spec=importlib.util.spec_from_file_location('positive_cluster_oracles',Path(__file__).with_name('test_xy_cluster.py'));oracle=importlib.util.module_from_spec(spec);spec.loader.exec_module(oracle)

@pytest.mark.parametrize('kappa',[-.1,-.5,-2.])
@pytest.mark.parametrize('epsilon',[1.7,2.4,np.pi])
def test_signed_energy_graph_oracle(kappa,epsilon):
 action=XY(Lattice(3,4),kappa,epsilon);update=signed.SignedCluster(action);rng=np.random.default_rng(913);phi=action.canonical(rng.uniform(0,.3,(4,4,4)))
 for _ in range(20):
  seed=int(rng.integers(64));alpha=rng.uniform(0,2*np.pi);tape=rng.uniform(size=(3,64));reference=np.asarray(phi).reshape(-1).copy()
  count,tests=reference_cluster(reference,update.plus,update.minus,seed,alpha,tape,kappa,epsilon)
  out=update.cluster(phi,seed,alpha,tape)
  np.testing.assert_array_equal(np.asarray(out).reshape(-1),reference);assert action.valid({'phi':out});assert update.last_cluster_size==count and update.last_bond_tests==tests;phi=out

@pytest.mark.parametrize('D,N,kappa,epsilon',[(1,2,-.1,np.pi),(1,3,-.5,1.8),(2,2,-.5,np.pi),(2,2,-2.,1.8)])
def test_exact_signed_transition_balance(D,N,kappa,epsilon,monkeypatch):
 monkeypatch.setattr(oracle,'ClusterUpdate',signed.SignedCluster)
 oracle.test_exact_embedded_state_transition_matrix_obeys_detailed_balance(D,N,kappa,epsilon)

def test_signed_checkpoint_replay_and_native_refusal(tmp_path):
 from supervillain.xy import run_chain as native
 kw=dict(D=3,N=4,kappa=-.25,epsilon=2.1,warmup=2,seed=667,backend='numba',parallel=False,winding=True,initial_winding=[0,0,0],checkpoint_every=5,update='hybrid',clusters_per_step=4,buffer_rows=5)
 a=tmp_path/'a.h5';b=tmp_path/'b.h5';signed.run_chain(a,sweeps=20,**kw);signed.run_chain(b,sweeps=10,**kw);signed.run_chain(b,sweeps=20,resume=True,**kw)
 with h5py.File(a) as aa,h5py.File(b) as bb:
  for key in aa:np.testing.assert_array_equal(aa[key][:],bb[key][:],err_msg=key)
 with pytest.raises((ValueError,RuntimeError)):native(b,sweeps=25,resume=True,**kw)

def test_signed_can_remove_initial_winding():
 action=XY(Lattice(2,8),-.25,np.pi/2);u=signed.SignedCluster(action);phi=action.canonical(np.broadcast_to(np.arange(8)[:,None]*2*np.pi/8,(8,8)))
 out=u.cluster(phi,0,0.,np.full((2,64),.99999));after=np.asarray(out)[0];w=oracle._principal(np.roll(after,-1,axis=0)-after).sum(0)/(2*np.pi)
 assert np.allclose(w,0);assert action.valid({'phi':out})
