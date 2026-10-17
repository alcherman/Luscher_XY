"""Replay and path acceptance checks for the isolated auxiliary transition."""
from pathlib import Path
import sys,h5py,numpy as np,pytest,json
P=Path(__file__).resolve().parents[1]/'experiments/2026-10-17-xy-3d-negative-kappa';sys.path.insert(0,str(P));import tempered_cluster as t

def test_path_reversal_and_palindrome():
 levels=t.ladder(-.25);values=np.random.default_rng(7).normal(size=len(levels)+1)
 assert t.path_log_ratio(-.25,levels,values)==pytest.approx(-t.path_log_ratio(-.25,levels,values[::-1]),abs=1e-14)
 with pytest.raises(ValueError):t.path_log_ratio(-.25,[0,1],[1,2,3])

def test_tempered_exact_replay_and_diagnostics(tmp_path):
 kw=dict(D=3,N=4,kappa=-.25,epsilon=2.1,warmup=2,seed=986,backend='numba',parallel=False,winding=True,initial_winding=[1,0,0],checkpoint_every=5,update='hybrid',clusters_per_step=4,buffer_rows=5,record_timing=True)
 a=tmp_path/'a.h5';b=tmp_path/'b.h5';t.run_chain(a,sweeps=80,**kw);t.run_chain(b,sweeps=35,**kw);t.run_chain(b,sweeps=80,resume=True,**kw)
 with h5py.File(a) as aa,h5py.File(b) as bb:
  for key in aa:
   if not key.endswith('_seconds'):np.testing.assert_array_equal(aa[key][:],bb[key][:],err_msg=key)
  assert aa['tempered_attempted'][:].sum()>0
  np.testing.assert_array_equal(aa['winding_accepted'][:],aa['native_twist_accepted'][:]+aa['tempered_accepted'][:])
  for i in np.flatnonzero(aa['tempered_attempted'][:]):
   actual=t.path_log_ratio(-.25,aa['tempered_kappa_path'][i],aa['tempered_C_path'][i]);assert actual==pytest.approx(aa['tempered_log_ratio'][i])
 from supervillain.xy import run_chain
 with pytest.raises((ValueError,RuntimeError)):run_chain(b,sweeps=90,resume=True,**kw)

def test_no_winding_cannot_create_stale_diagnostics(tmp_path):
 path=tmp_path/'no-winding.h5'
 with pytest.raises(ValueError,match='winding=True'):
  t.run_chain(path,D=3,N=4,kappa=-.25,epsilon=2.,sweeps=10,warmup=1,winding=False,update='hybrid')
 assert not path.exists()
