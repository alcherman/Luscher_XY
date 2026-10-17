from pathlib import Path
import sys,importlib.util
import numpy as np
import pytest
P=Path(__file__).resolve().parents[1]/'experiments/2026-10-17-xy-3d-negative-kappa'
sys.path.insert(0,str(P))
spec=importlib.util.spec_from_file_location('negative_scout',P/'scout.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
a_spec=importlib.util.spec_from_file_location('negative_audit',P/'audit.py');a=importlib.util.module_from_spec(a_spec);a_spec.loader.exec_module(a)

def test_finite_design():
 ps=m.make_points();assert len(ps)==48;assert len({p['id'] for p in ps})==48;assert len({p['seed'] for p in ps})==48
 assert {p['kappa'] for p in ps}=={0,-.1,-.25,-.5};assert {p['N'] for p in ps}=={8,12}
 assert all(p['update']=='metropolis' and p['steps']==4000 for p in ps)

def test_negative_full_audit_and_resume(tmp_path):
 p=dict(D=3,N=4,beta=-.25,kappa=-.25,epsilon=2.4,seed=135,initial_winding=[0,0,0],steps=20)
 kwargs=dict(D=3,N=4,kappa=-.25,epsilon=2.4,warmup=2,seed=135,backend='numba',parallel=False,interval_phi=np.pi,winding=True,initial_winding=[0,0,0],checkpoint_every=5,update='metropolis',record_timing=True,buffer_rows=100)
 f=tmp_path/'a.h5';m.run_chain(f,sweeps=10,**kwargs);m.run_chain(f,sweeps=20,resume=True,**kwargs)
 assert a.audit(f,p)['passed']
 import h5py
 with h5py.File(f,'r+') as h:h['magnetization_squared'][7]+=.1
 with pytest.raises(AssertionError):a.audit(f,p)

def test_root_knots_and_ambiguities():
 spec=importlib.util.spec_from_file_location('negative_analysis',P/'analyse.py');v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
 assert v.root([0,1,2],[.5,.25,0])==1
 assert v.root([0,1,2],[.5,0,.5]) is None
 assert v.root([0,1,2],[.25,.25,0]) is None
 assert v.root([0,1,2],[.5,np.nan,0]) is None
 assert v.root([0,1,2],[.5,.4,.3]) is None
