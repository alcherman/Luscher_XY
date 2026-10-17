"""Recover known finite-size roots and reject underidentified correction fits."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

PATH=Path(__file__).resolve().parents[1]/'experiments/2026-10-13-xy-infinite-volume/analyze_calibration.py'
SPEC=importlib.util.spec_from_file_location('xy_calibration_analysis',PATH)
analysis=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(analysis)


@pytest.mark.parametrize('D,critical,betas',[(3,.454,[.440,.455,.470]),(4,.3017,[.290,.302,.314])])
def test_finite_size_fit_recovers_root_and_bootstrap_interval(D,critical,betas):
    rows=[];samples=[];rng=np.random.default_rng(1801)
    for N in [6,8,12]:
        thermal=N**(1/.67169) if D==3 else N*N*np.log(N)**.1
        correction=N**(-.789) if D==3 else np.log(N)**(-.5)
        for beta in betas:
            value=(.35 if D==3 else .6)+.02*(beta-critical)*thermal+.1*correction
            rows.append(dict(D=D,N=N,beta=beta,summary=dict(values=[value]*3,errors=[.002]*3,flags=[])))
            samples.append(rng.normal(value,.002,size=(1000,3)))
    fit=analysis.fss_fit(rows,samples,D)
    assert fit['status']=='conditional_pilot_fit'
    assert fit['beta_c']==pytest.approx(critical,abs=1e-10)
    assert fit['bootstrap_95'][0]<critical<fit['bootstrap_95'][1]
    rejected=analysis.fss_fit(rows,samples,D,subleading=True)
    assert rejected['status']=='unidentifiable'
    rejected=analysis.fss_fit(rows,samples,D,minimum_N=8)
    assert rejected['status']=='unidentifiable'
