"""Independent quadrature and exact limiting cases for the supplied prediction."""
import sys
from pathlib import Path
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import i1e,i0e
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments/2026-10-13-xy-infinite-volume'))
from tree_level_prediction import v1,endpoint,beta_critical

@pytest.mark.parametrize('beta,epsilon',[(0,.9*np.pi),(.3,.8*np.pi),(.8,np.pi)])
def test_independent_integral(beta,epsilon):
    num=quad(lambda x:np.exp(beta*np.cos(x))*np.cos(x),-epsilon,epsilon,epsabs=1e-12)[0]
    den=quad(lambda x:np.exp(beta*np.cos(x)),-epsilon,epsilon,epsabs=1e-12)[0]
    assert v1(beta,epsilon)==pytest.approx(num/den,abs=2e-13)

@pytest.mark.parametrize('D',[2,3,4])
def test_endpoints_and_prediction_residual(D):
    q=endpoint(D);assert np.sin(np.pi*q)/(np.pi*q)==pytest.approx(1/(2*D-1),abs=2e-13)
    assert beta_critical(D,q)==0;assert np.isnan(beta_critical(D,q-.02))
    b=beta_critical(D,1);assert i1e(b)/i0e(b)==pytest.approx(1/(2*D-1),abs=2e-13)
    for x in np.linspace(q,1,9):assert (2*D-1)*v1(beta_critical(D,x),np.pi*x)==pytest.approx(1,abs=1e-12)
