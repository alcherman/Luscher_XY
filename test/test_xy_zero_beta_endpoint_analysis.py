"""Endpoint epsilon scans have the opposite ordering direction to beta scans."""
import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments/2026-10-13-xy-infinite-volume'))
from analyze_zero_beta_endpoints import epsilon_crossing

@pytest.mark.parametrize('method',['linear','pchip'])
def test_epsilon_scan_orientation_and_measured_bracket(method):
    x=np.array([.80,.806,.812]); decreasing=np.array([.45,.27,.09])
    root=epsilon_crossing(x,decreasing,.25,False,method=method)
    assert root==pytest.approx(.806+(.006*.02/.18))
    assert epsilon_crossing(x,-decreasing,-.25,True,method=method)==pytest.approx(root)
    assert np.isnan(epsilon_crossing(x,decreasing,.5,False,method=method))

def test_bootstrap_bracket_loss_is_retained():
    x=np.array([.558,.565,.572])
    draws=[np.array([.8,.7,.6]),np.array([.7,.65,.6])]
    roots=np.array([epsilon_crossing(x,v,.75,False) for v in draws])
    assert np.isfinite(roots[0]);assert np.isnan(roots[1])

def test_local_endpoint_window_preserves_neighbor_when_near_anchor_is_weak():
    from analyze_zero_beta_endpoints import local_indices
    x=np.array([.553,.558,.565]);y=np.array([.762,.75065,.7433])
    # Last point is less than3SD belowtarget, so retain it rather than stop
    # at the nearly-on-target central bracket's right edge.
    B=np.array([y+d for d in [-.004,.004,-.003,.003]])
    np.testing.assert_array_equal(local_indices(x,y,B,.7506912,-1),[0,1,2])
    # Actual sign reversal, rather than reversing epsilon coordinates.
    np.testing.assert_array_equal(local_indices(x,-y,-B,-.7506912,1),[0,1,2])
