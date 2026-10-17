"""Independent algebra/temporal-weight oracles for cross-fitted spin CV."""
from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).parents[1]/'experiments/2026-10-15-xy-3d-spin-dimension'))
import spin_control_variate as cv

def row(gq,gh,mcp,mdense):
    gq,gh,mcp,mdense=np.broadcast_arrays(gq,gh,mcp,mdense)
    return np.stack([gq,gh,mdense,mdense*100,mdense*20,np.zeros_like(mcp),mcp,gq*mcp,gh*mcp,mcp*mcp],axis=-1)

def test_coefficient_raw_regression_oracle():
    m=np.arange(1.,8.);a=row(2*m+3,4*m+1,m,m*.9)
    np.testing.assert_allclose(cv.coefficients(a.mean(0)),[2,4],atol=1e-13)

def test_opposite_start_coefficients_and_weighted_pool():
    m=np.arange(1.,8.);x=row(2*m+3,4*m+1,m,m-1).mean(0);y=row(3*m+7,5*m+4,m,m-2).mean(0)
    O,C,c=cv.crossfit(x,y,7,14)
    np.testing.assert_allclose(c,[[3,5],[2,4]],atol=1e-13)
    np.testing.assert_allclose(C[:2],(7*(x[:2]-[3,5])+14*(y[:2]-2*np.array([2,4])))/21)
    np.testing.assert_allclose(C[2:],O[2:]);np.testing.assert_allclose(O,(7*x[:6]+14*y[:6])/21)

def test_zero_mean_difference_leaves_G_unchanged():
    m=np.arange(1.,11.);x=row(m,2*m,m,m).mean(0);O,C,c=cv.crossfit(x,x,10,10)
    np.testing.assert_array_equal(O,C)

def test_training_variance_zero_explicit_nan():
    x=row(np.arange(4.),np.arange(4.)*2,np.ones(4),np.ones(4)).mean(0)
    assert np.isnan(cv.coefficients(x)).all()

def test_augmented_time_alignment_and_products():
    t=np.arange(0,401);series=np.ones((5,6));series[:,0]=np.arange(5);series[:,1]=np.arange(5)*2;series[:,2]=np.arange(5)*3
    d=dict(sweep=np.arange(0,401,100),series=series,dense_sweep=t,dense_m2=t*.1,dense_S0=t,dense_F=t*.2)
    times,a=cv.augmented_series(d,100)
    np.testing.assert_array_equal(times,[200,300,400]);np.testing.assert_array_equal(a[:,6],series[2:,2]);np.testing.assert_array_equal(a[:,7],series[2:,0]*series[2:,2]);assert a[:,2].mean()==pytest.approx(t[101:].mean()*.1)

def test_gain_keeps_negative_draws_before_any_log_mask():
    O=np.ones((5,6));C=O.copy();O[:,:2]=np.arange(5)[:,None]+1;C[:,:2]=np.array([-3,-1,1,3,5])[:,None]
    out=cv.raw_gain(O,C)
    np.testing.assert_allclose(out['SD_ratio'],[2,2]);assert out['finite_pair_fraction']==1
    np.testing.assert_array_equal(out['negative_mean_draw_counts']['CV'],[2,2])

def test_bootstrap_refits_training_coefficients_and_preserves_scalar_draws():
    rng=np.random.default_rng(9);m=np.linspace(.1,1,80);x=row(2*m+rng.normal(0,.2,80),3*m+rng.normal(0,.3,80),m,m*.95);y=row(2*m+rng.normal(0,.2,80),3*m+rng.normal(0,.3,80),m,m*1.05)
    a,_=cv.base.bootstrap_blocks(x,10,rng);b,_=cv.base.bootstrap_blocks(y,10,rng)
    O,C,c=cv.crossfit(a[0]/a[1][:,None],b[0]/b[1][:,None],a[1],b[1])
    assert (c.std(0)>.01).all();np.testing.assert_array_equal(O[:,2:],C[:,2:])
    # Training for first estimator depends only on second bank.
    np.testing.assert_allclose(c[:,0],cv.coefficients(b[0]/b[1][:,None]))

def test_insufficient_group_is_identifiable(monkeypatch):
    m=np.linspace(.1,1,80);r=row(m,2*m,m,.95*m)
    metas=[dict(id='start0'),dict(id='start1')]
    monkeypatch.setattr(cv,'read_group',lambda g,cut:([r,r],metas,[]))
    monkeypatch.setattr(cv.base,'bootstrap_blocks',lambda data,w,rng:(None,dict(blocks=1,rows=len(data),width_checkpoints=w)))
    g=dict(L=24,q=.8,flags=[],source_metadata=[],diagnostics=[dict(id=m['id'],blocking=dict(width_checkpoints=10)) for m in metas])
    variants,Z,O,C,gain,_=cv.summarize(g,1000,2)
    assert gain['L']==24 and gain['q']==.8 and gain['ids']==['start0','start1']
    assert gain['status']=='insufficient_blocks';assert Z is None
    assert all(v['status']=='insufficient_blocks' for v in variants)
