from pathlib import Path
import sys
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).parents[1]/'experiments/2026-10-15-xy-3d-spin-dimension'))
import joint_fss_existing as j

def data(model='leading',Q=.02):
    L=np.repeat([12,16,24,32,48,64],4);q=np.tile([.800625,.800750,.800875,.801000],6)
    X,_,names=j.design(L,q,Q,model)
    known=dict(Rstar=.59,R_x=-.12,R_x2=.02,R_correction=.01,log_amplitude=-3.4,Delta=.519,m_x=-.5,m_x2=.05,m_correction=.02,R_mixed=.007,m_mixed=-.01,m_analytic_amplitude=.03,R_second_correction=.003,m_second_correction=-.004,R_x3=.002,m_x3=.001)
    coef=np.array([known[n] for n in names]);y=np.einsum('nip,p->ni',X,coef);W=np.tile(np.eye(2),(len(L),1,1))
    return L,q,coef,y,W

@pytest.mark.parametrize('model',['leading','mixed','analytic','second','cubic'])
def test_noiseless_exact_model_recovery(model):
    L,q,c,y,W=data(model);r=j.solve_profile(L,q,y,W,model=model,grid=81,full=True)
    assert r['qc']==pytest.approx(j.QREF+.02/j.QSCALE,abs=1e-9)
    assert r['Delta']==pytest.approx(.519,abs=1e-7)
    assert r['chi2']<1e-15
    assert r['jacobian_rank']==len(c)+1

def test_block_whitening_matches_explicit_GLS():
    L,q,c,y,_=data();C=np.array([[.03,.004],[.004,.007]]);W=np.tile(np.linalg.inv(np.linalg.cholesky(C)),(len(L),1,1))
    X=j.design(L,q,.02)[0];rng=np.random.default_rng(42);yn=y+rng.normal(size=y.shape)*.001
    result=j.linear_solution(j.whiten(W,X),j.whiten(W,yn))[0]
    A=X.reshape(-1,len(c));precision=np.kron(np.eye(len(L)),np.linalg.inv(C));oracle=np.linalg.solve(A.T@precision@A,A.T@precision@yn.ravel())
    np.testing.assert_allclose(result,oracle,atol=2e-9)

def test_singular_sizes_not_identifiable():
    L,q,c,y,W=data();L[:]=32;X=j.design(L,q,.02)[0];coef,score,rank,condition=j.linear_solution(X,y)
    assert rank<len(c)

def test_profile_boundary_flagged():
    L,q,c,y,W=data(Q=2.);r=j.solve_profile(L,q,y,W,grid=81)
    assert r['boundary']

def test_raw_background_is_additive_and_signed():
    L,q,c,y,W=data();p=np.r_[.02,c,-.01];pred,J=j.raw_prediction(L,q,p,background=True)
    np.testing.assert_allclose(pred[:,1],np.exp(y[:,1])-.01*.033*(L/32.)**-3)
    p[-1]=0;pred,_=j.raw_prediction(L,q,p,background=True);np.testing.assert_allclose(pred[:,1],np.exp(y[:,1]))

@pytest.mark.parametrize('bg',[False,True])
def test_raw_analytic_jacobian(bg):
    L,q,c,y,W=data('mixed');p=np.r_[.02,c,[-.002] if bg else []];pred,J=j.raw_prediction(L,q,p,model='mixed',background=bg)
    for k in range(len(p)):
        e=np.zeros_like(p);e[k]=1e-6
        fd=(j.raw_prediction(L,q,p+e,model='mixed',background=bg)[0]-j.raw_prediction(L,q,p-e,model='mixed',background=bg)[0])/(2e-6)
        np.testing.assert_allclose(J[:,:,k],fd,rtol=2e-6,atol=2e-9)

def test_heldout_residual_uses_full_prediction_covariance():
    rng=np.random.default_rng(7);draws=rng.normal(size=(2000,3,2));draws[:,:,0]+=rng.normal(size=(2000,1))*2
    obs=np.array([[.1,.2],[-.1,.3],[.2,-.2]]);pred=np.zeros_like(obs)
    out=j.predictive_summary(obs,pred,draws);C=np.cov(draws.reshape(2000,-1),rowvar=False)
    assert out['chi2']==pytest.approx(obs.ravel()@np.linalg.solve(C,obs.ravel()))
    assert abs(C[0,2])>1

def test_holdout_never_includes_training_size():
    rows=[dict(L=L,q=.80075) for L in [12,16,24,32,48,64]];c=dict(window='inner',minimum=12,maximum=48)
    assert j.select(rows,c).tolist()==[0,1,2,3,4]
    assert j.select(rows,c,True).tolist()==[5]

def test_every_bootstrap_draw_checks_joint_rank(monkeypatch):
    L,q,c,y,W=data();rng=np.random.default_rng(50);B=y[None,:,:]+rng.normal(size=(12,len(L),2))*.001
    rows=[dict(L=int(n),q=float(qq),mean=yy,ids=[str(i)],flags=[]) for i,(n,qq,yy) in enumerate(zip(L,q,y))]
    original=j.solve_profile;calls=0
    def bad_draw(*args,**kwargs):
        nonlocal calls
        calls+=1;r=original(*args,**kwargs)
        if calls>1:r['jacobian_rank']=r['jacobian_rank']-1
        return r
    monkeypatch.setattr(j,'solve_profile',bad_draw)
    config=dict(name='test',window='inner',minimum=12,maximum=64,model='leading',nu=j.NU,omega=j.OMEGA,holdout=None)
    result,_=j.fit_one(config,rows,B,draws=12)
    assert result['successful_draws']==0
    assert result['bootstrap_failure_reasons']['joint_rank_failure']==12
    assert 'bootstrap_failure_unresolved' in result['flags']

def test_existing_output_folder_is_never_relabelled(monkeypatch,tmp_path):
    from types import SimpleNamespace
    folder=tmp_path/'prior';folder.mkdir();(folder/'design.json').write_text('original provenance')
    monkeypatch.setattr(j,'load_bank',lambda p:([],np.zeros((2,0,2)),{}))
    with pytest.raises(FileExistsError):j.main(SimpleNamespace(output=str(folder),draws=2000,central_only=False))
    assert (folder/'design.json').read_text()=='original provenance'
