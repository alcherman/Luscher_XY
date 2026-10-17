from pathlib import Path
import sys
import numpy as np
import h5py,pytest
sys.path.insert(0,str(Path(__file__).parents[1]/'experiments/2026-10-15-xy-3d-spin-dimension'))
import measure_spin_correlations as c

def test_uniform_spin():
    m=c.spin_measure(np.zeros((1,8,8,8)));np.testing.assert_allclose(m['axis_real'],1);assert m['m2']==1

def test_unit_helix_and_half_sine():
    L=8;phi=np.broadcast_to(2*np.pi*np.arange(L)[:,None,None]/L,(L,L,L));m=c.spin_measure(phi)
    np.testing.assert_allclose(m['axis_real'],[[0,1,1],[-1,1,1]],atol=1e-14)
    np.testing.assert_allclose(m['axis_real'].mean(1),[2/3,1/3]);np.testing.assert_allclose(m['axis_imag'][1],0,atol=1e-14)
    assert abs(m['axis_imag'][0,0])==pytest.approx(1)

def test_rotation_and_translation_invariance():
    rng=np.random.default_rng(8);phi=rng.uniform(0,2*np.pi,(8,8,8));a=c.spin_measure(phi);b=c.spin_measure(np.roll(phi+.739,2,axis=1))
    for k in ['axis_real','axis_imag','m2','S0','F']:np.testing.assert_allclose(a[k],b[k],atol=1e-13)

def test_against_explicit_pair_sum():
    L=4;phi=np.random.default_rng(9).uniform(0,2*np.pi,(L,L,L));m=c.spin_measure(phi)
    for k,r in enumerate([L//4,L//2]):
        for axis in range(3):
            vals=[]
            for x in np.ndindex(phi.shape):
                y=list(x);y[axis]=(y[axis]+r)%L;vals.append(np.exp(1j*(phi[tuple(y)]-phi[x])))
            np.testing.assert_allclose(m['axis_real'][k,axis]+1j*m['axis_imag'][k,axis],np.mean(vals),atol=1e-14)

def test_structure_min_matches_full_fft():
    L=8;phi=np.random.default_rng(5).normal(size=(L,L,L));m=c.spin_measure(phi);f=np.fft.fftn(np.exp(1j*phi));expected=[abs(f[1,0,0])**2/L**3,abs(f[0,1,0])**2/L**3,abs(f[0,0,1])**2/L**3]
    np.testing.assert_allclose(m['F'],expected,atol=1e-12)

def test_missing_field_rejected():
    a=np.zeros((4,4,4));a[0]=np.nan
    with pytest.raises(ValueError,match='Missing'):c.spin_measure(a)

def test_checkpoint_missing_or_misnumbered_rejected(tmp_path):
    with h5py.File(tmp_path/'a.h5','w') as h:
        h.create_dataset('phi',shape=(201,1,4,4,4),dtype='f8');h.create_dataset('sweep',data=np.arange(201));h.attrs['committed_rows']=201;h['phi'].attrs['retained_indices']=[0,200]
        with pytest.raises(ValueError,match='absent'):c.retained_checkpoint_rows(h,200)
        h['phi'].attrs['retained_indices']=[0,100,200];h['sweep'][100]=99
        with pytest.raises(ValueError,match='cycle'):c.retained_checkpoint_rows(h,200)

def test_prescribed_checkpoint_count(tmp_path):
    with h5py.File(tmp_path/'a.h5','w') as h:
        h.create_dataset('phi',shape=(12001,1,4,4,4),dtype='f8');h.create_dataset('sweep',data=np.arange(12001));h.attrs['committed_rows']=12001;h['phi'].attrs['retained_indices']=np.arange(0,12001,100)
        rows=c.retained_checkpoint_rows(h,12000);assert np.sum(rows>1000)==110

import analyze_spin_correlations as a

def test_dense_scalar_alignment_and_central_mean():
    t=np.arange(401);series=np.ones((5,6));series[:,0]=np.arange(5)+.2
    d=dict(sweep=np.arange(0,401,100),series=series,dense_sweep=t,dense_m2=t**2,dense_S0=t*3.,dense_F=t+7.)
    times,b=a.aligned_series(d,100,'dense_bins')
    np.testing.assert_array_equal(times,[200,300,400]);np.testing.assert_array_equal(b[:,0],series[2:,0])
    np.testing.assert_allclose(b[:,2:5].mean(0),[np.mean(t[101:]**2),np.mean(3*t[101:]),np.mean(t[101:]+7)])
    d['dense_sweep']=np.delete(t,200)
    with pytest.raises(ValueError,match='align'):a.aligned_series(d,100,'dense_bins')

def test_physical_block_units():
    cp,cycles=a.block_width_cycles({'G':.5},{'m2':12.},factor=2)
    assert cp==10 and cycles==1000
    assert a.block_width_cycles({'G':.5},{'m2':12.},factor=4)==(20,2000)

def test_joint_block_covariance_weighted_pooling():
    x=np.arange(11.)+1;raw=np.column_stack([x,2*x,3*x,4*x,5*x,6*x])
    (s,n),info=a.bootstrap_blocks(raw,3,np.random.default_rng(7),draws=100)
    means=s/n[:,None]
    np.testing.assert_allclose(means[:,1],2*means[:,0]);assert info['actual_block_lengths_checkpoints']==[3,4,4]
    assert set(n).issubset(set(range(9,13)))

def test_logs_of_ensemble_means_signed_configurations():
    raw=np.array([[-1,-2,.1,2.,1.,0],[3,4,.3,4.,1.,0]])
    z=a.derived(raw.mean(0),24)
    np.testing.assert_allclose(z[1:3],[0.,0.]);assert np.isfinite(z).all()

import spin_correlation_quotients as sq

@pytest.mark.parametrize('kind',['linear','pchip'])
def test_common_R_spin_quotient_oracle(kind):
    q=np.array([.7,.8,.9]);small=np.zeros((3,4));large=np.zeros((3,4))
    small[:,0]=1+2*(q-.8);large[:,0]=1+3*(q-.8)
    for j in range(1,4):
        small[:,j]=j+q;large[:,j]=j+q-2*.519*np.log(2)
    out,n=sq.quotient(q,small,large,kind)
    assert n==1;np.testing.assert_allclose(out,[.8,1,.519,.519,.519],atol=1e-13)

def test_no_clipping_and_ambiguous_crossings():
    q=np.array([.7,.8,.9]);small=np.ones((3,4));large=np.ones((3,4));large[:,0]=2
    assert sq.quotient(q,small,large)[0] is None
    large[:,0]=[2,0,2]
    out,n=sq.quotient(q,small,large);assert out is None and n==2

@pytest.mark.parametrize('kind',['linear','pchip'])
def test_crossing_exact_support_endpoint_once(kind):
    q=np.array([.7,.8,.9]);s=np.zeros((3,4));l=np.zeros((3,4));s[:,0]=q;l[:,0]=2*q-.7
    v,n=sq.quotient(q,s,l,kind);assert n==1;assert v[0]==pytest.approx(.7)

def test_nonmultiple_cut_rejected():
    d=dict(sweep=np.arange(0,401,100),series=np.ones((5,6)),dense_sweep=np.arange(401),dense_m2=np.ones(401),dense_S0=np.ones(401),dense_F=np.ones(401))
    with pytest.raises(ValueError,match='align'):a.aligned_series(d,101,'dense_bins')

def test_fixed_q_GLS_shared_slope_analytic_oracle():
    rng=np.random.default_rng(19);N=np.array([24,32,48,64]);delta=.519
    y=np.column_stack([-.5-2*delta*np.log(N/32),-.9-2*delta*np.log(N/32)])
    C=np.array([[.0004,.00024],[.00024,.0009]])
    B=np.zeros((5000,4,4));B[:,:,1:3]=y+rng.multivariate_normal([0,0],C,size=(5000,4))
    groups=[dict(L=int(n),q=.8,derived=[.59,*row,-3.]) for n,row in zip(N,y)]
    result=a.fixed_q_slope(groups,B)
    assert result['Delta']==pytest.approx(delta,abs=1e-12)
    X=np.zeros((8,3));X[::2,0]=1;X[1::2,1]=1;X[:,2]=np.repeat(-2*np.log(N/32),2)
    precision=np.zeros((8,8))
    for i in range(4):precision[2*i:2*i+2,2*i:2*i+2]=np.linalg.inv(np.cov(B[:,i,1:3].T))
    T=np.linalg.solve(X.T@precision@X,X.T@precision)
    expected=B[:,:,1:3].reshape(5000,8)@T.T
    assert result['SD']==pytest.approx(expected[:,-1].std(ddof=1),rel=1e-12)
    np.testing.assert_allclose(result['log_amplitude_at_L32'],[-.5,-.9],atol=1e-12)

def test_unresolved_groups_report_without_intervals(tmp_path,monkeypatch):
    import json,argparse
    import report_spin_correlations as report
    source=tmp_path/'source';source.mkdir();groups=[]
    for L in [24,32,48,64]:
        for q in [.8005,.80075,.801]:
            groups.append(dict(L=L,q=q,status='insufficient_blocks',flags=['insufficient_blocks'],samples_total=220,scalar_cycles_total=22000,raw_mean=[.1,.08,.2,2.,1.,0.],derived=[1+(q-.80075)*L,np.log(.1),np.log(.08),np.log(.2)]))
    (source/'analysis.json').write_text(json.dumps(dict(groups=groups,settings={},fixed_q_slope=dict(status='incomplete_sizes'))))
    monkeypatch.setattr(report,'save',lambda fig,path:report.plt.close(fig))
    out=tmp_path/'report';report.main(argparse.Namespace(summary=str(source),output=str(out)))
    result=json.loads((out/'analysis.json').read_text())
    assert all(p['status']=='unresolved_crossing_tails' and p['valid_fraction']==0 for p in result['pairs'])
    assert all(p['MC95'] is None and p['missing_draw_groups'] for p in result['pairs'])
