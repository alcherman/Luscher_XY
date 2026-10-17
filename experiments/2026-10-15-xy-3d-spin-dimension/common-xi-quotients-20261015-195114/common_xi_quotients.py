"""Same-epsilon xi/L crossings and joint magnetization quotients; existing MC only."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path
import numpy as np
import h5py
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.stats import chi2
import measure as m

VERSION = 1


def ratio(v, N):
    # Columns m2,m4,|m|,S(0),direction-averaged S(k_min).
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.sqrt(v[..., 3] / v[..., 4] - 1) / (2*N*np.sin(np.pi/N))


def summary(points, cut, factor):
    ids = sorted(p['id'] for p in points)
    key = hashlib.sha256(json.dumps([ids,cut,factor,VERSION,'common_xi']).encode()).hexdigest()
    cache = m.OUT/'xi-cache'; cache.mkdir(exist_ok=True)
    meta, data = cache/(key+'.json'), cache/(key+'.npz')
    if meta.exists():
        return json.loads(meta.read_text()), np.load(data)['samples'] if data.exists() else None
    N = points[0]['N']; rng = np.random.default_rng(int(key[:16],16))
    total = np.zeros((m.DRAWS,5)); lengths = np.zeros(m.DRAWS)
    central_sum = np.zeros(5); count = 0; diagnostics = []; flags = []; individual = []
    for i in ids:
        with h5py.File(m.CAMPAIGN/'chains'/f'{i}.h5') as h:
            a = np.column_stack([h[k][cut:] for k in ['magnetization_squared','magnetization_fourth','magnetization_abs','structure_zero']]+[h['structure_min'][cut:].mean(1)])
        if len(a)<4:
            return dict(status='insufficient_rows',ids=ids),None
        mu = a.mean(0); centered = a-mu
        monitor = {k:a[:,j] for j,k in enumerate(['m2','m4','mabs','S0','F'])}
        monitor['R_influence'] = centered[:,3]/mu[4]-mu[3]*centered[:,4]/mu[4]**2
        monitor['binder_influence'] = mu[1]*centered[:,0]/mu[0]**3-centered[:,1]/(2*mu[0]**2)
        taus = {k:m.helper.tau_int(v) for k,v in monitor.items()}
        width = factor*max(20,int(np.ceil(10*max([1.]+[float(t) for t in taus.values() if t is not None and np.isfinite(t)]))))
        nb = len(a)//width
        diag = dict(id=i,rows=len(a),width=width,blocks=nb,taus=taus)
        diagnostics.append(diag)
        if nb<2:
            return dict(status='insufficient_blocks',ids=ids,diagnostics=diagnostics),None
        if nb<50: flags.append('fewer_than_50_blocks_per_start')
        edges = np.linspace(0,len(a),nb+1,dtype=int)
        sums = np.add.reduceat(a,edges[:-1],axis=0); sizes = np.diff(edges)
        ix = rng.integers(nb,size=(m.DRAWS,nb))
        bs = sums[ix].sum(1); bn = sizes[ix].sum(1); draw = bs/bn[:,None]
        total += bs; lengths += bn; central_sum += a.sum(0); count += len(a)
        individual.append(dict(id=i,mean=mu.tolist(),m2_SD=float(draw[:,0].std(ddof=1)),R=float(ratio(mu,N)),R_SD=float(np.nanstd(ratio(draw,N),ddof=1))))
        halves = np.array_split(a,2); half_stats = []
        for half in halves:
            nbh = len(half)//width
            if nbh<2:
                flags.append('insufficient_half_blocks'); continue
            bm = np.array([v.mean(0) for v in np.array_split(half,nbh)])
            half_stats.append((np.array([half[:,0].mean(),ratio(half.mean(0),N)]),np.array([bm[:,0].std(ddof=1),np.nanstd(ratio(bm,N),ddof=1)])/np.sqrt(nbh)))
        if len(half_stats)==2:
            z = (half_stats[1][0]-half_stats[0][0])/np.maximum(np.hypot(half_stats[0][1],half_stats[1][1]),1e-15)
            diag['half_z_m2_R'] = z.tolist()
            if np.max(abs(z))>4: flags.append('half_shift_m2_or_R')
    B = total/lengths[:,None]; mu = central_sum/count
    for name in ['m2','R']:
        vals = np.array([v['mean'][0] if name=='m2' else v['R'] for v in individual]); sd = np.array([v[name+'_SD'] for v in individual]); good = np.isfinite(vals)&np.isfinite(sd)&(sd>0)
        if good.sum()>1:
            w = 1/sd[good]**2; center = np.average(vals[good],weights=w)
            if chi2.sf(np.sum(w*(vals[good]-center)**2),good.sum()-1)<.001: flags.append('start_disagreement_'+name)
    s = dict(status='conditional_joint_block_summary',ids=ids,mean=mu.tolist(),R=float(ratio(mu,N)),R_SD=float(np.nanstd(ratio(B,N),ddof=1)),flags=sorted(set(flags)),diagnostics=diagnostics,individual=individual)
    # A finite analysis process publishes data first; no parallel writer uses this cache.
    np.savez_compressed(data,samples=B); m.run.atomic_json(meta,s)
    return s,B


def interpolator(q, values, method):
    if method=='pchip': return PchipInterpolator(q,values,axis=0,extrapolate=False)
    return lambda x: np.array([np.interp(x,q,values[:,j]) for j in range(values.shape[1])])


def crossing(q1,v1,N1,q2,v2,N2,lo,hi,method='linear'):
    """Reject absent/multiple roots; interpolate moments then form ratios."""
    f1,f2 = interpolator(q1,v1,method),interpolator(q2,v2,method)
    knots = np.unique(np.r_[lo,hi,q1[(q1>lo)&(q1<hi)],q2[(q2>lo)&(q2<hi)]])
    # Subdivide to detect multiple PCHIP crossings, including within a knot cell.
    if method=='pchip': knots = np.unique(np.concatenate([np.linspace(a,b,5) for a,b in zip(knots[:-1],knots[1:])]))
    def difference(x): return float(ratio(f1(x),N1)-ratio(f2(x),N2))
    d = np.array([difference(x) for x in knots]); roots = []
    if not np.isfinite(d).all(): return np.full(5,np.nan)
    for k in range(len(knots)-1):
        if d[k]==0: roots.append(knots[k])
        elif d[k]*d[k+1]<0: roots.append(brentq(difference,knots[k],knots[k+1],xtol=1e-13))
    if d[-1]==0: roots.append(knots[-1])
    roots = sorted(set(round(float(x),13) for x in roots))
    if len(roots)!=1: return np.full(5,np.nan)
    root = roots[0]; a,b = f1(root),f2(root)
    if min(a[0],b[0])<=0: return np.full(5,np.nan)
    Delta = -np.log(b[0]/a[0])/(2*np.log(N2/N1))
    return np.array([root,a[0],b[0],ratio(a,N1),Delta])


def paired_record(q1,v1,B1,N1,q2,v2,B2,N2,method='linear',window='local'):
    lo,hi = max(q1[0],q2[0]),min(q1[-1],q2[-1])
    if lo>=hi: return dict(status='no_common_epsilon_support'),None
    center = crossing(q1,v1,N1,q2,v2,N2,lo,hi,method)
    if not np.isfinite(center).all(): return dict(status='central_crossing_unresolved',common_support=[lo,hi]),None
    union = np.unique(np.r_[lo,hi,q1[(q1>lo)&(q1<hi)],q2[(q2>lo)&(q2<hi)]])
    k = np.searchsorted(union,center[0],side='right')-1; k = max(0,min(k,len(union)-2))
    if window=='local': lo,hi = union[max(0,k-1)],union[min(len(union)-1,k+2)]
    B = np.array([crossing(q1,a,N1,q2,b,N2,lo,hi,method) for a,b in zip(B1,B2)])
    valid = np.isfinite(B).all(1); f = float(valid.mean())
    record = dict(status='conditional_quotient' if f>=.99 else 'unresolved_bootstrap_crossing',N=N1,N2=N2,method=method,window=window,epsilon_over_pi=float(center[0]),xi_over_L=float(center[3]),Delta_eff=float(center[4]),SD=float(B[valid,4].std(ddof=1)) if valid.sum()>1 else None,MC95=np.quantile(B[valid],[.025,.975],axis=0).T.tolist() if valid.sum()>1 else None,valid_fraction=f,search_window=[float(lo),float(hi)],central_m2_pair=center[1:3].tolist())
    return record,B


def correction_fit(records,draws,minimum=6,omega=.789):
    keep = [i for i,r in enumerate(records) if r['N']>=minimum and r.get('status')=='conditional_quotient']
    if len(keep)<2: return dict(status='insufficient_pairs',minimum=minimum)
    N = np.array([records[i]['N'] for i in keep]); y = np.array([records[i]['Delta_eff'] for i in keep]); B = np.column_stack([draws[i][:,4] for i in keep]); valid = np.isfinite(B).all(1)
    if valid.mean()<.99: return dict(status='unresolved_joint_bootstrap',minimum=minimum,valid_fraction=float(valid.mean()),sizes=N.tolist())
    C = np.cov(B[valid],rowvar=False); e,Q = np.linalg.eigh(C)
    if e.min()<=e.max()*1e-10: return dict(status='rank_unresolved',minimum=minimum)
    W = (Q/np.sqrt(e))@Q.T; X = np.column_stack([np.ones(len(N)),N.astype(float)**(-omega)]); T = np.linalg.pinv(W@X)@W; coef = T@y; samples = B[valid]@T.T; res = W@(y-X@coef); stat = float(res@res); dof = len(N)-2; p = float(chi2.sf(stat,dof)) if dof else None
    return dict(status='zero_dof_diagnostic' if not dof else 'rejected_correction_fit' if p<.01 else 'conditional_correction_fit',minimum=minimum,sizes=N.tolist(),omega=omega,Delta=float(coef[0]),SD=float(samples[:,0].std(ddof=1)),MC95=np.quantile(samples[:,0],[.025,.975]).tolist(),chi2=stat,dof=dof,p=p,valid_fraction=float(valid.mean()),scope='WF-informed fixed omega assumption; no calibrated systematic error.')


def main(cut=1001,factor=2):
    folder = m.OUT/('common-xi-quotients-'+m.run.now()[:19].replace('-','').replace(':','').replace('T','-'));folder.mkdir()
    shutil.copyfile(__file__,folder/'common_xi_quotients.py'); groups = {}; sources = []
    for f in sorted((m.CAMPAIGN/'jobs').glob('*.json')):
        j = json.loads(f.read_text()); p = j['point']
        if p['D']!=3 or p['beta']!=0 or j['state']!='complete' or not j.get('audit',{}).get('passed'): continue
        q = float(p.get('epsilon',np.pi*p['epsilon_over_pi']))/np.pi
        # Outside this interval the preliminary scan is massive and irrelevant to these crossings.
        if not .79<=q<=.817: continue
        groups.setdefault((p['N'],q),[]).append(p); sources.append(dict(id=p['id'],point=p,audit=j['audit']))
    curves = {}; group_records = []
    for (N,q),points in sorted(groups.items()):
        s,B = summary(points,cut,factor); group_records.append(dict(N=N,q=q,summary=s))
        if B is not None: curves.setdefault(N,[]).append((q,s,B))
    records = []; samples = []; fits = []
    for method,window in [('linear','local'),('linear','full'),('pchip','local')]:
        rr=[];bb=[]
        for N in sorted(curves):
            if 2*N not in curves: continue
            args = []
            for size in [N,2*N]:
                gs=curves[size];args += [np.array([v[0] for v in gs]),np.array([v[1]['mean'] for v in gs]),np.stack([v[2] for v in gs],axis=1),size]
            r,B = paired_record(*args,method,window)
            r.update(N=N,N2=2*N,method=method,window=window)
            # Only flags from knots which support the retained search window (plus interpolation neighbors).
            flags=[]; contributing=[]
            if 'search_window' in r:
                for size in [N,2*N]:
                    gs=curves[size];q=np.array([v[0] for v in gs]);lo,hi=r['search_window'];left=max(0,np.searchsorted(q,lo,side='right')-1);right=min(len(q),np.searchsorted(q,hi,side='left')+1)
                    if method=='pchip': left=max(0,left-1);right=min(len(q),right+1)
                    for v in gs[left:right]: flags+=v[1]['flags'];contributing+=v[1]['ids']
            r['input_flags']=sorted(set(flags));r['contributing_source_ids']=sorted(set(contributing))
            records.append(r)
            if B is not None:
                np.savez_compressed(folder/f'pair-{N}-{2*N}-{method}-{window}.npz',samples=B)
                rr.append(r);bb.append(B)
        for minimum in [6,8,12,16,24]: fits.append(dict(method=method,window=window,**correction_fit(rr,bb,minimum)))
    report = dict(created=m.run.now(),settings=dict(cut=cut,factor=factor,draws=m.DRAWS,columns=['m2','m4','mabs','S0','F'],bootstrap='Independent within-start block resamples, pooled by rows. Same coordinate draws reused across all pairs.',omega=.789),source_jobs=sources,groups=group_records,pairs=records,fits=fits,scope='Existing beta=0 3D MC only. Common-epsilon xi/L crossings and m2 quotients. Interpolation and correction fits remain conditional; no numerical sigma_sys or total precision certified. No new simulations. Old anchors and new scans pooled at exact coordinates, not rounded.')
    m.run.atomic_json(folder/'analysis.json',report)
    plot(folder,curves,records)
    m.run.atomic_json(m.OUT/'latest-common-xi-quotients.json',dict(folder=folder.name,status='pending_independent_review'))
    print(folder,flush=True)
    for r in records: print(json.dumps(r),flush=True)
    for f in fits: print(json.dumps(f),flush=True)


def plot(folder,curves,records):
    plt=m.plt;fig,(ax,bx)=plt.subplots(1,2,figsize=(12,5))
    for N,gs in sorted(curves.items()):
        q=np.array([v[0] for v in gs]);v=np.array([v[1]['R'] for v in gs]);sd=np.array([v[1]['R_SD'] for v in gs]);keep=(q>=.8)&(q<=.804)
        ax.errorbar(q[keep],v[keep],yerr=sd[keep],fmt='o-',ms=3,lw=.8,capsize=2,label=f'L={N}')
    ax.set(xlabel=r'$\epsilon/\pi$ at $\beta=0$',ylabel=r'$\xi/L$ (bars: one bootstrap SD)',xlim=(.79995,.80405),ylim=(.15,.9));ax.legend(fontsize=8,ncol=2);ax.grid(alpha=.2)
    for method,window,offset,color in [('linear','local',-.15,'#2369a5'),('pchip','local',.15,'#c76422')]:
        rs=[r for r in records if r.get('MC95') and r['method']==method and r['window']==window]
        x=np.array([r['N'] for r in rs]);y=np.array([r['Delta_eff'] for r in rs]);ci=np.array([r['MC95'][4] for r in rs]);bx.errorbar(x+offset,y,yerr=np.maximum(0,np.array([y-ci[:,0],ci[:,1]-y])),fmt='o',capsize=3,label=method)
    bx.set(xlabel=r'$L$ in the pair $(L,2L)$',ylabel=r'$\Delta_{\rm eff}$ at common-$\epsilon$ crossing');bx.legend();bx.grid(alpha=.2)
    fig.suptitle('3D spin dimension: common-coupling correlation-length quotients')
    fig.text(.5,.01,'Quotient bars: conditional MC 95%; interpolation/correction sensitivity is not included as a systematic error.',ha='center',fontsize=9);fig.tight_layout(rect=(0,.04,1,.94))
    for ext in ['png','pdf','svg']:fig.savefig(folder/f'common-xi-quotients.{ext}',dpi=180)
    plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cut',type=int,default=1001);p.add_argument('--factor',type=int,default=2);a=p.parse_args();main(a.cut,a.factor)
