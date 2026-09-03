"""Disk-only scouting diagnostics. Inspect histories before choosing --cut."""
import argparse
import csv
import json
from pathlib import Path
import sys
import h5py
import numpy as np
from run import OUT, SIZES, EPSILON_MILLI, KAPPA_MILLI, atomic_json, sha256, now
sys.path.insert(0,str(OUT.parent/'2026-09-02-xy-binder-bkt'))
from experiment import tau_int
import matplotlib.pyplot as plt

COLORS={16:'#0072B2',32:'#D55E00',64:'#009E73'}


def completed():
    return sorted([j['point'] for p in (OUT/'jobs').glob('*.json')
                   if (j:=json.loads(p.read_text())).get('state')=='complete'],
                   key=lambda p:(p['epsilon_milli'],p['kappa_milli'],p['L']))


def read(p):
    with h5py.File(OUT/'chains'/f"{p['id']}.h5",'r') as f:
        assert f.attrs['committed_rows']==p['steps']+1
        return {k:v[:] for k,v in f.items() if k not in ('phi','rng_state')}


def finite(value):
    if isinstance(value,dict):return {k:finite(v) for k,v in value.items()}
    if isinstance(value,list):return [finite(v) for v in value]
    if isinstance(value,(float,np.floating)):return float(value) if np.isfinite(value) else None
    if isinstance(value,np.integer):return int(value)
    return value


def derived(d):
    return dict(m2=d['magnetization_squared'],m4=d['magnetization_fourth'],
        m_real=d['magnetization_real'],m_imag=d['magnetization_imag'],m_abs=d['magnetization_abs'],
        F0=d['structure_zero'],Fx=d['structure_min'][:,0],Fy=d['structure_min'][:,1],
        link_cosine=d['mean_link_cosine'],action_density=d['action_density'],
        vortex_density=d['vortex_density'],max_link_angle=d['max_link_angle'],
        winding_mean=d['winding'].mean((1,2)),winding_squared=(d['winding'].astype(float)**2).mean((1,2)),
        winding_range=np.ptp(d['winding'],axis=(1,2)),local_accepted=d['local_accepted'],
        cluster_size=d['cluster_sizes'].sum(1),cluster_bond_tests=d['cluster_bond_tests'].sum(1),
        winding_accepted=d['winding_accepted'])


def histories(selected,focus=False):
    folder=OUT/'histories';folder.mkdir(exist_ok=True)
    if focus:
        selected=[p for p in selected if p['kappa_milli']==0 and p['epsilon_milli']<=620]
        eps=sorted({p['epsilon_milli'] for p in selected})
        datasets={p['id']:read(p) for p in selected}
        for stop,label in [(1001,'cold'),(max(p['steps'] for p in selected)+1,'full')]:
            for key in ('magnetization_squared','mean_link_cosine','vortex_density'):
                fig,axes=plt.subplots(len(eps),3,figsize=(13,2.1*len(eps)),squeeze=False)
                for row,e in enumerate(eps):
                    for col,L in enumerate(SIZES):
                        name=f'e{e:04d}_k0000_L{L:03d}';ax=axes[row,col]
                        if name not in datasets:continue
                        y=datasets[name][key][:stop];x=np.arange(len(y));ax.plot(x,y,lw=.35,alpha=.5,color=COLORS[L])
                        bins=(len(y)-1)//50
                        ax.plot(x[1:1+50*bins].reshape(-1,50).mean(1),y[1:1+50*bins].reshape(-1,50).mean(1),'k-',lw=.8)
                        ax.plot(0,y[0],'ko',ms=2);ax.set_title(f'ε/π={e/1000:g}, L={L}',fontsize=9)
                        ax.grid(alpha=.15)
                fig.suptitle(f'κ=0: {key}; seed at cycle 0, black = 50-cycle means')
                fig.tight_layout(rect=(0,0,1,.985));fig.savefig(folder/f'focus_{key}_{label}.png',dpi=100);plt.close(fig)
        return
    # Each page includes every physical observable and the update/mixing diagnostics.
    # Raw points are retained on disk; plotting all 10,001 samples does not thin them.
    for e in sorted({p['epsilon_milli'] for p in selected}):
        group=[p for p in selected if p['epsilon_milli']==e]
        if not group:continue
        keys=list(derived(read(group[0])))
        for stop,label in [(1001,'cold'),(max(p['steps'] for p in selected)+1,'full')]:
            path=folder/f'e{e:04d}_n{len(group):02d}_{label}.png'
            if path.exists():continue
            fig,axes=plt.subplots(len(keys),3,figsize=(14,2*len(keys)),squeeze=False)
            for p in group:
                d=derived(read(p));col=SIZES.index(p['L']);color=plt.cm.viridis(p['kappa']/1.2)
                for row,key in enumerate(keys):
                    y=d[key][:stop];x=np.arange(len(y));ax=axes[row,col]
                    bins=(len(y)-1)//50
                    xb=x[1:1+50*bins].reshape(-1,50).mean(1)
                    yb=y[1:1+50*bins].reshape(-1,50)
                    ax.fill_between(xb,yb.min(1),yb.max(1),alpha=.035,color=color,lw=0)
                    ax.plot(xb,yb.mean(1),
                            lw=.8,color=color,label=f"κ={p['kappa']:g}")
                    ax.plot(0,y[0],'.',ms=3,color=color)
                    ax.set_title(f'{key}, L={p["L"]}',fontsize=8);ax.grid(alpha=.12)
            axes[0,0].legend(fontsize=6,ncol=3)
            fig.suptitle(f'ε/π={e/1000:g}; 50-cycle means and min/max ranges of all samples; cold seed at 0',y=.997)
            fig.tight_layout(rect=(0,0,1,.985));fig.savefig(path,dpi=85);plt.close(fig)
            print('history',path.name,flush=True)


def quantities(means,L):
    a,b,f0,fx,fy,cosine,vortex,w2=means.T
    U=1-b/(2*a*a)
    ratio=f0/((fx+fy)/2)-1
    with np.errstate(invalid='ignore'):
        R=np.sqrt(ratio)/(2*L*np.sin(np.pi/L))
    return np.array([U,R,a,b,cosine,vortex,w2]).T


def summarize(d,p,cut,draws=2000,width_factor=1,seed_offset=0):
    a=d['magnetization_squared'][cut:];b=d['magnetization_fourth'][cut:]
    f0=d['structure_zero'][cut:];fx,fy=d['structure_min'][cut:].T;fm=(fx+fy)/2
    cosine=d['mean_link_cosine'][cut:];vortex=d['vortex_density'][cut:]
    w2=(d['winding'][cut:].astype(float)**2).mean((1,2))
    iu=b.mean()/a.mean()**3*(a-a.mean())-(b-b.mean())/(2*a.mean()**2)
    # Correlation-length influence differs from the ratio influence by a constant.
    ir=(f0-f0.mean())/fm.mean()-f0.mean()/fm.mean()**2*(fm-fm.mean())
    monitored=dict(m2=a,m4=b,binder=iu,xi_ratio=ir,link_cosine=cosine,vortex=vortex)
    taus={k:tau_int(v) for k,v in monitored.items()}
    width=width_factor*max(20,int(np.ceil(10*max(t for t in taus.values() if t is not None))))
    n=len(a);nb=n//width
    values=np.column_stack((a,b,f0,fx,fy,cosine,vortex,w2))
    central=quantities(values.mean(0)[None,:],p['L'])[0]
    samples=np.full((draws,7),np.nan)
    if nb>=2:
        edges=np.linspace(0,n,nb+1,dtype=int);lengths=np.diff(edges)
        sums=np.add.reduceat(values,edges[:-1],axis=0)
        rng=np.random.default_rng(np.random.SeedSequence([p['seed'],cut,width,seed_offset]))
        for start in range(0,draws,100):
            stop=min(start+100,draws);ix=rng.integers(nb,size=(stop-start,nb))
            means=sums[ix].sum(1)/lengths[ix].sum(1)[:,None]
            samples[start:stop]=quantities(means,p['L'])
    keys=['binder','xi_over_L','m2','m4','link_cosine','vortex_density','winding_squared']
    valid=np.isfinite(samples)
    errors={k:float(np.std(samples[valid[:,i],i],ddof=1)) if valid[:,i].sum()>1 else None for i,k in enumerate(keys)}
    invalid_R=1-float(valid[:,1].mean())
    result=dict(cut=cut,retained=n,blocks=nb,block_width=width,taus=taus,
        enough_blocks=nb>=20,values=dict(zip(keys,map(float,central))),errors=errors,
        invalid_xi_bootstrap_fraction=invalid_R,
        ess={k:n/(2*t) if t else None for k,t in taus.items()},
        winding_changes=int(np.count_nonzero(np.any(np.diff(d['winding'][cut:],axis=0),axis=(1,2)))),
        local_acceptance=float(d['local_accepted'][cut:].mean()/p['L']**2),
        mean_cluster_fraction=float(d['cluster_sizes'][cut:].mean()/p['L']**2),
        winding_acceptance=float(d['winding_accepted'][cut:].mean()),
        vortex_present_cycles=int(np.count_nonzero(vortex)),
        winding_tau=tau_int(w2))
    return finite(result),samples


def analyze(selected,cut):
    folder=OUT/'analysis';folder.mkdir(exist_ok=True)
    rows=[]
    for p in selected:
        path=folder/f"{p['id']}.json"
        if path.exists():
            row=json.loads(path.read_text())
            if row['summary']['cut']==cut and row['analysis_source_sha256']==sha256(Path(__file__)):
                rows.append(row);continue
        d=read(p);summary,samples=summarize(d,p,cut)
        np.savez_compressed(folder/f"{p['id']}_bootstrap.npz",samples=samples)
        checks=[]
        for c in sorted(set([1001,cut,4001])):
            for factor in ([1,2] if c==cut else [1]):
                check,_=summarize(d,p,c,draws=500,width_factor=factor,seed_offset=1)
                checks.append(check)
        flags=[]
        if not summary['enough_blocks']:flags.append('fewer_than_20_blocks')
        for check in checks:
            if check['cut']==cut:continue
            for key in ['binder','xi_over_L','m2','link_cosine','vortex_density']:
                a,b=summary['values'][key],check['values'][key]
                ea,eb=summary['errors'][key],check['errors'][key]
                if None not in (a,b,ea,eb) and abs(a-b)>2*np.hypot(ea,eb):
                    flags.append(f'cut_sensitive_{key}')
        # Disjoint retained halves provide a less conservative drift check than
        # comparing strongly overlapping cuts with independent-error arithmetic.
        middle=cut+(len(d['sweep'])-cut)//2
        first,_=summarize({k:v[:middle] for k,v in d.items()},p,cut,draws=500,seed_offset=2)
        second,_=summarize(d,p,middle,draws=500,seed_offset=3)
        for key in ['binder','xi_over_L','m2','link_cosine','vortex_density']:
            a,b=first['values'][key],second['values'][key]
            ea,eb=first['errors'][key],second['errors'][key]
            if None not in (a,b,ea,eb) and abs(a-b)>3*np.hypot(ea,eb):
                flags.append(f'half_sensitive_{key}')
        row=dict(point=p,summary=summary,stability=checks,halves=[first,second],flags=sorted(set(flags)),
            analysis_source_sha256=sha256(Path(__file__)),status='independent-start refinement trajectory; pooling and phase inference separate')
        atomic_json(path,row);rows.append(row)
        print('analyzed',p['id'],summary['blocks'],'blocks',row['flags'],flush=True)
    atomic_json(OUT/'results.json',rows)
    with (OUT/'observables.csv').open('w',newline='') as f:
        fields=['id','L','epsilon_over_pi','kappa','binder','binder_error','xi_over_L','xi_error',
                'm2','m2_error','vortex_density','vortex_error','blocks','flags']
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for r in rows:
            p,s=r['point'],r['summary'];v,se=s['values'],s['errors']
            writer.writerow({k:p[k] for k in fields[:4]}|dict(binder=v['binder'],binder_error=se['binder'],
                xi_over_L=v['xi_over_L'],xi_error=se['xi_over_L'],m2=v['m2'],m2_error=se['m2'],
                vortex_density=v['vortex_density'],vortex_error=se['vortex_density'],blocks=s['blocks'],
                flags=';'.join(r['flags'])))
    atomic_json(OUT/'analysis-provenance.json',dict(updated_utc=now(),chains=len(rows),cut=cut,
        source_sha256=sha256(Path(__file__)),draws=2000,
        block_rule='max(20,ceil(10*largest tau among m2,m4,U influence,xi ratio influence,cosine,vortex))',
        estimator='joint balanced contiguous block bootstrap, all retained rows; ratios of means',
        scope='per-trajectory refinement diagnostics; no thermodynamic certification'))
    return rows



if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cut',type=int)
    parser.add_argument('--histories',action='store_true')
    parser.add_argument('--batch')
    args=parser.parse_args()
    selected=completed()
    if args.batch:
        ids={p['id'] for p in json.loads((OUT/'batches'/f'{args.batch}.json').read_text())['points']}
        selected=[p for p in selected if p['id'] in ids]
    if args.histories:histories(selected)
    if args.cut is not None:analyze(selected,args.cut)
