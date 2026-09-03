"""Disk-only comparison: inspect histories before selecting a retained interval."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

from run import OUT, ARMS, LABELS, COLORS, STEPS, atomic_json
PREVIOUS=OUT.parent/'2026-09-02-xy-binder-bkt'
sys.path.insert(0,str(PREVIOUS))
import experiment as previous
from analyze import bootstrap
import h5py
import matplotlib.pyplot as plt
import numpy as np


def read(arm):
    with h5py.File(OUT/'chains'/f"{arm['name']}.h5",'r') as f:
        assert f.attrs['committed_rows']==STEPS+1
        assert all(len(ds)==STEPS+1 for ds in f.values())
        return {k:ds[:] for k,ds in f.items() if k not in ('phi','rng_state')}


def histories(names):
    previous.style()
    arms=[a for a in ARMS if a['name'] in names]
    directory=OUT/'histories'; directory.mkdir(exist_ok=True)
    for arm in arms:
        data=read(arm); keys=[k for k in data if k!='sweep']
        fig,axes=plt.subplots((len(keys)+2)//3,3,figsize=(15,22),squeeze=False)
        for ax,key in zip(axes.flat,keys):
            ax.plot(data['sweep'],data[key].reshape(STEPS+1,-1),lw=.4,alpha=.5)
            ax.set_title(key,fontsize=9);ax.set_xlabel('cycle (cold seed at 0)',fontsize=8)
        for ax in list(axes.flat)[len(keys):]: ax.axis('off')
        fig.suptitle(f"{LABELS[arm['name']]}: all recorded observables, L48 κ=.80 ε=π")
        fig.tight_layout();fig.savefig(directory/f"{arm['name']}.png");plt.close(fig)
    for stop,name in ((1001,'thermalization'),(STEPS+1,'overview')):
        fig,axes=plt.subplots(len(arms),3,figsize=(13,3.5*len(arms)),squeeze=False)
        for row,arm in enumerate(arms):
            data=read(arm)
            for col,(key,label) in enumerate((('magnetization_squared',r'$|m|^2$'),
                    ('action_density',r'$S/L^2$'),('vortex_density','vortex density'))):
                ax=axes[row,col]; t=data['sweep'][:stop];y=data[key][:stop]
                ax.plot(t,y,lw=.45,alpha=.6,color=COLORS[arm['name']])
                bins=(len(t)-1)//50
                ax.plot(t[1:1+50*bins].reshape(-1,50).mean(1),
                        y[1:1+50*bins].reshape(-1,50).mean(1),color='black',lw=.9)
                ax.plot(0,y[0],'ko',ms=3)
                if col==0: ax.set_yscale('log')
                ax.set_title(f"{LABELS[arm['name']]}: {label}",fontsize=9)
                ax.grid(alpha=.15)
        for ax in axes[-1]:ax.set_xlabel('fixed update cycle; all arms add one winding proposal')
        fig.suptitle('D=2, L=48, κ=0.80, ε=π; cold seed retained, black line = 50-cycle mean')
        fig.tight_layout(rect=(0,0,1,.965));fig.savefig(OUT/f'{name}.png');plt.close(fig)


def series(data,cut):
    a,b=data['magnetization_squared'][cut:],data['magnetization_fourth'][cut:]
    influence=b.mean()/a.mean()**3*(a-a.mean())-(b-b.mean())/(2*a.mean()**2)
    return dict(m2=a,m4=b,binder_influence=influence,energy=data['action_density'][cut:],
                vortex=data['vortex_density'][cut:])


def tau_and_width(data,cut):
    values=series(data,cut)
    taus={k:previous.tau_int(v) for k,v in values.items()}
    width=max(20,int(np.ceil(10*max(t for t in taus.values() if t is not None))))
    return values,taus,width


def acf(x):
    x=np.asarray(x,dtype=float)-np.mean(x);n=len(x)
    nfft=1<<(2*n-1).bit_length();f=np.fft.rfft(x,nfft)
    c=np.fft.irfft(f*f.conj(),nfft)[:n]/np.arange(n,0,-1)
    return c/c[0]


def analyze(cut):
    previous.style();rows=[];checks=[]
    status=json.loads((OUT/'status.json').read_text())
    cost={r['name']:r['seconds'] for r in status['completed']}
    for index,arm in enumerate(ARMS):
        data=read(arm);values,taus,width=tau_and_width(data,cut)
        result=bootstrap(values['m2'],values['m4'],width,seed=2026092001+index,include_tail=True)
        samples=result.pop('samples',None)
        if samples is None: raise ValueError('Not enough blocks to compare this arm')
        np.save(OUT/f"bootstrap_{arm['name']}.npy",samples)
        n=STEPS+1-cut;ess={k:n/(2*t) for k,t in taus.items() if t is not None}
        total_seconds=cost[arm['name']]
        sampler_seconds=float(data['update_seconds'][cut:].sum())
        means={k:float(v.mean()) for k,v in values.items() if k!='binder_influence'}
        # Equal-length batch means for scalar uncertainty diagnostics only;
        # Binder bootstrap above uses every retained sample with balanced blocks.
        blocks=n//width
        scalar_se={k:float(v[:blocks*width].reshape(blocks,width).mean(1).std(ddof=1)/np.sqrt(blocks))
                   for k,v in values.items() if k!='binder_influence'}
        row=arm|result|dict(first_retained_cycle=cut,retained=n,tau=taus,block_width=width,
                means=means,scalar_block_se=scalar_se,ess=ess,total_seconds=total_seconds,
                retained_update_seconds=sampler_seconds,
                ess_per_total_second={k:e/total_seconds for k,e in ess.items()},
                ess_per_update_second={k:e/sampler_seconds for k,e in ess.items()},
                binder_variance_times_total_seconds=result['error']**2*total_seconds,
                retained_timing={k:float(data[k][cut:].sum()) for k in data if k.endswith('_seconds')},
                total_timing={k:float(data[k].sum()) for k in data if k.endswith('_seconds')},
                local_acceptance=(float(data['local_accepted'][cut:].mean()/48**2)
                                  if 'local_accepted' in data else None),
                mean_cluster_fraction=(float(data['cluster_sizes'][cut:].mean()/48**2)
                                       if 'cluster_sizes' in data else None),
                winding_acceptance=float(data['winding_accepted'][cut:].mean()))
        rows.append(row)
        for other in sorted(set([cut,501,1001,2001,4001])):
            v,t,w=tau_and_width(data,other)
            for factor in (1,2,4):
                check=bootstrap(v['m2'],v['m4'],w*factor,seed=2026093000+index*100+other+factor,
                                draws=2000,include_tail=True)
                check.pop('samples',None)
                checks.append(dict(arm=arm['name'],cut=other,factor=factor,tau=t,block_width=w*factor)|check)
        for quarter,positions in enumerate(np.array_split(np.arange(n),4),1):
            check=bootstrap(values['m2'][positions],values['m4'][positions],width,
                            seed=2026096000+index*10+quarter,draws=2000,include_tail=True)
            check.pop('samples',None)
            checks.append(dict(arm=arm['name'],cut=cut,quarter=quarter)|check)
        print(arm['name'], 'U',row['binder'],'+/-',row['error'],'tau',taus,'ESS/s',row['ess_per_total_second'],flush=True)
    atomic_json(OUT/'results.json',rows);atomic_json(OUT/'stability.json',checks)
    with (OUT/'comparison.csv').open('w',newline='') as f:
        names=['name','binder','error','blocks','block_width','total_seconds',
               'tau_energy','tau_m2','tau_binder','ess_s_energy','ess_s_m2','ess_s_binder']
        writer=csv.DictWriter(f,fieldnames=names,extrasaction='ignore');writer.writeheader()
        for r in rows:
            writer.writerow(r|dict(tau_energy=r['tau']['energy'],tau_m2=r['tau']['m2'],
                tau_binder=r['tau']['binder_influence'],ess_s_energy=r['ess_per_total_second']['energy'],
                ess_s_m2=r['ess_per_total_second']['m2'],ess_s_binder=r['ess_per_total_second']['binder_influence']))
    plot(rows,cut,checks)
    atomic_json(OUT/'analysis-provenance.json',dict(cut=cut,
        analysis_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        bootstrap_source_sha256=hashlib.sha256((PREVIOUS/'analyze.py').read_bytes()).hexdigest(),
        tau_source_sha256=hashlib.sha256((PREVIOUS/'experiment.py').read_bytes()).hexdigest(),
        block_rule='max(20,ceil(10*max(tau_m2,tau_m4,tau_Binder_influence,tau_energy,tau_vortex)))',
        bootstrap='4000 joint block resamples, ratio recomputed, balanced blocks keep every retained row',
        efficiency='retained n/(2*tau) divided by full generation wall time including cold start, measurements and I/O',
        limits='one seed per arm, single-thread sequential runs, no uncertainty on tau or timing ratios; fixed tuning'))


def plot(rows,cut,checks):
    names=[r['name'] for r in rows];short=['Metropolis','Hybrid','100 clusters']
    fig,axes=plt.subplots(1,3,figsize=(15,4.9))
    for i,r in enumerate(rows):
        axes[0].errorbar(i,r['binder'],yerr=r['error'],fmt='o',capsize=4,color=COLORS[r['name']])
    axes[0].set_xticks(range(3),short);axes[0].set_ylabel('Binder cumulant U (1σ block bootstrap)')
    axes[0].set_title('Same equilibrium target')
    for i,key in enumerate(['energy','m2','binder_influence']):
        axes[1].bar(np.arange(3)+(i-1)*.25,[r['ess_per_total_second'][key] for r in rows],
                    width=.24,label={'energy':'Energy','m2':r'$|m|^2$','binder_influence':'Binder influence'}[key])
    axes[1].set_xticks(range(3),short);axes[1].set_ylabel('Estimated effective samples / second')
    axes[1].set_title('Full runtime, including warmup and I/O');axes[1].legend(frameon=False,fontsize=8)
    bottom=np.zeros(3)
    for key,label in [('metropolis_seconds','Metropolis'),('cluster_seconds','Cluster'),
                      ('winding_seconds','Winding'),('measurement_seconds','Measurement')]:
        value=np.array([r['total_timing'][key] for r in rows])
        axes[2].bar(np.arange(3),value,bottom=bottom,label=label);bottom+=value
    overhead=np.array([r['total_seconds'] for r in rows])-bottom
    axes[2].bar(np.arange(3),overhead,bottom=bottom,label='I/O and other overhead')
    axes[2].set_xticks(range(3),short);axes[2].set_ylabel('Total wall time (s)');axes[2].set_title('Measured runtime components')
    axes[2].legend(frameon=False,fontsize=8)
    for ax in axes:ax.grid(axis='y',alpha=.15)
    fig.suptitle('D=2, L=48, κ=0.80, ε=π — one cold-start chain per method')
    fig.text(.5,.012,f'10,000 fixed cycles per method; discard first {cut-1}; all methods add one winding proposal per cycle',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.045,1,.95))
    for extension in ('png','pdf','svg'):fig.savefig(OUT/f'comparison.{extension}')
    plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(13,4.3))
    for arm in ARMS:
        data=read(arm);v=series(data,cut)
        for ax,key in zip(axes,['energy','m2','binder_influence']):
            rho=acf(v[key]); ax.plot(np.arange(101),rho[:101],label=LABELS[arm['name']],color=COLORS[arm['name']])
            ax.set_title({'energy':'Energy','m2':r'$|m|^2$','binder_influence':'Binder influence'}[key])
            ax.set_xlabel('Lag in method-specific cycles');ax.grid(alpha=.15)
    axes[0].set_ylabel('Normalized autocorrelation');axes[1].legend(frameon=False,fontsize=8)
    fig.tight_layout();fig.savefig(OUT/'autocorrelations.png');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.7))
    for r in rows:
        g=[c for c in checks if c['arm']==r['name'] and c.get('factor')==1]
        axes[0].errorbar([c['cut']-1 for c in g],[c['binder'] for c in g],
                        yerr=[c['error'] for c in g],fmt='o-',capsize=3,label=LABELS[r['name']],color=COLORS[r['name']])
        g=[c for c in checks if c['arm']==r['name'] and c.get('cut')==cut and 'factor' in c]
        axes[1].plot([c['factor'] for c in g],[c['error'] for c in g],'o-',label=LABELS[r['name']],color=COLORS[r['name']])
    axes[0].set(xlabel='Discarded cycles',ylabel='Binder cumulant U');axes[1].set(xlabel='Block-width multiplier',ylabel='Bootstrap standard error')
    for ax in axes:ax.grid(alpha=.15);ax.legend(frameon=False,fontsize=8)
    fig.suptitle('Cut and block sensitivity: correlated estimates from the same chains')
    fig.tight_layout();fig.savefig(OUT/'stability.png');plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--histories',action='store_true');parser.add_argument('--cut',type=int)
    parser.add_argument('--arms',nargs='+',default=[a['name'] for a in ARMS])
    args=parser.parse_args()
    if args.histories:histories(args.arms)
    if args.cut is not None:
        if not 1<=args.cut<STEPS:raise ValueError('Invalid analysis cut')
        analyze(args.cut)
