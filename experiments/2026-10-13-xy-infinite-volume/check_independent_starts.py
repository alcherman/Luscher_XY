"""Retain history plots and cut/block checks for targeted start discrepancies."""
import os
for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[name]='1'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
os.environ['MPLBACKEND']='Agg'
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import run
from analyze_calibration import read,summarize
from scout import json_finite

OUT=Path(__file__).resolve().parent


def main():
    stamp=run.now()[:19].replace('-','').replace(':','').replace('T','-')
    folder=OUT/('independent-start-check-'+stamp);folder.mkdir()
    results=[]
    for D,q,N,beta in [(3,92,16,293),(4,98,10,261)]:
        points=[];pending=[]
        expected=[p for f in sorted((OUT/'batches').glob('*.json')) for p in json.loads(f.read_text())['points']
                  if (p['D'],round(p['epsilon_over_pi']*100),p['N'],round(p['beta']*1000))==(D,q,N,beta)]
        for point in expected:
            f=OUT/'jobs'/f"{point['id']}.json"
            j=json.loads(f.read_text()) if f.exists() else dict(state='pending',point=point)
            if j['state']=='complete' and j['audit']['passed']:points.append(j['point'])
            else:pending.append(dict(id=j['point']['id'],state=j['state']))
        points=sorted(points,key=lambda p:p['replica'])
        arrays=[read(p) for p in points]
        if len(arrays)<2:continue
        stem=f'd{D}_e{q:03d}_N{N:03d}_b{beta:04d}'
        cuts=[]
        for cut,factor in [(501,1),(1001,1),(2001,1),(1001,2)]:
            summary,boot=summarize(arrays,points,cut,factor=factor)
            if boot is not None:np.savez_compressed(folder/f'{stem}_cut{cut}_factor{factor}.npz',samples=boot)
            cuts.append(dict(cut=cut,block_factor=factor,summary=summary))
        starts=[]
        fig,axes=plt.subplots(4,1,figsize=(10,9),sharex=True)
        for a,p in zip(arrays,points):
            estimates=[]
            for cut in [1001,2001,5001,8001]:
                if len(a)<=cut+1000:continue
                m2,m4=a[cut:,:2].mean(0)
                estimates.append(dict(cut=cut,retained=len(a)-cut,U=float(1-m4/(2*m2*m2))))
            starts.append(dict(id=p['id'],replica=p['replica'],initial_winding=p['initial_winding'],estimates=estimates))
            length=(len(a)-1)//50
            for ax,index in zip(axes,[0,4,6]):
                values=a[1:1+50*length,index].reshape(length,50).mean(1)
                ax.plot(np.arange(length)*50+25,values,lw=.8,label=f"r{p['replica']}: winding {p['initial_winding']}")
            retained=a[1001:];counts=np.arange(100,len(retained)+1,100)
            if len(counts):
                sum_m2=np.cumsum(retained[:,0])[counts-1];sum_m4=np.cumsum(retained[:,1])[counts-1]
                running_U=1-counts*sum_m4/(2*sum_m2**2)
                axes[3].plot(1000+counts,running_U,lw=.9,label=f"r{p['replica']}")
        for ax,label in zip(axes,['Magnetization squared','Mean link cosine','Mean squared winding']):
            ax.set_ylabel(label);ax.grid(alpha=.2)
            ax.axvline(1001,color='black',ls='--',lw=.8);ax.axvline(2001,color='gray',ls=':',lw=.8)
        axes[0].legend(fontsize=8);axes[-1].set_xlabel('Cycle (upper panels: 50-cycle means; bottom: cumulative Binder estimate)')
        axes[3].set_ylabel('Running Binder U\n(after cut1001)');axes[3].grid(alpha=.2)
        axes[3].axhline(.25,color='gray',ls=':',lw=.8)
        fig.suptitle(f'Independent-start diagnostic: d={D}, epsilon/pi={q/100}, N={N}, beta={beta/1000}')
        fig.tight_layout()
        for extension in ['png','pdf','svg']:fig.savefig(folder/f'{stem}_histories.{extension}',dpi=150)
        plt.close(fig)
        results.append(dict(D=D,epsilon_over_pi=q/100,N=N,beta=beta/1000,starts=starts,pending=pending,cut_block_checks=cuts))
    run.atomic_json(folder/'diagnostics.json',json_finite(dict(created=run.now(),results=results,source_sha256=run.sha256(Path(__file__)),scope='Targeted below-crossing diagnostic; start/cut/block flags are not by themselves proof of root bias or equilibration. No raw histories thinned. New independent replicas may still be pending.')))
    run.atomic_json(OUT/'latest-independent-start-check.json',dict(created=run.now(),folder=folder.name))
    print('independent_start_snapshot',folder,flush=True)
    for r in results:print(r['D'],r['epsilon_over_pi'],'closed_starts',len(r['starts']),'pending',len(r['pending']),flush=True)


if __name__=='__main__':main()
