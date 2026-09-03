"""All-observable histories with explicit cold/prepared-start labels."""
from run import OUT,SIZES,atomic_json,sha256
from analyze import completed,read,derived
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def main():
    selected=[p for p in completed() if p['replica']<100];folder=OUT/'histories_final';folder.mkdir(exist_ok=True)
    for e in sorted({p['epsilon_milli'] for p in selected}):
        group=[p for p in selected if p['epsilon_milli']==e];keys=list(derived(read(group[0])));ks=sorted({p['kappa'] for p in group})
        signature=dict(source_sha256=sha256(Path(__file__)),chains=sorted(p['id'] for p in group))
        record=folder/f'e{e:04d}_inputs.json'
        if record.exists() and json.loads(record.read_text())==signature and all((folder/f'e{e:04d}_{label}.png').exists() for label in ['initial','full','retained']):
            continue
        colors={k:plt.cm.viridis(i/max(1,len(ks)-1)) for i,k in enumerate(ks)}
        for stop,label in [(1001,'initial'),(12001,'full'),(12001,'retained')]:
            fig,axes=plt.subplots(len(keys),3,figsize=(14,2*len(keys)),squeeze=False)
            for p in group:
                d=derived(read(p));col=SIZES.index(p['L']);color=colors[p['kappa']]
                for row,key in enumerate(keys):
                    start=2001 if label=='retained' else 0
                    y=d[key][start:stop];x=np.arange(start,start+len(y));ax=axes[row,col]
                    offset=0 if label=='retained' else 1;nb=(len(y)-offset)//50
                    xx=x[offset:offset+50*nb].reshape(-1,50).mean(1);yy=y[offset:offset+50*nb].reshape(-1,50)
                    ax.fill_between(xx,yy.min(1),yy.max(1),alpha=.035,color=color,lw=0)
                    ax.plot(xx,yy.mean(1),lw=.65,color=color,ls=['-','--',':','-.'][p['replica']%4],
                        label=f"κ={p['kappa']:g}" if p['replica']==0 else None)
                    if label!='retained':ax.plot(0,y[0],'.',ms=3,color=color)
                    if label=='full':ax.axvline(2000,color='.4',lw=.4,ls=':')
                    ax.set_title(f"{key}, L={p['L']}",fontsize=8);ax.grid(alpha=.12)
            axes[0,0].legend(fontsize=6,ncol=3)
            caption='retained cycles 2001–12000' if label=='retained' else 'row 0: cold and prepared winding starts'
            fig.suptitle(f'ε/π={e/1000:g}; {caption}; 50-cycle means and full min/max ranges',y=.997,fontsize=11)
            fig.tight_layout(rect=(0,0,1,.985));fig.savefig(folder/f'e{e:04d}_{label}.png',dpi=90);plt.close(fig)
            print('saved',e,label,flush=True)
        atomic_json(record,signature)

if __name__=='__main__':main()
