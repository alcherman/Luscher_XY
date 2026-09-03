"""Compare measured finite-size observables across the complete scouting grid."""
import json
from run import OUT, SIZES, EPSILON_MILLI, KAPPA_MILLI
import numpy as np
import matplotlib.pyplot as plt


def edges(x):
    x=np.asarray(x,dtype=float)
    return np.r_[x[0]-(x[1]-x[0])/2,(x[1:]+x[:-1])/2,x[-1]+(x[-1]-x[-2])/2]


def main():
    rows=json.loads((OUT/'results.json').read_text())
    assert len(rows)==588
    table={(r['point']['L'],r['point']['epsilon_milli'],r['point']['kappa_milli']):r['summary']['values'] for r in rows}
    ex=edges(np.array(EPSILON_MILLI)/1000);ky=edges(np.array(KAPPA_MILLI)/1000)
    for key,title,limits in [('binder','Binder cumulant',(-.05,.5)),('xi_over_L',r'$\xi_{2nd}/L$',(0,1.2))]:
        fig,axes=plt.subplots(1,3,figsize=(13,4.6),layout='constrained',sharey=True)
        for ax,L in zip(axes,SIZES):
            z=np.array([[table[L,e,k][key] for e in EPSILON_MILLI] for k in KAPPA_MILLI],dtype=float)
            mesh=ax.pcolormesh(ex,ky,z,cmap='viridis',vmin=limits[0],vmax=limits[1],rasterized=True)
            ax.set(xlabel=r'$\epsilon/\pi$',title=f'L={L}',xlim=(.495,1.01),ylim=(-.05,1.25))
        axes[0].set_ylabel(r'$\kappa$')
        fig.colorbar(mesh,ax=axes,label=title,shrink=.85)
        fig.suptitle(f'{title}: measured grid values, 8,000 retained hybrid cycles per point')
        for ext in ('png','pdf'):fig.savefig(OUT/f'{key}_heatmap.{ext}',dpi=160)
        plt.close(fig)


if __name__=='__main__':main()
