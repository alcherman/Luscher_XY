"""Display finite-box quarter crossings and saved leading size fits."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
os.environ['MPLBACKEND']='Agg'
import argparse,json,shutil
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import run
from threshold_calibration import coordinate
parser=argparse.ArgumentParser();parser.add_argument('snapshot',type=Path);args=parser.parse_args()
a=json.loads((args.snapshot/'analysis.json').read_text())
folder=args.snapshot/'finite-size-panels';folder.mkdir(exist_ok=True)
shutil.copyfile(__file__,folder/Path(__file__).name)
for D in [3,4]:
    qs=sorted({r['q'] for r in a['locations'] if r['D']==D})
    fig,axes=plt.subplots((len(qs)+2)//3,3,figsize=(12,2.8*((len(qs)+2)//3)),squeeze=False)
    for ax,q in zip(axes.flat,qs):
        rows=sorted([r for r in a['locations'] if (r['D'],r['q'],r['threshold'],r['interpolation'])==(D,q,.25,'linear') and r['beta_location'] is not None],key=lambda r:r['N'])
        ax.set_title(f'epsilon/pi={q/100:.2f}',fontsize=10)
        if rows:
            ns=np.array([r['N'] for r in rows]);x=coordinate(D,ns);y=np.array([r['beta_location'] for r in rows]);ci=np.array([r['bootstrap_95'] for r in rows])
            ax.errorbar(x,y,yerr=np.maximum(0,np.array([y-ci[:,0],ci[:,1]-y])),fmt='o',ms=3,capsize=2,color='black')
            for xi,yi,N in zip(x,y,ns):
                offset=-13 if N==20 and 24 in ns else 5
                ax.annotate(str(N),(xi,yi),xytext=(3,offset),textcoords='offset points',fontsize=7)
            fit=next((f for f in a['fits'] if (f['D'],f['q'],f['threshold'],f['interpolation'],f['minimum_N'],f['correction'])==(D,q,.25,'linear',min(r['N'] for r in rows if r['status']=='conditional_crossing'),False)),None) if any(r['status']=='conditional_crossing' for r in rows) else None
            if fit:
                r=fit['result'];xx=np.linspace(0,max(x)*1.05,100);rejected=r['p'] is not None and r['p']<.01;color='firebrick' if rejected else 'tab:blue'
                ax.plot(xx,r['coefficients'][0]+r['coefficients'][1]*xx,'--' if rejected else '-',color=color,lw=1)
                lo,hi=r['bootstrap_95'];ax.errorbar([0],[r['beta_infinite']],yerr=[[r['beta_infinite']-lo],[hi-r['beta_infinite']]],fmt='x' if rejected else 's',color=color,capsize=2,ms=4)
                ax.text(.98,.96,f"chi2/dof={r['chi2']:.1f}/{r['df']}"+('\nrejected diagnostic' if rejected else '\nconditional fit'),transform=ax.transAxes,ha='right',va='top',fontsize=7,color=color,bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=2))
        else:ax.text(.5,.5,'No bracketed U=1/4 crossing\nNo infinite-volume estimate',ha='center',va='center',transform=ax.transAxes,fontsize=9)
        ax.set_xlabel('N^(-1/nu)' if D==3 else 'N^(-2) [ln N]^(-1/10)',fontsize=8);ax.set_ylabel('beta at U=1/4',fontsize=8);ax.grid(alpha=.2);ax.tick_params(labelsize=8)
    for ax in list(axes.flat)[len(qs):]:ax.set_visible(False)
    fig.suptitle(f'd={D}: finite-box Binder-quarter locations and leading size extrapolations\nBars: conditional block-bootstrap 95% intervals; labels: N. Red dashed fits are rejected diagnostics.',fontsize=11)
    fig.tight_layout(rect=(0,0,1,.95))
    for ext in ['png','pdf','svg']:fig.savefig(folder/f'd{D}-quarter-size-fits.{ext}',dpi=160)
    plt.close(fig)
run.atomic_json(folder/'provenance.json',dict(created=run.now(),source=str(args.snapshot/'analysis.json'),source_sha256=run.sha256(args.snapshot/'analysis.json'),plot_sha256=run.sha256(Path(__file__)),scope='Saved local linear quarter roots, conditional MC intervals and existing leading size fits. Red diagnostic fits are not accepted critical couplings; figure is not a new uncertainty analysis.'))
print(folder,flush=True)
