"""Render saved rejected calibration fits, verifying their residual statistics."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
os.environ['MPLBACKEND']='Agg'
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT=Path(__file__).resolve().parent


def predict(beta,N,origin,D,parameters):
    beta=np.asarray(beta)
    thermal=N**(1/.67169) if D==3 else N*N*np.log(N)**.1
    correction=N**(-.789) if D==3 else np.log(N)**(-.5)
    c=parameters
    result=c[0]+c[1]*(beta-origin)*thermal+c[2]*thermal+c[3]*correction
    if len(c)==5:
        result+=c[4]*(N**(-2*.789) if D==3 else np.log(np.log(N))/np.log(N))
    return result


def main():
    target=ROOT/'rejected-fit-overlays-2026-10-13'
    target.mkdir(exist_ok=True)
    source=ROOT/'calibration-analysis/diagnostics.json'
    raw=source.read_bytes();data=json.loads(raw)
    # Freeze the measurements and coefficients that these figures display.
    snapshot=target/'diagnostics-snapshot.json'
    if snapshot.exists() and snapshot.read_bytes()!=raw:
        raise RuntimeError('Existing plot snapshot differs: choose a new output directory')
    shutil.copyfile(source,snapshot)
    checks=[]
    plt.rcParams.update({'font.size':11,'axes.titlesize':12,'axes.labelsize':12})
    with PdfPages(target/'rejected-calibration-fit-overlays.pdf') as pdf:
        for D in [3,4]:
            allrows=[r for r in data['groups'] if r['D']==D]
            sizes=sorted({r['N'] for r in allrows})
            colors={N:plt.get_cmap('tab10')(i) for i,N in enumerate(sizes)}
            index=0 if D==3 else 2
            for minimum in [6,8]:
                fig,axes=plt.subplots(2,2,figsize=(12,7.5),sharex='col',sharey='row',
                                     gridspec_kw={'height_ratios':[2.3,1]})
                for column,subleading in enumerate([False,True]):
                    entry=next(f for f in data['fits'] if f['D']==D and f['minimum_N']==minimum and f['subleading']==subleading)
                    result=entry['result'];top,bottom=axes[:,column]
                    selected=[r for r in allrows if r['N']>=minimum]
                    origin=np.mean([r['beta'] for r in selected])
                    total=0.
                    identifiable=isinstance(result.get('parameters'),list)
                    for N in sizes:
                        rows=sorted([r for r in selected if r['N']==N],key=lambda r:r['beta'])
                        if not rows:continue
                        beta=np.array([r['beta'] for r in rows])
                        values=np.array([r['summary']['values'][index] for r in rows])
                        sigma=np.array([r['summary']['errors'][index] for r in rows])
                        top.errorbar(beta,values,yerr=sigma,fmt='o',ms=5,color=colors[N],
                                     capsize=3,lw=1,label=f'N={N}',zorder=3)
                        flagged=np.array([bool(r['summary']['flags']) for r in rows])
                        if flagged.any():
                            top.scatter(beta[flagged],values[flagged],s=90,facecolors='none',
                                        edgecolors='black',linewidths=.9,zorder=4)
                        if identifiable:
                            grid=np.linspace(beta.min(),beta.max(),151)
                            top.plot(grid,predict(grid,N,origin,D,result['parameters']),
                                     color=colors[N],lw=1.6,ls='--')
                            residual=(values-predict(beta,N,origin,D,result['parameters']))/sigma
                            bottom.plot(beta,residual,'o',ms=5,color=colors[N])
                            total+=float(residual@residual)
                    if identifiable:
                        np.testing.assert_allclose(total,result['chi2'],rtol=1e-10,atol=1e-9)
                        checks.append(dict(D=D,minimum_N=minimum,subleading=subleading,
                                           chi2_recomputed=total,chi2_saved=result['chi2'],df=result['df']))
                        subtitle=rf"REJECTED: $\chi^2/\mathrm{{dof}}={total:.1f}/{result['df']}$"
                    else:
                        subtitle='Unidentifiable: no unique fit to overlay'
                        top.text(.5,.12,f"Rank {result['rank']} for {result['parameters']} parameters",
                                 ha='center',transform=top.transAxes,fontsize=10)
                    top.set_title(('Leading correction' if not subleading else 'Leading + subleading correction')+'\n'+subtitle)
                    top.legend(ncol=2,fontsize=9,loc='best')
                    top.grid(alpha=.2);bottom.grid(alpha=.2)
                    bottom.axhspan(-2,2,color='gray',alpha=.12)
                    bottom.axhline(0,color='black',lw=.8)
                    bottom.set_xlabel(r'$\beta$')
                    if not identifiable:
                        bottom.text(.5,.5,'Residuals undefined',ha='center',transform=bottom.transAxes)
                axes[0,0].set_ylabel(r'$U=1-\langle m^4\rangle/(2\langle m^2\rangle^2)$' if D==3 else r'$Q=\langle m^2\rangle^2/\langle m^4\rangle$')
                axes[1,0].set_ylabel(r'$(\mathrm{data}-\mathrm{fit})/\sigma$')
                fig.suptitle(rf'Ordinary XY calibration: $d={D}$, $\epsilon/\pi=1$, $N\geq {minimum}$',fontsize=15)
                fig.text(.5,.035,'Dots: paired-start measurements; vertical bars: 1 bootstrap SD. Dashed curves: saved joint fits.\n'
                         'Black rings: sampling-diagnostic flags. Residual shading: ±2 SD. These fits give no accepted critical coupling.',
                         ha='center',fontsize=9)
                fig.tight_layout(rect=[0,.09,1,.94])
                stem=target/f'd{D}-rejected-fits-minN{minimum}'
                fig.savefig(str(stem)+'.png',dpi=170)
                fig.savefig(str(stem)+'.svg')
                fig.savefig(str(stem)+'.pdf')
                pdf.savefig(fig)
                plt.close(fig)
    (target/'provenance-and-checks.json').write_text(json.dumps(dict(
        diagnostics_created=data['created'],diagnostics_sha256=hashlib.sha256(raw).hexdigest(),
        measurements='54 initial endpoint trajectories; 27 paired-start groups. Newly completed narrower-window trajectories are not included.',
        error_bars='One standard deviation of stratified contiguous-block bootstrap; yerr only.',
        fit_scope='Saved rejected fits; curves reconstructed from saved coefficients. No new fits or simulation.',
        chi_square_checks=checks),indent=2)+'\n')
    print('verified',len(checks),'chi-square reconstructions; plots in',target)


if __name__=='__main__':main()
