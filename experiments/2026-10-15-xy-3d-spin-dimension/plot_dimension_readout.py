"""Read-off plots for the spin dimension; reviewed MC only, no Delta imposed."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[k]='1'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json
import shutil
from pathlib import Path
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FixedLocator,FormatStrFormatter
import local_xi_response as local

x=local.x;m=x.m;plt=m.plt


def main():
    source=m.OUT/'local-xi-response-20261015-213702';a=json.loads((source/'analysis.json').read_text());rawpath=Path(a['sources']['primary']['path']);assert m.run.sha256(rawpath)==a['sources']['primary']['sha256'];raw=json.loads(rawpath.read_text());points={v['id']:v['point'] for v in raw['source_jobs']}
    folder=m.OUT/('dimension-readout-plots-'+m.run.now()[:19].replace('-','').replace(':','').replace('T','-'));folder.mkdir();shutil.copyfile(__file__,folder/Path(__file__).name)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    def group(N,q):
        rows=[g for g in raw['groups'] if g['N']==N and abs(g['q']-q)<1e-14];assert len(rows)==1
        return x.summary([points[i] for i in rows[0]['summary']['ids']],raw['settings']['cut'],raw['settings']['factor'])
    figures=[];data={}
    def save(fig,name):
        for ext in ['png','pdf','svg']:fig.savefig(folder/f'{name}.{ext}',dpi=190,facecolor='white')
        figures.append(fig)
    # RMS magnetization at an existing, fixed representative coupling.
    N=np.array([24,32,48,64]);q=.80075;means=[];Bs=[]
    for size in N:
        s,B=group(int(size),q);means.append(s['mean'][0]);Bs.append(B[:,0])
    rms=np.sqrt(means);B=np.sqrt(np.column_stack(Bs));logB=np.log(B);C=np.cov(logB,rowvar=False);X=np.column_stack([np.ones(4),-np.log(N)])
    T=np.linalg.solve(X.T@np.linalg.solve(C,X),np.linalg.solve(C,X).T);coef=T@np.log(rms);coeff_draws=logB@T.T
    oracle=m.fit_dimension(N,rms,B,1,24);np.testing.assert_allclose(coef[1],oracle['Delta'],atol=1e-12,rtol=0);np.testing.assert_allclose(np.quantile(coeff_draws[:,1],[.025,.975]),oracle['MC95'],atol=1e-12,rtol=0)
    Delta=float(coef[1]);SD=float(coeff_draws[:,1].std(ddof=1));ci=np.quantile(B,[.025,.975],axis=0)
    dense=np.geomspace(23,67,301);V=np.column_stack([np.ones(len(dense)),-np.log(dense)]);fit=np.exp(V@coef);bounds=np.quantile(np.exp(coeff_draws@V.T),[.025,.975],axis=0)
    fig,ax=plt.subplots(figsize=(9.6,6));ax.set_xscale('log');ax.set_yscale('log')
    ax.fill_between(dense,bounds[0],bounds[1],color='#1768ac',alpha=.13);ax.plot(dense,fit,color='#1768ac',lw=1.8,label='Unconstrained power-law fit')
    ax.vlines(N,ci[0],ci[1],color='#222222',lw=1.4);ax.plot(N,rms,'o',color='#222222',ms=6,label='Measured RMS magnetization')
    ax.xaxis.set_major_locator(FixedLocator(N));ax.xaxis.set_major_formatter(FormatStrFormatter('%d'));ax.xaxis.set_minor_locator(FixedLocator([]));ax.yaxis.set_major_locator(FixedLocator([.12,.14,.16,.18,.20,.22]));ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'));ax.yaxis.set_minor_locator(FixedLocator([]))
    ax.set(xlabel=r'Box size $L$ (logarithmic axis)',ylabel=r'$m_{\rm rms}=\sqrt{\langle |m|^2\rangle}$ (logarithmic axis)',title=r'3D, $\beta=0$: read the dimension from the slope')
    ax.text(.97,.96,rf'$m_{{\rm rms}}\propto L^{{-\Delta_{{\rm eff}}}}$'+'\n'+rf'Log–log slope $=-{Delta:.3f}$'+'\n'+rf'$\Delta_{{\rm eff}}={Delta:.3f}\pm{SD:.3f}$ (stat SD)',transform=ax.transAxes,ha='right',va='top',fontsize=12,bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
    ax.legend(loc='lower left',fontsize=10);ax.grid(which='major',alpha=.15)
    fig.text(.5,.025,r'Fixed measured $\epsilon/\pi=0.80075$, near the $\xi/L$ crossings; four sizes fitted together.'+'\nBars/band: conditional MC 95%. Coupling-location and finite-size corrections are not included in this fit error.',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.105,1,.98));save(fig,'01-rms-magnetization-slope')
    oracle.pop('delta_samples',None);data['fixed_q_rms']=dict(q=q,N=N.tolist(),RMS=rms.tolist(),MC95=ci.T.tolist(),fit=oracle,Delta_SD=SD,log_amplitude=float(coef[0]),scope='Fixed representative measured q, not an estimate of epsilon_c. Effective finite-size slope; no certified thermodynamic Delta or total uncertainty.')
    # Crossing-style effective dimensions, evaluated at the SAME q in each pair.
    fig,ax=plt.subplots(figsize=(10.2,6));pairs=[]
    r24=next(r for r in a['records'] if r['variant']=='primary' and r['N']==24 and r['window']=='broad' and r['degree']==2);r32=next(r for r in a['records'] if r['variant']=='primary' and r['N']==32 and r['window']=='broad' and r['degree']==2)
    qlo=min(r24['MC95'][0][0],r32['MC95'][0][0]);qhi=max(r24['MC95'][0][1],r32['MC95'][0][1]);ax.axvspan(qlo,qhi,color='#777777',alpha=.12,label=r'$\xi/L$ pair-crossing region')
    for size,color,marker in [(12,'#1768ac','o'),(16,'#17835b','s'),(24,'#ca6416','^'),(32,'#954aa3','D')]:
        r=next(r for r in a['records'] if r['variant']=='primary' and r['N']==size and r['window']=='broad' and r['degree']==2);q=np.array(r['knots']);lo,hi=r['range'];mu=[];draws=[];fits=[]
        for n in [size,2*size]:
            ss=[];bb=[]
            for knot in q:
                s,B=group(n,knot);ss.append(s['mean']);bb.append(B)
            ss=np.array(ss);bb=np.stack(bb,axis=1);ff,diag=local.response_fit(q,ss,bb,n,lo,hi,2);ref=next(v for v in r['responses'] if v['N']==n);np.testing.assert_allclose(ff[0][1],ref['response_diagnostics'][1]['coefficients'],atol=1e-12,rtol=0);mu.append(ss[:,0]);draws.append(bb[:,:,0]);fits.append(ff)
        rawD=-np.log(mu[1]/mu[0])/(2*np.log(2));rawB=-np.log(draws[1]/draws[0])/(2*np.log(2));interval=np.quantile(rawB,[.025,.975],axis=0)
        dense=np.linspace(lo,hi,401);V=np.polynomial.polynomial.polyvander((dense-(lo+hi)/2)/((hi-lo)/2),2);m1=V@fits[0][0][1];m2=V@fits[1][0][1];b1=fits[0][1][:,1]@V.T;b2=fits[1][1][:,1]@V.T
        with np.errstate(invalid='ignore',divide='ignore'):curve=-np.log(m2/m1)/(2*np.log(2));boot=-np.log(b2/b1)/(2*np.log(2))
        valid=(b1>0)&(b2>0)&np.isfinite(boot);boot[~valid]=np.nan;fraction=valid.mean(0);bands=np.nanquantile(boot,[.025,.975],axis=0);curve[(m1<=0)|(m2<=0)|(fraction<.99)]=np.nan;bands[:,fraction<.99]=np.nan
        ax.fill_between(dense,bands[0],bands[1],color=color,alpha=.09);ax.plot(dense,curve,color=color,lw=1.7);ax.vlines(q,interval[0],interval[1],color=color,lw=1);ax.plot(q,rawD,marker,color=color,ms=5,label=f'Pair ({size}, {2*size})')
        pairs.append(dict(N=size,N2=2*size,q=q.tolist(),Delta_eff=rawD.tolist(),MC95=interval.T.tolist(),dense_q=dense.tolist(),fitted_Delta=curve.tolist(),fit_band_MC95=bands.T.tolist(),fitted_positive_fraction=fraction.tolist()))
    ax.set(title=r'Effective spin dimension at a common coupling in each size pair',xlabel=r'$\epsilon/\pi$',ylabel=r'$\Delta_{\rm eff}=-\log[\langle|m|^2\rangle_{2L}/\langle|m|^2\rangle_L]/(2\log2)$')
    ax.xaxis.set_major_locator(FixedLocator([.8005,.80075,.801,.80125,.8015]));ax.xaxis.set_major_formatter(FormatStrFormatter('%.5f'));ax.grid(alpha=.15);ax.legend(loc='upper left',fontsize=10)
    fig.text(.5,.025,'Points: direct measured ratios. Curves/bands: independent quadratic magnetization responses; conditional MC 95%.'+'\nGray strip: descriptive union of two finite-size crossing intervals, not an infinite-volume critical-coupling interval.\nFinite-size corrections remain; these effective-dimension curves need not cross at one exact point.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.125,1,.98));save(fig,'02-effective-dimension-curves')
    data['effective_dimension_curves']=dict(pairs=pairs,xi_crossing_strip=[qlo,qhi])
    with PdfPages(folder/'dimension-readout.pdf') as pdf:
        for fig in figures:pdf.savefig(fig)
    for fig in figures:plt.close(fig)
    m.run.atomic_json(folder/'plotted-data.json',dict(data=data,source_analysis=str(source/'analysis.json'),source_sha256=m.run.sha256(source/'analysis.json'),raw_analysis=str(rawpath),raw_sha256=m.run.sha256(rawpath),scope='Existing MC replot and illustrative fixed-q power-law fit only. No prescribed dimension, no common crossing imposed, no new MC, no calibrated total uncertainty.'))
    m.run.atomic_json(m.OUT/'latest-dimension-readout-plots.json',dict(folder=folder.name,status='pending_independent_review'));print(folder,flush=True);print(json.dumps(data['fixed_q_rms']['fit']),flush=True)


if __name__=='__main__':main()
