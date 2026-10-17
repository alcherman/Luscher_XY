"""Representative figures from reviewed 3D beta=0 MC; no new simulations."""
import csv
import json
import os
import shutil
from pathlib import Path
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[key]='1'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FormatStrFormatter
import local_xi_response as local

x=local.x; m=x.m; plt=m.plt
LOCAL=m.OUT/'local-xi-response-20261015-213702'
POINTWISE=m.OUT/'common-xi-refined-primary-20261015-213346'
BENCHMARK=.519050


def gls_line(records,minimum):
    rs=[r for r in records if r['variant']=='primary' and r['window']=='broad' and r['degree']==2 and r['N']>=minimum and r['status']=='conditional_local_quotient']
    N=np.array([r['N'] for r in rs]);y=np.array([r['Delta_eff'] for r in rs]);B=np.column_stack([np.load(LOCAL/f"primary-{r['N']}-{r['N2']}-broad-degree2.npz")['samples'][:,4] for r in rs]);valid=np.isfinite(B).all(1);B=B[valid]
    C=np.cov(B,rowvar=False);X=np.column_stack([np.ones(len(N)),N.astype(float)**(-.789)])
    # Full covariance includes overlapping size pairs; reproduce the reviewed fit.
    T=np.linalg.solve(X.T@np.linalg.solve(C,X),np.linalg.solve(C,X).T);coef=T@y;bs=B@T.T
    ledger=json.loads((LOCAL/'finite-size-diagnostics.json').read_text())['fits'];ref=next(r for r in ledger if r['variant']=='primary' and r['window']=='broad' and r['degree']==2 and r['minimum']==minimum)
    np.testing.assert_allclose(coef[0],ref['Delta'],atol=1e-11,rtol=0)
    np.testing.assert_allclose(bs[:,0].std(ddof=1),ref['SD'],atol=1e-11,rtol=0)
    return N,coef,bs,ref


def main():
    a=json.loads((LOCAL/'analysis.json').read_text());raw=json.loads((POINTWISE/'analysis.json').read_text());point={s['id']:s['point'] for s in raw['source_jobs']}
    folder=m.OUT/('representative-analysis-plots-'+m.run.now()[:19].replace('-','').replace(':','').replace('T','-'));folder.mkdir()
    shutil.copyfile(__file__,folder/Path(__file__).name)
    plt.rcParams.update({'font.size':10.5,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white'})
    figures=[];figure_data={};checks=[]
    def save(fig,name):
        for ext in ['png','pdf','svg']:fig.savefig(folder/f'{name}.{ext}',dpi=180,facecolor='white')
        figures.append((fig,name))
    # 1. Measured scalar histories and joint-bootstrap local response bands.
    fig,axes=plt.subplots(1,2,figsize=(12,4.9),sharey=False)
    for ax,N in zip(axes,[24,32]):
        record=next(r for r in a['records'] if r['variant']=='primary' and r['N']==N and r['window']=='broad' and r['degree']==2);q=np.array(record['knots']);lo,hi=record['range'];dense=np.linspace(lo,hi,301);V=np.polynomial.polynomial.polyvander((dense-(lo+hi)/2)/((hi-lo)/2),2)
        for size,color in [(N,'#1768ac'),(2*N,'#ca6416')]:
            means=[];draws=[]
            for knot in q:
                matches=[g for g in raw['groups'] if g['N']==size and abs(g['q']-knot)<1e-14];assert len(matches)==1
                group=matches[0];s,B=x.summary([point[i] for i in group['summary']['ids']],raw['settings']['cut'],raw['settings']['factor']);means.append(s['mean']);draws.append(B)
            means=np.array(means);B=np.stack(draws,axis=1);fit,diag=local.response_fit(q,means,B,size,lo,hi,2)
            ref=next(v for v in record['responses'] if v['N']==size)
            np.testing.assert_allclose(fit[0][0],ref['response_diagnostics'][0]['coefficients'],atol=1e-12,rtol=0)
            checks.append(f'R{size} coefficients reproduce reviewed pair {N}:{2*N}')
            R=x.ratio(means,size);bounds=np.quantile(x.ratio(B,size),[.025,.975],axis=0)
            ax.vlines(q,bounds[0],bounds[1],color=color,lw=1);ax.plot(q,R,'o',ms=4.5,color=color,label=f'MC, $L={size}$')
            center=V@fit[0][0];bands=np.quantile(fit[1][:,0]@V.T,[.025,.975],axis=0)
            ax.plot(dense,center,color=color,lw=1.7);ax.fill_between(dense,bands[0],bands[1],color=color,alpha=.13)
            figure_data[f'R-pair{N}-size{size}']=dict(q=q.tolist(),R=R.tolist(),MC95=bounds.T.tolist(),fit_coefficients=fit[0][0].tolist())
        qcross=record['epsilon_over_pi'];qlo,qhi=record['MC95'][0]
        ax.axvspan(qlo,qhi,color='#222222',alpha=.07);ax.axvline(qcross,color='#333333',ls=':',lw=1);ax.plot(qcross,record['xi_over_L'],'x',ms=7,color='black')
        ax.text(.97,.95,rf'$\epsilon_\times/\pi={qcross:.7f}$'+'\n2000/2000 crossing draws',transform=ax.transAxes,ha='right',va='top',fontsize=9,bbox=dict(facecolor='white',edgecolor='none',alpha=.85))
        ax.set(title=f'Pair $({N},{2*N})$',xlabel=r'$\epsilon/\pi$',ylabel=r'$\xi/L$');ax.xaxis.set_major_formatter(FormatStrFormatter('%.5f'));ax.set_xticks([.8005,.80075,.801,.80125,.8015]);ax.tick_params(axis='x',labelrotation=20);ax.legend(loc='lower left',fontsize=9);ax.grid(alpha=.17)
    fig.suptitle(r'3D, $\beta=0$: measured correlation-length crossings',fontsize=14)
    fig.text(.5,.015,'Points and response bands: conditional MC 95%. Curves: broad-window quadratic fits; shading is not total model uncertainty.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.065,1,.92));save(fig,'01-measured-xi-crossings')
    # 2. Correlated finite-size extrapolations, two explicit size cuts.
    fig,ax=plt.subplots(figsize=(10.5,5.2));rs=[r for r in a['records'] if r['variant']=='primary' and r['window']=='broad' and r['degree']==2 and r['status']=='conditional_local_quotient']
    for r in rs:
        xx=r['N']**(-.789);lo,hi=r['MC95'][4];ax.vlines(xx,lo,hi,color='#333333');ax.plot(xx,r['Delta_eff'],'o',color='#333333',ms=5);ax.annotate(f"{r['N']}:{r['N2']}",(xx,r['Delta_eff']),xytext=(5,7),textcoords='offset points',fontsize=9)
    for minimum,color,style in [(8,'#1768ac','-'),(12,'#9a479b','--')]:
        N,coef,bs,ref=gls_line(a['records'],minimum);dense=np.linspace(0,float(N.min()**(-.789)),251);V=np.column_stack([np.ones(len(dense)),dense]);bands=np.quantile(bs@V.T,[.025,.975],axis=0)
        ax.plot(dense,V@coef,color=color,ls=style,lw=1.7,label=f"Leading-correction fit, $L\geq{minimum}$: $\Delta_\infty={ref['Delta']:.4f}\pm{ref['SD']:.4f}$ (stat SD)")
        ax.fill_between(dense,bands[0],bands[1],color=color,alpha=.10);ci=np.quantile(bs[:,0],[.025,.975]);ax.vlines(0,*ci,color=color,lw=1.6);ax.plot(0,coef[0],'s',color=color,ms=5)
        figure_data[f'extrapolation-min{minimum}']=dict(reference=ref,coefficients=coef.tolist())
    ax.axhline(BENCHMARK,color='#777777',ls=':',lw=1.3,label=r'XY Monte Carlo benchmark: $0.519050$ (Hasenbusch, 2019)');ax.axvline(0,color='#aaaaaa',lw=.7)
    ax.set(xlabel=r'$L^{-\omega}$, $L$ is the smaller box in $(L,2L)$; $0$ is infinite volume',ylabel=r'$\Delta_{\rm eff}$',title=r'Finite-size extrapolation: shared-size covariance retained');ax.set_xlim(-.008,.218);ax.grid(alpha=.17);ax.legend(loc='upper right',fontsize=9)
    fig.text(.5,.015,r'Bars/bands: conditional MC 95%; fixed $\omega=0.789$ is WF-informed. Size cuts are alternatives, not independent measurements. $6{:}12$ is unresolved.',ha='center',fontsize=8.5)
    fig.tight_layout(rect=(0,.075,1,.97));save(fig,'02-correlated-size-extrapolation')
    # 3. All relevant accepted/unresolved/rejected sensitivity examples.
    extras=json.loads((m.OUT/'common-xi-refinement-status.json').read_text())['snapshots'];raw_variants={k:json.loads(Path(v['path']).read_text()) for k,v in extras.items()}
    rows=[('Broad quadratic; pooled','local','primary','broad',2),('Broad quadratic; longer blocks','local','block4','broad',2),('Broad quadratic; later cut','local','latecut','broad',2),('Broad quadratic; even starts','local','even','broad',2),('Broad quadratic; odd starts','local','odd','broad',2),('Inner linear','local','primary','inner',1),('Inner quadratic','local','primary','inner',2),('Pointwise linear; pooled','raw','primary','linear',0),('Pointwise PCHIP; pooled','raw','primary','pchip',0),('Pointwise linear; even starts','raw','even','linear',0),('Pointwise linear; odd starts','raw','odd','linear',0),('Broad linear; inadequate curvature','local','primary','broad',1)]
    fig,axes=plt.subplots(1,2,figsize=(13,7),sharey=True)
    sensitivity=[]
    for ax,N in zip(axes,[24,32]):
        for y,(label,kind,variant,window,degree) in enumerate(rows):
            if kind=='local':r=next(r for r in a['records'] if r['N']==N and r['variant']==variant and r['window']==window and r['degree']==degree);reject=r['response_status']=='rejected_local_response';good=r['status']=='conditional_local_quotient'
            else:r=next(r for r in raw_variants[variant]['pairs'] if r['N']==N and r['window']=='local' and r['method']==window);reject=False;good=r['status']=='conditional_quotient'
            sensitivity.append(dict(pair=[N,2*N],label=label,record=r))
            if 'Delta_eff' not in r or not r.get('MC95'):
                ax.text(.52,y,'no unique central root',ha='center',va='center',fontsize=8);continue
            color='#bf3830' if reject else '#1768ac' if good else '#666666';lo,hi=r['MC95'][4];ax.hlines(y,lo,hi,color=color,lw=1.2);ax.plot(r['Delta_eff'],y,'x' if reject else 'o',color=color,markerfacecolor=color if good else 'white',ms=5)
        ax.axvline(BENCHMARK,color='#999999',ls=':',lw=1);ax.set(title=f'Pair $({N},{2*N})$',xlabel=r'$\Delta_{\rm eff}$');ax.grid(axis='x',alpha=.17);ax.set_yticks(np.arange(len(rows)),[v[0] for v in rows])
    axes[0].invert_yaxis()
    axes[1].tick_params(labelleft=False)
    fig.suptitle('Sensitivity to response shape, fit window, warmup cut and independent starts',fontsize=13)
    fig.text(.5,.015,'Bars: conditional MC 95%. Open circles: unresolved roots (intervals use surviving draws); red crosses: rejected response shape. Rows reuse correlated data.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.055,1,.94));save(fig,'03-analysis-sensitivity')
    figure_data['sensitivity']=sensitivity
    with PdfPages(folder/'representative-analysis.pdf') as pdf:
        for fig,name in figures:pdf.savefig(fig)
    for fig,name in figures:plt.close(fig)
    m.run.atomic_json(folder/'plotted-data.json',figure_data)
    provenance=dict(created=m.run.now(),status='pending_independent_visual_review',local_source=str(LOCAL),local_source_sha256=m.run.sha256(LOCAL/'analysis.json'),pointwise_source=str(POINTWISE),pointwise_source_sha256=m.run.sha256(POINTWISE/'analysis.json'),checks=checks+['Both GLS intercepts/bootstrap SDs reproduced to1e-11 from full shared-pair covariance.'],benchmark=dict(value=BENCHMARK,MC_SD=.000040,source='https://arxiv.org/abs/1910.05916',purpose='Comparison only; no dimension fixed by fitting.'),scope='Replot reviewed MC analysis only. Conditional MC intervals; no certified total CI or numerical systematic error. Unresolved/rejected alternatives retained.')
    m.run.atomic_json(folder/'provenance.json',provenance);m.run.atomic_json(m.OUT/'latest-representative-analysis-plots.json',dict(folder=folder.name,status=provenance['status']))
    print(folder,flush=True)


if __name__=='__main__':main()
