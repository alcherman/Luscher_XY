"""Overlay four independently fitted xi/L curves on one common set of axes."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:os.environ[key]='1'
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json
import shutil
from pathlib import Path
import numpy as np
from matplotlib.ticker import FormatStrFormatter
import local_xi_response as local

x=local.x;m=x.m;plt=m.plt


def main():
    source=m.OUT/'local-xi-response-20261015-213702';a=json.loads((source/'analysis.json').read_text());rawpath=Path(a['sources']['primary']['path'])
    assert m.run.sha256(rawpath)==a['sources']['primary']['sha256']
    raw=json.loads(rawpath.read_text());points={v['id']:v['point'] for v in raw['source_jobs']}
    folder=m.OUT/('four-size-xi-crossing-'+m.run.now()[:19].replace('-','').replace(':','').replace('T','-'));folder.mkdir();shutil.copyfile(__file__,folder/Path(__file__).name)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(figsize=(9.6,6));data=[];knots=None
    for N,color,marker in [(24,'#1768ac','o'),(32,'#17835b','s'),(48,'#ca6416','^'),(64,'#954aa3','D')]:
        pair=24 if N in [24,48] else 32
        record=next(r for r in a['records'] if r['variant']=='primary' and r['N']==pair and r['window']=='broad' and r['degree']==2)
        q=np.array(record['knots']);lo,hi=record['range']
        if knots is None:knots=q
        else:np.testing.assert_allclose(q,knots,atol=1e-14,rtol=0)
        mu=[];draws=[]
        for knot in q:
            groups=[g for g in raw['groups'] if g['N']==N and abs(g['q']-knot)<1e-14];assert len(groups)==1
            s,B=x.summary([points[i] for i in groups[0]['summary']['ids']],raw['settings']['cut'],raw['settings']['factor']);mu.append(s['mean']);draws.append(B)
        mu=np.array(mu);B=np.stack(draws,axis=1);fit,diagnostics=local.response_fit(q,mu,B,N,lo,hi,2)
        ref=next(r for r in record['responses'] if r['N']==N)
        np.testing.assert_allclose(fit[0][0],ref['response_diagnostics'][0]['coefficients'],atol=1e-12,rtol=0)
        R=x.ratio(mu,N);errors=np.quantile(x.ratio(B,N),[.025,.975],axis=0)
        dense=np.linspace(lo,hi,401);V=np.polynomial.polynomial.polyvander((dense-(lo+hi)/2)/((hi-lo)/2),2);curve=V@fit[0][0];band=np.quantile(fit[1][:,0]@V.T,[.025,.975],axis=0)
        ax.fill_between(dense,band[0],band[1],color=color,alpha=.09,zorder=1)
        ax.plot(dense,curve,color=color,lw=1.7,zorder=2)
        ax.vlines(q,errors[0],errors[1],color=color,lw=1.1,zorder=3)
        ax.plot(q,R,marker,ms=5.2,color=color,label=rf'$L={N}$',zorder=4)
        data.append(dict(N=N,q=q.tolist(),R=R.tolist(),MC95=errors.T.tolist(),fit_coefficients=fit[0][0].tolist(),response_diagnostics=diagnostics,dense_q=dense.tolist(),fitted_R=curve.tolist(),band_MC95=band.T.tolist()))
    ax.set(title=r'3D, $\beta=0$: correlation-length ratio for all four sizes',xlabel=r'$\epsilon/\pi$',ylabel=r'$\xi/L$')
    ax.set_xticks([.8005,.80075,.801,.80125,.8015]);ax.xaxis.set_major_formatter(FormatStrFormatter('%.5f'));ax.grid(alpha=.16)
    ax.legend(title='Measured MC points',loc='lower left',ncol=2,fontsize=11,title_fontsize=10,framealpha=.95)
    fig.text(.5,.025,'Four independent quadratic fits; no common crossing imposed.\nPoints and shaded bands: conditional MC 95%, excluding fit-model uncertainty.',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.095,1,.98))
    for ext in ['png','pdf','svg']:fig.savefig(folder/f'four-size-xi-crossing.{ext}',dpi=190,facecolor='white')
    plt.close(fig)
    m.run.atomic_json(folder/'plotted-data.json',dict(sizes=data,source_analysis=str(source/'analysis.json'),source_sha256=m.run.sha256(source/'analysis.json'),raw_analysis=str(rawpath),raw_sha256=m.run.sha256(rawpath),scope='Only a presentation change to reviewed curves. Identical physical knots/common axes; each size fitted independently; no common crossing, Delta or critical R imposed. Conditional MC intervals, not certified total confidence.'))
    m.run.atomic_json(m.OUT/'latest-four-size-xi-plot.json',dict(folder=folder.name,status='pending_visual_review'));print(folder,flush=True)


if __name__=='__main__':main()
