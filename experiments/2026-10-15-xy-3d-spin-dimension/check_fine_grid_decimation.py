"""Same-source grid sensitivity, including shared-draw difference covariance."""
import json
import numpy as np
import measure as m
import zero_beta as z

def main():
    root=m.OUT/'zero-beta';source=root/json.loads((root/'latest-fine-primary.json').read_text())['folder']/'analysis.json'
    a=json.loads(source.read_text());points={j['id']:j['point'] for j in a['source_jobs']};records=[]
    for parity in [0,1]:
        for method in ['linear','pchip','moment_linear']:
            roots=[];samples=[]
            for n in sorted({g['N'] for g in a['groups']}):
                gs=sorted([g for g in a['groups'] if g['N']==n],key=lambda g:g['q'])[parity::2]
                B=np.stack([m.group_summary([points[i] for i in g['summary']['ids']],1001,1)[1] for g in gs],axis=1)
                r,b=z.endpoint_contour(np.array([g['q'] for g in gs]),np.array([g['summary']['mean'] for g in gs]),B,.25,method,'full');r.update(N=n);roots.append(r)
                if r['status']=='conditional_contour':samples.append((r,b))
            fit=m.fit_dimension([r['N'] for r,b in samples],[r['location'][1] for r,b in samples],np.stack([b[:,1] for r,b in samples],axis=1),2,8)
            same_sizes=fit.get('sizes')==[8,12,16,24,32];paired=None
            if same_sizes and fit.get('valid_fraction')==1:
                rs=sorted([r for r in a['contours'] if r['method']==method and r['status']=='conditional_contour'],key=lambda r:r['N'])
                B=np.stack([np.load(source.parent/f'root_N{r["N"]}_{method}.npz')['samples'][:,1] for r in rs],axis=1)
                base=m.fit_dimension([r['N'] for r in rs],[r['location'][1] for r in rs],B,2,8)
                if base.get('valid_fraction')==1:
                    shifts=fit['delta_samples']-base['delta_samples']
                    paired=dict(central_shift=fit['Delta']-base['Delta'],conditional_shared_draw_MC95=np.quantile(shifts,[.025,.975]).tolist(),scope='Same coordinate bootstrap draws reused; captures correlated Monte Carlo variation, not an interpolation bias bound.')
            fit.pop('delta_samples',None);records.append(dict(parity=parity,method=method,roots=roots,fit=fit,same_complete_sizes=same_sizes,paired_same_size_difference=paired))
    path=m.OUT/'fine-grid-decimation-diagnostic.json'
    m.run.atomic_json(path,dict(source_sha256=m.run.sha256(source),scope='Same-source every-other-knot sensitivity; correlated, data-selected diagnostic, not independent bias estimate or total coverage guarantee. Missing-size variants are incomplete comparisons.',records=records))
    for r in records:print(r['parity'],r['method'],r['fit'].get('Delta'),r['paired_same_size_difference'],flush=True)

if __name__=='__main__':main()
