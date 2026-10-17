"""Joint Binder-contour/magnetization analysis at exactly beta=0 in 3D."""
from pathlib import Path
import json, argparse, shutil
import numpy as np
import measure as m


def endpoint_contour(q, means, draws, target=.25, method='moment_linear', window='local'):
    # U decreases with epsilon. x=-epsilon/pi increases toward the ordered phase.
    order=np.argsort(-np.asarray(q)); x=-np.asarray(q)[order]
    r,b=m.interpolate_record(x,means[order],draws[:,order],target,method,window)
    r['scan_coordinate']='minus_epsilon_over_pi'
    if r.get('location') is not None:
        r['epsilon_over_pi']=-r['location'][0]
        r['epsilon_MC95']=[-r['MC95'][0][1],-r['MC95'][0][0]]
    r['selected_original_indices']=[int(order[i]) for i in r.get('selected_indices',[])]
    return r,b


def analyze(tag='primary',cut=1001,factor=1,target=.25,starts='all',window='local'):
    out=m.OUT/'zero-beta';out.mkdir(exist_ok=True)
    folder=out/(tag+'-'+m.run.now()[:19].replace('-','').replace(':','').replace('T','-'));folder.mkdir()
    shutil.copyfile(__file__,folder/'zero_beta.py');shutil.copyfile(m.OUT/'measure.py',folder/'measure.py')
    groups={};sources=[]
    for f in sorted((m.CAMPAIGN/'jobs').glob('*.json')):
        j=json.loads(f.read_text());p=j['point']
        if p['D']!=3 or p['beta']!=0 or j['state']!='complete' or not j.get('audit',{}).get('passed'):continue
        if starts!='all' and p['replica']%2!=int(starts=='odd'):continue
        q=float(p.get('epsilon',np.pi*p['epsilon_over_pi']))/np.pi
        groups.setdefault((p['N'],q),[]).append(p);sources.append(dict(id=p['id'],point=p,audit=j['audit']))
    data=[];boots={}
    for (N,q),points in sorted(groups.items()):
        s,b=m.group_summary(points,cut,factor);data.append(dict(N=N,q=q,summary=s))
        if b is not None:boots[N,q]=b
    contours=[]
    for N in sorted({v['N'] for v in data}):
        gs=[v for v in data if v['N']==N and (N,v['q']) in boots]
        if len(gs)<2:continue
        q=np.array([v['q'] for v in gs]);means=np.array([v['summary']['mean'] for v in gs]);B=np.stack([boots[N,v['q']] for v in gs],axis=1)
        for method in ['linear','pchip','moment_linear']:
            r,b=endpoint_contour(q,means,B,target,method,window)
            r.update(N=N,method=method,flags=sorted({f for i in r.get('selected_original_indices',[]) for f in gs[i]['summary']['flags']}))
            contours.append(r)
            if b is not None:np.savez_compressed(folder/f'root_N{N}_{method}.npz',samples=b)
    fits=[]
    for method in ['linear','pchip','moment_linear']:
        rs=sorted([r for r in contours if r['method']==method and r['status']=='conditional_contour'],key=lambda r:r['N'])
        if not rs:continue
        B=np.stack([np.load(folder/f'root_N{r["N"]}_{method}.npz')['samples'] for r in rs],axis=1)
        for col,moment,power in [(1,'m2',2),(2,'mabs',1)]:
            for minimum,omega in [(6,None),(8,None),(12,None),(16,None),(6,.5),(6,.8),(6,1.)]:
                f=m.fit_dimension([r['N'] for r in rs],[r['location'][col] for r in rs],B[:,:,col],power,minimum,omega)
                f.pop('delta_samples',None)
                fits.append(dict(method=method,moment=moment,minimum=minimum,omega=omega,result=f,input_flags=sorted({flag for r in rs if r['N']>=minimum for flag in r['flags']}),omitted_sizes=[r['N'] for r in contours if r['method']==method and r['N']>=minimum and r['status']!='conditional_contour']))
    m.run.atomic_json(folder/'analysis.json',dict(created=m.run.now(),settings=dict(tag=tag,beta=0,D=3,cut=cut,factor=factor,target=target,starts=starts,window=window),source_jobs=sources,groups=data,contours=contours,fits=fits,scope='Conditional existing-MC effective spin dimension at beta=0. No nu/Delta fixed; no new simulations. Size/contour/correction and mixing sensitivity are not total confidence coverage.'))
    m.run.atomic_json(out/f'latest-{tag}.json',dict(folder=folder.name));print(folder,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--tag',default='primary');p.add_argument('--cut',type=int,default=1001);p.add_argument('--factor',type=int,default=1);p.add_argument('--target',type=float,default=.25);p.add_argument('--starts',choices=['all','even','odd'],default='all');p.add_argument('--window',choices=['local','full'],default='local');a=p.parse_args();analyze(a.tag,a.cut,a.factor,a.target,a.starts,a.window)
