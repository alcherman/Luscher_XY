"""Freeze conservative field-only reduction eligibility from measured size brackets."""
import json, math, hashlib
import run
from pathlib import Path
from collections import defaultdict
OUT=Path(__file__).resolve().parent

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def key(D,q,N,endpoint):return (D,round(float(q),12),N,endpoint)
def protected_interval(grid,window):
    grid=sorted(set(grid));lo,hi=min(window),max(window)
    left=[x for x in grid if x<lo];right=[x for x in grid if x>hi]
    width=hi-lo
    return [min(lo-width,left[-1] if left else lo),max(hi+width,right[0] if right else hi)]

def material_flags(row):
    flags=list(row.get('flags',[]))
    for values in row.get('individual_flags',{}).values():flags.extend(values)
    return [f for f in flags if f!='fewer_than_50_blocks_per_start']

def classify(point,curves):
    endpoint=point['beta']==0
    k=key(point['D'],0 if endpoint else point['epsilon_over_pi'],point['N'],endpoint)
    matching=curves.get(k,[]);x=point['epsilon_over_pi'] if endpoint else point['beta']
    if not matching:return False,'No measured same-size bracket coverage'
    for c in matching:
        if not c['resolved']:return False,'At least one contributing observable/variant is unresolved'
        if not any(abs(x-g)<1e-10 for g in c['grid']):return False,'Coordinate absent from a contributing measured curve'
        if c['protected'][0]-1e-10<=x<=c['protected'][1]+1e-10:return False,'Within protected bracket, neighbor or width margin'
    return True,'Outside every contributing measured bracket and conservative margin'

def input_signature():
    names=['latest-d2-eta-correlation.json','latest-d2-current-chain-sensitivity.json','latest-d3-d4-delivery.json','latest-zero-beta-endpoints.json','latest-endpoint-sensitivity.json']
    pointers={n:digest(OUT/n) for n in names}
    jobs=[]
    for f in (OUT/'jobs').glob('*.json'):
        j=json.loads(f.read_text())
        if j['state']=='complete' and j.get('audit',{}).get('passed'):jobs.append((j['point']['id'],j['audit']['sha256']))
    return hashlib.sha256(json.dumps([pointers,sorted(jobs)],sort_keys=True).encode()).hexdigest()

def refresh():
    path=OUT/'outside-bracket-retention-design-20261014.json'
    if path.exists() and json.loads(path.read_text()).get('input_signature')==input_signature():return
    build()

def build():
    signature=input_signature()
    curves=defaultdict(list);sources=[]
    def load(folder,name='analysis.json'):
        path=Path(folder)/name;sources.append(dict(path=str(path.resolve()),sha256=digest(path)))
        return json.loads(path.read_text())
    def pointer(name):return OUT/json.loads((OUT/name).read_text())['folder']
    def add(D,q,N,endpoint,grid,window,resolved,label):
        curves[key(D,q,N,endpoint)].append(dict(grid=grid,protected=protected_interval(grid,window) if resolved else None,resolved=resolved,source=label))
    d2folders={str(pointer('latest-d2-eta-correlation.json'))}
    sens=load(pointer('latest-d2-current-chain-sensitivity.json'))
    d2folders.update(v['folder'] for v in sens['sources'].values())
    for folder in sorted(d2folders):
        a=load(folder);obs={(r['epsilon_over_pi'],r['observable'],r['N']):r for r in a['observations']}
        for r in a['locations']:
            o=obs[r['epsilon_over_pi'],r['observable'],r['N']]
            for N in ([r['N'],2*r['N']] if r['observable']=='eta_eff' else [r['N']]):
                add(2,r['epsilon_over_pi'],N,False,o['beta'],r['bootstrap_window_beta'],r['status']=='conditional_crossing' and r['valid_fraction']>=.99 and not material_flags(r),folder)
    delivered=load(pointer('latest-d3-d4-delivery.json'),'phase-boundary-analysis.json')
    for folder in sorted({v['source'] for v in delivered['variants']}):
        a=load(folder)
        for r in a['locations']:
            grid=[g['beta'] for g in a['groups'] if (g['D'],g['q'],g['N'])==(r['D'],r['q'],r['N'])]
            add(r['D'],r['q']/100,r['N'],False,grid,r.get('bracket',grid),r['status']=='conditional_crossing' and r.get('valid_fraction',0)>=.99 and not material_flags(r),folder)
    es=load(pointer('latest-endpoint-sensitivity.json'))
    endpointfolders=set(es['snapshots'].values());endpointfolders.add(str(pointer('latest-zero-beta-endpoints.json')))
    for folder in sorted(endpointfolders):
        a=load(folder)
        for r in a['locations']:
            Ns=[r['N'],2*r['N']] if r['observable']=='eta_eff' else [r['N']]
            grid=[g['epsilon_over_pi'] for g in a['groups'] if (g['D'],g['N'])==(r['D'],r['N'])]
            for N in Ns:add(r['D'],0,N,True,grid,r['measured_epsilon_over_pi'],r['status']=='conditional_crossing' and r.get('valid_fraction',0)>=.99 and not material_flags(r),folder)
    eligible={};counts=defaultdict(int);physical=0
    for f in (OUT/'jobs').glob('*.json'):
        j=json.loads(f.read_text())
        if j['state']!='complete' or not j.get('audit',{}).get('passed'):continue
        ok,reason=classify(j['point'],curves);counts[reason]+=1
        if ok:
            p=j['point'];eligible[p['id']]=dict(point=p,audit_sha256=j['audit']['sha256'],reason=reason,bracket_evidence=curves[key(p['D'],0 if p['beta']==0 else p['epsilon_over_pi'],p['N'],p['beta']==0)])
            physical+=(OUT/'chains'/f"{p['id']}.h5").stat().st_size
    result=dict(input_signature=signature,created=run.now(),version=1,policy='checkpoints_only_outside_brackets',scope='Current campaign complete audited job IDs only; no historical external archives; future jobs require a new evidence-based manifest.',rule='Every contributing size/observable and tested cut/block/threshold/interpolation variant must be resolved >=99%; protect union of local measured windows plus one neighbor each side and one window-width margin. Unknown coordinates, unresolved curves, and material start/half/pooling flags remain adaptive; low-block-count alone does not authorize discarding observable histories, which remain complete.',sources=sources,eligible=eligible,classification_counts=dict(counts),eligible_current_bytes=physical)
    if input_signature()!=signature:
        print('Bracket inputs changed during classification; preserve prior manifest until next pass',flush=True)
        return
    path=OUT/'outside-bracket-retention-design-20261014.json'
    archive=OUT/'outside-bracket-retention-manifests';archive.mkdir(exist_ok=True)
    if path.exists():
        old=json.loads(path.read_text());run.atomic_json(archive/(digest(path)+'.json'),old)
    run.atomic_json(path,result)
    run.atomic_json(archive/(digest(path)+'.json'),result)
    print(json.dumps(dict(eligible=len(eligible),counts=dict(counts),current_GiB=physical/2**30)))
if __name__=='__main__':build()
