"""Descriptive finite-box response bands; not equilibrium confidence intervals."""
from pathlib import Path
import json,hashlib
import numpy as np
P=Path(__file__).resolve().parent
s=P/'analysis.json';a=json.loads(s.read_text());groups={}
for d in a['rows']:
 if d['point'].get('role') in ['production','fine_transport_validation']:groups.setdefault((d['point']['N'],d['point']['epsilon_over_pi'],tuple(d['point']['kappas'])),[]).append(d)
bands=[]
for (n,q,kt),rr in sorted(groups.items()):
 if len(rr)!=2:continue
 ks=np.asarray(rr[0]['point']['kappas']);u=np.array([half for r in rr for half in r['half_U']]);m=np.array([half for r in rr for half in r['half_m2']]);mean=m.mean(0)
 stable=np.ptp(m,axis=0)<=.1*mean+.002
 zero_dominant=np.array([all(r['spectrum_peaks'][i]==[0,0,0] for r in rr) for i in range(len(ks))])
 ordered=stable & (u.min(0)>.45) & (m.min(0)>.10) & zero_dominant
 low=stable & (u.max(0)<.15) & (m.max(0)<.03)
 # First measured low-uniform-moment anchor on descent from kappa=0,
 # with a higher-kappa stable ordered anchor. No fitted central coupling.
 pairs=[(i,j) for i in np.where(low)[0] for j in np.where(ordered)[0] if j>i]
 if not pairs:
  bands.append(dict(role=rr[0]['point']['role'],N=n,epsilon_over_pi=q,resolved_response_band=False,reason='No common stable four-half response anchors'));continue
 i=max(v[0] for v in pairs);j=min(v[1] for v in pairs if v[0]==i)
 bands.append(dict(role=rr[0]['point']['role'],N=n,epsilon_over_pi=q,resolved_response_band=True,kappa_low=float(ks[i]),kappa_high=float(ks[j]),anchor_indices=[int(i),int(j)],anchor_half_m2=m[:,[i,j]].tolist(),anchor_half_U=u[:,[i,j]].tolist(),independent_ladders=2,postcut_cycles_per_ladder=rr[0]['postcut_cycles'],source_ids=[r['id'] for r in rr],global_sampling_flags=sorted(set(f for r in rr for f in r['flags'])),round_trips=[r['walker_round_trips']['round_trips'] for r in rr],nonzero_warning_fractions=[[r['coherent_nonzero_warning_fraction'][k] for k in [i,j]] for r in rr]))
result=dict(bands=bands,definition='N8 observed response band. Equilibrium transition location unresolved. NOT a critical-coupling estimate or confidence interval. Four half-series (two independent ladders) must have m2 range <=0.1*mean+0.002; U<0.15 and m2<0.03 at lower anchor; U>0.45, m2>0.10 and dominant averaged zero mode at upper anchor. Thresholds selected descriptively after inspecting q0.65; do not certify stationarity or absence of shared metastability.',analysis_sha256=hashlib.sha256(s.read_bytes()).hexdigest(),analysis_source=str(s))
(P/'response-bands.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(bands,indent=2))
