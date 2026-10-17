"""Zero-beta epsilon scans, with independent namespace and endpoint role."""
import os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/luscher-mpl')
import json,numpy as np
from pathlib import Path
import run
OUT=Path(__file__).resolve().parent
choices={2:dict(epsilon_milli=[558,565,572],sizes=[64,128,256],steps=12000),3:dict(epsilon_milli=[806,812],sizes=[12,16,24,32],steps=8000),4:dict(epsilon_milli=[864,868,872],sizes=[10,12,16,20,24],steps=4000)}
points=[]
for D,c in choices.items():
 for e in c['epsilon_milli']:
  for N in c['sizes']:
   for replica in [4,5]:
    v=run.point(D,e//10,0,N,replica,c['steps']);v.update(id=f'zero_d{D}_e{e:04d}_N{N:03d}_r{replica:02d}_s{c["steps"]:06d}',epsilon_over_pi=e/1000,epsilon=e*np.pi/1000,seed=int(np.random.SeedSequence([2026101403,D,e,N,replica]).generate_state(1,dtype=np.uint64)[0]),campaign_role='zero_beta_endpoint',epsilon_sampling_milli=e,primary_observables=['eta_eff','xi_over_L'] if D==2 else ['binder_quarter'])
    points.append(v)
existing=[v for f in (OUT/'batches').glob('*.json') for v in json.loads(f.read_text())['points']];assert not {v['id'] for v in points}&{v['id'] for v in existing};assert not {v['seed'] for v in points}&{v['seed'] for v in existing}
m=run.freeze('zero-beta-endpoint-pilot-20261014',points,'User-requested beta=0 transitionlocations in epsilon/pi forD2/D3/D4 with horizontal uncertaintybars. D2useseta_eff andxi/L with matchingepsilon acrosssizes. D3/D4usefixedBinderU=.25 anddimensionappropriateFSS. Reuseoldzerocouplingdata; fresh epsilon scans nearknownendpointbrackets. Dedicated endpointrole prevents off-grid epsilon points frombeingpooled into fixed-epsilonpositivebeta curves. Independentseedsnamespace2026101403,replicas4/5,completefieldsRNGcheckpointslogsandall-fieldaudits. All316oldD2Binder-quarter jobsremain deferred;D4Nmax24,no64^4.')
run.atomic_json(OUT/'zero-beta-endpoint-pilot-20261014-design.json',dict(created=run.now(),choices=choices,trajectories=len(points),estimated_GiB=m['estimated_bytes']/2**30,source_sha256=run.sha256(Path(__file__)),plot_requirement='Atbeta0 plot inferred infinite-volume epsilon/pi location with horizontalMC andsensitivity bars only; beta coordinateexact0 has noverticalbar. No fixed-epsilon-grid horizontalbars elsewhere.',warmup_scope='Provisional labels; inspectcuts/blocks/starts andfinite-size sensitivity,reportwiderorunresolvedendpointifneeded.'))
print('zero_beta_pilot',len(points),'estimated_GiB',m['estimated_bytes']/2**30,flush=True)
for D in choices:print('D',D,'trajectories',sum(v['D']==D for v in points),flush=True)
