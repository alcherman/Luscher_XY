"""Attained-source summary and exact energy-event screen, without pooling starts."""
from pathlib import Path
import json,math,importlib.util
import numpy as np,h5py
P=Path(__file__).resolve().parent;S=P.parent/'scaled-asymptote'
s=importlib.util.spec_from_file_location('energy',S/'energy_sector_check.py');energy=importlib.util.module_from_spec(s);s.loader.exec_module(energy)
def start(d):return int(d['id'].split('_r')[1].split('_')[0])%10
def peak(d):return float(np.mean([x['power'] for x in d['spectrum_samples']]))
if __name__=='__main__':
 sources=[S/'scaled-N6-analysis.json',S/'size-validation-analysis.json',P/'analysis.json'];rows=[];events=[]
 for source in sources:
  a=json.loads(source.read_text())
  for d in a['rows']:
   p=d['point'];rows.append(d);eta=1/(12*abs(p['target_kappa']));c=math.cos(p['epsilon']-4*eta)+(math.log(math.pi/eta)+1)/(3*abs(p['target_kappa']));bound=p['N']**3*energy.log_bound_per_site(3,p['target_kappa'],p['epsilon'],eta,c)
   with h5py.File(source.parent/'chains'/f"{d['id']}.h5") as h:x=h['mean_link_cosine'][p['warmup']+1:,0]
   events.append(dict(id=d['id'],N=p['N'],kappa=p['target_kappa'],q=p['epsilon_over_pi'],initial_label=p['initial_label'],threshold=c,eta=eta,log_probability_bound=bound,frames_above=int((x>=c).sum()),frames=len(x),source_analysis=str(source),original_audit_sha256=d['full_file_sha256']))
 summary=[]
 for k in [-4,-8,-16,-32,-64]:
  for u in [.1,.15,.4,1.,1.6]:
   rr=[d for d in rows if d['point']['N']==8 and d['point']['target_kappa']==k and abs(d['point']['scaled_u']-u)<1e-12 and d['point']['role']!='off_scaling_energy_sector_control'];regular=[d for d in rr if start(d)!=1];unit=[d for d in rr if start(d)==1]
   summary.append(dict(kappa=k,u=u,q=.5+u/abs(k)/np.pi,three_start_vector_range=[min(map(peak,regular)),max(map(peak,regular))],unit_vector=peak(unit[0]),three_start_nematic_range=[min(d['mean']['staggered_nematic_squared'] for d in regular),max(d['mean']['staggered_nematic_squared'] for d in regular)],three_start_vortex_range=[min(d['mean']['vortex_density'] for d in regular),max(d['mean']['vortex_density'] for d in regular)]))
 result=dict(source_analyses=[str(x) for x in sources],trajectories=len(rows),low_block_count=sum('low_block_count' in d['flags'] for d in rows),insufficient_blocks=sum('insufficient_blocks' in d['flags'] for d in rows),generation_plus_audit_worker_seconds=sum(d['generation_seconds']+d['audit_seconds'] for d in rows),actual_bytes=sum(d['actual_bytes'] for d in rows),N8_scaled_comparison=summary,energy_events=events,scope='Finite-box response alignment and nonequilibrium diagnostics. Parameter-only energy-event threshold has equilibrium probability at most exp(-N^3). Event frame frequencies are descriptive correlated samples, not independent binomial events; complementary sectors are not certified equilibrated. No accepted thermodynamic roots, precision intervals or fitted asymptote.')
 (P/'asymptote-assessment.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['trajectories','low_block_count','insufficient_blocks','generation_plus_audit_worker_seconds','actual_bytes']}))
