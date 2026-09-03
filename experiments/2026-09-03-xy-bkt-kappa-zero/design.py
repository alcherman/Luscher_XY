"""Frozen small kappa=0 scan, retaining every generated field and measurement."""
import json
from pathlib import Path
from run import OUT,point,sha256,atomic_json,now,_provenance


def main():
    path=OUT/'batches/endpoint.json';path.parent.mkdir(exist_ok=True)
    if path.exists():raise FileExistsError(path)
    points=[point(e,0,L,rep,7000) for e in [550,570,575] for L in [16,32,64] for rep in [0,1]]
    record=dict(batch='endpoint',created_utc=now(),D=2,points=points,runner_sha256=sha256(OUT/'run.py'),
        model_source_sha256=_provenance()['source_sha256'],design_sha256=sha256(Path(__file__)),
        selection=dict(previous_endpoint=[.54,.58],epsilon_over_pi=[.55,.57,.575],kappa=0,
            method='Add a quasi-ordered-side candidate and two massive-side candidates around the existing unresolved epsilon/pi=0.565 point. Select phase labels only from direct measurements with the existing U/R/eta two-SE rule.',
            source_sha256=sha256(OUT.parent/'2026-09-03-xy-bkt-low-kappa/conservative-brackets.json')),
        update='1 full Metropolis sweep + 4 fixed Wolff clusters + 1 winding proposal per cycle',
        starts=[[0,0],[1,0]],replicas=2,analysis_cut_provisional=2001,
        retention='Initial field and every cycle field, measurement, RNG state, cluster work and timing. No thinning or deletion.',
        raw_bytes_upper_estimate=sum((p['steps']+1)*(8*p['L']**2+2600) for p in points))
    atomic_json(path,record);print(len(points),'chains;',record['raw_bytes_upper_estimate']/2**30,'GiB upper storage estimate')


if __name__=='__main__':main()
