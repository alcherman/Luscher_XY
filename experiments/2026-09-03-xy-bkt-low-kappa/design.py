"""Freeze a bounded low-kappa batch before simulation; preserve all raw records."""
import json
from pathlib import Path
from run import OUT,point,atomic_json,sha256,now,_provenance


def freeze(name,pairs,steps=32000):
    points=[point(e,k,L,rep,steps) for e,k in pairs for L in [16,32,64] for rep in range(4)]
    path=OUT/'batches'/f'{name}.json';path.parent.mkdir(exist_ok=True)
    if path.exists():raise FileExistsError(path)
    for previous in path.parent.glob('*.json'):
        owned={p['id'] for p in json.loads(previous.read_text())['points']}
        assert not owned.intersection(p['id'] for p in points)
    record=dict(schema='xy-low-kappa.v1',batch=name,created_utc=now(),D=2,points=points,
        runner_sha256=sha256(OUT/'run.py'),model_source_sha256=_provenance()['source_sha256'],
        design_sha256=sha256(Path(__file__)),selection=dict(
            method='Direct edge refinement selected from previous complete U/R/eta observations; fixed couplings, no interpolated phase labels.',
            pairs=pairs,source_sha256=sha256(OUT.parent/'2026-09-02-xy-bkt-refine/conservative-brackets.json'),
            current_pooled_sha256=sha256(OUT/'pooled.json') if (OUT/'pooled.json').exists() else None,
            target_vertical_halfwidth=.05,analysis_cut_provisional=2001),
        update='1 Metropolis sweep + 4 fixed Wolff clusters + 1 winding proposal per cycle',
        retention='Initial field and every cycle field, measurement, RNG state, cluster work and timing; no thinning or deletion',
        raw_bytes_upper_estimate=sum((p['steps']+1)*(8*p['L']**2+2600) for p in points))
    atomic_json(path,record)
    (path.parent/f'{name}-design.py').write_bytes(Path(__file__).read_bytes())
    print(name,len(points),'chains',record['raw_bytes_upper_estimate']/2**30,'GiB upper estimate')


if __name__=='__main__':
    freeze('initial',[(575,40),(575,140),(580,100),(580,180),(590,180),(590,260)])
