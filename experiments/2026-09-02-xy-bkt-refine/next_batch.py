"""Freeze the next disjoint adaptive batch from fully generated/analyzed slices."""
from run import OUT,sha256
from design import freeze
from supplement import propose
import argparse
import json


def main():
    parser=argparse.ArgumentParser();parser.add_argument('name');parser.add_argument('--freeze',action='store_true');args=parser.parse_args()
    manifests=[json.loads(p.read_text()) for p in (OUT/'batches').glob('*.json')]
    planned=[p for b in manifests for p in b['points']];owned={p['id'] for p in planned}
    pooled=json.loads((OUT/'pooled.json').read_text());analyzed={name for r in pooled for name in r['chains']}
    mature=set()
    for e in {p['epsilon_milli'] for p in planned if p['replica']<100}:
        positive=[p for p in planned if p['epsilon_milli']==e and p['kappa']>0 and p['replica']<100]
        if positive and all(p['id'] in analyzed and json.loads((OUT/'jobs'/f"{p['id']}.json").read_text())['state']=='complete' for p in positive):mature.add(e)
    points,decisions=propose();points=[p for p in points if p['epsilon_milli'] in mature and p['id'] not in owned]
    eps={p['epsilon_milli'] for p in points};decisions=[r for r in decisions if r['epsilon_milli'] in eps]
    assert len(points)%12==0 and len({p['id'] for p in points})==len(points)
    print(json.dumps(dict(mature_epsilon=sorted(mature),parameter_points=len(points)//12,chains=len(points),decisions=decisions),indent=2),flush=True)
    if args.freeze:
        if not points:raise RuntimeError('No mature additions proposed')
        snapshot=OUT/'batches'/f'{args.name}-pooled-input.txt';snapshot.write_bytes((OUT/'pooled.json').read_bytes())
        freeze(args.name,points,dict(method=__doc__,decisions=decisions,input_sha256=sha256(snapshot),
            selection_source_sha256=sha256(OUT/'supplement.py'),concurrent_runner_sha256=sha256(OUT/'run_edges.py')))
        (OUT/'batches'/f'{args.name}-design.py').write_bytes((OUT/'design.py').read_bytes())
        (OUT/'batches'/f'{args.name}-selection.py').write_bytes((OUT/'supplement.py').read_bytes())
        (OUT/'batches'/f'{args.name}-maturity.py').write_bytes((OUT/'next_batch.py').read_bytes())

if __name__=='__main__':main()
