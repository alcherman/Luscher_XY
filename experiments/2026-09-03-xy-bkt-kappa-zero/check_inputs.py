"""Record immutable analysis inputs and reproduce the preceding phase labels."""
import json
from pathlib import Path
from run import OUT,atomic_json,sha256
from conservative import data,classify


def main():
    files=[]
    for name in ['2026-09-02-xy-bkt-scout','2026-09-02-xy-bkt-refine','2026-09-03-xy-bkt-low-kappa']:
        folder=OUT.parent/name
        result=folder/('results.json' if name.endswith('scout') else 'pooled.json')
        files.append(result)
        for r in json.loads(result.read_text()):
            files.append(folder/'analysis'/f"{r['point']['id']}_bootstrap.npz" if name.endswith('scout') else folder/'pooled'/f"{r['id']}_bootstrap.npz")
        if (folder/'conservative-brackets.json').exists():files.append(folder/'conservative-brackets.json')
    target=OUT/'input-provenance.json'
    digests={str(p):sha256(p) for p in files}
    if target.exists():assert json.loads(target.read_text())==digests
    else:atomic_json(target,digests)
    lookup=data()
    previous=json.loads((OUT.parent/'2026-09-03-xy-bkt-low-kappa/conservative-brackets.json').read_text())['labels']
    for r in previous:
        e,k=r['epsilon_milli'],r['kappa_milli']
        assert classify({L:lookup[(e,k,L)] for L in [16,32,64]})=={key:v for key,v in r.items() if key not in ['epsilon_milli','kappa_milli']},(e,k)
    atomic_json(OUT/'classifier-parity.json',dict(passed=True,previous_labels=len(previous),classifier_sha256=sha256(OUT/'conservative.py')))
    print('Passed',len(previous),'prior classification comparisons;',len(files),'input hashes')


if __name__=='__main__':main()
