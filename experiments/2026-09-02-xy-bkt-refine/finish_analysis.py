"""Analyze closed chains as they arrive, then build the complete deliverables."""
import json
from time import sleep
from run import OUT,now
from analyze import completed,analyze
from pool import main as pool


def main():
    manifests=sorted((OUT/'batches').glob('*.json'))
    expected={p['id'] for path in manifests for p in json.loads(path.read_text())['points']}
    last=None
    while True:
        selected=[p for p in completed() if p['replica']<100]
        ids={p['id'] for p in selected}
        if ids!=last:
            print(now(),'analyzing',len(selected),'closed production chains',flush=True)
            analyze(selected,2001);pool();last=ids
        states=[json.loads((OUT/f'{p.stem}-status.json').read_text())['state'] for p in manifests]
        if all(state=='complete' for state in states):
            selected=[p for p in completed() if p['replica']<100]
            assert {p['id'] for p in completed()}==expected
            analyze(selected,2001);pool();break
        if any(state not in ['running','complete'] for state in states):
            raise RuntimeError(f'A batch needs attention: {states}')
        sleep(20)
    from conservative import main as conservative
    from boundary import main as boundary
    from plot import main as plot
    from histories_final import main as histories
    from finalize import main as finalize
    for stage in [conservative,boundary,plot,histories,finalize]:
        print(now(),'stage',stage.__module__,flush=True);stage()
    print(now(),'analysis and artifact reconciliation complete',flush=True)


if __name__=='__main__':main()
