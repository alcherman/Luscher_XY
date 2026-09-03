"""Analyze completed files without launching simulations; then build artifacts."""
import json
from time import sleep
from run import OUT,now
from analyze import completed,analyze
from pool import main as pool


def main():
    last=None
    while True:
        selected=completed();ids={p['id'] for p in selected}
        if selected and ids!=last:
            print(now(),'analyzing',len(selected),'closed chains',flush=True)
            analyze(selected,2001);pool();last=ids
        manifests=sorted((OUT/'batches').glob('*.json'))
        states=[json.loads(s.read_text())['state'] if (s:=OUT/f'{p.stem}-status.json').exists() else 'pending' for p in manifests]
        if states and all(s=='complete' for s in states):
            expected={p['id'] for path in manifests for p in json.loads(path.read_text())['points']}
            selected=completed();assert {p['id'] for p in selected}==expected
            analyze(selected,2001);pool();break
        if any(s not in ['complete','pending','running'] for s in states):raise RuntimeError(states)
        sleep(20)
    from check_inputs import main as check_inputs
    from conservative import main as conservative
    from plot import main as plot
    from histories_final import main as histories
    from finalize import main as finalize
    for stage in [check_inputs,conservative,plot,histories,finalize]:
        print(now(),'stage',stage.__module__,flush=True);stage()
    print(now(),'analysis complete',flush=True)


if __name__=='__main__':main()
