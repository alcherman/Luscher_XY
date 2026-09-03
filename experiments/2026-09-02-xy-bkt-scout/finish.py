"""Finish disk-only analysis after the scouting sampler completes."""
import json
from pathlib import Path
import subprocess
import sys
import time

OUT=Path(__file__).resolve().parent


def main():
    while True:
        state=json.loads((OUT/'status.json').read_text())
        if 'finished_utc' in state:
            if state.get('state')!='complete' or state['failures']:
                raise RuntimeError('Campaign incomplete; inspect status and logs')
            break
        time.sleep(5)
    # Preserve the venv invocation path instead of depending on base-interpreter packages.
    python=sys.executable
    with (OUT/'histories.log').open('w') as log:
        history=subprocess.Popen([python,str(OUT/'analyze.py'),'--histories'],stdout=log,stderr=subprocess.STDOUT)
        try:
            for name,args in [('analysis',['analyze.py','--cut','2001','--plots']),
                              ('threshold-plot',['plot_threshold.py']),('phase',['phase.py']),('maps',['plot_maps.py'])]:
                with (OUT/f'{name}.log').open('w') as output:
                    subprocess.run([python,str(OUT/args[0]),*args[1:]],stdout=output,stderr=subprocess.STDOUT,check=True)
                print(name,'complete',flush=True)
        finally:
            if history.wait():raise RuntimeError('History rendering failed')
        print('histories complete',flush=True)
    with (OUT/'finalize.log').open('w') as output:
        subprocess.run([python,str(OUT/'finalize.py')],stdout=output,stderr=subprocess.STDOUT,check=True)
    print('final audit and report complete',flush=True)


if __name__=='__main__':main()
