"""Analyze only closed, fully audited production files; keep full diagnostic logs."""
import json
from time import sleep
from analyze import OUT,completed,analyze,histories
from pool import main as pool

def main():
    while True:
        selected=[p for p in completed() if p['replica']<100]
        if selected:
            analyze(selected,2001)
            pool()
            # Complete epsilon neighborhoods get all-observable cold and full histories.
            for e in sorted({p['epsilon_milli'] for p in selected}):
                group=[p for p in selected if p['epsilon_milli']==e]
                if len(group)>=36:histories(group)
        status=json.loads((OUT/'production-status.json').read_text())
        if status['state']!='running':break
        sleep(20)
    if status['state']!='complete':raise RuntimeError('Generation did not complete')
    histories(selected)

if __name__=='__main__':main()
