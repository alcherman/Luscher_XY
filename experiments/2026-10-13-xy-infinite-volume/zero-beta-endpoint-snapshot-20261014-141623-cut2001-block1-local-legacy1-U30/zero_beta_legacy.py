"""Read-only adapters for archived zero-coupling bootstrap summaries."""
import json
from pathlib import Path
import numpy as np
import run

def load():
    exp=run.OUT.parent; entries={};provenance=[]
    # Later pooled archives supersede earlier summaries at identical coordinates.
    for name in ['2026-09-02-xy-bkt-refine','2026-09-03-xy-bkt-low-kappa','2026-09-03-xy-bkt-kappa-zero']:
        p=exp/name;f=p/'pooled.json';audit=p/'execution-audit.json'
        assert json.loads(audit.read_text())['passed']
        provenance.append(dict(path=str(f),sha256=run.sha256(f),audit_sha256=run.sha256(audit)))
        for r in json.loads(f.read_text()):
            if r['kappa_milli']!=0 or not .53<=r['epsilon_milli']/1000<=.59:continue
            q=r['epsilon_milli']/1000;N=r['L'];npz=p/'pooled'/f"{r['id']}_bootstrap.npz"
            b=np.load(npz)['samples'];assert b.shape==(2000,7)
            entries[2,N,q]=(dict(D=2,N=N,epsilon_over_pi=q,summary=dict(values=r['values'],errors=r['errors'],flags=r['flags'],individual_flags=r.get('individual_flags',{}),retained=r['retained'],archived_cut=2001),source_ids=r.get('chains',[]),legacy=True,source=str(f),samples_path=str(npz),samples_sha256=run.sha256(npz)),b)
    for D in [3,4]:
        p=exp/f'2026-09-03-xy-{D}d-kappa-zero';f=p/'results.json';prov=p/'analysis-provenance.json'
        meta=json.loads(prov.read_text());assert meta['cut']==2001 and meta['bootstrap_draws']==4000
        provenance.append(dict(path=str(f),sha256=run.sha256(f),analysis_provenance_sha256=run.sha256(prov),archive_audit_scope='Archived per-run closure audits; not rerun by this adapter.'))
        for r in json.loads(f.read_text()):
            q=r['epsilon']/np.pi if D==3 else r['epsilon_over_pi'];N=r['L']
            token=f"e{r['epsilon_centi']:03d}" if D==3 else f'q{round(q*1000):04d}'
            npz=p/f'bootstrap_{token}_L{N:03d}.npz';b=np.load(npz)['samples']
            assert b.shape==(4000,)
            # Fixed first 2000 draws align Monte Carlo array lengths; originals retained.
            entries[D,N,q]=(dict(D=D,N=N,epsilon_over_pi=q,summary=dict(values=[r['binder']],errors=[r['error']],flags=r['flags'],diagnostics=r['diagnostics'],archived_cut=2001),source_ids=[],legacy=True,source=str(f),samples_path=str(npz),samples_sha256=run.sha256(npz)),b[:2000,None])
    return entries,provenance
