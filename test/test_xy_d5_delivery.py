"""Attained reporting must preserve rejected/omitted-size qualifications."""
import importlib.util
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
CAMPAIGN=ROOT/'experiments/2026-10-13-xy-infinite-volume'
sys.path.insert(0,str(CAMPAIGN))
spec=importlib.util.spec_from_file_location('d5_delivery',CAMPAIGN/'build_attained_d5_delivery.py')
delivery=importlib.util.module_from_spec(spec);spec.loader.exec_module(delivery)


def test_rejected_fit_remains_stress_not_accepted():
    f=dict(included_root_sizes=[6,8,10,12],result=dict(status='rejected_size_fit',dof=2,location_infinite=.23,MC95=[.22,.24]))
    assert delivery.eligible_stress(f,[6,8,10,12])==(True,'rejected_fit_stress_only')
    assert delivery.eligible_stress(f,[8,10,12])==(False,'missing_or_different_sizes')


def test_zero_dof_cannot_enter_envelope():
    f=dict(included_root_sizes=[8,10,12],result=dict(status='zero_dof_diagnostic',dof=0,location_infinite=.23,MC95=[.22,.24]))
    assert not delivery.eligible_stress(f,[8,10,12])[0]


def test_uniform_baseline_preserves_rejection_at_q1():
    import json
    analysis=json.loads((CAMPAIGN/'d5-calibration-analysis-20261014-173811-cut1001-block1-U25-local-all/analysis.json').read_text())
    rows,ledger=delivery.assemble_d5({'primary':{'analysis':analysis}})
    assert len(rows)==7
    q1=next(r for r in rows if r['error_axis']=='vertical' and r['epsilon_over_pi']==1)
    assert not q1['accepted_extrapolation']
    assert q1['baseline_fit']['minimum_N']==6
    assert sum(r['accepted_extrapolation'] for r in rows)==6
    assert all(r['sensitivity_low']<=r['MC_low']<=r['MC_high']<=r['sensitivity_high'] for r in rows)
    assert any(x['included'] and x['reason']=='rejected_fit_stress_only' for x in ledger)
