import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments/2026-10-13-xy-infinite-volume'))
import run_d5_final_sensitivity as manager

def test_failed_audit_cannot_join_final_cohort(tmp_path,monkeypatch):
    monkeypatch.setattr(manager,'OUT',tmp_path);(tmp_path/'jobs').mkdir()
    (tmp_path/'jobs'/'bad.json').write_text(json.dumps(dict(state='complete',audit=dict(passed=False))))
    states,failed=manager.coverage({'bad':{},'unstarted':{}})
    assert failed==['bad'];assert states=={'complete':1,'pending':1}

def row(sizes,status,location):
    return dict(included_root_sizes=sizes,result=dict(status=status,location_infinite=location,MC95=[location-.001,location+.001]))

def test_common_size_comparison_does_not_hide_missing_root():
    base=row([6,8,10,12],'conditional_D5_size_fit',.2)
    comparison=manager.compare_fit(base,row([8,10,12],'conditional_D5_size_fit',.3))
    assert not comparison['same_sizes'];assert comparison['central_shift'] is None
    assert not comparison['both_conditional_fits']

def test_rejected_and_zero_dof_fits_never_count_as_validated_comparisons():
    base=row([6,8,10],'conditional_D5_size_fit',.2)
    for status in ['rejected_size_fit','zero_dof_diagnostic']:
        comparison=manager.compare_fit(base,row([6,8,10],status,.201))
        assert comparison['central_shift'] is not None
        assert not comparison['both_conditional_fits']
