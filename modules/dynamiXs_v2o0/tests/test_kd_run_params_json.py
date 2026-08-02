# ABOUTME: Tests that every Kd run writes an importable <prefix>_kd_params.json next to results.
# ABOUTME: Reproducibility: re-running an analysis needs no manual re-entry of binding params.

import sys
from pathlib import Path

import pytest

_KD_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_Kd"
sys.path.insert(0, str(_KD_DIR))

_CSV = """spectrum_name,assignment,ppm_x,ppm_y,height,volume
0,K14,7.500,110.00,200.0,2000.0
1,K14,7.510,110.05,150.0,1500.0
2,K14,7.530,110.20,90.0,900.0
0,A17,8.000,120.00,100.0,1000.0
1,A17,8.020,120.10,80.0,800.0
2,A17,8.050,120.30,50.0,500.0
"""


# A 3-point fit is deliberately degenerate (covariance can't be estimated); that
# OptimizeWarning is expected here and must not pollute test output.
@pytest.mark.filterwarnings("ignore:Covariance of the parameters could not be estimated")
def test_run_writes_importable_params_json(tmp_path):
    from kd_fit import run_kd_analysis_with_params
    from kd_params import load_params

    csv = tmp_path / "series_analysis_tidy.csv"
    csv.write_text(_CSV)
    out = tmp_path / "kd_analysis"

    params = dict(input_csv_file=str(csv), output_dir=str(out), output_prefix="TESTSER",
                  concentrations=[0.0, 6.25, 12.5], intensity_scales=[1.2, 1.0, 0.5],
                  protein_conc=50.0, alpha=0.14, observables=['csp', 'intensity'],
                  intensity_value='height', n_bootstrap=0)
    result = run_kd_analysis_with_params(params)

    pf = result.get('params_file')
    assert pf and Path(pf).exists()
    assert Path(pf) == out / "TESTSER_kd_params.json"

    loaded = load_params(pf)
    assert loaded['concentrations'] == [0.0, 6.25, 12.5]
    assert loaded['intensity_scales'] == [1.2, 1.0, 0.5]
    assert loaded['protein_conc'] == 50.0
    assert loaded['observables'] == ['csp', 'intensity']
    assert loaded['intensity_value'] == 'height'
