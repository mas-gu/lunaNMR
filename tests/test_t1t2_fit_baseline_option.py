# ABOUTME: A T1/T2 run asked to fit the baseline fits it, and says so in its summary and JSON.
# ABOUTME: Real JA1 600 MHz T1 series matrix.

"""
run_analysis_with_params read params['fit_baseline'] and never passed it to the fitter,
and reported 'baseline_fixed' from the module default - so a caller asking for a fitted
baseline got C = 0 and a summary saying what it had asked for was not done, or the reverse.
"""

import json
import sys
from pathlib import Path

import pytest

RUN = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
       / "run_20260825_1405")
MATRIX = RUN / "series" / "600_T1" / "peak_intensity_matrix.csv"

pytestmark = pytest.mark.skipif(not MATRIX.exists(), reason="data_test not present")


def _run(tmp_path, fit_baseline):
    import matplotlib
    matplotlib.use("Agg")
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "modules" / "dynamiXs_v2o0"
                           / "dynamiXs_T1_T2"))
    from fit_Tx_NMRRE import run_analysis_with_params
    out = tmp_path / ("fitted" if fit_baseline else "fixed")
    out.mkdir()
    summary = run_analysis_with_params({
        'input_csv_file': str(MATRIX), 'output_prefix': str(out / "t1"),
        'results_txt_file': str(out / "t1_fit_results.txt"), 'experiment_type': 'T1',
        'time_units': 'ms', 'error_method': 'analytical', 'n_bootstrap': 0,
        'field_name': 'f600', 'field_freq': 600.0, 'json_folder': str(out),
        'fit_baseline': fit_baseline})
    fits = json.loads(Path(summary['json_file']).read_text())['fits']
    return summary, fits


def test_a_fitted_baseline_is_fitted_and_reported(tmp_path):
    summary, fits = _run(tmp_path, True)
    assert summary['baseline_fixed'] is False
    assert all(f['baseline_fixed'] is False for f in fits)
    assert any(f['C'] != 0.0 for f in fits)


def test_the_default_keeps_the_baseline_at_zero(tmp_path):
    summary, fits = _run(tmp_path, False)
    assert summary['baseline_fixed'] is True
    assert all(f['C'] == 0.0 for f in fits)
