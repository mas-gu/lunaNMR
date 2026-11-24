# ABOUTME: The T1/T2 page's account of a finished fit names the residues whose fit raised an error,
# ABOUTME: with the error's text, and `dynamixs t1t2` carries the same list. Real 700 MHz T1 matrix.

import os
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

REPO_ROOT = Path(__file__).resolve().parents[1]
MATRIX = (REPO_ROOT.parent / "data_test" / "data_JA1_relaxation" / "run_20260825_1405"
          / "series" / "700_T1" / "peak_intensity_matrix.csv")

pytestmark = pytest.mark.skipif(not MATRIX.exists(), reason="data_test not present")


@pytest.fixture
def page():
    import lunaNMR.gui.dialogs.dynamixs_dialog  # noqa: F401 - puts the module on sys.path
    from dynamiXs_gui import T1T2FittingPage
    return T1T2FittingPage


def _summary(tmp_path, **params):
    """The summary of the page's fitter on the real matrix."""
    from dynamiXs_T1_T2 import fitmulti__Tx_NMRRE
    summary = fitmulti__Tx_NMRRE.run_analysis_with_params({
        'input_csv_file': str(MATRIX), 'output_prefix': str(tmp_path / 't1'),
        'results_txt_file': str(tmp_path / 't1_fit_results.txt'), 'experiment_type': 'T1',
        'time_units': 'ms', 'error_method': 'analytical', **params})
    return {**summary, 'experiment_type': 'T1'}


def test_a_fit_that_raised_is_named_with_its_error(page, tmp_path):
    """A starting amplitude that is not a number makes lmfit raise for every residue."""
    summary = _summary(tmp_path, initial_A=float('nan'))
    assert len(summary['fit_errors']) == 60

    lines = page._fit_report(summary)

    assert "ERROR: 60 fit(s) raised an error:" in lines
    for error in summary['fit_errors']:
        assert f"  {error['residue']}: {error['error']}" in lines


def test_a_run_without_errors_reports_none(page, tmp_path):
    lines = page._fit_report(_summary(tmp_path))
    assert not [line for line in lines if "raised an error" in line]
    assert "Not fitted (2): 10ValH, 43GluH" in lines


def test_the_t1rho_command_records_the_fit_errors_of_both_its_fits(tmp_path):
    import json
    from lunaNMR.cli import main
    peaks = next((MATRIX.parents[3] / "600_WT" / "T1").glob("*assignment*.txt"))
    # Any decay fits; only the record is under test, so one matrix serves as both.
    assert main(["dynamixs", "t1rho", "--input", str(MATRIX), "--t1", str(MATRIX),
                 "--peaks", str(peaks), "--omega1", "2000", "--carrier", "118",
                 "--theta", "45", "--field-freq", "700", "--out", str(tmp_path),
                 "--prefix", "f"]) == 0
    record = json.loads((tmp_path / "f_dynamixs.json").read_text())
    assert record['summary']['fit_errors'] == {'T1': [], 'T1rho': []}
