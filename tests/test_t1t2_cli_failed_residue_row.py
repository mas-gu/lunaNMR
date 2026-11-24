# ABOUTME: `dynamixs t1t2` writes a residue whose fit failed as a NaN row with Success=No, the row
# ABOUTME: the T1/T2 page's fitter writes. Real 700 MHz T1 matrix: two residues do not decay in the window.

import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
MATRIX = (REPO_ROOT.parent / "data_test" / "data_JA1_relaxation" / "run_20260825_1405"
          / "series" / "700_T1" / "peak_intensity_matrix.csv")
NO_DECAY = ['10ValH', '43GluH']

pytestmark = pytest.mark.skipif(not MATRIX.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def cli_run(tmp_path_factory):
    """The envelope and the results file of `dynamixs t1t2` on the real matrix."""
    from lunaNMR.cli import main
    out = tmp_path_factory.mktemp("cli_t1")
    envelope = out / "envelope.json"
    stdout, sys.stdout = sys.stdout, envelope.open("w")
    try:
        code = main(["dynamixs", "t1t2", "--input", str(MATRIX), "--out", str(out / "fits"),
                     "--exp", "T1", "--field-freq", "700", "--time-units", "ms",
                     "--format", "json"])
    finally:
        sys.stdout.close()
        sys.stdout = stdout
    assert code == 0
    summary = json.loads(envelope.read_text())
    return summary, Path(summary['results_file'])


def _rows(results_file):
    lines = Path(results_file).read_text().splitlines()
    return {line.split("\t")[0]: line.split("\t")[1:] for line in lines[1:]}


def test_a_failed_residue_is_a_nan_row(cli_run):
    summary, results_file = cli_run
    rows = _rows(results_file)
    for residue in NO_DECAY:
        assert rows[residue] == ['NaN'] * 6 + ['No', 'NaN']
    assert summary['n_excluded'] == len(NO_DECAY)
    assert summary['n_fitted'] == len(rows) - len(NO_DECAY)
    assert summary['fit_errors'] == []


def test_the_marginal_windows_counted_are_those_the_file_reports(cli_run):
    """A failed residue has no fitted time constant to compare the delay window with."""
    summary, results_file = cli_run
    marginal = [residue for residue, cells in _rows(results_file).items() if float(cells[7]) < 2]
    assert marginal and not set(marginal) & set(NO_DECAY)
    assert summary['n_window_marginal'] == len(marginal)


def test_every_other_residue_keeps_its_fit(cli_run):
    _, results_file = cli_run
    fitted = {residue: cells for residue, cells in _rows(results_file).items()
              if residue not in NO_DECAY}
    assert all(cells[6] == 'Yes' for cells in fitted.values())
    assert all(float(cells[1]) > 0 for cells in fitted.values())


def test_the_results_file_is_the_one_the_multicore_fitter_writes(cli_run, tmp_path):
    _, results_file = cli_run
    sys.path.insert(0, str(REPO_ROOT / "modules" / "dynamiXs_v2o0" / "dynamiXs_T1_T2"))
    from fitmulti__Tx_NMRRE import run_analysis_with_params
    multicore = run_analysis_with_params({
        'input_csv_file': str(MATRIX), 'output_prefix': str(tmp_path / 't1'),
        'results_txt_file': str(tmp_path / 't1_fit_results.txt'), 'experiment_type': 'T1',
        'time_units': 'ms', 'error_method': 'analytical'})
    assert Path(multicore['results_file']).read_text() == results_file.read_text()
