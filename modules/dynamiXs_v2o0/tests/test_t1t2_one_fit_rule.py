# ABOUTME: The single-core and the multicore T1/T2 fitters fit a residue with one function and fail
# ABOUTME: it by one rule, and a refit follows that rule. Real 700 MHz T1 matrix, baseline fixed and free.

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_T1_T2"
sys.path.insert(0, str(_DIR))
REPO_ROOT = Path(__file__).resolve().parents[3]

REAL_T1 = (REPO_ROOT.parent / "data_test" / "data_JA1_relaxation" / "run_20260825_1405"
           / "series" / "700_T1" / "peak_intensity_matrix.csv")

pytestmark = pytest.mark.skipif(not REAL_T1.exists(), reason="data_test not present")

# Two residues of this matrix fit to a time constant far beyond the longest delay.
NO_DECAY = ['10ValH', '43GluH']
FAILED_CELLS = ['NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'No', 'NaN']


def _residues():
    """(name, delays, intensities) of every residue of the real matrix, and of its first
    residue measured at its two shortest delays only."""
    from delay_parser import parse_delay_columns
    matrix = pd.read_csv(REAL_T1)
    x, keep, _ = parse_delay_columns(list(matrix.columns[4:]))
    rows = matrix.iloc[:, 4:].astype(float).to_numpy()[:, keep]
    residues = [(name, x, y) for name, y in zip(matrix['Assignment'], rows)]
    two_points = rows[0].copy()
    two_points[np.argsort(x)[2:]] = np.nan
    return residues + [('two delays', x, two_points)]


def _differences(single, multi):
    """The keys on which two per-residue fits disagree."""
    if set(single) != set(multi):
        return sorted(set(single) ^ set(multi))
    differing = []
    for key, value in single.items():
        other = multi[key]
        if key == 'result':
            same = (value is None) == (other is None)
        elif key in ('x', 'y'):
            same = np.array_equal(value, other, equal_nan=True)
        elif isinstance(value, float) and np.isnan(value):
            same = isinstance(other, float) and np.isnan(other)
        else:
            same = value == other
        if not same:
            differing.append(key)
    return differing


@pytest.mark.parametrize("fit_baseline", [False, True], ids=["fixed", "free"])
def test_both_fitters_return_the_same_fit_of_a_residue(fit_baseline):
    from fit_Tx_NMRRE import fit_single_residue
    from fitmulti__Tx_NMRRE import fit_single_residue_parallel
    residues = _residues()
    disagreements = {}
    outcomes = set()
    for idx, (name, x, y) in enumerate(residues):
        single = fit_single_residue(x, y, name, fit_baseline=fit_baseline)
        multi = fit_single_residue_parallel((x, y, name, None, None, None, 1000, 'analytical',
                                             idx, len(residues), fit_baseline))
        assert multi.pop('idx') == idx
        outcomes.add(single['success'])
        differing = _differences(single, multi)
        if differing:
            disagreements[name] = differing
    assert outcomes == {True, False}
    assert disagreements == {}


def test_a_residue_with_too_few_points_keeps_its_data_and_carries_no_fit():
    from fit_Tx_NMRRE import fit_single_residue
    name, x, y = _residues()[-1]
    fit = fit_single_residue(x, y, name)
    assert fit['success'] is False and fit['result'] is None
    assert all(np.isnan(fit[key]) for key in ('A', 't2', 'C', 'A_err', 't2_err', 'C_err'))
    assert np.array_equal(fit['y'], y, equal_nan=True)


@pytest.fixture
def single_core_run(tmp_path):
    """A single-core T1 run of the real matrix: its fit JSON and its results file."""
    import matplotlib
    matplotlib.use("Agg")
    from fit_Tx_NMRRE import run_analysis_with_params
    out = tmp_path
    summary = run_analysis_with_params({
        'input_csv_file': str(REAL_T1), 'output_prefix': str(out / 't1'),
        'results_txt_file': str(out / 't1_fit_results.txt'), 'experiment_type': 'T1',
        'time_units': 'ms', 'error_method': 'analytical', 'field_name': 'f700',
        'field_freq': 700.0, 'json_folder': str(out)})
    return Path(summary['json_file']), Path(summary['results_file'])


def _cells(results_file, residue):
    for line in Path(results_file).read_text().splitlines()[1:]:
        cells = line.split("\t")
        if cells[0] == residue:
            return cells[1:]
    raise KeyError(residue)


def test_a_refit_that_finds_no_decay_writes_the_failed_row(single_core_run):
    from refit import refit_residue, update_tsv_row
    json_file, results_file = single_core_run
    payload = json.loads(json_file.read_text())
    entry = next(f for f in payload['fits'] if f['residue'] == NO_DECAY[0])
    last = int(np.argmax(payload['metadata']['time_points']))

    refitted = refit_residue(entry, payload['metadata'], [last])
    update_tsv_row(results_file, NO_DECAY[0], refitted, 'T1')

    assert refitted['success'] is False
    assert _cells(results_file, NO_DECAY[0]) == FAILED_CELLS


def test_a_refit_left_with_too_few_points_writes_the_failed_row(single_core_run):
    from refit import refit_residue, update_tsv_row
    json_file, results_file = single_core_run
    payload = json.loads(json_file.read_text())
    entry = next(f for f in payload['fits'] if f['residue'] not in NO_DECAY)
    assert _cells(results_file, entry['residue'])[6] == 'Yes'
    all_but_two = list(range(2, len(payload['metadata']['time_points'])))

    refitted = refit_residue(entry, payload['metadata'], all_but_two)
    update_tsv_row(results_file, entry['residue'], refitted, 'T1')

    assert refitted['success'] is False
    assert _cells(results_file, entry['residue']) == FAILED_CELLS


def test_a_residue_whose_refit_failed_is_fitted_again_from_its_data(single_core_run):
    """A failed fit leaves no value to start the next one from."""
    from refit import refit_residue
    json_file, _ = single_core_run
    payload = json.loads(json_file.read_text())
    entry = next(f for f in payload['fits'] if f['residue'] not in NO_DECAY)
    all_but_two = list(range(2, len(payload['metadata']['time_points'])))
    failed = refit_residue(entry, payload['metadata'], all_but_two)
    assert failed['success'] is False and np.isnan(failed['t2'])

    recovered = refit_residue(failed, payload['metadata'], [])

    assert recovered['success'] is True
    assert recovered['t2'] == pytest.approx(entry['t2'], rel=1e-6)
