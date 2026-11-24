# ABOUTME: The single-core and the multicore T1/T2 fitters are one analysis run: the same fit JSON
# ABOUTME: with every residue flagged, the same summary, and a fit that raises kept as a failed residue.

import json
import sys
from pathlib import Path

import numpy as np
import pytest

_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_T1_T2"
sys.path.insert(0, str(_DIR))
REPO_ROOT = Path(__file__).resolve().parents[3]

REAL_T1 = (REPO_ROOT.parent / "data_test" / "data_JA1_relaxation" / "run_20260825_1405"
           / "series" / "700_T1" / "peak_intensity_matrix.csv")

pytestmark = pytest.mark.skipif(not REAL_T1.exists(), reason="data_test not present")

# Two residues of this matrix fit to a time constant far beyond the longest delay.
NO_DECAY = ['10ValH', '43GluH']
# What differs between a run in this process and one over a pool, and nothing else may.
RUN_FACTS = ('n_cores_used', 'processing_time', 'results_file', 'json_file', 'plots_prefix')


def _run(module_name, folder, **params):
    import importlib
    import matplotlib
    matplotlib.use("Agg")
    module = importlib.import_module(module_name)
    folder.mkdir()
    return module.run_analysis_with_params({
        'input_csv_file': str(REAL_T1), 'output_prefix': str(folder / 't1'),
        'results_txt_file': str(folder / 't1_fit_results.txt'), 'experiment_type': 'T1',
        'time_units': 'ms', 'error_method': 'analytical', 'field_name': 'f700',
        'field_freq': 700.0, 'json_folder': str(folder), **params})


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    """The real matrix fitted by each fitter: {fitter: (summary, folder)}."""
    base = tmp_path_factory.mktemp("both_fitters")
    return {name: (_run(module, base / name), base / name)
            for name, module in (('single', 'fit_Tx_NMRRE'), ('multi', 'fitmulti__Tx_NMRRE'))}


@pytest.mark.parametrize("fitter", ["single", "multi"])
def test_the_fit_json_holds_every_residue_with_its_flag(runs, fitter):
    summary, folder = runs[fitter]
    payload = json.loads((folder / 'f700_T1_fit_data.json').read_text())
    flags = {fit['residue']: fit['success'] for fit in payload['fits']}
    assert len(flags) == summary['n_total'] == payload['metadata']['n_residues']
    assert [residue for residue, fitted in flags.items() if not fitted] == NO_DECAY


def test_both_fitters_write_the_same_files(runs):
    single, multi = runs['single'][1], runs['multi'][1]
    for name in ('f700_T1_fit_data.json', 't1_fit_results.txt'):
        assert (single / name).read_text() == (multi / name).read_text()


def test_both_fitters_return_the_same_summary(runs):
    single, multi = runs['single'][0], runs['multi'][0]
    assert set(single) == set(multi)
    assert single['n_cores_used'] == 1 and multi['n_cores_used'] >= 1
    for key in set(single) - set(RUN_FACTS):
        assert single[key] == multi[key], key
    assert single['failed_residues'] == NO_DECAY
    assert single['n_excluded'] == len(NO_DECAY)
    assert single['n_fitted'] + single['n_excluded'] == single['n_total']
    assert single['fit_errors'] == []


@pytest.mark.parametrize("module", ["fit_Tx_NMRRE", "fitmulti__Tx_NMRRE"])
def test_a_fit_that_raises_is_a_failed_residue_whose_error_is_kept(module, tmp_path, capsys):
    """A starting amplitude that is not a number makes the model undefined, and lmfit
    raises for every residue."""
    summary = _run(module, tmp_path / 'out', initial_A=float('nan'))

    assert summary['n_fitted'] == 0 and summary['n_excluded'] == summary['n_total'] == 60
    assert [error['residue'] for error in summary['fit_errors']] == summary['failed_residues']
    assert all("NaN" in error['error'] for error in summary['fit_errors'])
    assert np.isnan(summary['mean_t2']) and np.isnan(summary['mean_error'])
    assert np.isnan(summary['t2_range']).all()
    rows = (tmp_path / 'out' / 't1_fit_results.txt').read_text().splitlines()[1:]
    assert all(row.split("\t")[1:] == ['NaN'] * 6 + ['No', 'NaN'] for row in rows)
    printed = capsys.readouterr().out
    first = summary['fit_errors'][0]
    assert f"{first['residue']}: {first['error']}" in printed
