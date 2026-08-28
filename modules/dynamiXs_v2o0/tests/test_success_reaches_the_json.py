# ABOUTME: The fitter's own reliability verdict must survive into the JSON it writes.
# ABOUTME: It was computed, printed as "unreliable ... excluded", then dropped.

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_T1_T2"
sys.path.insert(0, str(_DIR))


def _matrix(tmp_path, decaying=True):
    """A T2 matrix with one real decay and one flat residue.

    A flat series has no measurable decay, so its time constant is meaningless -- the
    fitter says so, and the summary counts it as excluded.
    """
    times = [8.0, 51.0, 102.0, 204.0, 400.0]
    rows = [["Peak_Number", "Assignment", "Reference_X", "Reference_Y"] + [str(t) for t in times]]
    rows.append([1, "K3", 8.0, 120.0] + list(1000.0 * np.exp(-np.array(times) / 120.0)))
    rows.append([2, "E4", 8.1, 121.0] + [900.0] * len(times))     # flat: no decay
    path = tmp_path / "peak_intensity_matrix.csv"
    pd.DataFrame(rows).to_csv(path, header=False, index=False)
    return path


def _run(tmp_path):
    from fit_Tx_NMRRE import run_analysis_with_params
    return run_analysis_with_params({
        'input_csv_file': str(_matrix(tmp_path)),
        'output_prefix': str(tmp_path / 'f'),
        'results_txt_file': str(tmp_path / 'f.txt'),
        'experiment_type': 'T2', 'time_units': 'ms',
        'error_method': 'analytical', 'n_bootstrap': 0,
        'field_name': 'field1', 'field_freq': 600.0,
        'json_folder': str(tmp_path),
    })


class TestSuccessIsSerialised:
    """save_fit_data_json already carries baseline_fixed, window_ratio and
    window_marginal, with a comment saying they were dropped here and thereby made
    unreachable from the CLI. `success` was the one still missing -- so a residue the
    fitter itself printed "unreliable ... excluded" for still entered the R2 table,
    where the only downstream filter is finite-and-positive.
    """

    def test_every_fit_carries_its_verdict(self, tmp_path):
        result = _run(tmp_path)
        fits = json.loads(Path(result['json_file']).read_text())['fits']
        assert fits, "no fits written"
        for fit in fits:
            assert 'success' in fit, f"{fit['residue']} has no success flag"
            assert isinstance(fit['success'], bool)

    def test_the_verdict_matches_the_summary_count(self, tmp_path):
        """n_fitted counts reliable residues only; the JSON must agree with it."""
        result = _run(tmp_path)
        fits = json.loads(Path(result['json_file']).read_text())['fits']
        assert sum(1 for f in fits if f['success']) == result['n_fitted']
        assert sum(1 for f in fits if not f['success']) == result['n_excluded']
