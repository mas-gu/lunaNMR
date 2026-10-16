# ABOUTME: A relaxation fit skips a point with no measurement and fits the residue on the rest.
# ABOUTME: An empty matrix cell must cost one point, never the residue or the whole run.

"""
A series writes an empty cell where a fit stopped on a bound: the number there came from
the bound, not the data. lmfit raises on NaN by default, so before this every fitter
either aborted the run (single-core T1/T2) or dropped the whole residue (multicore,
methyl). Only Kd skipped the point.
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_T1_T2"
sys.path.insert(0, str(_DIR))
REPO_ROOT = Path(__file__).resolve().parents[3]

X = np.array([10.0, 30.0, 50.0, 70.0, 90.0, 110.0, 150.0])
T2 = 60.0
Y = 1.0e9 * np.exp(-X / T2)


def _with_gap(y, index):
    y = y.copy()
    y[index] = np.nan
    return y


class TestMonoExponential:

    def test_single_core_fits_the_measured_points(self):
        from fit_Tx_NMRRE import fit_single_residue
        r = fit_single_residue(X, _with_gap(Y, 3), "A1")
        assert r['success']
        assert r['t2'] == pytest.approx(T2, rel=1e-3)

    def test_multicore_fits_the_measured_points(self):
        from fitmulti__Tx_NMRRE import fit_single_residue_parallel
        r = fit_single_residue_parallel(
            (X, _with_gap(Y, 3), "A1", None, None, None, 0, 'analytical', 0, 1))
        assert r['success']
        assert r['t2'] == pytest.approx(T2, rel=1e-3)

    def test_the_result_keeps_the_full_delay_grid(self):
        """The JSON writes every residue's intensities against ONE shared time axis, so
        a residue that returned only its measured points would misalign with it."""
        from fit_Tx_NMRRE import fit_single_residue
        r = fit_single_residue(X, _with_gap(Y, 3), "A1")
        assert len(r['x']) == len(X) and len(r['y']) == len(X)
        assert np.isnan(r['y'][3])

    def test_too_few_measured_points_is_a_failed_fit_not_an_exception(self):
        from fit_Tx_NMRRE import fit_single_residue
        y = np.full_like(Y, np.nan)
        y[0] = Y[0]
        r = fit_single_residue(X, y, "A1")
        assert r['success'] is False
        assert np.isnan(r['t2'])


class TestMethyl:

    def test_methyl_fits_the_measured_points(self):
        from fit_methyl_T2 import fit_single_residue_methyl
        y = 0.5 * 1.0e9 * (np.exp(-X / 80.0) + np.exp(-X / 20.0))
        r = fit_single_residue_methyl(X, _with_gap(y, 2), "M1")
        assert np.isfinite(r['T2_avg'])
        assert len(r['y']) == len(X)


def test_t1t2_cli_survives_an_empty_cell(tmp_path):
    """One empty cell used to abort the whole run with no output written."""
    matrix = pd.DataFrame({'Peak_Number': [1, 2], 'Assignment': ['A1', 'A2'],
                           'Reference_X': [8.0, 8.5], 'Reference_Y': [120.0, 121.0]})
    for i, delay in enumerate(X):
        matrix[f'{int(delay)}'] = [Y[i], 0.8 * Y[i]]
    matrix.loc[0, '70'] = np.nan
    path = tmp_path / 'peak_intensity_matrix.csv'
    matrix.to_csv(path, index=False)

    out = subprocess.run(
        [sys.executable, '-m', 'lunaNMR', 'dynamixs', 't1t2', '--input', str(path),
         '--out', str(tmp_path / 'fits'), '--exp', 'T2', '--field-freq', '600',
         '--time-units', 'ms', '--format', 'json'],
        cwd=str(REPO_ROOT), capture_output=True, text=True)
    assert out.returncode == 0, out.stderr[-2000:]
    summary = json.loads(out.stdout)
    assert summary['n_fitted'] == 2, summary
