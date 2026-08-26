# ABOUTME: T1/T2 fitters must converge on real-scale data via data-driven initials.
# ABOUTME: Guards against the fixed A=5/t2=100 non-convergence and the arg-order bug.

import sys
from pathlib import Path

import numpy as np
import pytest

_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_T1_T2"
sys.path.insert(0, str(_DIR))

# Real-scale decay: amplitudes ~1e6 over a 0-1.5 delay range, decay constant ~0.5.
X = np.array([0.0, 0.1, 0.3, 0.6, 0.9, 0.15, 1.2, 1.5])
Y = np.array([3691303.2, 3030295.3, 2026739.6, 1614479.0,
              1172109.4, 2938805.1, 860203.4, 699484.9])


def test_single_core_converges_with_data_driven_initials():
    from fit_Tx_NMRRE import fit_single_residue
    r = fit_single_residue(X, Y, "3.0")          # no initial_A/t2 -> data-driven
    # Band widened when the baseline was fixed at zero: this trace only decays to 19%
    # of its start (t_max/T ~ 1.9), where a free C used to absorb part of the decay
    # and shorten T. The test guards convergence, not the convention.
    assert 0.3 < r['t2'] < 1.5                    # sensible, not a degenerate ~1e9
    assert np.isfinite(r['t2_err'])


def test_multicore_worker_converges_with_correct_arg_order():
    from fitmulti__Tx_NMRRE import fit_single_residue_parallel
    # 10-tuple: (x, y, name, initial_A, initial_t2, initial_C, n_bootstrap,
    #            error_method, idx, total)
    r = fit_single_residue_parallel((X, Y, "3.0", None, None, None, 0, 'analytical', 0, 1))
    assert r['success']
    assert 0.3 < r['t2'] < 1.5


def test_degenerate_residues_flagged_unreliable():
    from fit_Tx_NMRRE import fit_single_residue
    from fitmulti__Tx_NMRRE import fit_single_residue_parallel
    zeros = np.zeros_like(Y)
    assert fit_single_residue(X, Y, "3.0")['success'] is True          # clean decay
    assert fit_single_residue(X, zeros, "zero")['success'] is False    # no signal
    r = fit_single_residue_parallel((X, zeros, "zero", None, None, None, 0, 'analytical', 0, 1))
    assert r['success'] is False


@pytest.mark.filterwarnings("ignore:FigureCanvasAgg is non-interactive")
@pytest.mark.filterwarnings("ignore:Attempting to set identical low and high ylims")
def test_summary_excludes_degenerate_residues(tmp_path):
    from fit_Tx_NMRRE import run_analysis_with_params
    import pandas as pd
    labels = [f"600_T1_{str(t).replace('.', 'o')}" for t in X]
    header = ["Peak_Number", "Assignment", "Reference_X", "Reference_Y"] + list(labels)
    rows = [
        [1, "3.0", 8.3, 128.3] + list(Y),          # clean
        [2, "4.0", 8.6, 122.8] + list(Y),          # clean
        [3, "9.0", 9.1, 121.6] + [0.0] * len(X),   # all-zeros -> excluded
    ]
    csv = tmp_path / "m.csv"
    pd.DataFrame(rows, columns=header).to_csv(csv, index=False)
    result = run_analysis_with_params({
        'input_csv_file': str(csv), 'output_prefix': str(tmp_path / 'f'),
        'results_txt_file': str(tmp_path / 'f.txt'), 'experiment_type': 'T1',
    })
    assert result['n_fitted'] == 2 and result['n_excluded'] == 1
    assert 0.3 < result['mean_t2'] < 1.5          # not polluted by the 1e11 outlier
