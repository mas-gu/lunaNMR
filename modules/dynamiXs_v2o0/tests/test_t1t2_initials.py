# ABOUTME: T1/T2 fitters must converge on real-scale data via data-driven initials.
# ABOUTME: Guards against the fixed A=5/t2=100 non-convergence and the arg-order bug.

import sys
from pathlib import Path

import numpy as np

_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_T1_T2"
sys.path.insert(0, str(_DIR))

# Real-scale decay: amplitudes ~1e6 over a 0-1.5 delay range, decay constant ~0.5.
X = np.array([0.0, 0.1, 0.3, 0.6, 0.9, 0.15, 1.2, 1.5])
Y = np.array([3691303.2, 3030295.3, 2026739.6, 1614479.0,
              1172109.4, 2938805.1, 860203.4, 699484.9])


def test_single_core_converges_with_data_driven_initials():
    from fit_Tx_NMRRE import fit_single_residue
    r = fit_single_residue(X, Y, "3.0")          # no initial_A/t2 -> data-driven
    assert 0.3 < r['t2'] < 0.7                    # sensible, not a degenerate ~1e9
    assert np.isfinite(r['t2_err'])


def test_multicore_worker_converges_with_correct_arg_order():
    from fitmulti__Tx_NMRRE import fit_single_residue_parallel
    # 10-tuple: (x, y, name, initial_A, initial_t2, initial_C, n_bootstrap,
    #            error_method, idx, total)
    r = fit_single_residue_parallel((X, Y, "3.0", None, None, None, 0, 'analytical', 0, 1))
    assert r['success']
    assert 0.3 < r['t2'] < 0.7
