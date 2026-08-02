# ABOUTME: The Kd fit viewer must render non-finite/None fit values as 'n/a', not crash.
# ABOUTME: A degenerate (few-point) fit yields an unbounded error -> json null -> None.

import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')

_DIR = Path(__file__).resolve().parent.parent / "visualization"
sys.path.insert(0, str(_DIR))


def test_fmt_num_handles_none_and_non_finite():
    from kd_titration_fit_viewer import _fmt_num
    assert _fmt_num(None, '.2g') == 'n/a'          # json_safe wrote null for a non-finite error
    assert _fmt_num(float('inf'), '.2g') == 'n/a'  # singular covariance
    assert _fmt_num(float('nan'), '.3g') == 'n/a'
    assert _fmt_num(48.5, '.3g') == '48.5'
    assert _fmt_num(0.058, '.2g') == '0.058'


def test_kd_with_err_omits_missing_error():
    from kd_titration_fit_viewer import _kd_with_err
    assert _kd_with_err(44.53, 11.4) == '44.53 ± 11.4'
    assert _kd_with_err(44.53, None) == '44.53'          # old fit, no stored error
    assert _kd_with_err(44.53, float('inf')) == '44.53'  # singular covariance -> no '± n/a'
    assert _kd_with_err(None, None) == 'n/a'
