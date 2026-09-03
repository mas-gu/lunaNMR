# ABOUTME: The Kd loader must tell a decayed peak from a fit that stopped on a bound.
# ABOUTME: Zero is a measurement for I/I0 and never one for CSP; a bound is neither.

"""
Two states look identical in a series and need opposite handling.

A residue whose peak has DECAYED reports a near-zero intensity, and that zero IS the
measurement - blanking it from an I/I0 series invents a plateau that is not there. It
has no position, though, so it cannot contribute a chemical-shift perturbation.

A residue whose FIT FAILED - pinned on a position or an intensity bound - reports a
number taken from the bound. The intensity floor is cluster-relative (0.1% of the
brightest member's seed), so that number is a product of its NEIGHBOURS' seeds. It is
not a measurement of anything.

Neither is visible in `success` or R-squared: R-squared is per CLUSTER, so a member
with nothing at its position rides the one its healthy neighbours earn.
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

_DIR = Path(__file__).resolve().parent.parent / 'dynamiXs_Kd'
sys.path.insert(0, str(_DIR))

from kd_input import load_titration_tidy       # noqa: E402


def _row(point, assignment, **kw):
    base = {'spectrum_name': str(point), 'assignment': assignment,
            'ppm_x': 8.2, 'ppm_y': 121.5, 'height': 1.0e9, 'volume': 2.0e9,
            'quality': 'Excellent', 'r_squared': 0.99, 'detection_sigma': 40.0,
            'position_pinned': False, 'intensity_pinned': False}
    base.update(kw)
    return base


def _csv(tmp_path, rows):
    path = tmp_path / 'series_analysis_tidy.csv'
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def test_a_decayed_peak_keeps_its_intensity_and_loses_its_position(tmp_path):
    """Zero IS the measurement for I/I0; blanking it would invent a plateau. But
    there is no position to measure a shift from, so CSP must not see it."""
    path = _csv(tmp_path, [_row(0.0, 'A1'),
                           _row(1.0, 'A1', height=3.0e6, volume=5.0e6,
                                detection_sigma=0.8)])
    _, residues = load_titration_tidy(path)
    a1 = residues['A1']
    assert a1['height'][1] == pytest.approx(3.0e6), 'the decay is the measurement'
    assert pd.isna(a1['ppm_x'][1]) and pd.isna(a1['ppm_y'][1]), (
        'a peak that is gone has no position')


def test_a_bound_pinned_point_is_blanked_everywhere(tmp_path):
    """Its intensity came from a cluster-relative floor - its neighbours' seeds."""
    path = _csv(tmp_path, [_row(0.0, 'A1'),
                           _row(1.0, 'A1', height=4.0e6, volume=8.0e6,
                                detection_sigma=0.4, intensity_pinned=True)])
    _, residues = load_titration_tidy(path)
    a1 = residues['A1']
    for key in ('height', 'volume', 'ppm_x', 'ppm_y'):
        assert pd.isna(a1[key][1]), f'{key} survived a bound-pinned point'


def test_a_legacy_csv_without_the_flags_still_loads(tmp_path):
    """A series run before the flags existed has no such columns at all. Refusing it
    would be worse than the risk; the flags are absent, not negative."""
    rows = []
    for point in (0.0, 1.0):
        row = _row(point, 'A1')
        for key in ('detection_sigma', 'position_pinned', 'intensity_pinned'):
            row.pop(key)
        rows.append(row)
    _, residues = load_titration_tidy(_csv(tmp_path, rows))
    assert all(not pd.isna(v) for v in residues['A1']['height'])
