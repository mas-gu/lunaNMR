# ABOUTME: Tests the wide vs-sequence table builder used for the Kd CSV/JSON export.
# ABOUTME: rows = residue (sequence-ordered), cols = titration points, blank = peak absent.

import sys
import math
from pathlib import Path

_KD_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_Kd"
sys.path.insert(0, str(_KD_DIR))


def _fits():
    # A17 listed before K14 to prove sequence-number ordering (K14 must come first).
    # Point 2 of K14 is undetected (0.0 sentinel) -> that cell must be blank.
    return [
        {'residue': 'A17', 'series': {
            'ppm_x': [8.00, 8.02, 8.05], 'ppm_y': [120.0, 120.1, 120.3],
            'height': [100.0, 80.0, 50.0]}},
        {'residue': 'K14', 'series': {
            'ppm_x': [7.50, 7.51, 0.0], 'ppm_y': [110.0, 110.05, 0.0],
            'height': [200.0, 150.0, 0.0]}},
    ]


def test_header_and_row_order():
    from kd_models import ref_vs_point_table
    header, rows = ref_vs_point_table(_fits(), [0.0, 1.0, 2.0], 'csp')
    assert header == ['Residue', '0.0', '1.0', '2.0']
    assert [r[0] for r in rows] == ['K14', 'A17']       # sequence number, not letter


def test_reference_column_and_blank_for_absent():
    from kd_models import ref_vs_point_table
    _, rows = ref_vs_point_table(_fits(), [0.0, 1.0, 2.0], 'csp')
    k14 = next(r for r in rows if r[0] == 'K14')
    assert k14[1] == 0.0            # CSP of the reference point vs itself is 0
    assert k14[3] == ''             # point 2 absent -> blank


def test_intensity_ratio_reference_is_one():
    from kd_models import ref_vs_point_table
    _, rows = ref_vs_point_table(_fits(), [0.0, 1.0, 2.0], 'intensity')
    a17 = next(r for r in rows if r[0] == 'A17')
    assert a17[1] == 1.0                       # I(0)/I(0)
    assert math.isclose(a17[2], 80.0 / 100.0)  # I(1)/I(0)


def test_residue_sort_key_number_before_letter():
    from kd_models import residue_sort_key
    assert sorted(['A17', 'K14', 'G2'], key=residue_sort_key) == ['G2', 'K14', 'A17']
