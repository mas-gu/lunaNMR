# ABOUTME: The shared delay-column parser must handle letter-suffix duplicate measurements.
# ABOUTME: e.g. '600_T1_sample_0o0b' (second acquisition at delay 0.0) -> 0.0.

import sys
from pathlib import Path

_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_T1_T2"
sys.path.insert(0, str(_DIR))


def test_letter_suffix_duplicate_measurements():
    from delay_parser import parse_delay_column as p
    # 'b' marks a repeat acquisition at the same delay -> same value as without it.
    assert p("600_T1_sample_0o0b") == 0.0
    assert p("600_T1_sample_0o3b") == 0.3
    assert p("600_T2_sample_51b") == 51.0
    assert p("600_T2_sample_102b") == 102.0


def test_unit_with_duplicate_suffix():
    from delay_parser import parse_delay_column as p
    # 600_WT naming: explicit ms unit plus a 'b' repeat-acquisition marker.
    assert p("600_T1_sample_0ms") == 0.0
    assert p("600_T1_sample_0msb") == 0.0
    assert p("600_T1_sample_100ms") == 100.0
    assert p("600_T1_sample_300msb") == 300.0
    assert p("600_T1_sample_600msb") == 600.0
    assert p("600_T1_sample_2400ms") == 2400.0
    assert p("600_T2_sample_8ms") == 8.0
    assert p("600_T2_sample_51msb") == 51.0
    assert p("600_T2_sample_102msb") == 102.0
    # unit conversion still applies with a duplicate marker
    assert p("..._5sb") == 5000.0          # 5 s (repeat) -> 5000 ms
    assert p("..._500usb") == 0.5          # 500 us (repeat) -> 0.5 ms


def test_non_delay_headers_return_none():
    from delay_parser import parse_delay_column as p
    # hetNOE spectra are labelled sat/unsat, not delays.
    assert p("600_hetnoe_sample_sat") is None
    assert p("600_hetnoe_sample_unsat") is None


def test_existing_formats_still_parse():
    from delay_parser import parse_delay_column as p
    assert p("300") == 300.0
    assert p("300_2") == 300.0            # underscore-number duplicate
    assert p("600_T1_0o3") == 0.3          # o-decimal
    assert p("003_T2_ADDA_3ms") == 3.0     # embedded delay + unit
    assert p("..._5s") == 5000.0           # unit conversion
    assert p("no_number") is None


def test_an_acquisition_index_is_not_a_delay():
    """A trailing _<digits> with no unit, no decimal separator and no repeat marker is
    an acquisition index, not a time. Parsing it produced a complete relaxation table
    from delays nobody measured: three real spectra named
    03_2D_sample_ref_001..003 read as 1, 2 and 3 ms, and the fit reported success.
    """
    from delay_parser import parse_delay_column as p
    assert p("03_2D_sample_ref_001") is None
    assert p("03_2D_sample_ref_002") is None
    assert p("sample_2") is None
    assert p("experiment_002") is None


def test_a_descriptive_name_needs_a_unit_separator_or_repeat_marker():
    """The three things that distinguish a delay from an index, kept explicit so the
    rule survives a regex edit."""
    from delay_parser import parse_delay_column as p
    assert p("600_T2_sample_102ms") == 102.0    # unit
    assert p("600_T1_sample_0o3") == 0.3       # decimal separator ('o' for '.')
    assert p("600_T2_sample_51b") == 51.0      # single-letter repeat marker
    assert p("600_T2_sample_51") is None       # none of the three


def test_a_series_produced_matrix_is_untouched():
    """`series` normalises delays to bare numbers and marks repeats with _2, which take
    the bare-number path. Tightening the descriptive-name branch must not reach them —
    a real 600_T2 matrix is labelled 8, 17, ..., 271, 51_2, 102_2."""
    from delay_parser import parse_delay_column as p
    for col, expected in (("8", 8.0), ("271", 271.0), ("51_2", 51.0),
                          ("102_2", 102.0), ("2400", 2400.0), ("0", 0.0)):
        assert p(col) == expected, f"{col!r} -> {p(col)}"
