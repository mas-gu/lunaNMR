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
