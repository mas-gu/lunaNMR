# ABOUTME: Tests for DelayExtractor time-series and titration value extraction.
# ABOUTME: Covers the o-for-dot titration suffix convention and value-agnostic sequencing.

import pytest

from lunaNMR.utils.delay_extractor import DelayExtractor


class TestTimeModeUnchanged:
    """Default mode is 'time'; existing ms/s behavior must be byte-for-byte preserved."""

    def test_default_mode_is_time(self):
        assert DelayExtractor().mode == "time"

    def test_ms_suffix(self):
        assert DelayExtractor().extract_value("T1_50ms.ft") == 50.0

    def test_seconds_converted_to_ms(self):
        assert DelayExtractor().extract_value("T1_1s.ft") == 1000.0

    def test_time_mode_rejects_bare_titration_suffix(self):
        # A dimensionless titration-style suffix is not a valid delay in time mode.
        assert DelayExtractor(mode="time").extract_value("titr_1o0.ft") is None

    def test_extract_value_matches_legacy_extract_delay_ms(self):
        e = DelayExtractor()
        for name in ("T1_50ms.ft", "T2_2.5s.ft", "no_delay_here.ft"):
            assert e.extract_value(name) == e.extract_delay_ms(name)


class TestTitrationExtraction:
    """Titration mode parses the _<value> suffix with the o-for-dot convention."""

    @pytest.mark.parametrize("name,expected", [
        ("sample_1o0.ft", 1.0),
        ("sample_0o5.ft", 0.5),
        ("sample_1.0.ft", 1.0),     # literal dot also accepted
        ("sample_2o0.ft", 2.0),     # whole-number point written with the o separator
        ("sample_0o25.ft", 0.25),
        ("titr_10o0.ft", 10.0),
    ])
    def test_titration_suffix(self, name, expected):
        assert DelayExtractor(mode="titration").extract_value(name) == expected

    def test_titration_without_extension(self):
        assert DelayExtractor(mode="titration").extract_value("sample_0o5") == 0.5

    @pytest.mark.parametrize("name", [
        "reference_spectrum.ft",   # no value suffix
        "sample_2.ft",             # bare integer: needs a separator (2o0)
        "experiment_002.ft",       # index-style file must not become a titration point
        "ref_0.ft",
    ])
    def test_titration_requires_decimal_separator(self, name):
        assert DelayExtractor(mode="titration").extract_value(name) is None


class TestTitrationSequencing:
    """The sort / sequence / column-name machinery is value-agnostic across modes."""

    def test_sort_orders_by_titration_value(self):
        files = ["a_1o0.ft", "a_0o5.ft", "a_2o0.ft"]
        ordered = DelayExtractor(mode="titration").sort_files_with_sequence(files)
        values = [v for _, v, _ in ordered]
        assert values == [0.5, 1.0, 2.0]

    def test_column_mapping_titration(self):
        files = ["a_1o0.ft", "a_0o5.ft"]
        mapping = DelayExtractor(mode="titration").build_column_mapping(files)
        assert mapping["a_0o5.ft"] == "0.5"
        assert mapping["a_1o0.ft"] == "1"   # whole number formats as int, plots as 1.0

    def test_duplicate_titration_points_get_sequence(self):
        files = ["a_0o5.ft", "b_0o5.ft"]
        mapping = DelayExtractor(mode="titration").build_column_mapping(files)
        cols = sorted(mapping.values())
        assert cols == ["0.5", "0.5_2"]

    def test_two_decimal_titration_point_keeps_precision(self):
        # 0o15 -> 0.15 must NOT be truncated to "0.1" in the column name,
        # which would also collide with a real 0.1 point.
        ext = DelayExtractor(mode="titration")
        assert ext.get_column_name(0.15, 1) == "0.15"
        mapping = ext.build_column_mapping(["x_0o15.ft", "x_0o1.ft"])
        assert mapping["x_0o15.ft"] == "0.15"
        assert mapping["x_0o1.ft"] == "0.1"

    def test_time_mode_column_names_unchanged(self):
        ext = DelayExtractor()  # time
        assert ext.get_column_name(50.0, 1) == "50"
        assert ext.get_column_name(100.0, 2) == "100_2"
        assert ext.get_column_name(12.5, 1) == "12.5"


class TestParserConsistency:
    """The DelayExtractor and decay-widget titration parsers must agree.

    The CSV/column path uses DelayExtractor; the plot path uses the widget
    function. They must accept exactly the same filename forms or points
    appear in one but not the other.
    """

    @pytest.mark.parametrize("name", [
        "sample_1o0.ft",
        "sample_1.0.ft",
        "sample_0o5.ft",
        "sample_0o25.ft",
        "titr_10o0.ft",
        "sample_2.ft",          # bare integer: both reject
        "experiment_002.ft",    # index: both reject
        "reference.ft",
        "T1_50ms.ft",           # time form: both reject as titration
    ])
    def test_both_parsers_agree(self, name):
        from lunaNMR.gui.components.intensity_decay_widget import (
            extract_titration_from_spectrum_name,
        )
        de = DelayExtractor(mode="titration").extract_value(name)
        widget = extract_titration_from_spectrum_name(name)
        assert de == widget, f"{name}: DelayExtractor={de} widget={widget}"
