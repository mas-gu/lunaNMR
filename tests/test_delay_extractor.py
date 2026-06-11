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
        ("sample_2.ft", 2.0),       # bare integer
        ("sample_0o25.ft", 0.25),
        ("titr_10o0.ft", 10.0),
    ])
    def test_titration_suffix(self, name, expected):
        assert DelayExtractor(mode="titration").extract_value(name) == expected

    def test_titration_without_extension(self):
        assert DelayExtractor(mode="titration").extract_value("sample_0o5") == 0.5

    def test_titration_no_suffix_returns_none(self):
        assert DelayExtractor(mode="titration").extract_value("reference_spectrum.ft") is None


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
