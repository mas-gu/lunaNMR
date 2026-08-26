# ABOUTME: A series run must record which spectrum each matrix column came from.
# ABOUTME: The column labels are delays, so without this the filename link is lost for good.

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pytest

from lunaNMR.processors.multi_spectrum_processor import (
    build_series_metadata, measure_repeat_scale)
from lunaNMR.utils.delay_extractor import DelayExtractor


def _mapping(names, mode="time"):
    return DelayExtractor(mode=mode).build_column_mapping(names)


class TestSeriesMetadata:
    def test_every_column_maps_back_to_its_spectrum(self):
        names = ["T1_50ms.ft", "T1_300ms.ft", "T1_600ms.ft"]
        meta = build_series_metadata(names, _mapping(names))
        pairs = {c["column"]: c["spectrum"] for c in meta["columns"]}
        assert pairs == {"50": "T1_50ms.ft", "300": "T1_300ms.ft", "600": "T1_600ms.ft"}

    def test_repeat_acquisitions_record_which_file_owns_the_suffixed_column(self):
        """Which duplicate gets `_2` follows a case-sensitive sort, so it must be recorded,
        not inferred: `T2_8ms` wins the bare label over `bT2_8ms` only because T < b."""
        names = ["T2_8ms.ft", "bT2_8ms.ft"]
        meta = build_series_metadata(names, _mapping(names))
        pairs = {c["column"]: c["spectrum"] for c in meta["columns"]}
        assert pairs == {"8": "T2_8ms.ft", "8_2": "bT2_8ms.ft"}

    def test_delays_are_recorded_in_ms_whatever_the_filename_unit(self):
        """`series` normalises to ms — _2s becomes 2000. The sidecar states the unit
        so no downstream step has to guess."""
        names = ["T1_2s.ft", "T1_500ms.ft"]
        meta = build_series_metadata(names, _mapping(names))
        assert meta["value_units"] == "ms"
        by_file = {c["spectrum"]: c["value"] for c in meta["columns"]}
        assert by_file["T1_2s.ft"] == 2000.0
        assert by_file["T1_500ms.ft"] == 500.0

    def test_an_unparseable_delay_is_null_not_a_boolean_flag(self):
        """A filename the delay parser cannot read falls back to a stem column. `value`
        is null there — a separate boolean beside it invites reading the flag as the
        number, and True == 1 in arithmetic, so the mistake would be silent."""
        names = ["T1_50ms.ft", "T1_2400.ft"]
        meta = build_series_metadata(names, _mapping(names))
        by_file = {c["spectrum"]: c for c in meta["columns"]}
        assert by_file["T1_50ms.ft"]["value"] == 50.0
        assert by_file["T1_2400.ft"]["value"] is None      # bare number: ms or s unknowable
        assert by_file["T1_2400.ft"]["column"] == "T1_2400"
        assert meta["n_value_unparsed"] == 1

    def test_each_column_entry_has_exactly_the_documented_keys(self):
        names = ["T1_50ms.ft"]
        meta = build_series_metadata(names, _mapping(names))
        assert set(meta["columns"][0]) == {"column", "spectrum", "value"}

    def test_titration_mode_records_points_not_milliseconds(self):
        names = ["titr_0.ft", "titr_1.ft"]
        meta = build_series_metadata(names, _mapping(names, mode="titration"),
                                     series_mode="titration")
        assert meta["value_units"] == "point"

    def test_counts_are_reported(self):
        names = ["T1_50ms.ft", "T1_300ms.ft"]
        meta = build_series_metadata(names, _mapping(names))
        assert meta["n_spectra"] == 2
        assert meta["series_mode"] == "time"


class TestRepeatScale:
    """A series with repeat acquisitions must report how well they agree, measured on the
    fitted intensities it just produced — the pre-flight estimate is box maxima."""

    @staticmethod
    def _matrix(columns_and_scales, n_peaks=20):
        import pandas as pd
        base = np.linspace(1e6, 1e4, n_peaks)
        data = {'Peak_Number': range(1, n_peaks + 1),
                'Assignment': [f'{i}ValH' for i in range(1, n_peaks + 1)],
                'Reference_X': [8.0] * n_peaks, 'Reference_Y': [120.0] * n_peaks}
        for col, scale in columns_and_scales.items():
            data[col] = base * scale
        return pd.DataFrame(data)

    def _meta_columns(self, names):
        return build_series_metadata(names, _mapping(names))['columns']

    def test_matched_repeats_report_a_ratio_of_one(self):
        names = ['T2_8ms.ft', 'T2_102ms.ft', 'T2_b_8ms.ft', 'T2_b_102ms.ft']
        m = self._matrix({'8': 1.0, '102': 0.5, '8_2': 1.0, '102_2': 0.5})
        got = measure_repeat_scale(m, self._meta_columns(names))
        assert got['ratio'] == pytest.approx(1.0, abs=1e-6)
        assert got['needs_normalisation'] is False

    def test_a_scaled_sub_series_is_measured_and_the_correction_reported(self):
        names = ['T2_8ms.ft', 'T2_102ms.ft', 'T2_b_8ms.ft', 'T2_b_102ms.ft']
        m = self._matrix({'8': 1.0, '102': 0.5, '8_2': 0.88, '102_2': 0.44})
        got = measure_repeat_scale(m, self._meta_columns(names))
        assert got['ratio'] == pytest.approx(0.88, abs=1e-6)
        assert got['scale'] == pytest.approx(1 / 0.88, abs=1e-6)
        assert got['needs_normalisation'] is True
        assert sorted(got['per_value']) == [8.0, 102.0]

    def test_the_second_column_is_the_one_carrying_the_suffix(self):
        """The ratio is second-relative-to-first, so which is which must be unambiguous."""
        names = ['T2_8ms.ft', 'T2_b_8ms.ft']
        m = self._matrix({'8': 1.0, '8_2': 0.5})
        got = measure_repeat_scale(m, self._meta_columns(names))
        rec = got['per_value'][8.0]
        assert rec['first'] == '8' and rec['second'] == '8_2'
        assert rec['ratio'] == pytest.approx(0.5, abs=1e-6)

    def test_a_series_without_repeats_reports_nothing(self):
        names = ['T2_8ms.ft', 'T2_17ms.ft']
        m = self._matrix({'8': 1.0, '17': 0.9})
        assert measure_repeat_scale(m, self._meta_columns(names)) is None

    def test_it_measures_but_never_applies(self):
        """Reporting is in scope; rescaling the data is a decision for the user."""
        names = ['T2_8ms.ft', 'T2_b_8ms.ft']
        m = self._matrix({'8': 1.0, '8_2': 0.88})
        before = m.copy()
        measure_repeat_scale(m, self._meta_columns(names))
        assert m.equals(before)


class TestMetadataFileWriting:
    """The sidecar must be written even when nothing fitted, and carry the repeat scale
    when something did."""

    def _processor(self, tmp_path, names):
        from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor
        p = MultiSpectrumProcessor.__new__(MultiSpectrumProcessor)
        p.output_folder = str(tmp_path)
        p.delay_mapping = _mapping(names)
        p.series_mode = 'time'
        return p

    def test_written_with_no_intensity_matrix(self, tmp_path):
        import json
        names = ['T2_8ms.ft', 'T2_17ms.ft']
        self._processor(tmp_path, names)._create_series_metadata_file(
            {'results': {n: {} for n in names}}, None)
        meta = json.loads((tmp_path / 'series_metadata.json').read_text())
        assert meta['n_spectra'] == 2
        assert 'repeat_scale' not in meta

    def test_repeat_scale_lands_in_the_file(self, tmp_path):
        import json
        names = ['T2_8ms.ft', 'T2_b_8ms.ft']
        matrix = TestRepeatScale._matrix({'8': 1.0, '8_2': 0.88})
        self._processor(tmp_path, names)._create_series_metadata_file(
            {'results': {n: {} for n in names}}, matrix)
        meta = json.loads((tmp_path / 'series_metadata.json').read_text())
        assert meta['repeat_scale']['ratio'] == pytest.approx(0.88, abs=1e-6)
        assert meta['repeat_scale']['needs_normalisation'] is True
