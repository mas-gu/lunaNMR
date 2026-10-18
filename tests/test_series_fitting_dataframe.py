# ABOUTME: Tests that BatchResults.to_fitting_dataframe honors the series mode.
# ABOUTME: Guards the titration column headers and series_mode propagation path.

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from lunaNMR.processors.multi_spectrum_processor import (
    BatchResults,
    MultiSpectrumProcessor,
)


def _titration_results():
    return {
        '600_T1_0o5.ft': {'success': True, 'integration_results': [
            {'assignment': 'R1', 'volume': 100.0}]},
        '600_T1_1o0.ft': {'success': True, 'integration_results': [
            {'assignment': 'R1', 'volume': 50.0}]},
    }


def test_batch_results_defaults_to_time():
    assert BatchResults().series_mode == 'time'


def test_to_fitting_dataframe_uses_titration_columns():
    br = BatchResults()
    br.series_mode = 'titration'
    br.results = _titration_results()

    df = br.to_fitting_dataframe()
    cols = [c for c in df.columns if c != 'Assignment']
    assert cols == ['0.5', '1']  # titration points, sorted ascending


def test_time_mode_leaves_titration_names_unparsed():
    # Without titration mode the o-suffix is not a delay, so headers fall back
    # to the raw spectrum name (the bug was this happening for titration runs).
    br = BatchResults()  # series_mode == 'time'
    br.results = _titration_results()

    df = br.to_fitting_dataframe()
    cols = [c for c in df.columns if c != 'Assignment']
    assert cols == ['600_T1_0o5', '600_T1_1o0']


def test_convert_propagates_series_mode():
    proc = MultiSpectrumProcessor.__new__(MultiSpectrumProcessor)
    proc.series_mode = 'titration'
    batch = proc._convert_to_batch_results({'results': {}, 'summary': {}}, None, None)
    assert batch.series_mode == 'titration'


# A real T1 series whose fits include Failed points and points on their intensity floor.
def _newest_flagged_t1_series():
    found = sorted((Path(__file__).resolve().parents[2] / "data_test" / "data_asyn").glob(
        "*/T1_*/series_results_*/series_analysis_tidy.csv"))
    flagged = [p for p in found if 'intensity_pinned' in p.read_text().splitlines()[0]]
    return flagged[-1] if flagged else None


_SERIES_WITH_UNMEASURED_POINTS = _newest_flagged_t1_series()


@pytest.mark.skipif(_SERIES_WITH_UNMEASURED_POINTS is None, reason="data_test not present")
def test_the_fitting_frame_leaves_failed_and_floor_pinned_points_empty(series_from_tidy):
    """As the intensity and volume matrices do: those numbers are a placeholder or a
    bound, not a measurement."""
    br, tidy = series_from_tidy(_SERIES_WITH_UNMEASURED_POINTS)
    frame = br.to_fitting_dataframe(value_column='volume').set_index('Assignment')
    unmeasured = (tidy.quality == 'Failed') | (tidy.intensity_pinned == True)  # noqa: E712
    assert unmeasured.sum() > 0
    for row in tidy.itertuples():
        cell = frame.loc[row.assignment, row.spectrum_name]
        if unmeasured[row.Index]:
            assert np.isnan(cell), (row.assignment, row.spectrum_name)
        else:
            assert cell == pytest.approx(row.volume), (row.assignment, row.spectrum_name)
