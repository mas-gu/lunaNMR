# ABOUTME: Tests that BatchResults.to_fitting_dataframe honors the series mode.
# ABOUTME: Guards the titration column headers and series_mode propagation path.

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
