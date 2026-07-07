# ABOUTME: Guards that the never-wired ML prediction hook does not run by default.
# ABOUTME: With no ml_config, _apply_ml_predictions must no-op silently (no warning).

import pandas as pd

from lunaNMR.processors.single_spectrum_processor import SingleSpectrumProcessor


def _proc():
    # _apply_ml_predictions gates before touching instance state, so skip __init__.
    return SingleSpectrumProcessor.__new__(SingleSpectrumProcessor)


def test_ml_prediction_off_by_default_no_warning(capsys):
    peaks = pd.DataFrame({'Position_X': [8.0], 'Position_Y': [120.0]})
    result = _proc()._apply_ml_predictions(peaks, None)
    assert result is None
    captured = capsys.readouterr()
    assert "ML prediction failed" not in (captured.out + captured.err)


def test_ml_prediction_off_with_empty_processing_options(capsys):
    peaks = pd.DataFrame({'Position_X': [8.0], 'Position_Y': [120.0]})
    result = _proc()._apply_ml_predictions(peaks, {})
    assert result is None
    captured = capsys.readouterr()
    assert "ML prediction failed" not in (captured.out + captured.err)
