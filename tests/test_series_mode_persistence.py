# ABOUTME: Tests that a series' series_mode survives project save/load.
# ABOUTME: A titration series must reload as 'titration' so the Kd drag list keeps it.

import types
from pathlib import Path

import pytest

from lunaNMR.processors.multi_spectrum_processor import BatchResults
from lunaNMR.utils.project_manager import ProjectManager


def _save_and_reload(tmp_path, series_mode):
    """Round-trip one series through _save_series_results / _load_series_results."""
    batch = BatchResults()
    batch.series_mode = series_mode
    batch.metadata['output_folder'] = str(tmp_path / "series_results_run")
    batch.add_result('spectrum_0', {'status': 'success', 'height': 1.0})

    mw_save = types.SimpleNamespace(saved_series={'ERDJ6_HSPA1A': batch})
    bundle = tmp_path / "proj.lunaNMR"
    bundle.mkdir()
    ProjectManager(mw_save)._save_series_results(bundle)

    mw_load = types.SimpleNamespace()
    ProjectManager(mw_load)._load_series_results(bundle)
    return mw_load.saved_series['ERDJ6_HSPA1A']


def test_titration_series_mode_survives_roundtrip(tmp_path):
    reloaded = _save_and_reload(tmp_path, 'titration')
    assert reloaded.series_mode == 'titration'


def test_time_series_mode_survives_roundtrip(tmp_path):
    reloaded = _save_and_reload(tmp_path, 'time')
    assert reloaded.series_mode == 'time'
