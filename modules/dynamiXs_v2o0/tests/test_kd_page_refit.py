# ABOUTME: Tests the Kd page "Load for refit" — restores params + input data into the
# ABOUTME: panel, rebuilding the input from embedded series when the original CSV is gone.

import os
import sys
import types
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_MOD = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_MOD))
sys.path.insert(0, str(_MOD / "dynamiXs_Kd"))


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _page(monkeypatch):
    import kd_titration_page as mod
    monkeypatch.setattr(mod, "show_info", lambda *a, **k: None)
    monkeypatch.setattr(mod, "show_warning", lambda *a, **k: None)
    from kd_titration_page import KdTitrationPage
    return KdTitrationPage(types.SimpleNamespace(show_main_menu=lambda: None, current_dir=""))


_META = {
    'detected_points': ['0', '1', '2'],
    'output_dir': '/tmp',
    'concentrations': [0.0, 10.0, 25.0],
    'intensity_scales': [1.0, 0.5, 0.5],
    'p0': 40.0, 'alpha': 0.2, 'observable': 1,
    'intensity_value': 'height', 'bootstrap_iter': 100,
}

_FIT_DATA = {
    'metadata': {'points': ['0', '1', '2'], 'concentrations': [0.0, 10.0, 25.0]},
    'fits': [
        {'residue': 'A1', 'series': {'ppm_x': [8.0, 8.01, 8.02], 'ppm_y': [120, 120, 120],
                                     'height': [100, 80, 50], 'volume': [200, 160, 100]}},
        {'residue': 'K2', 'series': {'ppm_x': [7.5, 7.5, 7.5], 'ppm_y': [110, 110, 110],
                                     'height': [90, 70, 40], 'volume': [180, 140, 80]}},
    ],
}


def test_refit_restores_params_and_keeps_existing_csv(app, tmp_path, monkeypatch):
    csv = tmp_path / "series_analysis_tidy.csv"
    csv.write_text("spectrum_name,assignment,ppm_x,ppm_y,height,volume\n"
                   "0,A1,8.0,120,100,200\n1,A1,8.01,120,80,160\n2,A1,8.02,120,50,100\n")
    meta = dict(_META, input_file=str(csv))
    page = _page(monkeypatch)
    page.load_saved_for_refit({'meta': meta, 'fit_data': _FIT_DATA}, name="MYFIT")
    assert page.input_file == str(csv)                       # original CSV kept
    assert [s.value() for s in page.conc_spins] == [0.0, 10.0, 25.0]
    assert [s.value() for s in page.scale_spins] == [1.0, 0.5, 0.5]   # scales preserved
    assert page.p0_spin.value() == 40.0
    assert page.obs_combo.currentIndex() == 1
    assert page._current_analysis_name == "MYFIT"


def test_refit_rebuilds_input_when_csv_missing(app, tmp_path, monkeypatch):
    meta = dict(_META, input_file=str(tmp_path / "gone.csv"))   # does not exist
    page = _page(monkeypatch)
    page.load_saved_for_refit({'meta': meta, 'fit_data': _FIT_DATA}, name="REBUILT")
    assert page.input_file and os.path.exists(page.input_file)   # rebuilt temp CSV
    assert page.input_file.endswith("_series_analysis_tidy.csv")
    assert [s.value() for s in page.scale_spins] == [1.0, 1.0, 1.0]  # reset to avoid double-scaling
    # the rebuilt CSV is a valid titration the fitter can read
    from kd_input import load_titration
    points, residues = load_titration(page.input_file)
    assert len(points) == 3
    assert set(residues) == {'A1', 'K2'}
