# ABOUTME: Tests the Kd page auto-detects a sibling params JSON and populates the panel.
# ABOUTME: Also that a positionless (intensity-only) CSV keeps Observable = intensity.

import os
import sys
import json
import types
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_MOD = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_MOD))


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _page():
    from kd_titration_page import KdTitrationPage
    return KdTitrationPage(types.SimpleNamespace(show_main_menu=lambda: None, current_dir=""))


_CSV = ("spectrum_name,assignment,ppm_x,ppm_y,height,volume\n"
        "0,K14,7.50,110.00,200,2000\n1,K14,7.51,110.05,150,1500\n2,K14,7.53,110.20,90,900\n")


def test_autodetect_populates_panel(app, tmp_path):
    csv = tmp_path / "series_analysis_tidy.csv"
    csv.write_text(_CSV)
    (tmp_path / "kd_params.json").write_text(json.dumps({
        "concentrations": [0.0, 6.25, 12.5], "intensity_scales": [1.2, 1.0, 0.5],
        "protein_conc": 33.0, "alpha": 0.2, "observable": "csp",
        "intensity_value": "volume", "n_bootstrap": 250}))

    page = _page()
    page._set_input_file(str(csv))

    assert [s.value() for s in page.conc_spins] == [0.0, 6.25, 12.5]
    assert [s.value() for s in page.scale_spins] == [1.2, 1.0, 0.5]
    assert page.p0_spin.value() == 33.0
    assert page.alpha_spin.value() == 0.2
    assert page.obs_combo.currentIndex() == 1          # CSP only
    assert page.intvalue_combo.currentText() == "volume"
    assert page.boot_spin.value() == 250


def test_no_params_leaves_defaults(app, tmp_path):
    csv = tmp_path / "series_analysis_tidy.csv"
    csv.write_text(_CSV)
    page = _page()
    page._set_input_file(str(csv))
    # [L] fields default to the point labels; params-driven fields keep defaults.
    assert page.p0_spin.value() == 50.0
    assert page.boot_spin.value() == 0


def test_positionless_csv_keeps_intensity_even_if_params_say_csp(app, tmp_path):
    csv = tmp_path / "series_analysis_tidy.csv"
    csv.write_text("spectrum_name,assignment,ppm_x,ppm_y,height,volume\n"
                   "0,K14,,,200,2000\n1,K14,,,150,1500\n2,K14,,,90,900\n")
    (tmp_path / "kd_params.json").write_text(json.dumps({"observable": "csp"}))
    page = _page()
    page._set_input_file(str(csv))
    assert page.obs_combo.currentIndex() == 2          # intensity-only wins
