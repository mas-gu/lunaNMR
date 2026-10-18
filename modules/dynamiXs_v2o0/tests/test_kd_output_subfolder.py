# ABOUTME: Tests that loading a titration input defaults the Kd output to a
# ABOUTME: 'kd_analysis' subfolder inside the peak-integration (series_results) folder.

import sys
import types
from pathlib import Path

import pandas as pd
import pytest

_MOD = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_MOD))

from PySide6.QtWidgets import QApplication


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication(sys.argv)


def _fake_main_window():
    return types.SimpleNamespace(show_main_menu=lambda: None, current_dir="")


def _tidy_csv(folder):
    """Minimal titration tidy CSV that load_titration can read."""
    pts = [0.0, 10.0, 25.0]
    rows = []
    for name in ("A17", "K14"):
        for p in pts:
            rows.append((str(p), name, 8.0 + p * 0.001, 120.0, 1000.0 - 2 * p, 2000.0))
    df = pd.DataFrame(rows, columns=["spectrum_name", "assignment",
                                     "ppm_x", "ppm_y", "height", "volume"])
    csv = folder / "series_analysis_tidy.csv"
    df.to_csv(csv, index=False)
    return csv


def test_input_defaults_output_to_kd_analysis_subfolder(app, tmp_path):
    from kd_titration_page import KdTitrationPage
    series_dir = tmp_path / "series_results_ABC"
    series_dir.mkdir()
    csv = _tidy_csv(series_dir)

    page = KdTitrationPage(_fake_main_window())
    page._set_input_file(str(csv))

    assert Path(page.output_dir) == series_dir / "kd_analysis"


DNAJA1_TIDY = (Path(__file__).resolve().parents[4] / "data_test" / "data_binding" / "tidy"
               / "DNAJA1_HSPA8_tidy.csv")


@pytest.mark.skipif(not DNAJA1_TIDY.exists(), reason="data_test not present")
def test_fit_output_prefix_is_series_name(app, tmp_path, monkeypatch):
    """Every series writes the same generic tidy CSV name, so the fit's files must carry
    the series' name. The real worker fits the real DNAJA1 titration; only the message
    box that reports the result is answered, since nothing closes a modal offscreen."""
    import shutil
    import kd_titration_page as mod
    series_dir = tmp_path / "series_results_ABC"
    series_dir.mkdir()
    csv = series_dir / "series_analysis_tidy.csv"
    shutil.copy(DNAJA1_TIDY, csv)

    reported = []
    monkeypatch.setattr(mod, "show_info", lambda parent, title, msg: reported.append(title))
    monkeypatch.setattr(mod, "show_warning", lambda parent, title, msg: reported.append(title))

    page = mod.KdTitrationPage(_fake_main_window())
    page._set_input_file(str(csv))
    page.series_name = "HSPA8"                       # as if a named series was dropped
    page._start_analysis()
    assert page.worker.wait(300_000), "the fit did not finish"
    app.processEvents()
    assert reported, "the page never reported the finished fit"

    written = [f.name for f in (series_dir / "kd_analysis").rglob("*") if f.is_file()]
    assert written
    assert all(name.startswith("HSPA8") for name in written), written
