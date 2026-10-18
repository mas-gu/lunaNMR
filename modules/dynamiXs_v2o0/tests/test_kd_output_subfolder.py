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


def _fit(app, tmp_path, equivalents):
    """Fit the real DNAJA1 titration through the page with the real worker; return the
    page and the series folder."""
    import shutil
    import kd_titration_page as mod
    series_dir = tmp_path / "series_results_ABC"
    series_dir.mkdir()
    csv = series_dir / "series_analysis_tidy.csv"
    shutil.copy(DNAJA1_TIDY, csv)

    page = mod.KdTitrationPage(_fake_main_window())
    page._set_input_file(str(csv))
    page.series_name = "HSPA8"                       # as if a named series was dropped
    # DNAJA1's points are equivalents; absolute is chosen deliberately, never defaulted.
    page.conc_units_combo.setCurrentIndex(1 if equivalents else 0)
    page._start_analysis()
    assert page.worker.wait(300_000), "the fit did not finish"
    app.processEvents()
    return page, series_dir


@pytest.mark.skipif(not DNAJA1_TIDY.exists(), reason="data_test not present")
def test_fit_output_prefix_is_series_name(app, tmp_path):
    """Every series writes the same generic tidy CSV name, so the fit's files must carry
    the series' name."""
    _, series_dir = _fit(app, tmp_path, equivalents=True)
    written = [f.name for f in (series_dir / "kd_analysis").rglob("*") if f.is_file()]
    assert written
    assert all(name.startswith("HSPA8") for name in written), written


@pytest.mark.skipif(not DNAJA1_TIDY.exists(), reason="data_test not present")
def test_the_page_states_a_fit_it_can_report(app, tmp_path):
    """On the page, not in a modal box: a modal waits for a click nothing gives here."""
    page, _ = _fit(app, tmp_path, equivalents=True)
    assert page.result_label.isVisibleTo(page)
    assert page.result_label.text().startswith("Fitted")


@pytest.mark.skipif(not DNAJA1_TIDY.exists(), reason="data_test not present")
def test_the_page_states_a_fit_that_is_not_usable(app, tmp_path):
    """The point labels (0-3 equivalents) read as micromolar: the verdict says the Kd is
    outside what the titration resolves, and the page says so."""
    page, _ = _fit(app, tmp_path, equivalents=False)
    assert "not usable" in page.result_label.text()
