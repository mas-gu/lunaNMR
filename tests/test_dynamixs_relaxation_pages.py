# ABOUTME: The DynamiXs relaxation pages fit what the CLI fits and say on the page what it reports.
# ABOUTME: Real JA1 600 MHz T1 series and the sampleNBD methyl series, driven through the real pages.

import os
import sys
import types
from pathlib import Path

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

RUN = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
       / "run_20260825_1405")
T1_SERIES = RUN / "series" / "600_T1"
T1_MATRIX = T1_SERIES / "peak_intensity_matrix.csv"

pytestmark = pytest.mark.skipif(not T1_MATRIX.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


def _dialog(window=None):
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    return DynamiXsDialog(parent=None, main_window=window or types.SimpleNamespace())


@pytest.fixture
def t1t2_page(app, tmp_path):
    d = _dialog()
    page = d.t1t2_page
    page.output_dir = str(tmp_path / "out")
    os.makedirs(page.output_dir)
    yield page
    d.deleteLater()


def _wait(app, page):
    worker = page.worker
    if worker is not None:
        worker.wait()
    app.processEvents()


def _fit_t1(app, page, matrix):
    page.field1_t1_file = str(matrix)
    page.analysis_combo.setCurrentIndex(0)          # Field 1 - T1
    page._start_analysis()
    _wait(app, page)
    return page.results_text.toPlainText()


def _matrix(tmp_path, columns=None, rename=None):
    df = pd.read_csv(T1_MATRIX)
    if columns is not None:
        df = df[["Peak_Number", "Assignment", "Reference_X", "Reference_Y", *columns]]
    if rename:
        df = df.rename(columns=rename)
    path = tmp_path / "matrix.csv"
    df.to_csv(path, index=False)
    return path


class TestTheT1T2PageReportsTheFit:

    def test_a_column_with_no_delay_is_named_on_the_page(self, app, t1t2_page, tmp_path):
        log = _fit_t1(app, t1t2_page, _matrix(tmp_path, rename={"2400": "last_plane"}))
        assert "last_plane" in log
        assert "no parseable delay" in log

    def test_the_fitted_count_is_given_out_of_the_total(self, app, t1t2_page, tmp_path):
        log = _fit_t1(app, t1t2_page, _matrix(tmp_path))
        n_total = len(pd.read_csv(T1_MATRIX))
        assert f"/{n_total} residues fitted" in log

    def test_a_fit_where_nothing_fitted_is_not_reported_as_a_success(self, app, t1t2_page,
                                                                      tmp_path):
        log = _fit_t1(app, t1t2_page, _matrix(tmp_path, columns=["0", "50"]))
        assert "completed successfully" not in log
        assert "No residue fitted" in log
        assert "field1_t1" not in t1t2_page.fitted_experiments


class TestAFileChosenForASlotReplacesTheSeriesThere:
    """A slot that once held a dropped series is fitted from the file chosen after it,
    not from the series still in memory."""

    @pytest.fixture
    def page(self, app, tmp_path, series_from_tidy):
        br, _ = series_from_tidy(T1_SERIES / "series_analysis_tidy.csv")
        d = _dialog(types.SimpleNamespace(saved_series={"t1_series": br},
                                          current_project_path=None))
        page = d.t1t2_page
        page.output_dir = str(tmp_path / "out")
        os.makedirs(page.output_dir)
        page._on_series_dropped("field1_t1", "t1_series",
                                str(T1_SERIES / "series_analysis_tidy.csv"))
        yield page
        d.deleteLater()

    def _run(self, app, page):
        page.analysis_combo.setCurrentIndex(0)
        page._start_analysis()
        _wait(app, page)
        return page.results_text.toPlainText()

    def test_the_dropped_series_is_fitted_from_memory(self, app, page):
        assert "Using data from series 't1_series' in memory" in self._run(app, page)

    def test_browsing_a_file_replaces_it(self, app, page, tmp_path, monkeypatch):
        import dynamiXs_gui
        matrix = _matrix(tmp_path)
        monkeypatch.setattr(dynamiXs_gui, "open_file_dialog", lambda *a, **k: str(matrix))
        page._browse_data_file("field1_t1")
        assert "in memory" not in self._run(app, page)
        assert page.field1_t1_file == str(matrix)

    def test_auto_loading_a_folder_replaces_it(self, app, page, tmp_path, monkeypatch):
        import dynamiXs_gui
        folder = tmp_path / "folder"
        folder.mkdir()
        pd.read_csv(T1_MATRIX).to_csv(folder / "T1_matrix.csv", index=False)
        monkeypatch.setattr(dynamiXs_gui, "open_directory_dialog", lambda *a, **k: str(folder))
        page._auto_load_folder("field1")
        assert "in memory" not in self._run(app, page)
