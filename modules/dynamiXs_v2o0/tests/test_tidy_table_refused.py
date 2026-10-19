# ABOUTME: A relaxation fit given a series' tidy table, not its intensity matrix, refuses it by name.
# ABOUTME: The shared check on real JA1 series files, and the T1/T2, methyl and model-free browse slots.

"""
A series writes series_analysis_tidy.csv (one row per spectrum and peak) beside its
peak_intensity_matrix.csv (one row per residue, one column per delay). The fitters read
the matrix; given the tidy table they parsed a residue name as an intensity and stopped
on "could not convert string to float". The check names the matrix to use instead.
"""

import os
import shutil
import sys
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DIR))
sys.path.insert(0, str(_DIR / "dynamiXs_T1_T2"))

SERIES = (Path(__file__).resolve().parents[4] / "data_test" / "data_JA1_relaxation"
          / "run_20260825_1405" / "series" / "600_T1")
TIDY = SERIES / "series_analysis_tidy.csv"
MATRIX = SERIES / "peak_intensity_matrix.csv"

pytestmark = pytest.mark.skipif(not TIDY.exists(), reason="data_test not present")


class TestTheCheck:

    def test_a_tidy_table_beside_its_matrix_names_that_matrix(self):
        from series_inputs import tidy_table_refusal
        refusal = tidy_table_refusal(str(TIDY))
        assert "series_analysis_tidy.csv is a series' tidy table" in refusal
        assert str(MATRIX) in refusal

    def test_the_matrix_named_is_the_one_asked_for(self):
        from series_inputs import tidy_table_refusal
        refusal = tidy_table_refusal(str(TIDY), "peak_volume_matrix.csv")
        assert str(SERIES / "peak_volume_matrix.csv") in refusal

    def test_a_tidy_table_alone_names_the_matrix_to_look_for(self, tmp_path):
        from series_inputs import tidy_table_refusal
        alone = tmp_path / "series_analysis_tidy.csv"
        shutil.copy(TIDY, alone)
        refusal = tidy_table_refusal(str(alone))
        assert "give the series' peak_intensity_matrix.csv instead" in refusal

    @pytest.mark.parametrize("path", [MATRIX, SERIES / "comprehensive_peak_tracking.csv",
                                      SERIES / "no_such_file.csv"])
    def test_anything_else_passes(self, path):
        from series_inputs import tidy_table_refusal
        assert tidy_table_refusal(str(path)) is None


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def window(app):
    from dynamiXs_gui import DynamiXsMainWindow
    w = DynamiXsMainWindow()
    w.saved_series = {}
    yield w
    w.deleteLater()


@pytest.fixture
def no_modal(monkeypatch):
    """Records any message box a page opens instead of blocking on it."""
    from PySide6.QtWidgets import QMessageBox
    shown = []
    for box in ("information", "warning", "critical", "question"):
        monkeypatch.setattr(QMessageBox, box, staticmethod(
            lambda *a, _box=box, **k: shown.append(_box) or QMessageBox.Ok))
    return shown


def _choose(monkeypatch, module, path):
    """The file dialog answers with `path`."""
    monkeypatch.setattr(module, "open_file_dialog", lambda *a, **k: str(path))


class TestTheBrowseSlots:

    def test_the_t1t2_page_refuses_it_on_the_page(self, window, monkeypatch, no_modal):
        import dynamiXs_gui
        page = window.t1t2_page
        _choose(monkeypatch, dynamiXs_gui, TIDY)
        page._browse_data_file("field1_t1")
        assert page.field1_t1_file is None
        assert str(MATRIX) in page.results_text.toPlainText()
        assert no_modal == []

    def test_the_t1t2_page_takes_the_matrix(self, window, monkeypatch, no_modal):
        import dynamiXs_gui
        page = window.t1t2_page
        _choose(monkeypatch, dynamiXs_gui, MATRIX)
        page._browse_data_file("field1_t1")
        assert page.field1_t1_file == str(MATRIX)

    def test_the_model_free_page_refuses_it_on_the_page(self, window, monkeypatch, no_modal):
        import dynamiXs_gui
        page = window.integrated_page
        _choose(monkeypatch, dynamiXs_gui, TIDY)
        page._browse_file("field1_t1", page.field1_t1_display)
        assert page.field1_t1_file is None
        assert str(MATRIX) in page.progress_text.toPlainText()
        assert no_modal == []

    def test_the_methyl_page_refuses_it_on_the_page(self, window, monkeypatch, no_modal):
        from PySide6.QtWidgets import QFileDialog
        page = window.methyl_t2_page
        monkeypatch.setattr(QFileDialog, "getOpenFileName",
                            staticmethod(lambda *a, **k: (str(TIDY), "")))
        page._browse_input_file()
        assert page.input_file is None
        assert str(MATRIX) in page.log_text.toPlainText()
        assert no_modal == []
