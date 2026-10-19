# ABOUTME: A DynamiXs page refuses a second Start while its analysis still runs, and says so on the page.
# ABOUTME: Real JA1 600 MHz T1 matrix through the real T1/T2 page; the running thread is never dropped.

"""
Every Start built a new worker into the same `self.worker`. A second click while a fit
ran dropped the only reference to the running QThread; Python collected it mid-run and
Qt aborted the application ("QThread: Destroyed while thread is still running"), after
both fits had written the same output files.
"""

import os
import sys
import types
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

T1_MATRIX = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
             / "run_20260825_1405" / "series" / "600_T1" / "peak_intensity_matrix.csv")

pytestmark = pytest.mark.skipif(not T1_MATRIX.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def page(app, tmp_path):
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    dialog = DynamiXsDialog(parent=None, main_window=types.SimpleNamespace())
    page = dialog.t1t2_page
    page.output_dir = str(tmp_path / "out")
    os.makedirs(page.output_dir)
    page.field1_t1_file = str(T1_MATRIX)
    page.analysis_combo.setCurrentIndex(0)          # Field 1 - T1
    yield page
    if page.worker is not None:
        page.worker.wait()
    app.processEvents()
    dialog.deleteLater()


def test_a_second_start_while_a_fit_runs_is_refused(app, page):
    page._start_analysis()
    first = page.worker
    assert first is not None and first.isRunning()

    page._start_analysis()

    assert page.worker is first
    assert "already running" in page.results_text.toPlainText()
    first.wait()
    app.processEvents()


def test_after_it_finishes_a_new_start_runs(app, page):
    page._start_analysis()
    page.worker.wait()
    app.processEvents()

    page._start_analysis()

    assert page.worker is not None and page.worker.isRunning()
