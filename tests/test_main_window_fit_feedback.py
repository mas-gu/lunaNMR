# ABOUTME: What the main window shows after Detect and Fit Spectrum: the statistics panel, the status line, the Voigt tabs.
# ABOUTME: Real LunaNMRMainWindow offscreen on the DNAJA1 spectrum, a 10-peak list cut from its peak list.

import os
import sys
import time
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
SPECTRUM = "2D_DNAJA1_0o0.ft"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def window(app, tmp_path, monkeypatch):
    """A window showing the spectrum, linked into a folder of its own so the fit's CSV
    lands there, and a 10-peak list."""
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.components import ProgressDialog
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    set_ps2d_config('15N')
    w = LunaNMRMainWindow()
    (tmp_path / SPECTRUM).symlink_to(DNAJA1 / SPECTRUM)
    lines = (DNAJA1 / "DNAJA1.txt").read_text().splitlines()
    peaks = tmp_path / "ten.txt"
    peaks.write_text("\n".join(lines[:11]) + "\n")
    w.load_data(str(tmp_path / SPECTRUM), str(peaks))
    monkeypatch.setattr(ProgressDialog, "exec", lambda self: 0)
    yield w
    for worker in (getattr(w, 'fitting_worker', None), getattr(w, 'detection_worker', None)):
        if worker is not None:
            worker.wait()
    app.processEvents()
    w.close()
    w.deleteLater()
    set_ps2d_config('15N')


def _wait_for(app, worker, timeout=600):
    start = time.time()
    while worker.isRunning() and time.time() - start < timeout:
        app.processEvents()
        time.sleep(0.05)
    assert not worker.isRunning()
    for _ in range(50):
        app.processEvents()


def _detect(app, window):
    window.on_detect_peaks()
    _wait_for(app, window.detection_worker)


def _fit(app, window):
    """Click Fit Spectrum and wait; returns the summary the fit handed the window."""
    summaries = []
    handle = window.on_fitting_complete

    def recording(summary, progress_dialog):
        summaries.append(summary)
        handle(summary, progress_dialog)

    window.on_fitting_complete = recording
    window.on_fit_spectrum()
    _wait_for(app, window.fitting_worker)
    return summaries[0]


def _suptitle(plotter):
    title = plotter.figure._suptitle
    return title.get_text() if title is not None else ''


def test_the_statistics_panel_starts_empty(window):
    assert window.stats_label.text() == "No processing performed yet"


def test_the_statistics_panel_reports_a_detection(app, window):
    _detect(app, window)

    stats = window.integrator.get_detection_statistics()
    assert stats['total_peaks'] == 10
    assert window.stats_label.text().splitlines() == [
        f"Detected: {stats['detected_peaks']}", "Total: 10", f"Rate: {stats['detection_rate']:.1f}%"]


def test_the_statistics_panel_reports_the_fit_that_followed_a_detection(app, window):
    _detect(app, window)

    summary = _fit(app, window)

    assert summary['successful_peaks'] > 0
    assert window.stats_label.text().splitlines() == [
        f"Fitted: {summary['successful_peaks']}/{summary['total_peaks']}",
        f"Failed fits: {summary['failed_peaks']}",
        f"Average R²: {summary['average_r_squared']:.3f}"]


def test_loading_new_data_empties_the_statistics_panel(app, window):
    _detect(app, window)

    window.load_data(window.current_nmr_file, window.current_peak_file)

    assert window.stats_label.text() == "No processing performed yet"


def test_the_status_line_gives_the_average_r_squared_of_the_fit(app, window):
    summary = _fit(app, window)

    assert summary['average_r_squared'] > 0
    assert (f"{summary['successful_peaks']}/{summary['total_peaks']} peaks "
            f"(avg R²: {summary['average_r_squared']:.3f})") in window.statusBar().currentMessage()


def test_the_voigt_tabs_show_the_first_successful_fit(app, window, recwarn):
    summary = _fit(app, window)

    assert [str(w.message) for w in recwarn if issubclass(w.category, UserWarning)] == []

    first = next(result for result in summary['fitted_results'] if result.get('success'))
    assert window.voigt_3d_plotter.current_result['assignment'] == first['assignment']
    assert str(first['assignment']) in _suptitle(window.voigt_3d_plotter)
    assert str(first['assignment']) in _suptitle(window.voigt_2d_plotter)
    assert window.current_voigt_result['assignment'] == first['assignment']


def test_reset_empties_the_statistics_panel(app, window, monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    _detect(app, window)
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.Yes))

    window.on_reset_results()

    assert window.stats_label.text() == "No processing performed yet"
