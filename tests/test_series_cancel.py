# ABOUTME: A cancelled series stops inside the spectrum being fitted and writes and saves nothing.
# ABOUTME: Real MultiSpectrumProcessor on two DNAJA1 titration points and a 10-peak list cut from DNAJA1.txt.

"""
Cancel stopped a series only between spectra: the spectrum being fitted ran to its end
first. The series now hands its stop event to every spectrum fit, which ends before its
next cluster and raises Cancelled out of the series, which then writes nothing: no
series_results folder, nothing in the session.
"""

import threading
import time
from pathlib import Path

import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
POINTS = ("2D_DNAJA1_0o0.ft", "2D_DNAJA1_HSPA8_0o125.ft")

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


@pytest.fixture
def series(tmp_path):
    """Two linked titration points and a 10-peak list."""
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.utils.file_manager import NMRFileManager
    set_ps2d_config('15N')
    spectra = tmp_path / "spectra"
    spectra.mkdir()
    for point in POINTS:
        (spectra / point).symlink_to(DNAJA1 / point)
    lines = (DNAJA1 / "DNAJA1.txt").read_text().splitlines()
    (tmp_path / "ten.txt").write_text("\n".join(lines[:11]) + "\n")
    yield {"files": [str(spectra / p) for p in POINTS],
           "peaks": NMRFileManager().load_peak_list(str(tmp_path / "ten.txt")),
           "out": tmp_path / "series_results_cancel"}
    set_ps2d_config('15N')


@pytest.mark.parametrize("parallel", [False, True], ids=["sequential", "parallel"])
def test_a_stop_during_the_first_spectrum_ends_the_series_there(series, parallel):
    from lunaNMR.cli import _default_series_params
    from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor
    from lunaNMR.utils.cancellation import Cancelled
    stop = threading.Event()
    messages = []
    set_at = []

    def progress(pct, task, log_msg=None, failed=False):
        messages.append(log_msg or task)
        if not stop.is_set() and "fitted (R²=" in (log_msg or ""):
            stop.set()
            set_at.append(time.monotonic())
        return True

    with pytest.raises(Cancelled):
        MultiSpectrumProcessor(_default_series_params(parallel)).process_nmr_series(
            series["files"], series["peaks"], str(series["out"]),
            peak_source_mode='cascade', progress_callback=progress,
            series_mode='titration', stop=stop)
    latency = time.monotonic() - set_at[0]

    assert not any(POINTS[1] in m for m in messages)
    assert latency < 1.0


@pytest.fixture(scope="module")
def app():
    import os
    import sys
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def dialog(app, series, monkeypatch):
    """The series dialog of a main window holding the whole DNAJA1 list, over the two points."""
    from PySide6.QtWidgets import QInputDialog, QMessageBox
    from lunaNMR.gui.dialogs.series_integration_dialog import SeriesIntegrationDialog
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    monkeypatch.setattr(QInputDialog, "getText", staticmethod(lambda *a, **k: ("cancel", True)))
    boxes = []
    for name in ("information", "warning", "critical"):
        monkeypatch.setattr(QMessageBox, name, staticmethod(
            lambda *a, _name=name, **k: boxes.append((_name, a[1:3]))))

    window = LunaNMRMainWindow()
    window.load_peak_file(str(DNAJA1 / "DNAJA1.txt"))
    window.current_nmr_folder = str(Path(series["files"][0]).parent)
    dlg = SeriesIntegrationDialog(parent=None, main_window=window)
    dlg.cascade_radio.setChecked(True)
    dlg.boxes = boxes
    yield dlg
    if dlg.thread is not None:
        dlg.worker.cancel()
        _wait_until_stopped(app, dlg)
    app.processEvents()
    dlg.deleteLater()
    window.close()
    window.deleteLater()


def _wait_until_stopped(app, dialog, timeout=300):
    start = time.time()
    while dialog.thread.isRunning() and time.time() - start < timeout:
        app.processEvents()
        time.sleep(0.01)
    for _ in range(50):
        app.processEvents()
    return not dialog.thread.isRunning()


def _cancel_once_a_peak_is_fitted(app, dialog, timeout=300):
    """Start the run, Cancel as the first peak is reported fitted; seconds Cancel took to land."""
    seen = []
    dialog._start_processing()
    dialog.worker.progress.connect(lambda current, total, message: seen.append(message or ""))
    start = time.time()
    while not any("fitted (R²=" in m for m in seen) and time.time() - start < timeout:
        app.processEvents()
        time.sleep(0.01)
    assert dialog.thread.isRunning()
    cancelled_at = time.monotonic()
    dialog._cancel_processing()
    assert _wait_until_stopped(app, dialog)
    return time.monotonic() - cancelled_at


def test_cancel_in_the_dialog_lands_inside_the_spectrum_being_fitted(app, dialog):
    latency = _cancel_once_a_peak_is_fitted(app, dialog)

    assert latency < 1.0
    assert "Cancelled" in dialog.status_label.text()
    assert dialog.boxes == []


def test_a_cancelled_series_in_the_dialog_writes_and_saves_nothing(app, dialog, series):
    window = dialog.main_window
    spectra = Path(series["files"][0]).parent

    _cancel_once_a_peak_is_fitted(app, dialog)

    assert list(spectra.glob("series_results_*")) == []
    assert getattr(window, 'batch_results', None) is None
    assert not getattr(window, 'saved_series', {})
    assert dialog.status_label.text() == "Cancelled - nothing written"
    assert "cancelled - nothing written" in dialog.log_text.toPlainText().lower()


def test_a_series_stopped_after_its_first_spectrum_removes_what_it_wrote(series):
    """With training data collection on, the series writes into its folder while it runs;
    a stop before the second spectrum leaves no folder behind."""
    from lunaNMR.cli import _default_series_params
    from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor
    from lunaNMR.utils.cancellation import Cancelled
    params = _default_series_params()
    params['gui_params']['collect_training_data'] = True
    stop = threading.Event()
    wrote_during_the_run = []

    def progress(pct, task, log_msg=None, failed=False):
        if POINTS[1] in (log_msg or ""):
            wrote_during_the_run.append(series["out"].exists())
            stop.set()
        return True

    with pytest.raises(Cancelled):
        MultiSpectrumProcessor(params).process_nmr_series(
            series["files"], series["peaks"], str(series["out"]),
            peak_source_mode='cascade', progress_callback=progress,
            series_mode='titration', stop=stop)

    assert wrote_during_the_run == [True]
    assert not series["out"].exists()


def test_a_folder_that_was_there_before_the_series_is_kept(series):
    """The CLI writes into a folder it made itself: a stopped series removes only its own."""
    from lunaNMR.cli import _default_series_params
    from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor
    from lunaNMR.utils.cancellation import Cancelled
    series["out"].mkdir()
    (series["out"] / "notes.txt").write_text("kept\n")
    stop = threading.Event()
    stop.set()

    with pytest.raises(Cancelled):
        MultiSpectrumProcessor(_default_series_params()).process_nmr_series(
            series["files"], series["peaks"], str(series["out"]),
            peak_source_mode='cascade', series_mode='titration', stop=stop)

    assert sorted(p.name for p in series["out"].iterdir()) == ["notes.txt"]
