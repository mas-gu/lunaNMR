# ABOUTME: The series dialog's run ends when it is cancelled, and the dialog stays open while it runs.
# ABOUTME: Real SeriesIntegrationDialog and LunaNMRMainWindow offscreen, two DNAJA1 points, a 10-peak list.

"""
Cancel set a flag and nothing else: the worker emitted neither finished nor error, the
only two signals that quit its thread, so the thread ran on and the dialog stayed on
'Cancelling...' with Start disabled for good. And the dialog could be closed, by its
Close button, Esc or the window's X, while the run went on fitting against the PS2D
config the main window was free to change again.
"""

import os
import sys
import time
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
POINTS = ("2D_DNAJA1_0o0.ft", "2D_DNAJA1_HSPA8_0o125.ft")

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def dialog(app, tmp_path, monkeypatch):
    """The series dialog of a main window holding a 10-peak list, over two linked points."""
    from PySide6.QtWidgets import QInputDialog, QMessageBox
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.dialogs.series_integration_dialog import SeriesIntegrationDialog
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    set_ps2d_config('15N')
    spectra = tmp_path / "spectra"
    spectra.mkdir()
    for point in POINTS:
        (spectra / point).symlink_to(DNAJA1 / point)
    lines = (DNAJA1 / "DNAJA1.txt").read_text().splitlines()
    (tmp_path / "ten.txt").write_text("\n".join(lines[:11]) + "\n")

    monkeypatch.setattr(QInputDialog, "getText", staticmethod(lambda *a, **k: ("lifetime", True)))
    boxes = []
    for name in ("information", "warning", "critical"):
        monkeypatch.setattr(QMessageBox, name, staticmethod(
            lambda *a, _name=name, **k: boxes.append((_name, a[1:3]))))

    window = LunaNMRMainWindow()
    window.load_peak_file(str(tmp_path / "ten.txt"))
    window.current_nmr_folder = str(spectra)
    dlg = SeriesIntegrationDialog(parent=None, main_window=window)
    dlg.cascade_radio.setChecked(True)
    dlg.boxes = boxes
    yield dlg
    if dlg.thread is not None:
        dlg.worker.cancel()
        _wait_until_stopped(app, dlg)
        dlg.thread.quit()
        dlg.thread.wait()
    app.processEvents()
    dlg.deleteLater()
    window.close()
    window.deleteLater()
    set_ps2d_config('15N')


def _wait_until_stopped(app, dialog, timeout=300):
    start = time.time()
    while dialog.thread.isRunning() and time.time() - start < timeout:
        app.processEvents()
        time.sleep(0.05)
    for _ in range(50):
        app.processEvents()
    return not dialog.thread.isRunning()


def test_cancel_ends_the_run_and_frees_the_dialog(app, dialog):
    dialog._start_processing()
    assert dialog.thread.isRunning()

    dialog._cancel_processing()

    assert _wait_until_stopped(app, dialog)
    assert dialog.start_btn.isEnabled()
    assert not dialog.cancel_btn.isEnabled()
    assert dialog.close_btn.isEnabled()
    assert "Cancelled" in dialog.status_label.text()
    assert dialog.boxes == []


def test_a_run_offers_cancel_and_no_pause(app, dialog):
    """Nothing in the series read a pause; the button paused nothing."""
    from PySide6.QtWidgets import QPushButton
    dialog._start_processing()

    labels = [b.text() for b in dialog.findChildren(QPushButton) if b.isVisibleTo(dialog)]
    assert dialog.cancel_btn.isEnabled()
    assert not any("Pause" in label or "Resume" in label for label in labels)
    assert not hasattr(dialog.worker, "set_paused")


def test_the_dialog_stays_open_while_its_run_goes_on(app, dialog):
    dialog.show()
    dialog._start_processing()
    assert not dialog.close_btn.isEnabled()

    dialog.close()
    dialog.reject()

    assert dialog.isVisible()
    assert dialog.thread.isRunning()

    dialog._cancel_processing()
    assert _wait_until_stopped(app, dialog)
    dialog.close()
    assert not dialog.isVisible()
