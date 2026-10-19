# ABOUTME: The Detect/Fit progress dialog ends on what the run produced: peaks detected, or fits
# ABOUTME: that succeeded and failed. Real LunaNMRMainWindow offscreen on DNAJA1 and a 10-peak list.

"""
Its Successful/Failed counters counted LOG LINES, not peaks: a detection logs none, so it
ended on "Successful: 0 | Failed: 0" whatever it found, and a fit counted the messages of
its second pass.
"""

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
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.dialogs import DataLoadingDialog
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    set_ps2d_config('15N')
    spectra = tmp_path / "spectra"
    spectra.mkdir()
    (spectra / SPECTRUM).symlink_to(DNAJA1 / SPECTRUM)
    lines = (DNAJA1 / "DNAJA1.txt").read_text().splitlines()
    peaks = tmp_path / "ten.txt"
    peaks.write_text("\n".join(lines[:11]) + "\n")
    w = LunaNMRMainWindow()
    monkeypatch.setattr(DataLoadingDialog, "exec", lambda self: 1)
    monkeypatch.setattr(DataLoadingDialog, "get_result", lambda self: {
        'nmr_file': str(spectra / SPECTRUM), 'peak_file': str(peaks),
        'nmr_folder': str(spectra), 'peak_folder': str(tmp_path)})
    w.on_load_data()
    yield w
    for worker in (getattr(w, 'fitting_worker', None), getattr(w, 'detection_worker', None)):
        if worker is not None:
            worker.wait()
    app.processEvents()
    w.close()
    w.deleteLater()
    set_ps2d_config('15N')


def _run(app, window, monkeypatch, start, worker_name):
    """Start a run from its button, keep the dialog it opens, and let it finish."""
    from lunaNMR.gui.components import ProgressDialog
    dialogs = []
    monkeypatch.setattr(ProgressDialog, "exec", lambda self: dialogs.append(self) or 0)
    start()
    worker = getattr(window, worker_name)
    deadline = time.time() + 600
    while worker.isRunning() and time.time() < deadline:
        app.processEvents()
        time.sleep(0.05)
    for _ in range(50):
        app.processEvents()
    return dialogs[0].stats_label.text()


def test_a_detection_ends_on_the_peaks_it_detected(app, window, monkeypatch):
    text = _run(app, window, monkeypatch, window.on_detect_peaks, "detection_worker")
    assert f"Detected: {len(window.integrator.fitted_peaks)} peaks" in text
    assert "Successful" not in text and "Failed" not in text


def test_a_fit_ends_on_the_fits_that_succeeded_and_failed(app, window, monkeypatch):
    text = _run(app, window, monkeypatch, window.on_fit_spectrum, "fitting_worker")
    results = window.last_fitting_results
    succeeded = sum(1 for r in results if r.get('success', True))
    assert f"Successful: {succeeded} | ❌ Failed: {len(results) - succeeded}" in text
    assert succeeded > 0
