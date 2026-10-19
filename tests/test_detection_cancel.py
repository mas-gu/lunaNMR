# ABOUTME: A stop event ends a peak detection between candidate peaks, in both detection modes.
# ABOUTME: Real DNAJA1 0-eq spectrum and its peak list; an unset event leaves every detected peak unchanged.

"""
Cancel on the detection's progress dialog only discarded the result once the detection
had run to its end. Detection now takes the same stop event a fit does, checked before
each candidate peak it examines and each reference peak it matches; a stopped detection
raises Cancelled and leaves the integrator's peaks as they were.
"""

import contextlib
import io
import threading
import time
from pathlib import Path

import numpy as np
import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
SPECTRUM = "2D_DNAJA1_0o0.ft"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")

MODES = pytest.mark.parametrize("detect", ["process_peaks", "detect_peaks_sn_native"],
                                ids=["peak-list", "sn-threshold"])


@pytest.fixture(scope="module")
def integrator():
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.utils.file_manager import NMRFileManager
    set_ps2d_config('15N')
    it = EnhancedVoigtIntegrator()
    with contextlib.redirect_stdout(io.StringIO()):
        it.load_nmr_file(str(DNAJA1 / SPECTRUM))
    it.peak_list = NMRFileManager().load_peak_list(str(DNAJA1 / "DNAJA1.txt"))
    it.set_search_window(0.070, 0.070)
    return it


def _fingerprint(peaks):
    return [(p.get('assignment'), p.get('ppm_x'), p.get('ppm_y'), p.get('intensity'),
             p.get('snr'), p.get('detected')) for p in peaks]


@MODES
def test_an_unset_stop_changes_no_detected_peak(integrator, detect):
    plain = getattr(integrator, detect)()
    watched = getattr(integrator, detect)(stop=threading.Event())

    assert len(plain) > 0
    assert _fingerprint(plain) == _fingerprint(watched)


@MODES
def test_a_stop_set_before_detection_detects_nothing(integrator, detect):
    from lunaNMR.utils.cancellation import Cancelled
    stop = threading.Event()
    stop.set()
    before = integrator.fitted_peaks = [{'assignment': 'kept'}]

    with pytest.raises(Cancelled):
        getattr(integrator, detect)(stop=stop)

    assert integrator.fitted_peaks is before


@MODES
def test_a_stop_set_during_detection_ends_it_at_once(integrator, detect):
    """The stop is set a third of the way into a detection that takes seconds."""
    from lunaNMR.utils.cancellation import Cancelled
    start = time.monotonic()
    getattr(integrator, detect)()
    duration = time.monotonic() - start

    stop = threading.Event()
    set_at = []

    def set_stop():
        set_at.append(time.monotonic())
        stop.set()

    timer = threading.Timer(duration / 3, set_stop)
    timer.start()
    try:
        with pytest.raises(Cancelled):
            getattr(integrator, detect)(stop=stop)
        landed = time.monotonic() - set_at[0]
    finally:
        timer.cancel()

    # Wall time is the point (it is what the user waits), so a loaded machine stretches it;
    # bounded by a share of a whole detection run under the same load as well. A Cancel
    # that waited for the detection to end lands at about 2/3 of it.
    bound = max(0.25, 0.2 * duration)
    assert landed < bound, f"Cancel landed {landed:.3f} s after the stop (detection {duration:.2f} s)"


class TestTheMainWindowsDetection:
    """Cancel on Detect Peaks' progress dialog sets the stop event its worker carries."""

    @pytest.fixture
    def window(self):
        import os
        import sys
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        from PySide6.QtWidgets import QApplication
        from lunaNMR.core.ps2d_config import set_ps2d_config
        from lunaNMR.gui.main_window import LunaNMRMainWindow
        app = QApplication.instance() or QApplication(sys.argv)
        set_ps2d_config('15N')
        with contextlib.redirect_stdout(io.StringIO()):
            w = LunaNMRMainWindow()
            w.load_nmr_file(str(DNAJA1 / SPECTRUM))
            w.load_peak_file(str(DNAJA1 / "DNAJA1.txt"))
        yield app, w
        worker = getattr(w, 'detection_worker', None)
        if worker is not None:
            worker.wait()
        app.processEvents()
        w.close()
        w.deleteLater()
        set_ps2d_config('15N')

    def test_cancel_ends_the_detection_and_keeps_the_previous_peaks(self, window, monkeypatch):
        from lunaNMR.gui.components import ProgressDialog
        app, w = window
        previous = w.integrator.fitted_peaks = [{'assignment': 'kept'}]
        cancelled_at = []

        def exec_then_cancel(dialog):
            time.sleep(0.3)   # the worker is detecting
            cancelled_at.append(time.monotonic())
            dialog._on_cancel_clicked()
            return 0
        monkeypatch.setattr(ProgressDialog, "exec", exec_then_cancel)

        w.on_detect_peaks()
        while w.detection_worker.isRunning():
            app.processEvents()
            w.detection_worker.wait(10)
        landed = time.monotonic() - cancelled_at[0]
        for _ in range(50):
            app.processEvents()

        assert landed < 0.5, f"Cancel landed {landed:.3f} s after the click"
        assert w.integrator.fitted_peaks is previous
        assert "cancelled" in w.statusBar().currentMessage().lower()
        assert (w.detect_peaks_button.isEnabled(), w.fit_all_peaks_button.isEnabled()) == (True, True)
