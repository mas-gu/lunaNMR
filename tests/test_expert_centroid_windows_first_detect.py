# ABOUTME: The Expert centroid windows reach the detection of a session's first Detect.
# ABOUTME: Real LunaNMRMainWindow offscreen on the DNAJA1 spectrum, a 10-peak list cut from its peak list.

"""Detection refines each maximum to a centroid inside a window it reads from the
integrator's gui_params. A window that has never fitted must hand its integrator those
parameters all the same, or the first detection runs on the integrator's own fallback."""

import os
import sys
import time
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
SPECTRUM = "2D_DNAJA1_0o0.ft"


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def window(app):
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    set_ps2d_config('15N')
    w = LunaNMRMainWindow()
    yield w
    worker = getattr(w, 'detection_worker', None)
    if worker is not None:
        worker.wait()
    app.processEvents()
    w.close()
    w.deleteLater()
    set_ps2d_config('15N')


def _set_centroid_windows(window, monkeypatch, x_ppm, y_ppm):
    """The user opens Expert Mode and types the two centroid windows."""
    from PySide6.QtWidgets import QDialog
    monkeypatch.setattr(QDialog, "exec", lambda self: 0)
    window.on_expert_mode()
    window.expert_centroid_x_spin.setValue(x_ppm)
    window.expert_centroid_y_spin.setValue(y_ppm)


def test_a_window_that_never_fitted_hands_the_integrator_the_centroid_windows(window, monkeypatch):
    _set_centroid_windows(window, monkeypatch, 0.03, 0.2)

    window.sync_parameters_to_integrator()

    gui_params = window.integrator.gui_params
    assert (gui_params['centroid_window_x_ppm'], gui_params['centroid_window_y_ppm']) == (0.03, 0.2)


@pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")
def test_the_first_detect_of_a_session_refines_inside_the_expert_windows(app, window, tmp_path,
                                                                        monkeypatch):
    from lunaNMR.gui.components import ProgressDialog
    lines = (DNAJA1 / "DNAJA1.txt").read_text().splitlines()
    peaks = tmp_path / "ten.txt"
    peaks.write_text("\n".join(lines[:11]) + "\n")
    window.load_data(str(DNAJA1 / SPECTRUM), str(peaks))
    _set_centroid_windows(window, monkeypatch, 0.03, 0.2)
    refined_inside = set()
    refine = window.integrator._calculate_top_contour_center

    def recording(y_idx, x_idx, **kwargs):
        refined_inside.add((kwargs['max_shift_x'], kwargs['max_shift_y']))
        return refine(y_idx, x_idx, **kwargs)

    monkeypatch.setattr(window.integrator, "_calculate_top_contour_center", recording)
    monkeypatch.setattr(ProgressDialog, "exec", lambda self: 0)

    window.on_detect_peaks()
    start = time.time()
    while window.detection_worker.isRunning() and time.time() - start < 600:
        app.processEvents()
        time.sleep(0.05)
    for _ in range(50):
        app.processEvents()

    assert refined_inside == {(0.03, 0.2)}
