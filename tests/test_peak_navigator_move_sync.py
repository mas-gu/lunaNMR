#!/usr/bin/env python3
# ABOUTME: A peak moved in the main window's plot moves in the Peak Navigator's list, and a save keeps it.
# ABOUTME: Real LunaNMRMainWindow offscreen on the DNAJA1 0-eq spectrum, its peak list and a real detection.

import os
import sys
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
SPECTRUM = "2D_DNAJA1_0o0.ft"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")

MOVE_X, MOVE_Y = 0.05, 0.4


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
    w.load_nmr_file(str(DNAJA1 / SPECTRUM))
    w.load_peak_file(str(DNAJA1 / "DNAJA1.txt"))
    yield w
    app.processEvents()
    w.close()
    w.deleteLater()
    set_ps2d_config('15N')


def _detect(window):
    """Run the detection worker in this thread and hand its result to the window's slot."""
    from lunaNMR.gui.main_window import PeakDetectionWorker
    window.sync_parameters_to_integrator()
    worker = PeakDetectionWorker(window.integrator, window.workflow_mode)
    worker.detection_complete.connect(lambda peaks: window.on_detection_complete(peaks, None))
    worker.run()


def _move(window, x, y, new_x, new_y):
    """Middle-click on the peak at (x, y), then middle-click where it belongs."""
    window._on_peak_select_for_move(x, y)
    window._on_peak_select_for_move(new_x, new_y)


def _entry(peaks, assignment):
    return next(p for p in peaks if str(p[0]) == assignment)


def test_a_moved_reference_peak_moves_in_the_navigator(window):
    nav = window.peak_navigator
    moved, other = nav.reference_peaks[1], nav.reference_peaks[0]
    name, x, y = moved[0], moved[1], moved[2]
    other_before = (other[1], other[2])

    _move(window, x, y, x + MOVE_X, y + MOVE_Y)

    assert _entry(nav.reference_peaks, name)[1:3] == [x + MOVE_X, y + MOVE_Y]
    assert (other[1], other[2]) == other_before
    row = window.integrator.peak_list[window.integrator.peak_list['Assignment'] == name].iloc[0]
    assert (row['Position_X'], row['Position_Y']) == (x + MOVE_X, y + MOVE_Y)


def test_a_moved_detected_peak_moves_in_the_navigator(window):
    _detect(window)
    nav = window.peak_navigator
    assert len(nav.detected_peaks) > 1
    moved, other = nav.detected_peaks[1], nav.detected_peaks[0]
    x, y = moved[1], moved[2]
    other_before = (other[1], other[2])
    # The middle-click finds the nearest of every peak, so the reference list must not
    # hold one closer to the click than the detected peak itself.
    window.integrator.peak_list = None

    _move(window, x, y, x + MOVE_X, y + MOVE_Y)

    assert (moved[1], moved[2]) == (x + MOVE_X, y + MOVE_Y)
    assert (other[1], other[2]) == other_before


def test_a_save_after_a_move_writes_the_moved_position(window, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QFileDialog, QMessageBox
    nav = window.peak_navigator
    nav.peak_type_combo.setCurrentText("Reference Peaks")
    name, x, y = nav.reference_peaks[1][:3]
    _move(window, x, y, x + MOVE_X, y + MOVE_Y)

    target = tmp_path / "refs.txt"
    monkeypatch.setattr(QFileDialog, "getSaveFileName",
                        lambda *a, **k: (str(target), "Text files (*.txt)"))
    shown = []
    for box in ("information", "warning", "critical"):
        monkeypatch.setattr(QMessageBox, box, lambda *a, _box=box, **k: shown.append(_box))
    nav._save_peak_changes()

    assert shown == ["information"]
    rows = {line.split(",")[0]: line for line in target.read_text().splitlines()[1:]}
    assert len(rows) == len(nav.reference_peaks)
    assert rows[name].startswith(f"{name}, {x + MOVE_X:.6f}, {y + MOVE_Y:.6f}")
