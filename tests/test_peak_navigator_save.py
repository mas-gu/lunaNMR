#!/usr/bin/env python3
# ABOUTME: The Peak Navigator's save button writes whichever list it shows, reference peaks included.
# ABOUTME: Real LunaNMRMainWindow offscreen on the DNAJA1 0-eq spectrum, its peak list and a real detection.

import os
import sys
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


def _save(nav, target, monkeypatch):
    """Press save, answering the file dialog with `target`; returns the message boxes shown."""
    from PySide6.QtWidgets import QFileDialog, QMessageBox
    monkeypatch.setattr(QFileDialog, "getSaveFileName",
                        lambda *a, **k: (str(target), "Text files (*.txt)"))
    shown = []
    for box in ("information", "warning", "critical"):
        monkeypatch.setattr(QMessageBox, box, lambda *a, _box=box, **k: shown.append(_box))
    nav._save_peak_changes()
    return shown


def _rows(target):
    lines = target.read_text().splitlines()
    return lines[0], {line.split(",")[0]: line for line in lines[1:]}


def test_the_save_button_is_enabled_for_a_loaded_reference_list(window):
    nav = window.peak_navigator
    nav.peak_type_combo.setCurrentText("Reference Peaks")

    assert nav.selected_peak_type == "reference"
    assert len(nav.reference_peaks) > 0
    assert nav.save_btn.isEnabled()


def test_saving_reference_peaks_writes_every_reference_peak(window, tmp_path, monkeypatch):
    nav = window.peak_navigator
    nav.peak_type_combo.setCurrentText("Reference Peaks")
    target = tmp_path / "refs.txt"

    shown = _save(nav, target, monkeypatch)

    assert shown == ["information"]
    header, rows = _rows(target)
    assert header == "Assignment, Position_X, Position_Y, Height"
    assert len(rows) == len(nav.reference_peaks)
    for name, x, y in (peak[:3] for peak in nav.reference_peaks):
        assert rows[str(name)].startswith(f"{name}, {x:.6f}, {y:.6f}")


def test_saving_detected_peaks_writes_every_detected_peak(window, tmp_path, monkeypatch):
    _detect(window)
    nav = window.peak_navigator
    assert nav.selected_peak_type == "detected"
    assert len(nav.detected_peaks) > 0
    target = tmp_path / "det.txt"

    shown = _save(nav, target, monkeypatch)

    assert shown == ["information"]
    header, rows = _rows(target)
    assert header == "Assignment, Position_X, Position_Y, Height"
    assert len(rows) == len(nav.detected_peaks)
    for name, x, y in (peak[:3] for peak in nav.detected_peaks):
        assert rows[str(name)].startswith(f"{name}, {x:.6f}, {y:.6f}")
