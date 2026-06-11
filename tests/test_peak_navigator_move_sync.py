#!/usr/bin/env python3
# ABOUTME: Regression test — Peak Navigator cache must track peak moves so saves persist them
# ABOUTME: Exercises PeakNavigator.update_peak_position for detected and reference peaks

import importlib.util
import os
import sys
import types
from pathlib import Path
from unittest.mock import patch

REPO_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def _install_namespace_packages():
    for dotted in ("lunaNMR", "lunaNMR.gui", "lunaNMR.gui.styles", "lunaNMR.gui.components"):
        if dotted in sys.modules:
            continue
        sub_path = REPO_ROOT / Path(*dotted.split("."))
        pkg = types.ModuleType(dotted)
        pkg.__path__ = [str(sub_path)]
        sys.modules[dotted] = pkg


def _load_module(dotted, file_path):
    spec = importlib.util.spec_from_file_location(dotted, file_path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[dotted] = mod
    spec.loader.exec_module(mod)
    return mod


_install_namespace_packages()
_load_module(
    "lunaNMR.gui.styles.design_system",
    REPO_ROOT / "lunaNMR/gui/styles/design_system.py",
)
peak_navigator = _load_module(
    "lunaNMR.gui.components.peak_navigator",
    REPO_ROOT / "lunaNMR/gui/components/peak_navigator.py",
)
PeakNavigator = peak_navigator.PeakNavigator

from PySide6.QtWidgets import QApplication, QFileDialog, QMessageBox


def _app():
    return QApplication.instance() or QApplication(sys.argv)


def _read(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


def test_update_peak_position_detected_mutates_cache():
    _app()
    nav = PeakNavigator()
    nav.load_detected_peaks([
        {"assignment": "D1", "ppm_x": 8.50, "ppm_y": 121.0, "height": 1.0e5, "r_squared": 0.92},
        {"assignment": "D2", "ppm_x": 7.30, "ppm_y": 118.0, "height": 2.0e5, "r_squared": 0.88},
    ])
    # peak_id for detected peaks is the enumerate index at load time
    peak_id = nav.detected_peaks[1][4]
    assert peak_id == 1

    nav.update_peak_position("detected", peak_id=peak_id, assignment="D2",
                             new_x=7.42, new_y=118.7)

    assert nav.detected_peaks[1][1] == 7.42
    assert nav.detected_peaks[1][2] == 118.7
    # untouched entry stays put
    assert nav.detected_peaks[0][1] == 8.50
    assert nav.detected_peaks[0][2] == 121.0


def test_update_peak_position_reference_mutates_cache():
    _app()
    nav = PeakNavigator()
    nav.load_reference_peaks([
        ["A1", 8.12, 120.5, 1.0e5],
        ["A2", 7.34, 118.2, 2.5e5],
    ])

    nav.update_peak_position("reference", peak_id=None, assignment="A2",
                             new_x=7.40, new_y=118.9)

    assert nav.reference_peaks[1][1] == 7.40
    assert nav.reference_peaks[1][2] == 118.9
    assert nav.reference_peaks[0][1] == 8.12
    assert nav.reference_peaks[0][2] == 120.5


def test_save_after_move_writes_updated_position(tmp_path):
    """End-to-end: move a reference peak, save, file should contain the new coords."""
    _app()
    nav = PeakNavigator()
    nav.load_reference_peaks([
        ["A1", 8.12, 120.5, 1.0e5],
        ["A2", 7.34, 118.2, 2.5e5],
    ])
    nav.peak_type_combo.setCurrentText("Reference Peaks")

    nav.update_peak_position("reference", peak_id=None, assignment="A2",
                             new_x=7.40, new_y=118.9)

    target = tmp_path / "refs.txt"
    with patch.object(
        QFileDialog, "getSaveFileName",
        return_value=(str(target), "Text files (*.txt)"),
    ), patch.object(QMessageBox, "information"), patch.object(QMessageBox, "warning"), \
            patch.object(QMessageBox, "critical"):
        nav._save_peak_changes()

    content = _read(target)
    assert "7.400000" in content, "moved X coordinate must be in saved file"
    assert "118.900000" in content, "moved Y coordinate must be in saved file"
    # Old coords must not leak back in for A2
    assert "7.340000" not in content
    assert "118.200000" not in content


if __name__ == "__main__":
    import tempfile
    test_update_peak_position_detected_mutates_cache()
    print("OK: test_update_peak_position_detected_mutates_cache")
    test_update_peak_position_reference_mutates_cache()
    print("OK: test_update_peak_position_reference_mutates_cache")
    with tempfile.TemporaryDirectory() as d:
        test_save_after_move_writes_updated_position(Path(d))
        print("OK: test_save_after_move_writes_updated_position")
    print("All tests passed")
