#!/usr/bin/env python3
# ABOUTME: Regression test for Peak Navigator save button — must save reference peaks too
# ABOUTME: Exercises _save_peak_changes for both reference and detected peak types

import importlib.util
import os
import sys
import types
from pathlib import Path
from unittest.mock import patch

REPO_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


def _install_namespace_packages():
    """Register lunaNMR.* as namespace packages to avoid heavy __init__ side effects."""
    for dotted in ("lunaNMR", "lunaNMR.gui", "lunaNMR.gui.styles", "lunaNMR.gui.components"):
        if dotted in sys.modules:
            continue
        sub_path = REPO_ROOT / Path(*dotted.split("."))
        pkg = types.ModuleType(dotted)
        pkg.__path__ = [str(sub_path)]
        sys.modules[dotted] = pkg


def _load_module(dotted: str, file_path: Path):
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


def test_save_reference_peaks_writes_file(tmp_path):
    _app()
    nav = PeakNavigator()
    nav.load_reference_peaks([
        ["A1", 8.12, 120.5, 1.0e5],
        ["A2", 7.34, 118.2, 2.5e5],
    ])
    nav.peak_type_combo.setCurrentText("Reference Peaks")
    assert nav.selected_peak_type == "reference"

    target = tmp_path / "refs.txt"
    with patch.object(
        QFileDialog, "getSaveFileName",
        return_value=(str(target), "Text files (*.txt)"),
    ), patch.object(QMessageBox, "information"), patch.object(QMessageBox, "warning"), \
            patch.object(QMessageBox, "critical"):
        nav._save_peak_changes()

    assert target.exists(), "reference peak save must produce a file"
    content = _read(target)
    assert "Assignment" in content
    assert "A1" in content and "A2" in content
    assert "8.120000" in content
    assert "120.500000" in content


def test_save_detected_peaks_still_works(tmp_path):
    _app()
    nav = PeakNavigator()
    nav.load_detected_peaks([
        {"assignment": "D1", "ppm_x": 8.5, "ppm_y": 121.0, "height": 5.0e5, "r_squared": 0.92},
    ])
    assert nav.selected_peak_type == "detected"

    target = tmp_path / "det.txt"
    with patch.object(
        QFileDialog, "getSaveFileName",
        return_value=(str(target), "Text files (*.txt)"),
    ), patch.object(QMessageBox, "information"), patch.object(QMessageBox, "warning"), \
            patch.object(QMessageBox, "critical"):
        nav._save_peak_changes()

    assert target.exists()
    content = _read(target)
    assert "D1" in content
    assert "8.500000" in content


def test_save_button_enabled_for_reference():
    _app()
    nav = PeakNavigator()
    nav.load_reference_peaks([["A1", 8.0, 120.0, 1.0]])
    nav.peak_type_combo.setCurrentText("Reference Peaks")
    nav._update_edit_button_states()
    assert nav.save_btn.isEnabled(), "save button must be enabled when reference peaks are loaded"


if __name__ == "__main__":
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        test_save_reference_peaks_writes_file(Path(d))
        print("OK: test_save_reference_peaks_writes_file")
    with tempfile.TemporaryDirectory() as d:
        test_save_detected_peaks_still_works(Path(d))
        print("OK: test_save_detected_peaks_still_works")
    test_save_button_enabled_for_reference()
    print("OK: test_save_button_enabled_for_reference")
    print("All tests passed")
