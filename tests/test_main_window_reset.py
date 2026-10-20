# ABOUTME: The main window shows the contour settings it draws with, and Reset leaves it as a newly opened one.
# ABOUTME: Real LunaNMRMainWindow offscreen; compares its widgets and plain state against a fresh window.

"""The window restores no preferences: nothing in the application saves any but the recent
projects list, so a loader of window, display and contour keys only ever read its own defaults."""

import os
import sys

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def make_window(app):
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    windows = []

    def make():
        set_ps2d_config('15N')
        windows.append(LunaNMRMainWindow())
        return windows[-1]

    yield make
    for w in windows:
        w.close()
        w.deleteLater()
    app.processEvents()
    set_ps2d_config('15N')


def _snapshot(window):
    """Every spin box, check box, radio button and combo box the window holds, and its plain state."""
    from PySide6.QtWidgets import QAbstractSpinBox, QCheckBox, QComboBox, QRadioButton
    out = {}
    for name, obj in vars(window).items():
        if isinstance(obj, QAbstractSpinBox) and hasattr(obj, 'value'):
            out[name] = obj.value()
        elif isinstance(obj, (QCheckBox, QRadioButton)):
            out[name] = obj.isChecked()
        elif isinstance(obj, QComboBox):
            out[name] = obj.currentText()
        elif isinstance(obj, (int, float, str, bool)) and not name.startswith('_'):
            out['state.' + name] = obj
    return out


def test_the_contour_spin_boxes_show_what_the_contours_are_drawn_with(make_window):
    window = make_window()
    assert (window.contour_levels_spin.value(), window.contour_min_spin.value(),
            window.contour_inc_spin.value()) == (window.contour_levels, window.contour_min,
                                                 window.contour_increment)


def test_reset_leaves_the_window_as_a_newly_opened_one(make_window, monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    fresh = _snapshot(make_window())
    window = make_window()
    window.contour_min_spin.setValue(0.2)
    window.contour_inc_spin.setValue(1.5)
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.Yes))

    window.on_reset_results()

    after = _snapshot(window)
    differ = {k: (fresh.get(k), after.get(k)) for k in fresh.keys() | after.keys()
              if fresh.get(k) != after.get(k)}
    assert differ == {}


def test_the_window_restores_no_preferences_nothing_saves():
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    assert not hasattr(LunaNMRMainWindow, 'apply_user_preferences')
