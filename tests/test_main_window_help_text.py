# ABOUTME: The User Guide and the series-complete message name only buttons, menus and shortcuts the window has.
# ABOUTME: Real LunaNMRMainWindow offscreen; the texts are read from the dialogs the window opens.

import os
import re
import sys

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


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
    w.close()
    w.deleteLater()
    set_ps2d_config('15N')


def _menus(window):
    """{menu title: [action texts]} of the menu bar, without mnemonics."""
    from PySide6.QtWidgets import QMenu
    return {menu.title().replace('&', ''): [action.text() for action in menu.actions() if action.text()]
            for menu in window.menuBar().findChildren(QMenu)}


def _user_guide(window, monkeypatch):
    """The text of Help > User Guide."""
    from PySide6.QtWidgets import QDialog, QMenu, QTextEdit
    shown = []
    monkeypatch.setattr(QDialog, "exec",
                        lambda self: shown.append(self.findChild(QTextEdit).toPlainText()))
    next(action for menu in window.menuBar().findChildren(QMenu)
         for action in menu.actions() if action.text() == "User Guide").trigger()
    return shown[0]


def test_the_guide_names_the_buttons_of_the_workflow(window, monkeypatch):
    guide = _user_guide(window, monkeypatch)
    for widget in (window.workflow_peak_list_radio, window.workflow_sn_threshold_radio,
                   window.detect_peaks_button, window.fit_all_peaks_button, window.series_button):
        assert widget.text() in guide, widget.text()
    assert "Load Data" in guide


def test_the_guide_names_no_mode_or_manager_the_window_lacks(window, monkeypatch):
    guide = _user_guide(window, monkeypatch)
    for absent in ("Full Detection", "In-Place", "Configuration", "preview"):
        assert absent not in guide, absent


def test_every_menu_the_guide_names_is_in_the_menu_bar(window, monkeypatch):
    guide = _user_guide(window, monkeypatch)
    named = re.findall(r"^• ([A-Za-z ]+) menu:", guide, flags=re.MULTILINE)
    assert named and set(named) <= set(_menus(window))


def test_every_shortcut_the_guide_names_is_one_the_window_handles(window, monkeypatch):
    guide = _user_guide(window, monkeypatch)
    action_shortcuts = {action.shortcut().toString()
                        for action in (window.undo_action, window.redo_action)}
    # Esc, Delete and Ctrl+D are handled in keyPressEvent
    handled = action_shortcuts | {"Esc", "Delete", "Ctrl+D"}
    named = set(re.findall(r"\b(?:Ctrl\+(?:Shift\+)?[A-Z]|F\d+|Esc|Delete)\b", guide))
    assert {"Ctrl+Z", "Ctrl+Shift+Z"} <= named
    assert named <= handled


def test_the_series_complete_message_names_a_menu_entry_that_exists(window, monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    from lunaNMR.processors.multi_spectrum_processor import BatchResults
    shown = []
    monkeypatch.setattr(QMessageBox, "information",
                        staticmethod(lambda parent, title, text, *a, **k: shown.append(text)))

    window._on_series_processing_complete(BatchResults())

    menu, entry = re.search(r"Use (.+?) > (.+?) to analyze results\.", shown[0]).groups()
    assert entry in _menus(window)[menu]
