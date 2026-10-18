# ABOUTME: The main window's Detect/Fit progress dialog offers Cancel and no Pause.
# ABOUTME: Real ProgressDialog offscreen; nothing in detection or fitting reads a pause.
import os
import sys

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def dialog(app):
    from lunaNMR.gui.components.progress_dialog import ProgressDialog
    dlg = ProgressDialog(title="Detecting peaks")
    dlg.show()
    yield dlg
    dlg.close()
    dlg.deleteLater()


def test_offers_cancel_and_no_pause(dialog):
    from PySide6.QtWidgets import QPushButton
    labels = [b.text() for b in dialog.findChildren(QPushButton) if b.isVisibleTo(dialog)]
    assert any("Cancel" in label for label in labels)
    assert not any("Pause" in label or "Resume" in label for label in labels)


def test_has_no_pause_state(dialog):
    assert not hasattr(dialog, "paused")
    assert not hasattr(dialog, "is_paused")


def test_cancel_still_cancels(dialog):
    received = []
    dialog.cancelled.connect(lambda: received.append(True))
    dialog.cancel_button.click()
    assert received == [True]
    assert dialog.is_cancelled


def test_completion_leaves_close_enabled(dialog):
    dialog.complete("done")
    assert not dialog.cancel_button.isEnabled()
    assert dialog.close_button.isEnabled()
