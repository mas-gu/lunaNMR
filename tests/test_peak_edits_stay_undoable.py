# ABOUTME: Peak edits stay undoable from the Edit menu after Reset Results, and a peak added in the navigator has an id.
# ABOUTME: Real LunaNMRMainWindow offscreen with a peak list read from a file; the Edit menu's own actions are triggered.

import os
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
    app.processEvents()
    set_ps2d_config('15N')


def _load_two_peaks(window, tmp_path):
    peaks = tmp_path / "peaks.csv"
    peaks.write_text("Assignment,Position_X,Position_Y\n1AlaH,8.10,120.5\n2GlyH,8.40,109.2\n")
    assert window.integrator.load_peak_list_file(str(peaks))


def test_the_edit_menu_undoes_a_peak_added_after_reset(window, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.Yes))
    window.on_reset_results()

    _load_two_peaks(window, tmp_path)
    window.add_new_peak(8.25, 115.0)
    assert len(window.integrator.peak_list) == 3
    assert window.undo_action.isEnabled()

    window.undo_action.trigger()

    assert list(window.integrator.peak_list['Assignment']) == ['1AlaH', '2GlyH']


def test_a_peak_added_in_the_navigator_has_an_id_of_its_own(window, tmp_path):
    _load_two_peaks(window, tmp_path)

    window._on_navigator_peak_added({'assignment': '3SerH', 'x_coord': 8.25, 'y_coord': 115.0})

    ids = window.integrator.peak_list['peak_id']
    assert ids.notna().all()
    assert ids.is_unique
