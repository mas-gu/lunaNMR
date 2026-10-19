# ABOUTME: A peak deleted in the spectrum view leaves the peak navigator too, and comes back on undo.
# ABOUTME: Real LunaNMRMainWindow offscreen on the DNAJA1 spectrum and a 10-peak list cut from its list.

"""
The deletion updated the integrator, but the navigator only redrew its own copy of the
peaks (`refresh_peak_list`), so the deleted peak stayed listed - and an undo restored the
integrator's peak without the navigator ever having lost it.
"""

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
def window(app, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.dialogs import DataLoadingDialog
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    set_ps2d_config('15N')
    spectra = tmp_path / "spectra"
    spectra.mkdir()
    (spectra / SPECTRUM).symlink_to(DNAJA1 / SPECTRUM)
    lines = (DNAJA1 / "DNAJA1.txt").read_text().splitlines()
    peaks = tmp_path / "ten.txt"
    peaks.write_text("\n".join(lines[:11]) + "\n")
    w = LunaNMRMainWindow()
    monkeypatch.setattr(DataLoadingDialog, "exec", lambda self: 1)
    monkeypatch.setattr(DataLoadingDialog, "get_result", lambda self: {
        'nmr_file': str(spectra / SPECTRUM), 'peak_file': str(peaks),
        'nmr_folder': str(spectra), 'peak_folder': str(tmp_path)})
    # The user confirms the deletion
    monkeypatch.setattr(QMessageBox, "question", lambda *a, **k: QMessageBox.Yes)
    w.on_load_data()
    yield w
    app.processEvents()
    w.close()
    w.deleteLater()
    set_ps2d_config('15N')


def _listed(navigator, kind):
    return [str(p[0]) for p in (navigator.reference_peaks if kind == "reference"
                                else navigator.detected_peaks)]


def _delete_nearest(window, kind):
    """Delete the peak a click lands on, as Ctrl+click does."""
    peaks = (window.integrator.peak_list.to_dict('records') if kind == "reference"
             else window.integrator.fitted_peaks)
    first = peaks[0]
    x = first.get('Position_X', first.get('ppm_x'))
    y = first.get('Position_Y', first.get('ppm_y'))
    # The user enables editing of that peak type, as the edit boxes do
    window.edit_reference_peaks = kind == "reference"
    window.edit_detected_peaks = kind == "detected"
    info = window.find_nearest_peak(float(x), float(y))
    assert info and info['type'] == kind
    window.delete_selected_peak(info)
    return str(info['assignment'])


def test_a_deleted_reference_peak_leaves_the_navigator(window):
    gone = _delete_nearest(window, "reference")
    assert gone not in _listed(window.peak_navigator, "reference")
    assert len(_listed(window.peak_navigator, "reference")) == 9


def test_undo_puts_it_back_in_the_navigator(window):
    gone = _delete_nearest(window, "reference")
    window.undo_manager.undo_stack.undo()
    assert gone in _listed(window.peak_navigator, "reference")
    assert len(_listed(window.peak_navigator, "reference")) == 10


def test_a_deleted_detected_peak_leaves_the_navigator(window):
    from lunaNMR.gui.main_window import PeakDetectionWorker
    window.sync_parameters_to_integrator()
    PeakDetectionWorker(window.integrator, window.workflow_mode).run()
    window.peak_navigator.load_detected_peaks(window.integrator.fitted_peaks)
    before = len(_listed(window.peak_navigator, "detected"))
    window.peak_navigator.selected_peak_type = "detected"
    gone = _delete_nearest(window, "detected")
    assert gone not in _listed(window.peak_navigator, "detected")
    assert len(_listed(window.peak_navigator, "detected")) == before - 1
