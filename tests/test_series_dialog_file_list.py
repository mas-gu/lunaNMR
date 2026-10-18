# ABOUTME: The series dialog's X-axis mode orders its spectra; it never drops or loses one.
# ABOUTME: A spectrum with no value in its name is kept after the rest, as the CLI orders it.

"""
Choosing Time series or Titration point replaced the dialog's file list with the spectra
whose names carry a value, so a reference named without one was left out of the series
without the '(N files without ...)' note written for it, and switching back to Off did
not bring it back.
"""

import os
import sys
import types
from pathlib import Path

import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")

OFF, TIME, TITRATION = 0, 1, 2
EXTRA = "2D_DNAJA1_HSPA8_extra.ft"


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def dialog(app, tmp_path):
    from lunaNMR.gui.dialogs.series_integration_dialog import SeriesIntegrationDialog
    for spectrum in DNAJA1.glob("*.ft"):
        (tmp_path / spectrum.name).symlink_to(spectrum)
    (tmp_path / EXTRA).symlink_to(DNAJA1 / "2D_DNAJA1_0o0.ft")
    dlg = SeriesIntegrationDialog(main_window=types.SimpleNamespace(current_nmr_folder=str(tmp_path)))
    yield dlg
    dlg.deleteLater()


def _names(dialog):
    return [os.path.basename(f) for f in dialog.nmr_files]


def test_a_titration_keeps_the_spectrum_with_no_point_last(dialog):
    dialog.xaxis_mode_combo.setCurrentIndex(TITRATION)
    assert len(dialog.nmr_files) == 8
    assert _names(dialog)[-1] == EXTRA
    assert _names(dialog)[:2] == ["2D_DNAJA1_0o0.ft", "2D_DNAJA1_HSPA8_0o125.ft"]
    assert dialog.file_list.count() == 8


def test_the_spectrum_with_no_point_is_counted_in_the_note(dialog):
    dialog.xaxis_mode_combo.setCurrentIndex(TITRATION)
    assert "(1 files without titration points)" in dialog.delay_info_label.text()


@pytest.mark.parametrize("modes", [(TITRATION, OFF), (TIME, OFF), (TIME, TITRATION, OFF)])
def test_switching_modes_loses_no_spectrum(dialog, modes):
    for index in modes:
        dialog.xaxis_mode_combo.setCurrentIndex(index)
        assert len(dialog.nmr_files) == 8, index
