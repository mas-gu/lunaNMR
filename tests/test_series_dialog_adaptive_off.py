# ABOUTME: A series run from the dialog never runs the adaptive radius search, whatever the Expert checkbox says.
# ABOUTME: Drives the real dialog and worker offscreen on DNAJA1 spectra linked into a scratch folder.

"""
The adaptive search chose a window that cuts the peak at half height on every real
titration it was tried on, so a series runs without it. The Expert checkbox still applies
to Fit All Peaks; the series dialog says in its log that the series does not use it, and
no tooltip promises it.
"""

import contextlib
import io
import sys
import types
import warnings
from pathlib import Path

import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")

SPECTRA = ("2D_DNAJA1_0o0.ft", "2D_DNAJA1_HSPA8_0o125.ft")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture(scope="module")
def fitted_peaks():
    """What Fit All Peaks leaves in integrator.fitted_peaks on the first spectrum."""
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.utils.file_manager import NMRFileManager
    with contextlib.redirect_stdout(io.StringIO()):
        integrator = EnhancedVoigtIntegrator()
        integrator.load_nmr_file(str(DNAJA1 / SPECTRA[0]))
        integrator.peak_list = NMRFileManager().load_peak_list(str(DNAJA1 / "DNAJA1.txt"))
        integrator.set_search_window(0.07, 0.07)
        return integrator._detect_peaks_reference_based()


def _main_window(folder, fitted_peaks):
    from PySide6.QtWidgets import QDoubleSpinBox
    contour = QDoubleSpinBox()
    contour.setDecimals(5)
    contour.setValue(0.05)
    return types.SimpleNamespace(
        current_nmr_folder=str(folder), search_window_x=0.07, search_window_y=0.07,
        use_parallel_processing=False, saved_series={}, contour_min_spin=contour,
        use_adaptive_optimization=True,
        integrator=types.SimpleNamespace(fitted_peaks=fitted_peaks))


@pytest.fixture(scope="module")
def adaptive_ticked_run(app, fitted_peaks, tmp_path_factory):
    """A series started with the Expert 'Adaptive optimization' box ticked."""
    from numba.core.errors import NumbaWarning
    from lunaNMR.gui.dialogs.series_integration_dialog import SeriesIntegrationDialog
    folder = tmp_path_factory.mktemp("adaptive_ticked")
    for name in SPECTRA:
        (folder / name).symlink_to(DNAJA1 / name)
    dialog = SeriesIntegrationDialog(main_window=_main_window(folder, fitted_peaks))
    dialog._prompt_series_name = lambda: "x"   # the user's answer to the name prompt
    dialog.detected_radio.setChecked(True)
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore", NumbaWarning)
        dialog._start_processing()
        while dialog.thread.isRunning():
            app.processEvents()
            dialog.thread.wait(50)
        app.processEvents()
    dialog.outcome_box.close()
    yield dialog
    dialog.deleteLater()


def test_the_series_parameters_turn_it_off(app, fitted_peaks, tmp_path):
    from lunaNMR.gui.dialogs.series_integration_dialog import SeriesIntegrationDialog
    dialog = SeriesIntegrationDialog(main_window=_main_window(tmp_path, fitted_peaks))
    for radio in (dialog.cascade_radio, dialog.detected_radio, dialog.independent_radio):
        radio.setChecked(True)
        assert dialog._get_voigt_params()['gui_params']['use_adaptive_optimization'] is False
    dialog.deleteLater()


def test_the_search_does_not_run(adaptive_ticked_run):
    processor = adaptive_ticked_run.worker.processor
    assert processor.integrator.enhanced_fitter.parallel_processor.optimal_params is None


def test_the_log_says_so(adaptive_ticked_run):
    assert "Adaptive radius search: off for a series" in adaptive_ticked_run.log_text.toPlainText()


def test_no_tooltip_promises_it(app, fitted_peaks, tmp_path):
    from lunaNMR.gui.dialogs.series_integration_dialog import SeriesIntegrationDialog
    dialog = SeriesIntegrationDialog(main_window=_main_window(tmp_path, fitted_peaks))
    for radio in (dialog.cascade_radio, dialog.detected_radio, dialog.independent_radio):
        assert "adaptive" not in radio.toolTip().lower()
    dialog.deleteLater()
