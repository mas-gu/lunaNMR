# ABOUTME: The series dialog reports a run as the CLI summary does: failures, failed spectra, fallbacks.
# ABOUTME: Drives the real dialog and worker offscreen on DNAJA1 spectra linked into a scratch folder.

"""
The dialog called every run 'complete' and saved it as a series, even one that fitted
nothing; it logged only 'Loading' lines, so the processors' warnings and errors, the
spectra that produced no fits and the fits that fell back to another algorithm never
reached the user. These tests run the dialog the way Start does and read what it shows.
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

FIRST = "2D_DNAJA1_0o0.ft"
SECOND = "2D_DNAJA1_HSPA8_0o125.ft"
TITRATION = 2


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
        integrator.load_nmr_file(str(DNAJA1 / FIRST))
        integrator.peak_list = NMRFileManager().load_peak_list(str(DNAJA1 / "DNAJA1.txt"))
        integrator.set_search_window(0.07, 0.07)
        return integrator._detect_peaks_reference_based()


def _main_window(app, folder, fitted_peaks, **settings):
    from PySide6.QtWidgets import QDoubleSpinBox
    contour = QDoubleSpinBox()
    contour.setDecimals(5)
    contour.setValue(0.05)
    return types.SimpleNamespace(
        current_nmr_folder=str(folder), search_window_x=0.07, search_window_y=0.07,
        use_parallel_processing=False, saved_series={}, contour_min_spin=contour,
        integrator=types.SimpleNamespace(fitted_peaks=fitted_peaks), **settings)


def _run(dialog, app):
    """Start the series as the Start button does and wait for its thread to finish."""
    from numba.core.errors import NumbaWarning
    dialog._prompt_series_name = lambda: "x"   # the user's answer to the name prompt
    dialog.detected_radio.setChecked(True)     # the Fit All Peaks result is the peak source
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore", NumbaWarning)
        dialog._start_processing()
        while dialog.thread.isRunning():
            app.processEvents()
            dialog.thread.wait(50)
        app.processEvents()
    box = dialog.outcome_box
    text = box.text()
    box.close()
    return text


def _dialog(app, main_window):
    from lunaNMR.gui.dialogs.series_integration_dialog import SeriesIntegrationDialog
    return SeriesIntegrationDialog(main_window=main_window)


@pytest.fixture(scope="module")
def failed_run(app, fitted_peaks, tmp_path_factory):
    """A run that fits nothing: neither of its two files is a spectrum."""
    folder = tmp_path_factory.mktemp("nothing_fitted")
    for name in (FIRST, SECOND):
        (folder / name).write_bytes(b"not a spectrum")
    main_window = _main_window(app, folder, fitted_peaks)
    dialog = _dialog(app, main_window)
    text = _run(dialog, app)
    yield dialog, main_window, text
    dialog.deleteLater()


class TestARunThatFitsNothing:

    def test_is_not_saved_as_a_series(self, failed_run):
        dialog, main_window, _ = failed_run
        assert main_window.saved_series == {}
        assert not hasattr(main_window, "batch_results")
        assert not dialog.overlay_viewer_btn.isEnabled()

    def test_is_reported_as_a_failure(self, failed_run):
        dialog, _, text = failed_run
        assert "complete" not in text.lower()
        assert "no spectrum" in text.lower()
        assert "fail" in dialog.status_label.text().lower()

    def test_the_log_carries_the_engines_error_and_names_each_spectrum(self, failed_run):
        log = failed_run[0].log_text.toPlainText()
        assert f"Failed to process {FIRST}: All loading methods failed" in log
        assert f"{FIRST} produced no fits" in log
        assert f"{SECOND} produced no fits" in log

    def test_an_empty_result_from_a_top_level_exception_is_a_failure(self, app, fitted_peaks, tmp_path):
        """process_nmr_series catches its own exception and returns it in .errors."""
        from lunaNMR.processors.multi_spectrum_processor import BatchResults
        main_window = _main_window(app, tmp_path, fitted_peaks)
        dialog = _dialog(app, main_window)
        result = BatchResults()
        result.errors.append({'error': 'boom'})
        dialog._on_finished(result)
        dialog.outcome_box.close()
        assert main_window.saved_series == {}
        assert "boom" in dialog.log_text.toPlainText()
        assert "fail" in dialog.status_label.text().lower()
        dialog.deleteLater()


@pytest.fixture(scope="module")
def partial_run(app, fitted_peaks, tmp_path_factory):
    """One spectrum fits, a second that is no spectrum at all fails."""
    folder = tmp_path_factory.mktemp("partial")
    (folder / FIRST).symlink_to(DNAJA1 / FIRST)
    (folder / "2D_unreadable.ft").write_bytes(b"not a spectrum")
    main_window = _main_window(app, folder, fitted_peaks)
    dialog = _dialog(app, main_window)
    text = _run(dialog, app)
    yield dialog, main_window, text
    dialog.deleteLater()


class TestARunWithAFailedSpectrum:

    def test_is_saved(self, partial_run):
        _, main_window, _ = partial_run
        assert list(main_window.saved_series) == ["x"]

    def test_names_the_spectrum_that_failed(self, partial_run):
        dialog, _, text = partial_run
        assert "2D_unreadable.ft produced no fits" in dialog.log_text.toPlainText()
        assert "2D_unreadable.ft" in text

    def test_states_the_peak_tracking_that_ran(self, partial_run):
        assert "Peak tracking: detected" in partial_run[0].log_text.toPlainText()


def test_a_titration_states_the_cascade_it_runs_in_place_of_detected(app, fitted_peaks, tmp_path):
    main_window = _main_window(app, tmp_path, fitted_peaks)
    dialog = _dialog(app, main_window)
    dialog.detected_radio.setChecked(True)
    dialog.xaxis_mode_combo.setCurrentIndex(TITRATION)
    assert dialog._tracking_statement("detected") == (
        "Peak tracking: cascade (requested 'detected'), drift limit off")
    dialog.deleteLater()


def test_a_spectrum_fitted_by_a_fallback_is_named_with_its_reason(app, fitted_peaks, tmp_path):
    from lunaNMR.processors.multi_spectrum_processor import BatchResults
    main_window = _main_window(app, tmp_path, fitted_peaks)
    dialog = _dialog(app, main_window)
    result = BatchResults()
    result.add_result(FIRST, {'status': 'success', 'fitting_fallbacks': []})
    result.add_result(SECOND, {'status': 'success',
                               'fitting_fallbacks': ['per-peak fit: cluster pipeline raised']})
    dialog._on_finished(result)
    text = dialog.outcome_box.text()
    dialog.outcome_box.close()
    log = dialog.log_text.toPlainText()
    assert f"{SECOND} was fitted by a fallback: per-peak fit: cluster pipeline raised" in log
    assert f"{FIRST} was fitted by a fallback" not in log
    assert SECOND in text
    dialog.deleteLater()


def test_a_warning_from_the_run_reaches_the_log(app, fitted_peaks, tmp_path):
    """Every log_warning / log_error in the series thread reaches the dialog's log."""
    from lunaNMR.gui.dialogs.series_integration_dialog import ProcessingWorker
    from lunaNMR.utils.output_manager import OutputManager, log_error, log_warning
    dialog = _dialog(app, _main_window(app, tmp_path, fitted_peaks))
    worker = ProcessingWorker([], None, "detected", {})
    worker.notice.connect(dialog._on_notice)
    worker.progress.connect(dialog._on_progress)
    OutputManager.set_callback(worker.progress_callback)
    try:
        log_warning("Falling back to per-peak fitting")
        log_error("Failed to process 2D_x.ft")
    finally:
        OutputManager.set_callback(None)
    log = dialog.log_text.toPlainText()
    assert "Warning: Falling back to per-peak fitting" in log
    assert "✕  Failed to process 2D_x.ft" in log
    dialog.deleteLater()


def test_the_log_names_no_noise_threshold(partial_run):
    """Neither detection nor the fit reads a noise threshold the dialog could state."""
    assert "Noise Threshold" not in partial_run[0].log_text.toPlainText()


def test_the_expert_panel_offers_no_noise_threshold(app):
    """Its spinner reached set_threshold_multiplier only, whose value no detection, fit or
    contour reads: Detect Peaks found the same peaks at 3 and at 9."""
    from PySide6.QtCore import QTimer
    from PySide6.QtWidgets import QLabel
    from lunaNMR.gui.main_window import LunaNMRMainWindow

    with contextlib.redirect_stdout(io.StringIO()):
        window = LunaNMRMainWindow()
        QTimer.singleShot(0, lambda: window._expert_dialog.reject())
        window.on_expert_mode()   # the Expert Mode dialog, answered by the timer above
    labels = [label.text() for label in window._expert_dialog.findChildren(QLabel)]
    assert not any("Noise Threshold" in text for text in labels)
    assert not hasattr(window, "expert_noise_spin")
    assert not hasattr(window, "noise_threshold")
    window.close()
    window.deleteLater()
