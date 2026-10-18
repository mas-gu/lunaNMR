# ABOUTME: The main window's state across runs: what a new load, a Reset and a series leave behind.
# ABOUTME: Real LunaNMRMainWindow offscreen on the DNAJA1 spectrum, a 10-peak list cut from its peak list.

"""
Fit Spectrum fits the detected peaks in preference to the loaded list, and a series locks
the clusters of the last fit. Loading new data kept both, so the next fit was of the
previous list. The PS2D config is one object for the whole process: a Reset left it on
the nucleus chosen before, and a series switches it to the nucleus it detected, behind
the window, which kept showing and fitting with its own.
"""

import os
import sys
import time
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


def _peak_list(path, rows):
    """A peak list of the given rows of the DNAJA1 list (its header kept)."""
    lines = (DNAJA1 / "DNAJA1.txt").read_text().splitlines()
    path.write_text("\n".join([lines[0]] + [lines[i] for i in rows]) + "\n")
    return path


@pytest.fixture
def data(tmp_path):
    """The spectrum linked into a folder of its own, so a fit's CSV lands there; two
    disjoint 10- and 5-peak lists cut from the DNAJA1 list."""
    spectra = tmp_path / "spectra"
    spectra.mkdir()
    (spectra / SPECTRUM).symlink_to(DNAJA1 / SPECTRUM)
    return {"spectrum": str(spectra / SPECTRUM), "folder": spectra,
            "ten": str(_peak_list(tmp_path / "ten.txt", range(1, 11))),
            "five": str(_peak_list(tmp_path / "five.txt", range(11, 16)))}


@pytest.fixture
def window(app):
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    set_ps2d_config('15N')
    w = LunaNMRMainWindow()
    yield w
    for worker in (getattr(w, 'fitting_worker', None), getattr(w, 'detection_worker', None)):
        if worker is not None:
            worker.wait()
    app.processEvents()
    w.close()
    w.deleteLater()
    set_ps2d_config('15N')


def _load(window, monkeypatch, nmr_file, peak_file):
    """The user picks the files in the Load Data dialog and clicks Load."""
    from lunaNMR.gui.dialogs import DataLoadingDialog
    folder = str(Path(nmr_file).parent)
    monkeypatch.setattr(DataLoadingDialog, "exec", lambda self: 1)
    monkeypatch.setattr(DataLoadingDialog, "get_result", lambda self: {
        'nmr_file': nmr_file, 'peak_file': peak_file,
        'nmr_folder': folder, 'peak_folder': str(Path(peak_file).parent)})
    window.on_load_data()


def _detect(window):
    from lunaNMR.gui.main_window import PeakDetectionWorker
    window.sync_parameters_to_integrator()
    PeakDetectionWorker(window.integrator, window.workflow_mode).run()


def _fit(window):
    """Run the Fit Spectrum worker to the end, in this thread."""
    from lunaNMR.gui.main_window import VoigtFittingWorker
    window.sync_parameters_to_integrator()
    VoigtFittingWorker(window.integrator, window.param_manager,
                       {'use_parallel': True, 'use_voigt_fitting': True,
                        'skip_series_params': True}).run()


def _wait_for(app, worker, timeout=600):
    start = time.time()
    while worker.isRunning() and time.time() - start < timeout:
        app.processEvents()
        time.sleep(0.05)
    assert not worker.isRunning()
    for _ in range(50):
        app.processEvents()


def _run_series(window, data):
    """Run the series dialog's worker, in this thread, over the 15N spectrum and the 10-peak list."""
    from lunaNMR.gui.dialogs.series_integration_dialog import ProcessingWorker, SeriesIntegrationDialog
    window.load_peak_file(data["ten"])
    dialog = SeriesIntegrationDialog(main_window=window)
    ProcessingWorker(None, [data["spectrum"]], window.integrator.peak_list.copy(), 'reference',
                     dialog._get_voigt_params(), series_name='nucleus').run()
    dialog.deleteLater()


def _assignments(peaks):
    return sorted(p.get('assignment') for p in peaks)


class TestLoadingNewDataStartsAFreshAnalysis:

    def test_the_previous_detection_and_clusters_are_forgotten(self, window, data, monkeypatch):
        _load(window, monkeypatch, data["spectrum"], data["ten"])
        _detect(window)
        _fit(window)
        assert len(window.integrator.fitted_peaks) == 10
        assert window.integrator.enhanced_fitter.series_params

        _load(window, monkeypatch, data["spectrum"], data["five"])

        assert window.integrator.fitted_peaks == []
        assert window.integrator.enhanced_fitter.series_params is None
        assert window.peak_navigator.detected_peaks == []

    def test_fit_spectrum_then_fits_the_list_just_loaded(self, app, window, data, monkeypatch):
        from lunaNMR.gui.components import ProgressDialog
        _load(window, monkeypatch, data["spectrum"], data["ten"])
        _detect(window)
        _load(window, monkeypatch, data["spectrum"], data["five"])
        monkeypatch.setattr(ProgressDialog, "exec", lambda self: 0)

        window.on_fit_spectrum()
        _wait_for(app, window.fitting_worker)

        five = sorted(window.integrator.peak_list['Assignment'])
        assert _assignments(window.integrator.fitted_peaks) == five

    def test_fit_spectrum_starts_without_the_last_fits_clusters(self, app, window, data, monkeypatch):
        from lunaNMR.gui.components import ProgressDialog
        _load(window, monkeypatch, data["spectrum"], data["ten"])
        _fit(window)
        assert window.integrator.enhanced_fitter.series_params
        at_start = []
        monkeypatch.setattr(ProgressDialog, "exec", lambda self: at_start.append(
            window.integrator.enhanced_fitter.series_params) or 0)

        window.on_fit_spectrum()
        _wait_for(app, window.fitting_worker)

        assert at_start == [None]


class TestThePs2dConfigFollowsTheWindow:

    def test_reset_puts_the_15n_defaults_back(self, window, monkeypatch):
        from PySide6.QtWidgets import QMessageBox
        from lunaNMR.core.ps2d_config import PS2DConfig, get_ps2d_config, restore_base_ps2d_config
        window._on_nucleus_type_change('13C')
        monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.Yes))

        window.on_reset_results()
        restore_base_ps2d_config()

        config = get_ps2d_config()
        assert (config.nucleus_type, config.radF1) == ('15N', PS2DConfig('15N').radF1)
        assert window.ps2d_radF1 == PS2DConfig('15N').radF1

    def test_loading_a_spectrum_applies_its_nucleus_when_the_config_holds_another(self, window, data,
                                                                                  monkeypatch):
        from lunaNMR.core.ps2d_config import PS2DConfig, get_ps2d_config, set_ps2d_config
        set_ps2d_config('13C')

        _load(window, monkeypatch, data["spectrum"], data["ten"])

        config = get_ps2d_config()
        assert (config.nucleus_type, config.radF1) == ('15N', PS2DConfig('15N').radF1)

    def test_reset_puts_the_15n_defaults_in_an_open_expert_panel(self, window, monkeypatch):
        from PySide6.QtWidgets import QDialog, QMessageBox
        from lunaNMR.core.ps2d_config import PS2DConfig, get_ps2d_config, restore_base_ps2d_config
        monkeypatch.setattr(QDialog, "exec", lambda self: 0)
        window.on_expert_mode()
        window._on_nucleus_type_change('13C')
        monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.Yes))

        window.on_reset_results()
        window.sync_parameters_to_integrator()
        restore_base_ps2d_config()

        config = get_ps2d_config()
        assert (config.nucleus_type, config.radF1) == ('15N', PS2DConfig('15N').radF1)
        assert window.expert_nucleus_combo.currentText() == '15N'

    def test_a_series_on_another_nucleus_does_not_change_what_the_window_fits(self, window, data):
        from lunaNMR.core.ps2d_config import PS2DConfig, get_ps2d_config, restore_base_ps2d_config
        window._on_nucleus_type_change('13C')
        _run_series(window, data)
        assert get_ps2d_config().nucleus_type == '15N'

        window.sync_parameters_to_integrator()
        restore_base_ps2d_config()

        config = get_ps2d_config()
        assert (config.nucleus_type, config.radF1) == ('13C', PS2DConfig('13C').radF1)


class TestTheExpertPanelIsWhatTheFitUses:

    @pytest.fixture
    def expert(self, window, monkeypatch):
        from PySide6.QtWidgets import QDialog
        monkeypatch.setattr(QDialog, "exec", lambda self: 0)
        window.on_expert_mode()
        return window

    def test_a_radius_typed_without_apply_survives_the_fresh_fit_restore(self, expert):
        from lunaNMR.core.ps2d_config import get_ps2d_config, restore_base_ps2d_config
        expert.expert_radF1_spin.setValue(0.25)

        expert.sync_parameters_to_integrator()
        restore_base_ps2d_config()

        assert get_ps2d_config().radF1 == 0.25

    def test_the_applied_13c_settings_come_back_after_a_series_on_15n(self, expert, data):
        from lunaNMR.core.ps2d_config import get_ps2d_config, restore_base_ps2d_config
        expert._on_nucleus_type_change('13C')
        expert.expert_radF1_spin.setValue(0.2)
        expert.expert_max_cluster_spin.setValue(5)
        expert._on_ps2d_params_apply()
        _run_series(expert, data)
        assert get_ps2d_config().nucleus_type == '15N'

        expert.sync_parameters_to_integrator()
        restore_base_ps2d_config()

        config = get_ps2d_config()
        assert (config.nucleus_type, config.radF1, config.max_cluster_size) == ('13C', 0.2, 5)


class TestCancelStopsTheRunFromCounting:
    """Cancel closes the progress dialog. A fit stops before its next cluster; a detection
    cannot be interrupted and runs on. Neither one's results may land, and no second run
    may start beside them."""

    @pytest.fixture
    def cancel_at_once(self, monkeypatch):
        from lunaNMR.gui.components import ProgressDialog

        def exec_then_cancel(self):
            self._on_cancel_clicked()
            return 0
        monkeypatch.setattr(ProgressDialog, "exec", exec_then_cancel)

    def _buttons(self, window):
        return (window.detect_peaks_button.isEnabled(), window.fit_all_peaks_button.isEnabled(),
                window.series_button.isEnabled())

    def test_a_cancelled_fit_leaves_the_results_it_found(self, app, window, data, monkeypatch,
                                                         cancel_at_once):
        _load(window, monkeypatch, data["spectrum"], data["ten"])
        _detect(window)
        detected = window.integrator.fitted_peaks

        window.on_fit_spectrum()
        _wait_for(app, window.fitting_worker)

        assert window.integrator.fitted_peaks is detected
        assert window.last_fitting_results is None
        assert not list(data["folder"].glob("*_fitting_results.csv"))

    @pytest.mark.parametrize("parallel", [False, True], ids=["sequential", "parallel"])
    def test_cancel_stops_the_fit_before_its_first_pass_ends(self, app, window, data, monkeypatch,
                                                             cancel_at_once, parallel):
        _load(window, monkeypatch, data["spectrum"], data["ten"])
        window.use_parallel_processing = parallel
        previous_clusters = window.integrator.enhanced_fitter.series_params

        window.on_fit_spectrum()
        _wait_for(app, window.fitting_worker)

        assert window.integrator.enhanced_fitter.parallel_processor.pass1_results is None
        assert window.integrator.enhanced_fitter.series_params is previous_clusters
        assert "cancelled" in window.statusBar().currentMessage().lower()
        assert self._buttons(window) == (True, True, True)

    def test_nothing_else_starts_while_a_cancelled_fit_runs(self, app, window, data, monkeypatch,
                                                            cancel_at_once):
        _load(window, monkeypatch, data["spectrum"], data["ten"])
        window.on_fit_spectrum()
        first = window.fitting_worker
        assert first.isRunning()
        assert self._buttons(window) == (False, False, False)

        window.on_fit_spectrum()
        window.on_detect_peaks()
        _load(window, monkeypatch, data["spectrum"], data["five"])

        assert window.fitting_worker is first
        assert getattr(window, 'detection_worker', None) is None
        assert len(window.integrator.peak_list) == 10
        _wait_for(app, first)
        assert self._buttons(window) == (True, True, True)

    def test_a_cancelled_detection_leaves_the_previous_one(self, app, window, data, monkeypatch,
                                                           cancel_at_once):
        _load(window, monkeypatch, data["spectrum"], data["ten"])
        _detect(window)
        previous = window.integrator.fitted_peaks

        window.on_detect_peaks()
        _wait_for(app, window.detection_worker)

        assert window.integrator.fitted_peaks is previous
        assert self._buttons(window) == (True, True, True)
