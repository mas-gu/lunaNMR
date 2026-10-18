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
