# ABOUTME: File > Export Results writes the Fit Spectrum results to the file the user names.
# ABOUTME: Real LunaNMRMainWindow offscreen on the DNAJA1 spectrum, a 10-peak list cut from its peak list.

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
def window(app, tmp_path):
    """A window showing the spectrum, linked into a folder of its own so the fit's CSV
    lands there, and a 10-peak list."""
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    set_ps2d_config('15N')
    w = LunaNMRMainWindow()
    spectra = tmp_path / "spectra"
    spectra.mkdir()
    (spectra / SPECTRUM).symlink_to(DNAJA1 / SPECTRUM)
    lines = (DNAJA1 / "DNAJA1.txt").read_text().splitlines()
    peaks = tmp_path / "ten.txt"
    peaks.write_text("\n".join(lines[:11]) + "\n")
    w.load_data(str(spectra / SPECTRUM), str(peaks))
    yield w
    w.close()
    w.deleteLater()
    set_ps2d_config('15N')


@pytest.fixture
def boxes(monkeypatch):
    """The message boxes the window opens, as (kind, title, text), each answered at once."""
    from PySide6.QtWidgets import QMessageBox
    shown = []
    for kind in ("critical", "information", "warning"):
        monkeypatch.setattr(QMessageBox, kind, staticmethod(
            lambda parent, title, text, *a, kind=kind, **k: shown.append((kind, title, text))))
    return shown


def _fit(window):
    """Run the Fit Spectrum worker to the end, in this thread, and hand the window its summary."""
    from lunaNMR.gui.main_window import VoigtFittingWorker
    window.sync_parameters_to_integrator()
    worker = VoigtFittingWorker(window.integrator, window.param_manager,
                                {'use_parallel': True, 'use_voigt_fitting': True,
                                 'skip_series_params': True})
    worker.fitting_complete.connect(lambda summary: window.on_fitting_complete(summary, None))
    worker.run()


def _export_results(window, monkeypatch, destination):
    """The user picks File > Export Results and names the file."""
    from PySide6.QtWidgets import QFileDialog, QMenu
    monkeypatch.setattr(QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(destination), "")))
    action = next(action for menu in window.menuBar().findChildren(QMenu)
                  for action in menu.actions() if action.text() == "Export Results...")
    action.trigger()


def test_export_results_writes_the_fit_results_to_the_named_file(window, tmp_path, monkeypatch, boxes):
    _fit(window)
    saved_by_the_fit = tmp_path / "spectra" / "2D_DNAJA1_0o0_fitting_results.csv"
    assert len(saved_by_the_fit.read_text().splitlines()) == 11

    _export_results(window, monkeypatch, tmp_path / "exported.csv")

    assert (tmp_path / "exported.csv").read_text() == saved_by_the_fit.read_text()
    assert [kind for kind, _, _ in boxes] == ["information"]
    assert "Results exported: 10 peaks" in window.statusBar().currentMessage()


def test_export_results_without_a_fit_says_there_is_none(window, tmp_path, monkeypatch, boxes):
    _export_results(window, monkeypatch, tmp_path / "exported.csv")

    assert boxes == [("critical", "Error", "No fit results to export. Run Fit Spectrum first.")]
    assert not (tmp_path / "exported.csv").exists()
