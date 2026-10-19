# ABOUTME: `lunaNMR batch` fits the peaks it detects, as the GUI's Fit Spectrum does, and writes
# ABOUTME: the GUI's fitting_results CSV per spectrum into --output. Real DNAJA1 spectrum, real main window.

"""
batch fitted each detected peak alone through `enhanced_peak_fitting`, which has nothing
to fit without the rest of the peak list: every fit came back empty, the spectrum still
counted as processed, the run exited 0, and nothing but a log was written.
"""

import contextlib
import io
import json
import os
import sys
import time
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
SPECTRUM = "2D_DNAJA1_0o0.ft"
CSV = "2D_DNAJA1_0o0_fitting_results.csv"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


def _spectrum_folder(root):
    spectra = root / "spectra"
    spectra.mkdir()
    (spectra / SPECTRUM).symlink_to(DNAJA1 / SPECTRUM)
    return spectra


@pytest.fixture(scope="module")
def batch_run(tmp_path_factory):
    """One `lunaNMR batch --preset 15N1H` over the spectrum: (exit code, envelope, --output)."""
    from lunaNMR.cli import main
    root = tmp_path_factory.mktemp("batch")
    out = root / "out"
    stdout = io.StringIO()
    with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(io.StringIO()):
        code = main(["batch", "--preset", "15N1H", str(_spectrum_folder(root)),
                     "--output", str(out), "--format", "json"])
    return code, json.loads(stdout.getvalue()), out


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


def _wait_for(app, worker, timeout=600):
    start = time.time()
    while worker.isRunning() and time.time() - start < timeout:
        app.processEvents()
        time.sleep(0.05)
    assert not worker.isRunning()
    for _ in range(50):
        app.processEvents()


@pytest.fixture(scope="module")
def gui_fit(app, batch_run, tmp_path_factory):
    """The main window on the same spectrum: Detect in S/N mode with the threshold and
    peak count the batch detected with, then Fit Spectrum. Returns (window, its CSV)."""
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.components import ProgressDialog
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    _, envelope, _ = batch_run
    spectrum_run = envelope["files"][0]
    spectra = _spectrum_folder(tmp_path_factory.mktemp("gui"))

    set_ps2d_config('15N')
    window = LunaNMRMainWindow()
    exec_ = ProgressDialog.exec
    ProgressDialog.exec = lambda self: 0
    try:
        window.load_data(str(spectra / SPECTRUM))
        window.workflow_sn_threshold_radio.setChecked(True)
        window.sn_threshold_spin.setValue(spectrum_run["sn_threshold"])
        window.expected_peaks_spin.setValue(spectrum_run["expected_peak_count"])
        window.on_detect_peaks()
        _wait_for(app, window.detection_worker)
        window.on_fit_spectrum()
        _wait_for(app, window.fitting_worker)
    finally:
        ProgressDialog.exec = exec_
    yield window, spectra / CSV
    window.close()
    window.deleteLater()
    set_ps2d_config('15N')


def test_the_spectrum_is_fitted_and_processed(batch_run):
    code, envelope, _ = batch_run
    assert code == 0
    assert envelope["ok"] is True
    assert envelope["processed_files"] == 1
    assert envelope["failed_files"] == 0
    spectrum_run = envelope["files"][0]
    assert spectrum_run["nucleus_type"] == "15N1H"
    assert spectrum_run["peaks_detected"] > 0
    assert spectrum_run["peaks_fitted"] > 0.8 * spectrum_run["peaks_detected"]


def test_the_fitting_results_csv_is_written_into_output(batch_run):
    import pandas as pd
    _, envelope, out = batch_run
    assert envelope["files"][0]["results_file"] == str(out / CSV)
    table = pd.read_csv(out / CSV)
    assert list(table.columns) == ['Assignment', 'Position_X', 'Position_Y', 'Height', 'Volume',
                                   'LW_X', 'LW_Y', 'R_Squared', 'Quality']
    assert len(table) == envelope["files"][0]["peaks_detected"]


def test_batch_writes_the_csv_the_gui_writes(batch_run, gui_fit):
    _, _, out = batch_run
    _, gui_csv = gui_fit
    assert (out / CSV).read_text() == gui_csv.read_text()


def test_batch_fits_each_peak_as_the_gui_does(batch_run, gui_fit, tmp_path):
    """At full precision: the batch processor's own results against the window's."""
    from lunaNMR.batch_processing.batch_processor import BatchProcessor
    _, envelope, _ = batch_run
    spectrum_run = envelope["files"][0]
    window, _ = gui_fit
    processor = BatchProcessor(output_dir=tmp_path)
    result = processor.process_single_spectrum(
        _spectrum_folder(tmp_path) / SPECTRUM,
        processor.get_optimal_parameters(spectrum_run["nucleus_type"], Path(SPECTRUM))
        | {'sn_threshold': spectrum_run["sn_threshold"],
           'expected_peak_count': spectrum_run["expected_peak_count"]})

    def per_peak(results):
        return {r['assignment']: (r['peak_x'], r['peak_y'], r['height'], r['volume'],
                                  r['fitting_quality'], r['success'])
                for r in results}

    gui = per_peak(window.last_fitting_results)
    assert len(gui) == spectrum_run["peaks_detected"]
    assert per_peak(result['fitted_results']) == gui


def test_a_spectrum_with_nothing_fitted_fails_the_run(tmp_path):
    """No peak above the threshold: nothing to fit, so the spectrum failed and so did the
    run - exit 1, BatchFailed naming the spectrum and why."""
    from lunaNMR.cli import main
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"processing": {"sn_thresholds": {"15N1H": 1000.0}},
                                  "optimization": {"enable_auto_optimization": False}}))
    stdout = io.StringIO()
    with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(io.StringIO()):
        code = main(["batch", "--config", str(config), str(_spectrum_folder(tmp_path)),
                     "--output", str(tmp_path / "out"), "--format", "json"])
    envelope = json.loads(stdout.getvalue())
    assert code == 1
    assert envelope["ok"] is False
    assert (envelope["processed_files"], envelope["failed_files"]) == (0, 1)
    assert envelope["error"]["type"] == "BatchFailed"
    assert SPECTRUM in envelope["error"]["message"]
    assert envelope["files"][0]["results_file"] is None
    assert not list((tmp_path / "out").glob("*.csv"))
