# ABOUTME: The three Results Analysis exports write, per spectrum, the counts a series result carries.
# ABOUTME: Real LunaNMRMainWindow offscreen; a real two-point DNAJA1 series of a 10-peak list.

"""A series result records, per spectrum, how many peaks it fitted and how many of the
fits succeeded. It records no detection rate, no detected-peak count and no processing
time, so an export column of that name could only ever hold 0."""

import json
import os
import sys
from pathlib import Path

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
SPECTRA = ("2D_DNAJA1_0o0.ft", "2D_DNAJA1_HSPA8_0o5.ft")

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")

ABSENT_FROM_A_SERIES_RESULT = {'detection_rate', 'detected_peaks', 'processing_time'}


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


def _new_window():
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    set_ps2d_config('15N')
    return LunaNMRMainWindow()


def _close(window):
    from lunaNMR.core.ps2d_config import set_ps2d_config
    window.close()
    window.deleteLater()
    set_ps2d_config('15N')


@pytest.fixture(scope="module")
def series(app, tmp_path_factory):
    """The BatchResults of the series dialog's worker over two spectra and ten peaks."""
    from lunaNMR.gui.dialogs.series_integration_dialog import ProcessingWorker, SeriesIntegrationDialog
    folder = tmp_path_factory.mktemp("series_exports")
    for name in SPECTRA:
        (folder / name).symlink_to(DNAJA1 / name)
    lines = (DNAJA1 / "DNAJA1.txt").read_text().splitlines()
    (folder / "ten.txt").write_text("\n".join(lines[:11]) + "\n")
    window = _new_window()
    window.load_data(str(folder / SPECTRA[0]), str(folder / "ten.txt"))
    dialog = SeriesIntegrationDialog(main_window=window)
    finished = []
    worker = ProcessingWorker([str(folder / name) for name in SPECTRA],
                              window.integrator.peak_list.copy(), 'reference',
                              dialog._get_voigt_params(), series_name='exports')
    worker.finished.connect(finished.append)
    worker.run()
    dialog.deleteLater()
    _close(window)
    assert len(finished[0].results) == 2
    return finished[0]


@pytest.fixture
def window(app, series, monkeypatch):
    """A window holding the series, its message boxes answered at once."""
    from PySide6.QtWidgets import QMessageBox
    w = _new_window()
    w.batch_results = series
    for kind in ("critical", "information", "warning"):
        monkeypatch.setattr(QMessageBox, kind, staticmethod(lambda *a, **k: None))
    yield w
    _close(w)


def _trigger(window, text):
    from PySide6.QtWidgets import QMenu
    next(action for menu in window.menuBar().findChildren(QMenu)
         for action in menu.actions() if action.text() == text).trigger()


def _counts(result):
    return {'total_peaks': result['total_peaks'], 'successful_fits': result['successful_fits'],
            'success_rate': result['success_rate']}


def test_a_series_result_counts_its_fits_and_no_detection(series):
    for result in series.results.values():
        assert result['total_peaks'] == 10 and result['successful_fits'] > 0
        assert not ABSENT_FROM_A_SERIES_RESULT & set(result)


def test_the_data_matrix_carries_each_spectrum_s_fit_counts(window, series, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QFileDialog
    monkeypatch.setattr(QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(tmp_path / "matrix.csv"), "")))

    _trigger(window, "Export Data Matrix...")

    matrix = pd.read_csv(tmp_path / "matrix.csv")
    assert not ABSENT_FROM_A_SERIES_RESULT & set(matrix.columns)
    assert len(matrix) == 20
    for name, result in series.results.items():
        rows = matrix[matrix['spectrum'] == name]
        assert len(rows) == 10
        for column, value in _counts(result).items():
            assert (rows[column] == value).all(), column


def test_the_analysis_report_carries_each_spectrum_s_fit_counts(window, series, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QFileDialog
    monkeypatch.setattr(QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(tmp_path / "report.json"), "")))

    _trigger(window, "Export Analysis Report...")

    report = json.loads((tmp_path / "report.json").read_text())
    averages = []
    for name, result in series.results.items():
        r_squared = [row['r_squared'] for row in result['integration_results']]
        averages.append(sum(r_squared) / len(r_squared))
        assert averages[-1] > 0
        assert report['spectra_details'][name] == dict(
            _counts(result), status='success', avg_r_squared=pytest.approx(averages[-1]))
    assert report['quality_statistics']['r_squared']['min'] == pytest.approx(min(averages))
    assert report['quality_statistics']['r_squared']['max'] == pytest.approx(max(averages))
    rates = [result['success_rate'] for result in series.results.values()]
    assert report['quality_statistics']['success_rate']['min'] == pytest.approx(min(rates))
    assert set(report['quality_statistics']) == {'r_squared', 'success_rate'}


def test_the_batch_export_summary_carries_each_spectrum_s_fit_counts(window, series, tmp_path,
                                                                    monkeypatch):
    from PySide6.QtWidgets import QDialog, QFileDialog
    monkeypatch.setattr(QDialog, "exec", lambda self: QDialog.Accepted)
    monkeypatch.setattr(QFileDialog, "getExistingDirectory",
                        staticmethod(lambda *a, **k: str(tmp_path)))

    _trigger(window, "Batch Export Results...")

    summary = pd.read_csv(next(tmp_path.glob("summary_*.csv")))
    assert list(summary.columns) == ['spectrum', 'status', 'total_peaks', 'successful_fits',
                                     'success_rate', 'avg_r_squared', 'avg_snr']
    for name, result in series.results.items():
        row = summary[summary['spectrum'] == name].iloc[0]
        for column, value in _counts(result).items():
            assert row[column] == pytest.approx(value), column
    report = json.loads(next(tmp_path.glob("analysis_report_*.json")).read_text())
    assert [set(row) for row in report['summary']] == [set(summary.columns)] * 2
