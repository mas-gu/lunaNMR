# ABOUTME: A series opened from a project bundle summarises itself as one run in the session does.
# ABOUTME: Real DNAJA1 series -> bundle -> the real main window: its summary and Export Analysis Report.

"""
A bundle stores a series' start and end times as text. Read back as text, the summary
of a loaded series could not subtract them and raised, and with it everything that asks
a series for its summary: the Export Analysis Report of the main window among them.
"""

import json
import os
import sys
from datetime import datetime

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

SERIES = "DNAJA1"


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture(scope="module")
def stored_times(dnaja1_series_bundle):
    """The series' start and end as the bundle stores them."""
    stored = json.loads((dnaja1_series_bundle / "series_results" / SERIES
                         / "batch_results.json").read_text())
    return stored['metadata']['start_time'], stored['metadata']['end_time']


@pytest.fixture
def window(app, dnaja1_series_bundle):
    """A main window that opened the project bundle."""
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    set_ps2d_config('15N')
    w = LunaNMRMainWindow()
    assert w._get_project_manager().load_project(dnaja1_series_bundle)[0]
    yield w
    w.close()
    w.deleteLater()
    set_ps2d_config('15N')


def _duration(stored_times):
    start, end = (datetime.fromisoformat(text) for text in stored_times)
    return str(end - start).split('.')[0]


def test_a_loaded_series_holds_its_times_as_a_series_run_in_the_session_does(window):
    metadata = window.batch_results.metadata
    assert isinstance(metadata['start_time'], datetime)
    assert isinstance(metadata['end_time'], datetime)


def test_a_loaded_series_summarises_itself(window, stored_times):
    summary = window.batch_results.get_summary()
    assert (summary['total_spectra'], summary['successful'], summary['failed']) == (3, 3, 0)
    assert summary['success_rate'] == 100
    assert summary['duration'] == _duration(stored_times)


def test_the_analysis_report_of_a_loaded_series_is_exported(window, stored_times, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QFileDialog, QMenu, QMessageBox
    said = []
    monkeypatch.setattr(QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(tmp_path / "report.json"), "")))
    monkeypatch.setattr(QMessageBox, "information",
                        staticmethod(lambda parent, title, text: said.append(title)))
    monkeypatch.setattr(QMessageBox, "critical",
                        staticmethod(lambda parent, title, text: said.append(f"{title}: {text}")))

    next(action for menu in window.menuBar().findChildren(QMenu)
         for action in menu.actions() if action.text() == "Export Analysis Report...").trigger()

    assert said == ["Export Complete"]
    report = json.loads((tmp_path / "report.json").read_text())
    assert report['summary']['total_spectra'] == 3
    assert report['summary']['successful'] == 3
    assert report['summary']['duration'] == _duration(stored_times)
    assert set(report['spectra_details']) == set(window.batch_results.results)
    for details in report['spectra_details'].values():
        assert details['status'] == 'success' and 0 < details['avg_r_squared'] <= 1


def test_saving_the_project_again_stores_the_same_times(window, stored_times, tmp_path):
    saved = tmp_path / "again.lunaNMR"
    assert window._get_project_manager().save_project(saved)[0]
    stored = json.loads((saved / "series_results" / SERIES / "batch_results.json").read_text())
    assert (stored['metadata']['start_time'], stored['metadata']['end_time']) == stored_times
