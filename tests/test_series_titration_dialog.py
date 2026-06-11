# ABOUTME: Smoke tests for the series integration dialog x-axis mode dropdown.
# ABOUTME: Verifies off/time/titration selection maps to the right extractor mode.

import sys

import pytest

from PySide6.QtWidgets import QApplication

from lunaNMR.gui.dialogs.series_integration_dialog import (
    SeriesIntegrationDialog,
    ProcessingWorker,
)


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def dialog(app):
    dlg = SeriesIntegrationDialog(main_window=None)
    yield dlg
    dlg.deleteLater()


class TestXAxisModeDropdown:
    def test_three_modes(self, dialog):
        assert dialog.xaxis_mode_combo.count() == 3

    @pytest.mark.parametrize("index,enabled,mode", [
        (0, False, "time"),     # off: disabled, mode falls back to time
        (1, True, "time"),
        (2, True, "titration"),
    ])
    def test_mode_selection(self, dialog, index, enabled, mode):
        dialog.xaxis_mode_combo.setCurrentIndex(index)
        assert dialog._xaxis_enabled() is enabled
        assert dialog._series_mode() == mode


class TestWorkerSeriesMode:
    def test_worker_stores_series_mode(self, app):
        worker = ProcessingWorker(
            None, [], None, "reference", {},
            extract_delays=True, series_mode="titration",
        )
        assert worker.series_mode == "titration"

    def test_worker_defaults_to_time(self, app):
        worker = ProcessingWorker(None, [], None, "reference", {})
        assert worker.series_mode == "time"
