# ABOUTME: Tests the manual position-edit (two-click move) in the result viewer.
# ABOUTME: Phase 1 — picking a peak and recording a corrected position; no refit yet.

import sys

import pytest

from PySide6.QtWidgets import QApplication

from lunaNMR.gui.dialogs.multi_spectrum_viewer_dialog import MultiSpectrumViewerDialog


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication(sys.argv)


def _results():
    # Two spectra, each with the same two assignments at slightly different spots.
    def spec(name, r1x, r1y):
        return {
            'spectrum_name': name,
            'spectrum_file': '',
            'integration_results': [
                {'assignment': 'R1', 'ppm_x': r1x, 'ppm_y': r1y, 'quality': 'Good'},
                {'assignment': 'R2', 'ppm_x': 7.0, 'ppm_y': 118.0, 'quality': 'Good'},
            ],
        }
    return [spec('p0', 8.00, 120.0), spec('p1', 8.05, 120.3)]


@pytest.fixture
def dialog(app):
    dlg = MultiSpectrumViewerDialog(all_results=_results())
    yield dlg
    dlg.deleteLater()


class TestTwoClickMove:
    def test_toggle_wires_callback(self, dialog):
        dialog._on_toggle_edit_mode(True)
        assert dialog.edit_mode is True
        assert dialog._nav_handler.on_peak_select == dialog._on_edit_peak_click
        dialog._on_toggle_edit_mode(False)
        assert dialog.edit_mode is False
        assert dialog._nav_handler.on_peak_select is None

    def test_first_click_picks_nearest_peak(self, dialog):
        dialog.overlay_selected_spectrum_idx = 0
        dialog._on_toggle_edit_mode(True)
        dialog._on_edit_peak_click(8.01, 120.1)   # near R1
        assert dialog._edit_selected_peak is not None
        assert dialog._edit_selected_peak['assignment'] == 'R1'

    def test_second_click_records_correction_and_moves_marker(self, dialog):
        dialog.overlay_selected_spectrum_idx = 0
        dialog._on_toggle_edit_mode(True)
        dialog._on_edit_peak_click(8.01, 120.1)   # pick R1
        dialog._on_edit_peak_click(8.40, 121.5)   # place
        assert dialog._pending_corrections[(0, 'R1')] == (8.40, 121.5)
        # selection cleared, ready for the next correction
        assert dialog._edit_selected_peak is None
        # in-memory peak moved so the marker follows
        cx, cy = dialog._peak_center(dialog._get_fitted_peaks(0)[0])
        assert (cx, cy) == (8.40, 121.5)

    def test_click_far_from_any_peak_selects_nothing(self, dialog):
        dialog.overlay_selected_spectrum_idx = 0
        dialog._on_toggle_edit_mode(True)
        dialog._on_edit_peak_click(5.0, 100.0)    # nowhere near a peak
        assert dialog._edit_selected_peak is None
