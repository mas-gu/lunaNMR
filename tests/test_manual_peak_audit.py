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


class TestRefitAndUndo:
    """Orchestration of the constrained re-fit (the fit itself is stubbed)."""

    @staticmethod
    def _stub_refit(dialog):
        # Replace the real fit with a deterministic write-back.
        def fake(spec_idx, assignment, x, y):
            peak = dialog._find_peak_obj(spec_idx, assignment)
            dialog._set_peak_position(peak, x, y)
            peak['r_squared'] = 0.95
            peak['quality'] = 'Excellent'
            return 0.95
        dialog._refit_peak_in_spectrum = fake

    def test_refit_marks_edited_and_clears_pending(self, dialog):
        dialog.overlay_selected_spectrum_idx = 0
        dialog._on_toggle_edit_mode(True)
        dialog._pending_corrections[(0, 'R1')] = (8.40, 121.5)
        self._stub_refit(dialog)
        dialog._on_refit_corrected()
        # R1 appears in both spectra, so both get re-fit (corrected anchor / current anchor).
        assert (0, 'R1') in dialog._edited
        assert (1, 'R1') in dialog._edited
        assert dialog._pending_corrections == {}
        assert dialog._find_peak_obj(0, 'R1')['r_squared'] == 0.95

    def test_undo_restores_pre_refit_state(self, dialog):
        dialog.overlay_selected_spectrum_idx = 0
        dialog._on_toggle_edit_mode(True)
        dialog._pending_corrections[(0, 'R1')] = (8.40, 121.5)
        self._stub_refit(dialog)
        dialog._on_refit_corrected()
        assert dialog._find_peak_obj(0, 'R1').get('r_squared') == 0.95
        dialog._on_undo_corrections()
        # snapshot was taken before re-fit, so r_squared is gone again
        assert 'r_squared' not in dialog._find_peak_obj(0, 'R1')
        assert dialog._edited == set()


class TestSaveCorrected:
    def test_rewrites_tracking_csv(self, dialog, tmp_path):
        import pandas as pd
        # Minimal tracking CSV with two points "0" and "1" matching the 2 spectra.
        cols = {'Peak_Number': [1], 'Assignment': ['R1'],
                'Reference_X': [8.0], 'Reference_Y': [120.0]}
        for label in ('0', '1'):
            cols[f'{label}_Position_X'] = [0.0]
            cols[f'{label}_Position_Y'] = [0.0]
            cols[f'{label}_R_Squared'] = [0.0]
        pd.DataFrame(cols).to_csv(tmp_path / "comprehensive_peak_tracking.csv", index=False)

        dialog.output_folder = str(tmp_path)
        dialog._edited = {(0, 'R1')}
        # corrected in-memory peak for spectrum 0
        peak = dialog._find_peak_obj(0, 'R1')
        dialog._set_peak_position(peak, 8.40, 121.5)
        peak['r_squared'] = 0.95

        dialog._on_save_corrected()

        out = pd.read_csv(tmp_path / "comprehensive_peak_tracking.csv")
        assert out.loc[0, '0_Position_X'] == 8.40
        assert out.loc[0, '0_Position_Y'] == 121.5
        assert out.loc[0, '0_R_Squared'] == 0.95
