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

    def test_undo_restores_original_position_and_state(self, dialog):
        dialog.overlay_selected_spectrum_idx = 0
        dialog._on_toggle_edit_mode(True)
        orig = dialog._peak_center(dialog._find_peak_obj(0, 'R1'))
        # drag then re-fit
        dialog._on_edit_peak_click(8.01, 120.1)
        dialog._on_edit_peak_click(8.40, 121.5)
        self._stub_refit(dialog)
        dialog._on_refit_corrected()
        assert dialog._find_peak_obj(0, 'R1').get('r_squared') == 0.95
        dialog._on_undo_corrections()
        # undo reverts the drag AND the re-fit, back to the original fit
        assert dialog._peak_center(dialog._find_peak_obj(0, 'R1')) == orig
        assert 'r_squared' not in dialog._find_peak_obj(0, 'R1')
        assert dialog._edited == set()
        assert dialog._pending_corrections == {}


class TestShowSyncsPeakNavigator:
    """Showing a spectrum focuses the Peak Navigator on it (user can still change)."""

    def test_visibility_show_sets_active_spectrum(self, dialog):
        name = dialog.spectra[1]['name']
        dialog._on_visibility_changed(name, True)
        assert dialog.overlay_selected_spectrum_idx == 1
        assert dialog.overlay_spectrum_list.currentRow() == 1

    def test_hide_does_not_change_active(self, dialog):
        dialog._set_overlay_active_spectrum(1)
        dialog._on_visibility_changed(dialog.spectra[0]['name'], False)  # hide a spectrum
        assert dialog.overlay_selected_spectrum_idx == 1  # unchanged


class TestRefitHelpers:
    def test_learned_linewidths_medians(self, dialog):
        for p, lwx, lwy in zip(dialog._get_fitted_peaks(0), (0.02, 0.04), (0.1, 0.3)):
            p['r_squared'], p['lw_x'], p['lw_y'] = 0.95, lwx, lwy
        lw = dialog._learned_linewidths(0)
        assert lw == {'lw_f2_median': 0.03, 'lw_f1_median': 0.2}

    def test_learned_linewidths_none_when_no_data(self, dialog):
        assert dialog._learned_linewidths(0) is None  # peaks have no lw keys

    def test_edit_window_scales_with_reference_lw(self, dialog):
        # Give the reference spectrum (idx 0) good fits with known linewidths.
        for p in dialog._get_fitted_peaks(0):
            p['r_squared'], p['lw_x'], p['lw_y'] = 0.99, 0.0227, 0.16
        dialog.__dict__.pop('_ref_lw_cache', None)  # recompute
        wf1, wf2 = dialog._edit_window()
        assert round(wf1, 4) == 0.48    # 3 x 0.16 (15N), under cap
        assert round(wf2, 4) == 0.0681  # 3 x 0.0227 (1H)

    def test_edit_window_caps_huge_lw(self, dialog):
        for p in dialog._get_fitted_peaks(0):
            p['r_squared'], p['lw_x'], p['lw_y'] = 0.99, 1.0, 1.0  # absurd LW
        dialog.__dict__.pop('_ref_lw_cache', None)
        wf1, wf2 = dialog._edit_window()
        assert wf1 == dialog._EDIT_MARGIN_CAP_F1
        assert wf2 == dialog._EDIT_MARGIN_CAP_F2

    def test_edit_window_falls_back_without_lw(self, dialog):
        dialog.__dict__.pop('_ref_lw_cache', None)  # peaks have no lw keys
        assert dialog._edit_window() == (dialog._EDIT_POS_MARGIN_F1, dialog._EDIT_POS_MARGIN_F2)

    def test_write_back_ignores_none_position(self, dialog):
        peak = dialog._find_peak_obj(0, 'R1')
        before = dialog._peak_center(peak)
        ok = dialog._write_back_fit(0, 'R1', {'pos_f2': None, 'pos_f1': None}, 0.9)
        assert ok is False
        assert dialog._peak_center(peak) == before  # not stamped with None


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
