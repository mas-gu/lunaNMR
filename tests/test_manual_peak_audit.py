# ABOUTME: Tests the manual position-edit (two-click move) in the result viewer: picking a peak,
# ABOUTME: recording its corrected position, and the real constrained re-fit and its undo.

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
    def test_toggle_wires_left_click_callback(self, dialog):
        dialog._on_toggle_edit_mode(True)
        assert dialog.edit_mode is True
        assert dialog._nav_handler.on_left_click == dialog._on_edit_left_click
        dialog._on_toggle_edit_mode(False)
        assert dialog.edit_mode is False
        assert dialog._nav_handler.on_left_click is None
        assert dialog._nav_handler.on_escape is None

    def test_first_click_selects_nearest_peak(self, dialog):
        dialog.overlay_selected_spectrum_idx = 0
        dialog._on_toggle_edit_mode(True)
        dialog._on_edit_left_click(8.01, 120.1)    # click near R1
        assert dialog._edit_selected_peak['assignment'] == 'R1'

    def test_second_click_moves_and_auto_deselects(self, dialog):
        dialog.overlay_selected_spectrum_idx = 0
        dialog._on_toggle_edit_mode(True)
        dialog._on_edit_left_click(8.01, 120.1)    # select R1
        dialog._on_edit_left_click(8.40, 121.5)    # move it here
        assert dialog._pending_corrections[(0, 'R1')] == (8.40, 121.5)
        assert dialog._edit_selected_peak is None  # auto-deselected
        cx, cy = dialog._peak_center(dialog._get_fitted_peaks(0)[0])
        assert (cx, cy) == (8.40, 121.5)           # marker followed

    def test_first_click_far_from_any_peak_selects_nothing(self, dialog):
        dialog.overlay_selected_spectrum_idx = 0
        dialog._on_toggle_edit_mode(True)
        dialog._on_edit_left_click(5.0, 100.0)     # nowhere near a peak
        assert dialog._edit_selected_peak is None

    def test_click_selects_peak_from_a_non_active_spectrum(self, dialog):
        # Active spectrum is 0, but the click is closest to spectrum 1's R1
        # (at 8.05, 120.3) — selection must follow the click across spectra.
        for s in dialog.spectra:                  # overlay both spectra
            s['visible'] = True
            s['show_peaks'] = True
        dialog._set_overlay_active_spectrum(0)
        dialog._on_toggle_edit_mode(True)
        dialog._on_edit_left_click(8.05, 120.3)
        assert dialog._edit_selected_peak['spec_idx'] == 1
        assert dialog.overlay_selected_spectrum_idx == 1   # switched to clicked spectrum
        # and the move applies to that spectrum
        dialog._on_edit_left_click(8.40, 121.5)
        assert dialog._pending_corrections[(1, 'R1')] == (8.40, 121.5)
        assert (0, 'R1') not in dialog._pending_corrections

    def test_middle_click_moves_without_edit_mode(self, dialog):
        # Middle-click is a quick select->move that needs no edit mode.
        for s in dialog.spectra:
            s['visible'] = True
            s['show_peaks'] = True
        dialog._set_overlay_active_spectrum(0)
        assert dialog.edit_mode is False
        dialog._on_overlay_middle_click(8.00, 120.0)   # select p0's R1
        assert dialog._edit_selected_peak['assignment'] == 'R1'
        dialog._on_overlay_middle_click(8.40, 121.5)   # move it
        assert dialog._pending_corrections[(0, 'R1')] == (8.40, 121.5)
        assert dialog._edit_selected_peak is None      # auto-deselected

    def test_deselect_clears_selection_without_moving(self, dialog):
        dialog.overlay_selected_spectrum_idx = 0
        dialog._on_toggle_edit_mode(True)
        assert dialog._nav_handler.on_escape == dialog._edit_cancel
        orig = dialog._peak_center(dialog._get_fitted_peaks(0)[0])
        dialog._on_edit_left_click(8.01, 120.1)    # select R1
        assert dialog._edit_selected_peak is not None
        dialog._edit_cancel()                      # ✖ Deselect / Esc
        assert dialog._edit_selected_peak is None
        assert (0, 'R1') not in dialog._pending_corrections
        assert dialog._peak_center(dialog._get_fitted_peaks(0)[0]) == orig


class TestRefitAndUndo:
    """The constrained re-fit, on a real 3-point DNAJA1 titration written by `lunaNMR series`."""

    ASSIGNMENT = '42GlyH'

    @pytest.fixture
    def series_viewer(self, app, dnaja1_titration_series, open_overlay_viewer):
        viewer = open_overlay_viewer(dnaja1_titration_series)
        yield viewer
        viewer.deleteLater()

    def test_refit_marks_edited_and_clears_pending(self, series_viewer):
        viewer = series_viewer
        viewer.overlay_selected_spectrum_idx = 0
        viewer._on_toggle_edit_mode(True)
        before = dict(viewer._find_peak_obj(0, self.ASSIGNMENT))
        cx, cy = viewer._peak_center(before)
        viewer._pending_corrections[(0, self.ASSIGNMENT)] = (cx + 0.002, cy)
        viewer._on_refit_corrected()
        # The peak is in every spectrum, so each is re-fit (corrected anchor / current anchor).
        n = len(viewer.spectra)
        assert viewer._edited == {(i, self.ASSIGNMENT) for i in range(n)}
        assert viewer._pending_corrections == {}
        assert viewer.edit_status_label.text().startswith(f"Re-fit {n} peak-fits (0 failed)")
        after = viewer._find_peak_obj(0, self.ASSIGNMENT)
        assert after['height'] != before['height']   # written back from the refit

    def test_undo_restores_original_position_and_state(self, series_viewer):
        import copy
        viewer = series_viewer
        viewer.overlay_selected_spectrum_idx = 0
        viewer._on_toggle_edit_mode(True)
        originals = [copy.deepcopy(viewer._find_peak_obj(i, self.ASSIGNMENT))
                     for i in range(len(viewer.spectra))]
        cx, cy = viewer._peak_center(originals[0])
        # select, move, then re-fit
        viewer._on_edit_left_click(cx, cy)
        viewer._on_edit_left_click(cx + 0.002, cy)
        viewer._on_refit_corrected()
        assert viewer._find_peak_obj(0, self.ASSIGNMENT) != originals[0]
        viewer._on_undo_corrections()
        # undo reverts the drag AND the re-fit, back to the original fit, in every spectrum
        assert [viewer._find_peak_obj(i, self.ASSIGNMENT)
                for i in range(len(viewer.spectra))] == originals
        assert viewer._edited == set()
        assert viewer._pending_corrections == {}


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
        import json
        (tmp_path / "series_metadata.json").write_text(json.dumps({'columns': [
            {'column': '0', 'spectrum': 'p0', 'value': 0.0},
            {'column': '1', 'spectrum': 'p1', 'value': 1.0}]}))

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
