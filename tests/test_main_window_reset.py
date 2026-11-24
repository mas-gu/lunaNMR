# ABOUTME: The main window shows the contour settings it draws with; Reset Results discards what Detect and Fit made and keeps the rest.
# ABOUTME: Real LunaNMRMainWindow offscreen, on the DNAJA1 spectrum and a 10-peak list cut from its peak list where a fit is needed.

"""The window restores no preferences: nothing in the application saves any but the recent
projects list, so a loader of window, display and contour keys only ever read its own defaults.

Reset Results returns the window to where it stood once the spectrum and the peak list were
loaded and the parameters set: what Detect and Fit Spectrum produced is discarded, what the
user loaded or set is kept. A project opened after another is what starts a clean session
(tests/test_project_session.py)."""

import os
import sys
from types import SimpleNamespace

import pandas as pd
import pytest

from test_main_window_fit_feedback import DNAJA1, SPECTRUM, _detect, _fit, _suptitle
from test_main_window_project_parameters import EXPERT_SETTINGS, MAIN_SETTINGS, _set, _showing

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

needs_data = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def make_window(app):
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    windows = []

    def make():
        set_ps2d_config('15N')
        windows.append(LunaNMRMainWindow())
        return windows[-1]

    yield make
    for w in windows:
        for worker in (getattr(w, 'fitting_worker', None), getattr(w, 'detection_worker', None)):
            if worker is not None:
                worker.wait()
        app.processEvents()
        w.close()
        w.deleteLater()
    app.processEvents()
    set_ps2d_config('15N')


@pytest.fixture
def confirmed(monkeypatch):
    """The confirmation is answered Yes; the texts it showed."""
    from PySide6.QtWidgets import QMessageBox
    shown = []

    def question(parent, title, text, *args, **kwargs):
        shown.append(text)
        return QMessageBox.Yes

    monkeypatch.setattr(QMessageBox, "question", staticmethod(question))
    return shown


def _snapshot(window):
    """Every spin box, check box, radio button and combo box the window holds, and its plain state."""
    from PySide6.QtWidgets import QAbstractSpinBox, QCheckBox, QComboBox, QRadioButton
    out = {}
    for name, obj in vars(window).items():
        if isinstance(obj, QAbstractSpinBox) and hasattr(obj, 'value'):
            out[name] = obj.value()
        elif isinstance(obj, (QCheckBox, QRadioButton)):
            out[name] = obj.isChecked()
        elif isinstance(obj, QComboBox):
            out[name] = obj.currentText()
        elif isinstance(obj, (int, float, str, bool)) and not name.startswith('_'):
            out['state.' + name] = obj
    return out


def _ps2d_config():
    from lunaNMR.core.ps2d_config import get_ps2d_config
    config = get_ps2d_config()
    return (config.nucleus_type, config.radF1, config.radF2, config.overlap_threshold_x,
            config.overlap_threshold_y, config.max_iterations, config.max_cluster_size)


def _navigator_rows(navigator):
    """The navigator's table as shown: each cell's text and background."""
    table = navigator.table
    return [[(table.item(row, column).text(), table.item(row, column).background().color().name())
             for column in range(table.columnCount())] for row in range(table.rowCount())]


def _result_displays(window):
    """Everything in the window that shows a detection or a fit."""
    navigator = window.peak_navigator
    return {
        'statistics': window.stats_label.text(),
        'navigator list': navigator.peak_type_combo.currentText(),
        'navigator rows': _navigator_rows(navigator),
        'navigator status': navigator.status_label.text(),
        'navigator selection': navigator.selected_peak_index,
        'peak parameters': window.peak_params_text.toPlainText(),
        'selected peak': window.selected_peak_label.text(),
        'markers': len(window.spectrum_plotter.peak_annotations),
        'voigt titles': (_suptitle(window.voigt_2d_plotter), _suptitle(window.voigt_3d_plotter)),
        'voigt plots': any(ax.has_data() for plotter in (window.voigt_2d_plotter, window.voigt_3d_plotter)
                           for ax in plotter.figure.axes),
        'buttons': [button.styleSheet() for button in (
            window.detect_peaks_button, window.fit_all_peaks_button, window.series_button)],
    }


def _view(window):
    axes = window.spectrum_plotter.ax
    return (*axes.get_xlim(), *axes.get_ylim())


def _peaks(window):
    """{assignment: every value of the peak} of the integrator's detected or fitted peaks, less
    the id the integrator numbers it with."""
    return {peak['assignment']: {key: value for key, value in peak.items()
                                 if isinstance(value, (int, float, str, bool)) and key != 'peak_id'}
            for peak in window.integrator.fitted_peaks}


@pytest.fixture
def analysed(app, make_window, tmp_path, monkeypatch):
    """A window in which the user loaded the spectrum and a 10-peak list, set parameters, moved
    a peak of the list and zoomed, then ran Detect and Fit Spectrum."""
    from lunaNMR.gui.components import ProgressDialog
    monkeypatch.setattr(ProgressDialog, "exec", lambda self: 0)
    window = make_window()
    (tmp_path / SPECTRUM).symlink_to(DNAJA1 / SPECTRUM)
    lines = (DNAJA1 / "DNAJA1.txt").read_text().splitlines()
    peaks = tmp_path / "ten.txt"
    peaks.write_text("\n".join(lines[:11]) + "\n")
    window.load_data(str(tmp_path / SPECTRUM), str(peaks))
    window.use_parallel_processing = False
    window.contour_levels_spin.setValue(12)
    moved = window.integrator.peak_list.iloc[0]
    window.undo_manager.move_peak(int(moved['peak_id']), 'reference',
                                  moved['Position_X'], moved['Position_Y'],
                                  moved['Position_X'] + 0.002, moved['Position_Y'] + 0.02)
    window.spectrum_plotter.set_zoom(8.2, 121.0, 1.0, 12.0)
    loaded = SimpleNamespace(displays=_result_displays(window), view=_view(window),
                             peak_list=window.integrator.peak_list.copy(),
                             spectrum=window.integrator.nmr_data, settings=_snapshot(window))
    _detect(app, window)
    detected = _peaks(window)
    _fit(app, window)
    return SimpleNamespace(window=window, loaded=loaded, detected=detected, fitted=_peaks(window))


def test_the_contour_spin_boxes_show_what_the_contours_are_drawn_with(make_window):
    window = make_window()
    assert (window.contour_levels_spin.value(), window.contour_min_spin.value(),
            window.contour_inc_spin.value()) == (window.contour_levels, window.contour_min,
                                                 window.contour_increment)


def test_the_window_restores_no_preferences_nothing_saves():
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    assert not hasattr(LunaNMRMainWindow, 'apply_user_preferences')


def test_reset_keeps_every_parameter_and_expert_setting_the_user_set(make_window, confirmed,
                                                                    monkeypatch):
    from PySide6.QtWidgets import QDialog, QPushButton
    monkeypatch.setattr(QDialog, "exec", lambda self: 0)
    window = make_window()
    for name, value in MAIN_SETTINGS.items():
        _set(getattr(window, name), value)
    window.on_expert_mode()
    for name, value in EXPERT_SETTINGS.items():
        _set(getattr(window, name), value)
    next(button for button in window._expert_dialog.findChildren(QPushButton)
         if button.text() == "Apply Changes").click()
    settings = {**MAIN_SETTINGS, **EXPERT_SETTINGS}
    assert _showing(window, settings) == settings
    before = _snapshot(window), _ps2d_config()

    window.on_reset_results()

    assert _showing(window, settings) == settings
    assert (_snapshot(window), _ps2d_config()) == before
    assert window.sn_params_widget.isVisibleTo(window)


@needs_data
def test_reset_keeps_the_spectrum_the_edited_peak_list_and_the_view(analysed, confirmed):
    window, loaded = analysed.window, analysed.loaded
    files = window.current_nmr_file, window.current_peak_file
    labels = window.nmr_file_label.text(), window.peak_file_label.text()

    window.on_reset_results()

    assert (window.current_nmr_file, window.current_peak_file) == files
    assert (window.nmr_file_label.text(), window.peak_file_label.text()) == labels
    assert window.integrator.nmr_data is loaded.spectrum
    pd.testing.assert_frame_equal(window.integrator.peak_list, loaded.peak_list)
    assert _view(window) == pytest.approx(loaded.view)
    assert window.spectrum_plotter.contour_plot is not None
    assert _snapshot(window) == loaded.settings


@needs_data
def test_reset_discards_what_detect_and_fit_produced(analysed, confirmed):
    window, loaded = analysed.window, analysed.loaded
    window.show_detected_peak_analysis_by_index(0, auto_switch_tab=False)
    window.peak_navigator.table.selectRow(1)
    first = window.integrator.peak_list.iloc[0]
    window._on_peak_select_for_move(first['Position_X'], first['Position_Y'])
    assert window.selected_peak_info is not None
    assert _result_displays(window) != loaded.displays

    window.on_reset_results()

    assert window.integrator.fitted_peaks == []
    assert window.last_fitting_results is None
    assert window.current_voigt_result is None
    assert window.integrator.enhanced_fitter.series_params is None
    assert window.voigt_3d_plotter.current_result is None
    assert window.peak_navigator.detected_peaks == []
    assert (window.selected_peak_info, window.selected_peaks_multi) == (None, [])
    assert window.spectrum_plotter._selection_highlight is None
    assert _result_displays(window) == loaded.displays
    assert [{text.get_text() for ax in plotter.figure.axes for text in ax.texts}
            for plotter in (window.voigt_2d_plotter, window.voigt_3d_plotter)] == [
        {"Select a peak to view Voigt analysis"}, {"Select a peak to view 3D Voigt analysis"}]


@needs_data
def test_detect_and_fit_after_reset_give_what_the_first_ones_gave(app, analysed, confirmed):
    window = analysed.window

    window.on_reset_results()
    assert window.integrator.nmr_data is not None and window.integrator.peak_list is not None
    _detect(app, window)
    detected = _peaks(window)
    _fit(app, window)

    assert len(analysed.fitted) == 10
    assert detected == analysed.detected
    assert _peaks(window) == analysed.fitted


@needs_data
def test_reset_discards_the_peak_list_a_detection_by_threshold_made(app, make_window, confirmed,
                                                                   tmp_path, monkeypatch):
    from lunaNMR.gui.components import ProgressDialog
    monkeypatch.setattr(ProgressDialog, "exec", lambda self: 0)
    window = make_window()
    window.workflow_sn_threshold_radio.setChecked(True)
    (tmp_path / SPECTRUM).symlink_to(DNAJA1 / SPECTRUM)
    window.load_data(str(tmp_path / SPECTRUM))
    assert window.integrator.peak_list is None
    _detect(app, window)
    made = window.integrator.peak_list.copy()
    assert len(made) > 0

    window.on_reset_results()

    assert window.integrator.peak_list is None
    assert window.integrator.fitted_peaks == []
    assert window.peak_navigator.reference_peaks == []
    assert len(window.spectrum_plotter.peak_annotations) == 0
    assert window.workflow_mode == "sn_threshold"
    assert window.integrator.nmr_data is not None
    _detect(app, window)
    pd.testing.assert_frame_equal(window.integrator.peak_list, made)


@needs_data
def test_a_list_loaded_after_a_detection_by_threshold_is_kept(app, make_window, confirmed,
                                                             tmp_path, monkeypatch):
    from lunaNMR.gui.components import ProgressDialog
    monkeypatch.setattr(ProgressDialog, "exec", lambda self: 0)
    window = make_window()
    window.workflow_sn_threshold_radio.setChecked(True)
    (tmp_path / SPECTRUM).symlink_to(DNAJA1 / SPECTRUM)
    window.load_data(str(tmp_path / SPECTRUM))
    _detect(app, window)
    lines = (DNAJA1 / "DNAJA1.txt").read_text().splitlines()
    peaks = tmp_path / "ten.txt"
    peaks.write_text("\n".join(lines[:11]) + "\n")
    window.load_peak_file(str(peaks))
    loaded = window.integrator.peak_list.copy()

    window.on_reset_results()

    pd.testing.assert_frame_equal(window.integrator.peak_list, loaded)


@needs_data
def test_the_list_an_opened_project_brought_is_kept(app, make_window, confirmed, tmp_path,
                                                    monkeypatch):
    from lunaNMR.gui.components import ProgressDialog
    monkeypatch.setattr(ProgressDialog, "exec", lambda self: 0)
    window = make_window()
    window.workflow_sn_threshold_radio.setChecked(True)
    (tmp_path / SPECTRUM).symlink_to(DNAJA1 / SPECTRUM)
    window.load_data(str(tmp_path / SPECTRUM))
    _detect(app, window)
    project = tmp_path / "threshold.lunaNMR"
    assert window._get_project_manager().save_project(project)[0]
    assert window._get_project_manager().load_project(project)[0]
    window._refresh_after_project_load()
    brought = window.integrator.peak_list.copy()
    assert len(brought) > 0

    window.on_reset_results()

    pd.testing.assert_frame_equal(window.integrator.peak_list, brought)


def test_reset_keeps_the_series_and_the_kd_analysis_of_the_session(make_window, confirmed):
    """They are deleted where they are listed: their own dialogs and the Project Browser."""
    from lunaNMR.processors.multi_spectrum_processor import BatchResults
    window = make_window()
    series = BatchResults()
    kd = {'titration': {'fit_data': {}, 'meta': {}}}
    window.saved_series = {'T1': series}
    window.batch_results = series
    window.kd_analyses = kd
    window.dynamixs_analyses = {'T1T2': {'state': {}, 'file_refs': {}}}
    window.current_project_path = "/somewhere/project.lunaNMR"

    window.on_reset_results()

    assert window.saved_series == {'T1': series} and window.batch_results is series
    assert window.kd_analyses is kd
    assert list(window.dynamixs_analyses) == ['T1T2']
    assert window.current_project_path == "/somewhere/project.lunaNMR"


def test_reset_leaves_the_kd_window_open(make_window, confirmed):
    window = make_window()
    window.launch_kd_titration()
    dialog = window.kd_titration_dialog
    assert dialog.isVisible()

    window.on_reset_results()

    assert dialog.isVisible()
    assert window.kd_titration_dialog is dialog


def _load_two_peaks(window, tmp_path):
    peaks = tmp_path / "peaks.csv"
    peaks.write_text("Assignment,Position_X,Position_Y\n1AlaH,8.10,120.5\n2GlyH,8.40,109.2\n")
    window.load_peak_file(str(peaks))


def test_reset_empties_the_undo_history_and_keeps_the_edit(make_window, confirmed, tmp_path):
    window = make_window()
    _load_two_peaks(window, tmp_path)
    window.add_new_peak(8.25, 115.0)
    assert window.undo_action.isEnabled()

    window.on_reset_results()

    assert not window.undo_action.isEnabled() and not window.redo_action.isEnabled()
    assert len(window.integrator.peak_list) == 3


def test_an_edit_made_after_reset_is_undone_and_redone(make_window, confirmed, tmp_path):
    window = make_window()
    _load_two_peaks(window, tmp_path)
    window.on_reset_results()

    window.add_new_peak(8.25, 115.0)
    window.undo_action.trigger()
    assert list(window.integrator.peak_list['Assignment']) == ['1AlaH', '2GlyH']
    window.redo_action.trigger()

    assert len(window.integrator.peak_list) == 3
    assert len(window.peak_navigator.reference_peaks) == 3


def test_the_confirmation_and_the_button_say_what_is_discarded_and_what_is_kept(make_window,
                                                                               confirmed):
    from PySide6.QtWidgets import QPushButton
    window = make_window()
    button = next(b for b in window.findChildren(QPushButton) if b.text() == "Reset Results")

    button.click()

    for text in (confirmed[0], button.toolTip()):
        discarded, kept = text.lower().split("kept")[0], text.lower()
        assert "detect" in discarded and "fit" in discarded, text
        assert "kept" in kept and "spectrum" in kept and "peak list" in kept, text
    assert "closing and reopening" not in confirmed[0]


def test_answering_no_discards_nothing(make_window, monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    window = make_window()
    detected = [{'assignment': '1AlaH', 'ppm_x': 8.10, 'ppm_y': 120.5}]
    window.integrator.fitted_peaks = detected
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.No))

    window.on_reset_results()

    assert window.integrator.fitted_peaks is detected
