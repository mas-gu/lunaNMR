# ABOUTME: A project carries the parameters its window showed, and a window that opens it shows and uses them.
# ABOUTME: Real LunaNMRMainWindow offscreen, saved and opened through the real ProjectManager.

"""The window keeps each parameter twice: as an attribute, which a project saves and
restores, and in a widget, which the user reads and which detection and fitting are
synchronised from. A load that restores only the attribute leaves the window showing, and
running with, the values of the session before."""

import json
import os
import sys
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
SPECTRUM = "2D_DNAJA1_0o0.ft"

# What the user sets in the main window: widget -> value
MAIN_SETTINGS = {
    'workflow_sn_threshold_radio': True,
    'sn_threshold_spin': 6.5,
    'expected_peaks_spin': 120,
    'contour_levels_spin': 12,
    'contour_min_spin': 0.1,
    'contour_inc_spin': 1.4,
    'shift_1h_spin': 0.012,
    'shift_15n_spin': 0.25,
    'show_exp_3d_checkbox': False,
    'show_fit_3d_checkbox': False,
    'show_individual_3d_checkbox': True,
    'show_resid_3d_checkbox': True,
    'limit_peak_display_3d_checkbox': False,
    'color_scheme_3d_dropdown': 'Dark',
}
# What the user sets in Expert Mode; the nucleus first, since choosing one resets the rest
EXPERT_SETTINGS = {
    'expert_nucleus_combo': '13C',
    'expert_search_x_spin': 0.05,
    'expert_search_y_spin': 0.2,
    'expert_centroid_x_spin': 0.03,
    'expert_centroid_y_spin': 0.2,
    'expert_auto_add_dummy_cb': False,
    'expert_show_ellipses_check': True,
    'expert_fix_lw_check': True,
    'expert_fix_pos_check': True,
    'expert_parallel_check': False,
    'expert_adaptive_check': True,
    'expert_lw_reuse_check': True,
    'expert_cascade_drift_check': False,
    'expert_radF1_spin': 0.5,
    'expert_radF2_spin': 0.05,
    'expert_overlap_x_spin': 0.05,
    'expert_overlap_y_spin': 0.9,
    'expert_max_iter_spin': 300,
    'expert_max_cluster_spin': 6,
}


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def make_window(app, monkeypatch):
    """Windows of one session each: the PS2D config is back at its startup state for each."""
    from PySide6.QtWidgets import QDialog
    from lunaNMR.core.ps2d_config import reset_ps2d_config
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    monkeypatch.setattr(QDialog, "exec", lambda self: 0)
    windows = []

    def make():
        reset_ps2d_config()
        windows.append(LunaNMRMainWindow())
        return windows[-1]

    yield make
    for w in windows:
        w.close()
        w.deleteLater()
    app.processEvents()
    reset_ps2d_config()


def _set(widget, value):
    from PySide6.QtWidgets import QAbstractButton, QComboBox
    if isinstance(widget, QAbstractButton):
        widget.setChecked(value)
    elif isinstance(widget, QComboBox):
        widget.setCurrentText(value)
    else:
        widget.setValue(value)


def _shown(widget):
    from PySide6.QtWidgets import QAbstractButton, QComboBox
    if isinstance(widget, QAbstractButton):
        return widget.isChecked()
    if isinstance(widget, QComboBox):
        return widget.currentText()
    return widget.value()


def _showing(window, settings):
    return {name: pytest.approx(_shown(getattr(window, name))) if not isinstance(value, (bool, str))
            else _shown(getattr(window, name)) for name, value in settings.items()}


@pytest.fixture
def project(make_window, tmp_path):
    """A project saved from a window in which the user set every parameter above."""
    from PySide6.QtWidgets import QPushButton
    window = make_window()
    for name, value in MAIN_SETTINGS.items():
        _set(getattr(window, name), value)
    window.on_expert_mode()
    for name, value in EXPERT_SETTINGS.items():
        _set(getattr(window, name), value)
    next(button for button in window._expert_dialog.findChildren(QPushButton)
         if button.text() == "Apply Changes").click()
    assert _showing(window, {**MAIN_SETTINGS, **EXPERT_SETTINGS}) == {**MAIN_SETTINGS, **EXPERT_SETTINGS}
    path = tmp_path / "parameters.lunaNMR"
    assert window._get_project_manager().save_project(path)[0]
    return path


def _open(window, project):
    assert window._get_project_manager().load_project(project)[0]
    window._refresh_after_project_load()


def test_the_project_stores_the_threshold_and_peak_count_the_window_showed(project):
    state = json.loads((project / "gui_state.json").read_text())
    assert (state['sn_threshold'], state['expected_peak_count']) == (6.5, 120)


def test_a_window_that_opens_the_project_shows_its_parameters(make_window, project):
    window = make_window()

    _open(window, project)

    assert _showing(window, MAIN_SETTINGS) == MAIN_SETTINGS
    assert window.sn_params_widget.isVisibleTo(window)


def test_expert_mode_opened_after_the_project_shows_its_parameters(make_window, project):
    window = make_window()

    _open(window, project)
    window.on_expert_mode()

    assert _showing(window, EXPERT_SETTINGS) == EXPERT_SETTINGS


def test_expert_mode_opened_before_the_project_shows_its_parameters(make_window, project):
    window = make_window()
    window.on_expert_mode()

    _open(window, project)
    window.on_expert_mode()

    assert _showing(window, EXPERT_SETTINGS) == EXPERT_SETTINGS


@pytest.mark.parametrize("expert_mode_was_open", [False, True])
def test_detection_and_fitting_run_with_the_project_s_parameters(make_window, project,
                                                                expert_mode_was_open):
    from lunaNMR.core.ps2d_config import get_ps2d_config
    window = make_window()
    if expert_mode_was_open:
        window.on_expert_mode()

    _open(window, project)
    window.sync_parameters_to_integrator()

    integrator = window.integrator
    assert (integrator.sn_threshold, integrator.expected_peak_count) == (6.5, 120)
    assert integrator.processing_mode == 'sn_native'
    assert (integrator.search_window_x, integrator.search_window_y) == pytest.approx((0.05, 0.2))
    assert {key: integrator.gui_params[key] for key in (
        'centroid_window_x_ppm', 'centroid_window_y_ppm', 'fix_linewidths', 'fix_positions',
        'use_parallel_processing', 'use_adaptive_optimization')} == {
        'centroid_window_x_ppm': pytest.approx(0.03), 'centroid_window_y_ppm': pytest.approx(0.2),
        'fix_linewidths': True, 'fix_positions': True, 'use_parallel_processing': False,
        'use_adaptive_optimization': True}
    config = get_ps2d_config()
    assert (config.nucleus_type, config.radF1, config.radF2, config.overlap_threshold_x,
            config.overlap_threshold_y, config.max_iterations, config.max_cluster_size) == (
        '13C', pytest.approx(0.5), pytest.approx(0.05), pytest.approx(0.05), pytest.approx(0.9), 300, 6)


@pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")
def test_the_voigt_tabs_show_the_first_successful_fit_of_an_opened_project(make_window, tmp_path):
    from lunaNMR.gui.main_window import VoigtFittingWorker
    (tmp_path / SPECTRUM).symlink_to(DNAJA1 / SPECTRUM)
    lines = (DNAJA1 / "DNAJA1.txt").read_text().splitlines()
    (tmp_path / "ten.txt").write_text("\n".join(lines[:11]) + "\n")
    saving = make_window()
    saving.load_data(str(tmp_path / SPECTRUM), str(tmp_path / "ten.txt"))
    saving.sync_parameters_to_integrator()
    worker = VoigtFittingWorker(saving.integrator, saving.param_manager,
                                {'use_parallel': True, 'use_voigt_fitting': True,
                                 'skip_series_params': True})
    worker.fitting_complete.connect(lambda summary: saving.on_fitting_complete(summary, None))
    worker.run()
    first = next(result for result in saving.last_fitting_results if result.get('success'))
    assert saving._get_project_manager().save_project(tmp_path / "fit.lunaNMR")[0]
    window = make_window()

    _open(window, tmp_path / "fit.lunaNMR")

    assert window.voigt_3d_plotter.current_result['assignment'] == first['assignment']
    assert str(first['assignment']) in window.voigt_2d_plotter.figure._suptitle.get_text()
