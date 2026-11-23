# ABOUTME: A project whose gui_state.json names parameters the window has no control for opens, and takes none of them.
# ABOUTME: Real LunaNMRMainWindow offscreen, saved and opened through the real ProjectManager.

"""The window holds no custom initial linewidths, no detection rectangle in pixels and none
of the 1D peak-picking and fit-quality settings: no engine code reads any of them. A
bundle can still name them in its gui_state.json. Opening it restores the parameters the
window has, and neither the window nor a series started from it carries the others."""

import json
import os
import sys

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

WITHOUT_A_CONTROL = {
    'use_custom_linewidths': True, 'lw_lorentz_1h': 0.002, 'lw_gauss_1h': 0.04,
    'lw_lorentz_15n': 0.0002, 'lw_gauss_15n': 0.4,
    'detection_square_size': 5, 'detection_rectangle_y': 2,
    'use_voigt_fitting': False, 'use_ps2d_multi_peak': False,
    'use_centroid_refinement': False, 'centroid_noise_multiplier': 4.5,
    'height_threshold': 0.5, 'distance_factor': 50.0, 'prominence_threshold': 0.3,
    'smoothing_sigma': 3.0, 'max_peaks_fit': 2, 'max_optimization_iterations': 2,
}


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def make_window(app, monkeypatch):
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


@pytest.fixture
def project(make_window, tmp_path):
    """A saved project whose gui_state.json also names the parameters above."""
    window = make_window()
    window.expected_peaks_spin.setValue(120)
    path = tmp_path / "parameters.lunaNMR"
    assert window._get_project_manager().save_project(path)[0]
    state_file = path / "gui_state.json"
    state = json.loads(state_file.read_text())
    assert not set(WITHOUT_A_CONTROL) & set(state)
    state_file.write_text(json.dumps({**state, **WITHOUT_A_CONTROL}, indent=2))
    return path


@pytest.fixture
def opened(make_window, project):
    window = make_window()
    assert window._get_project_manager().load_project(project)[0]
    window._refresh_after_project_load()
    return window


def test_the_project_opens_with_the_parameters_the_window_has(opened):
    assert opened.expected_peaks_spin.value() == 120


def test_the_window_takes_no_parameter_it_has_no_control_for(opened):
    for name in WITHOUT_A_CONTROL:
        assert not hasattr(opened, name), name


def test_a_series_started_from_it_is_given_none_of_them(opened):
    from lunaNMR.gui.dialogs.series_integration_dialog import SeriesIntegrationDialog
    dialog = SeriesIntegrationDialog(main_window=opened)
    try:
        params = dialog._get_voigt_params()
        assert not set(WITHOUT_A_CONTROL) & set(params['gui_params'])
        assert 'fitting_params' not in params
    finally:
        dialog.deleteLater()


def test_a_fit_started_from_it_is_given_none_of_them(opened):
    opened.sync_parameters_to_integrator()
    assert not set(WITHOUT_A_CONTROL) & set(opened.integrator.gui_params)


def test_saving_it_again_writes_none_of_them(opened, tmp_path):
    path = tmp_path / "again.lunaNMR"
    assert opened._get_project_manager().save_project(path)[0]
    assert not set(WITHOUT_A_CONTROL) & set(json.loads((path / "gui_state.json").read_text()))
