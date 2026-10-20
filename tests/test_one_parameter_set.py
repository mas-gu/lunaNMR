# ABOUTME: The parameter manager has one parameter set, which the main window hands the integrator as it is.
# ABOUTME: Pins that set and the main window's sync, and that no alternative parameter mode or its helpers remain.

import importlib
import os
import sys

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

GUI_PARAMS = {
    'centroid_noise_multiplier': 2.0, 'centroid_window_x_ppm': 0.02, 'centroid_window_y_ppm': 1.0,
    'distance_factor': 2.0, 'fix_linewidths': False, 'fix_positions': False,
    'height_threshold': 0.1, 'lw_gauss_15n': None, 'lw_gauss_1h': None,
    'lw_lorentz_15n': None, 'lw_lorentz_1h': None, 'max_optimization_iterations': 50,
    'max_peaks_fit': 50, 'prominence_threshold': 0.05, 'smoothing_sigma': 1.0,
    'use_adaptive_optimization': False, 'use_centroid_refinement': True,
    'use_parallel_processing': False, 'use_ps2d_linewidth_reuse': False,
    'use_ps2d_multi_peak': True,
}
FITTING_PARAMS = {'fitting_window_x': 0.15, 'fitting_window_y': 2.0,
                  'max_iterations': 1000, 'min_r_squared': 0.1}


@pytest.fixture
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


def test_the_manager_hands_the_integrator_its_default_parameter_set():
    from lunaNMR.core.ps2d_config import SEARCH_WINDOW_PPM, get_ps2d_config
    from lunaNMR.utils.parameter_manager import NMRParameterManager
    search_window_y = get_ps2d_config().search_window_f1
    assert NMRParameterManager().get_integrator_parameters() == {
        'detection_params': {'search_window_x': SEARCH_WINDOW_PPM, 'search_window_y': search_window_y},
        'fitting_params': FITTING_PARAMS,
        'gui_params': GUI_PARAMS,
        'processing_options': {'use_parallel_processing': False},
    }


def test_the_main_window_sync_hands_the_integrator_that_set(app):
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    set_ps2d_config('15N')
    window = LunaNMRMainWindow()
    try:
        window.integrator.gui_params = {}
        window.sync_parameters_to_integrator()
        assert window.integrator.gui_params == dict(
            GUI_PARAMS, centroid_window_x_ppm=0.01, centroid_window_y_ppm=0.1,
            use_parallel_processing=True)
        assert {k: window.integrator.fitting_parameters[k] for k in FITTING_PARAMS} == FITTING_PARAMS
        assert (window.integrator.search_window_x, window.integrator.search_window_y) == (0.07, 0.3)
    finally:
        window.close()
        window.deleteLater()
        set_ps2d_config('15N')


def test_there_is_no_simplified_parameter_module():
    with pytest.raises(ImportError):
        importlib.import_module('lunaNMR.utils.simplified_parameter_manager')


def test_the_parameter_manager_has_no_simplified_mode():
    from lunaNMR.utils.parameter_manager import NMRParameterManager
    manager = NMRParameterManager()
    for name in ('use_simplified_mode', 'simplified_manager', 'parameter_adapter',
                 'update_simplified_parameters', 'validate_simplified_parameters',
                 'get_effective_parameters'):
        assert not hasattr(manager, name), name


def test_the_fitter_has_no_simplified_parameter_manager():
    from lunaNMR.core.enhanced_voigt_fitter import EnhancedVoigtFitter
    fitter = EnhancedVoigtFitter()
    for name in ('simplified_param_manager', 'parameter_adapter'):
        assert not hasattr(fitter, name), name


def test_the_integrator_has_no_ppm_range_nucleus_detection():
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator, VoigtIntegrator
    assert not hasattr(VoigtIntegrator, 'detect_nucleus_type')
    assert not hasattr(EnhancedVoigtIntegrator, 'detect_nucleus_type')
