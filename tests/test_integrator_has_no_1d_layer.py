# ABOUTME: Every peak is fitted in 2D, so the integrator carries no 1D cross-section fitting layer.
# ABOUTME: Pins the absence of its 1D fits, 1D detection, cross-section extraction and 1D helpers.

import importlib

import pytest

ONE_D_LAYER = (
    'get_adaptive_window_parameters', 'voigt_profile_1d', 'gaussian_profile_1d',
    '_estimate_peak_width', '_get_fitting_bounds', 'fit_peak_1d', 'fit_peak_1d_ps2d_style',
    '_calculate_initial_guess_1d', 'detect_peaks_1d', 'fit_multi_peak_1d', 'select_target_peak',
    'iterative_parameter_optimization', 'extract_peak_region', 'calculate_r_squared',
)


@pytest.mark.parametrize('name', ONE_D_LAYER)
def test_the_integrator_has_no(name):
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    assert not hasattr(EnhancedVoigtIntegrator, name)


def test_there_is_no_1d_guess_clamp():
    module = importlib.import_module('lunaNMR.core.core_integrator')
    assert not hasattr(module, '_clamp_guess_to_bounds')
