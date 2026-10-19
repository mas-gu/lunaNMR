# ABOUTME: The base class of the live integrator holds only what a live path calls: loading,
# ABOUTME: ppm axes, nucleus detection and noise. Its detection, integration and plotting are gone.

"""
The live integrator (EnhancedVoigtIntegrator) overrides the base's detection and its
statistics, and nothing calls the base's own reference/traditional detection, box
integration, export or summary plot. A method that is defined but never reached reads as
the implementation of whatever its name says, so these are not kept.
"""

import importlib.util

import pytest

GONE = ['set_search_window', 'set_threshold_multiplier', 'set_noise_regions', 'ppm_to_point',
        'point_to_ppm', 'detect_peaks_reference_based', '_is_local_maximum',
        '_create_failed_detection', '_consolidate_y_dimension_peaks',
        '_refine_peaks_with_centroids', '_find_all_peaks_in_window',
        'detect_peaks_professional', 'detect_peaks_traditional', '_match_traditional',
        'integrate_peaks', 'get_detection_statistics', 'export_results', 'create_summary_plot']

KEPT = ['load_data', '_detect_nmr_format', '_calculate_ppm_axes', '_detect_nucleus_type',
        '_estimate_noise_level', '_estimate_noise_level_corners']


@pytest.mark.parametrize("name", GONE)
def test_the_base_does_not_define_an_unreached_method(name):
    from lunaNMR.integrators.inplace_advanced_nmr_integrator import InPlaceAdvancedNMRIntegrator
    assert name not in vars(InPlaceAdvancedNMRIntegrator)


@pytest.mark.parametrize("name", KEPT)
def test_the_base_keeps_what_the_live_integrator_calls(name):
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.integrators.inplace_advanced_nmr_integrator import InPlaceAdvancedNMRIntegrator
    assert getattr(EnhancedVoigtIntegrator, name) is vars(InPlaceAdvancedNMRIntegrator)[name]


def test_the_series_integrator_that_nothing_builds_is_gone():
    assert importlib.util.find_spec("lunaNMR.integrators.inplace_series_nmr_integrator") is None


def test_the_package_exports_the_base_class():
    import lunaNMR.integrators as integrators
    from lunaNMR.integrators.inplace_advanced_nmr_integrator import InPlaceAdvancedNMRIntegrator
    assert integrators.__all__ == ['InPlaceAdvancedNMRIntegrator']
    assert integrators.InPlaceAdvancedNMRIntegrator is InPlaceAdvancedNMRIntegrator


def test_the_main_window_has_no_integrate_action_without_a_menu():
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    assert 'integrate_peaks' not in vars(LunaNMRMainWindow)
