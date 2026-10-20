# ABOUTME: Every peak is fitted in 2D, so the Voigt fitter carries no 1D multi-peak fitting chain.
# ABOUTME: Pins the absence of its strategies, AIC peak counting, constraint coupling and parameter estimator.

import importlib

import pytest

ONE_D_CHAIN = (
    '_clamp_vector_to_bounds', '_evaluate_parameter_sanity', '_fallback_asymmetric_baseline',
    '_fit_multiple_peaks_simultaneously', 'apply_constraint_coupling',
    'assess_fit_quality_comprehensive', 'calculate_aic_for_peak_count',
    'calculate_parameter_uncertainties', 'enhanced_multi_peak_parameter_validation',
    'estimate_initial_parameters_from_resolved_peaks', 'fit_pairwise_then_combine',
    'fit_single_peak_with_expanded_context', 'fit_with_global_constraints',
    'fit_with_progressive_strategies', 'get_adaptive_bounds', 'monitored_curve_fit',
    'optimal_peak_count_by_aic', 'optimize_overlap_detection_iteratively', 'refine_peak_center',
    'sequential_peak_isolation', 'set_gui_parameters', 'validate_and_fix_initial_guess',
    'validate_initial_parameters_and_bounds', 'parameter_estimator',
    'detect_overlapping_peaks', 'get_adaptive_peak_detection_params',
    'comprehensive_quality_assessment',
)


@pytest.mark.parametrize('name', ONE_D_CHAIN)
def test_the_fitter_has_no(name):
    from lunaNMR.core.enhanced_voigt_fitter import EnhancedVoigtFitter
    assert not hasattr(EnhancedVoigtFitter(), name)


def test_there_is_no_parameter_estimator_class():
    module = importlib.import_module('lunaNMR.core.enhanced_voigt_fitter')
    assert not hasattr(module, 'RobustParameterEstimator')
