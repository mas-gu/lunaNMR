# ABOUTME: The 1D multi-peak detection settings went with the 1D fitting chain that read them.
# ABOUTME: No parameter set or fitter carries them, and a project saves none of them.

import pytest

KEYS = ('multi_peak_r2_threshold', 'multi_peak_improvement_threshold',
        'peak_detection_sensitivity', 'overlap_detection_factor', 'residual_analysis_threshold')


@pytest.mark.parametrize('key', KEYS)
def test_the_parameter_manager_does_not_define(key):
    from lunaNMR.utils.parameter_manager import NMRParameterManager
    manager = NMRParameterManager()
    assert key not in manager.parameter_definitions
    assert key not in manager.current_params


@pytest.mark.parametrize('key', KEYS)
def test_the_fitter_does_not_carry(key):
    from lunaNMR.core.enhanced_voigt_fitter import EnhancedVoigtFitter
    fitter = EnhancedVoigtFitter()
    assert not hasattr(fitter, key)
    assert not any(key in value for value in vars(fitter).values() if isinstance(value, dict))


@pytest.mark.parametrize('key', KEYS)
def test_a_project_does_not_save(key):
    from lunaNMR.utils.project_manager import ProjectManager
    assert key not in ProjectManager.GUI_PARAMS
