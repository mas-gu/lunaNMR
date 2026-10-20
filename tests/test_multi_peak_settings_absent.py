# ABOUTME: The 1D multi-peak detection settings went with the 1D fitting chain that read them.
# ABOUTME: No parameter set, fitter or default configuration carries them.

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
    assert key not in EnhancedVoigtFitter().fitting_parameters


@pytest.mark.parametrize('key', KEYS)
def test_the_default_configuration_does_not_carry(key, tmp_path, monkeypatch):
    from lunaNMR.utils.config_manager import ConfigurationManager, ProcessingParameters
    monkeypatch.setattr(ConfigurationManager, '_get_config_folder', lambda self: str(tmp_path))
    config = ConfigurationManager("lunaNMR_test")
    assert key not in config.default_config['processing_parameters']
    assert key not in ProcessingParameters(config).get_fitting_params()
