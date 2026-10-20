# ABOUTME: No parameter set carries a noise threshold: nothing on any detection or fit path reads one.
# ABOUTME: The CLI defaults, the parameter manager, the configuration file and the parallel workers.


def test_the_cli_series_defaults_carry_none():
    from lunaNMR.cli import _default_series_params
    assert 'noise_threshold' not in _default_series_params()['detection_params']


def test_the_parameter_manager_carries_none():
    from lunaNMR.utils.parameter_manager import NMRParameterManager
    manager = NMRParameterManager()
    assert 'noise_threshold' not in manager.parameter_definitions
    assert 'noise_threshold' not in manager.get_parameter_summary()
    assert 'noise_threshold' not in manager.get_integrator_parameters()['detection_params']


def test_the_configuration_file_defaults_carry_none():
    from lunaNMR.utils.config_manager import ConfigurationManager, ProcessingParameters
    config = ConfigurationManager("lunaNMR_v1")
    assert 'noise_threshold' not in config.default_config['processing_parameters']
    assert 'noise_threshold' not in ProcessingParameters(config).get_detection_params()


def test_the_parallel_workers_are_handed_none():
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.core.parallel_voigt_processor import ParallelVoigtProcessor
    integrator = EnhancedVoigtIntegrator()
    processor = ParallelVoigtProcessor(integrator.enhanced_fitter, in_process=True)
    integrator.noise_threshold = 3.0
    assert 'noise_threshold' not in processor._serialize_integrator_params()
