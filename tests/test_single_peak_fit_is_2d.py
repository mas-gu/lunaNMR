# ABOUTME: A single peak is fitted in 2D whatever the parameter mode, and nothing carries an overlap-resolution switch.
# ABOUTME: Real 600 MHz T1 reference plane from data_example; the integrator, its fitter, the workers and the config file.

import contextlib
import io
from pathlib import Path

import pytest

SPECTRA = Path(__file__).resolve().parents[1] / "data_example" / "2DHSQC"
REFERENCE_POINT = "600_T1_0o0.ft"

pytestmark = pytest.mark.skipif(not (SPECTRA / REFERENCE_POINT).exists(),
                                reason="data_example/2DHSQC not present")


@pytest.fixture(scope="module")
def integrator():
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.utils.file_manager import NMRFileManager
    set_ps2d_config('15N')
    it = EnhancedVoigtIntegrator()
    with contextlib.redirect_stdout(io.StringIO()):
        it.load_nmr_file(str(SPECTRA / REFERENCE_POINT))
    it.noise_level = it._estimate_noise_level_advanced()
    it.peak_list = NMRFileManager().load_peak_list(str(SPECTRA / "600_peaks.txt"))
    return it


def _fit_first_peak(integrator, gui_params):
    """(fit, captured stdout) of the first listed peak, fitted among the whole list."""
    context = [{'assignment': str(r['Assignment']), 'x_ppm': float(r['Position_X']),
                'y_ppm': float(r['Position_Y'])} for _, r in integrator.peak_list.iterrows()]
    row = integrator.peak_list.iloc[0]
    integrator.gui_params = gui_params
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        fit = integrator.enhanced_peak_fitting(
            float(row['Position_X']), float(row['Position_Y']), str(row['Assignment']),
            all_peaks_context=context)
    return fit, log.getvalue()


class TestASinglePeakWithSimplifiedParameters:
    BASE = {'use_ps2d_multi_peak': True, 'use_adaptive_optimization': False}

    def test_is_fitted_in_2d_exactly_as_without_them(self, integrator):
        plain, _ = _fit_first_peak(integrator, dict(self.BASE))
        simplified, _ = _fit_first_peak(integrator, dict(self.BASE, use_simplified_parameters=True))
        assert simplified['method'] == '2d_simultaneous_multi_peak'
        for key in ('peak_x', 'peak_y', 'height', 'volume', 'avg_r_squared'):
            assert simplified[key] == plain[key], key

    def test_logs_no_other_fitting_method(self, integrator):
        _, log = _fit_first_peak(integrator, dict(self.BASE, use_simplified_parameters=True))
        assert 'onsensus' not in log
        assert 'utomated fitting' not in log


def test_the_integrator_has_no_overlap_resolution_switch(integrator):
    for name in ('configure_overlap_resolution', 'enable_overlap_resolution_preset',
                 'overlap_resolution_enabled'):
        assert not hasattr(integrator, name), name


def test_the_fitter_has_no_consensus_or_overlap_routing(integrator):
    for name in ('fit_with_consensus', 'use_automated_fitting', '_convert_consensus_to_legacy',
                 '_legacy_fit_fallback', 'detect_potential_overlap'):
        assert not hasattr(integrator.enhanced_fitter, name), name


def test_a_worker_is_handed_no_overlap_resolution_settings(integrator):
    from lunaNMR.core.parallel_voigt_processor import (
        ParallelVoigtProcessor, _initialize_worker_fitter)
    fitter = integrator.enhanced_fitter
    fitter.nmr_data, fitter.ppm_x_axis, fitter.ppm_y_axis = (
        integrator.nmr_data, integrator.ppm_x_axis, integrator.ppm_y_axis)
    integrator.gui_params = {'use_ps2d_multi_peak': True, 'use_adaptive_optimization': False}
    processor = ParallelVoigtProcessor(fitter, enable_adaptive=False)
    try:
        shared = processor._prepare_shared_context([])
        worker = _initialize_worker_fitter(shared)
    finally:
        processor._cleanup_shared_memory()
    assert 'overlap_resolution_params' not in shared
    for name in ('overlap_detection_enabled', 'overlap_detection_threshold', 'overlap_config'):
        assert not hasattr(worker.enhanced_fitter, name), name


def test_the_configuration_file_has_no_overlap_resolution_section():
    from lunaNMR.utils.config_manager import ConfigurationManager
    config = ConfigurationManager("lunaNMR_v1")
    assert 'overlap_resolution' not in config.default_config
    assert not hasattr(config, 'load_overlap_config')
    assert not hasattr(config, 'enable_overlap_resolution')
