# ABOUTME: A single-spectrum fit has one strategy, the cluster pipeline: no processing option
# ABOUTME: diverts it, and the per-peak "global optimisation" survey it once could take is gone.

import importlib
from pathlib import Path

import pytest

HSQC = Path(__file__).resolve().parents[1] / "data_example" / "2DHSQC"


@pytest.fixture(scope="module")
def fit():
    """fit(processing_options) -> the results of fitting six 600 MHz HSQC peaks."""
    if not (HSQC / "600_T1_0o0.ft").exists():
        pytest.skip("data_example 2DHSQC not present")
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.processors.single_spectrum_processor import SingleSpectrumProcessor
    from lunaNMR.utils.file_manager import NMRFileManager
    from lunaNMR.utils.parameter_manager import NMRParameterManager

    peaks = NMRFileManager().load_peak_list(str(HSQC / "600_peaks.txt")).head(6)

    def run(processing_options):
        integrator = EnhancedVoigtIntegrator()
        integrator.load_nmr_file(str(HSQC / "600_T1_0o0.ft"))
        integrator.peak_list = peaks
        processor = SingleSpectrumProcessor(integrator, NMRParameterManager())
        results, _ = processor.process_peak_list(peaks, processing_options)
        return [(r['assignment'], r['height'], r['volume']) for r in results]

    return run


def test_a_global_optimisation_flag_fits_as_the_cluster_pipeline_does(fit):
    """A project's gui_state.json can still carry use_global_optimization: true. It must
    fit the spectrum, not route it to a per-peak survey that, having no peak list to fit
    each peak in 2D among, returned no results at all."""
    plain = fit({'use_parallel': False})
    flagged = fit({'use_parallel': False, 'use_global_optimization': True})
    assert len(plain) == 6
    assert flagged == plain


def test_the_global_optimisation_survey_is_gone():
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.processors.single_spectrum_processor import SingleSpectrumProcessor
    import lunaNMR.utils

    assert not hasattr(EnhancedVoigtIntegrator, 'optimize_peak_list_globally')
    assert not hasattr(SingleSpectrumProcessor, '_process_with_global_optimization')
    assert 'GlobalOptimizationManager' not in lunaNMR.utils.__all__
    with pytest.raises(ImportError):
        importlib.import_module('lunaNMR.utils.global_optimization_manager')


def test_a_peak_fit_takes_no_linewidth_constraints():
    """The 2D fit never read them; the survey was the only caller that passed any."""
    import inspect
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator

    for method in (EnhancedVoigtIntegrator.enhanced_peak_fitting,
                   EnhancedVoigtIntegrator.fit_peak_voigt_2d):
        assert 'linewidth_constraints' not in inspect.signature(method).parameters
