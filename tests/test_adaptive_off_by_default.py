# ABOUTME: The adaptive radius search is off unless asked for.
# ABOUTME: On real titrations it always chose the smallest window, which broke CSP positions.

"""
The search tries fitting windows of 1.0, 1.5 and 2.0 x the median linewidth and keeps the
one whose residual is closest to the noise. On DNAJA1 and ERDJ6 it chose 1.0x both times -
a window that cuts the peak at about half height - and ERDJ6's CSP Kd went from a
credible 282 uM to pinned at the 1500 uM bound. With the search off, both are credible.
"""

import contextlib
import io
from pathlib import Path

import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"


def test_the_parameter_manager_leaves_it_off():
    from lunaNMR.utils.parameter_manager import NMRParameterManager
    params = NMRParameterManager().get_integrator_parameters()
    assert params['gui_params']['use_adaptive_optimization'] is False


def test_a_processor_leaves_it_off():
    from lunaNMR.core.parallel_voigt_processor import ParallelVoigtProcessor
    assert ParallelVoigtProcessor(object()).enable_adaptive is False


@pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")
def test_a_fit_that_does_not_ask_for_it_does_not_run_it():
    """The CLI's parameters do not mention it at all."""
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.utils.file_manager import NMRFileManager
    it = EnhancedVoigtIntegrator()
    with contextlib.redirect_stdout(io.StringIO()):
        it.load_nmr_file(str(DNAJA1 / "2D_DNAJA1_0o0.ft"))
    it.noise_level = it._estimate_noise_level_advanced()
    it.gui_params = {'use_ps2d_multi_peak': True}
    peaks = NMRFileManager().load_peak_list(str(DNAJA1 / "DNAJA1.txt"))
    it.enhanced_fitter.enhanced_peak_fitting_parallel(peaks, parent_integrator=it,
                                                      in_process=True)
    assert it.enhanced_fitter.parallel_processor.optimal_params is None
