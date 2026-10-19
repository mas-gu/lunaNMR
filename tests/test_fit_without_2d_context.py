# ABOUTME: A peak the 2D fitter cannot place among its neighbours is an explicit failure that says why.
# ABOUTME: Real DNAJA1 data: the integrator called without a peak list, and the per-peak fallback of the cluster pipeline.

import contextlib
import io
from pathlib import Path

import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
REFERENCE_POINT = "2D_DNAJA1_0o0.ft"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def integrator():
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.utils.file_manager import NMRFileManager
    set_ps2d_config('15N')
    it = EnhancedVoigtIntegrator()
    with contextlib.redirect_stdout(io.StringIO()):
        it.load_nmr_file(str(DNAJA1 / REFERENCE_POINT))
    it.noise_level = it._estimate_noise_level_advanced()
    it.peak_list = NMRFileManager().load_peak_list(str(DNAJA1 / "DNAJA1.txt"))
    it.gui_params = {'use_ps2d_multi_peak': True, 'use_adaptive_optimization': False}
    return it


def _first_peak(integrator):
    row = integrator.peak_list.iloc[0]
    return float(row['Position_X']), float(row['Position_Y']), str(row['Assignment'])


class TestAFitWithoutAPeakContext:
    def test_is_an_explicit_failure_naming_the_missing_context(self, integrator):
        x, y, assignment = _first_peak(integrator)
        with contextlib.redirect_stdout(io.StringIO()):
            fit = integrator.fit_peak_voigt_2d(x, y, assignment)
        assert fit['success'] is False
        assert fit['assignment'] == assignment
        assert fit['error'].startswith("no 2D fit context: ")

    def test_through_enhanced_peak_fitting_is_not_fitted_and_the_log_says_why(self, integrator):
        x, y, assignment = _first_peak(integrator)
        log = io.StringIO()
        with contextlib.redirect_stdout(log):
            fit = integrator.enhanced_peak_fitting(x, y, assignment)
        assert fit is None
        assert f"{assignment} not fitted - no 2D fit context: " in log.getvalue()


class TestThePerPeakFallbackOfTheClusterPipeline:
    def test_fits_each_peak_in_2d_among_the_whole_peak_list(self, integrator):
        from lunaNMR.core.parallel_voigt_processor import ParallelVoigtProcessor
        peaks = integrator.peak_list.iloc[:4]
        processor = ParallelVoigtProcessor(integrator.enhanced_fitter, in_process=True)
        with contextlib.redirect_stdout(io.StringIO()):
            results = processor._sequential_fallback(peaks, None)
        assert [r['assignment'] for r in results] == [str(a) for a in peaks['Assignment']]
        assert all(r['success'] for r in results), [r.get('failure_reason') for r in results]
        assert {r['method'] for r in results} == {'2d_simultaneous_multi_peak'}
