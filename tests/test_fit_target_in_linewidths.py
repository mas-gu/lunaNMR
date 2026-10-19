# ABOUTME: A peak fitted in 2D among its neighbours takes the fit nearest it in linewidths, not in raw ppm.
# ABOUTME: The 49SerH/52TyrH pair of DNAJA1 0o0, fitted through fit_peak_voigt_2d and as its two fitted positions.

"""
fit_peak_voigt_2d fits the target's whole overlap group and reports the fit nearest the
target. 52TyrH sits 0.034 ppm from 49SerH in 1H and 0.143 ppm in 15N: nearly a 1H
linewidth, a third of a 15N one. A target a fifth of a 15N linewidth from 49SerH's fit,
toward 52TyrH, lay nearer 52TyrH's fit in raw ppm and was given it.
"""

import contextlib
import io
from pathlib import Path

import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"

# The two fits of the pair on DNAJA1 0o0, as fit_peak_voigt_2d returns them for 49SerH
SER49 = (7.999579, 116.560248)
TYR52 = (8.034129, 116.417319)


@pytest.fixture
def nucleus():
    from lunaNMR.core.ps2d_config import set_ps2d_config
    set_ps2d_config('15N')
    yield
    set_ps2d_config('15N')


def _fits(*positions):
    return [{'pos_f2': x, 'pos_f1': y} for x, y in positions]


def test_a_target_takes_the_fit_nearest_in_linewidths(nucleus):
    from lunaNMR.core.core_integrator import nearest_fit
    assert nearest_fit(_fits(SER49, TYR52), SER49[0], SER49[1] - 0.08) == 0
    assert nearest_fit(_fits(TYR52, SER49), SER49[0], SER49[1] - 0.08) == 1


def test_an_exact_target_takes_its_own_fit(nucleus):
    from lunaNMR.core.core_integrator import nearest_fit
    assert nearest_fit(_fits(SER49, TYR52), *TYR52) == 1


@pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")
def test_each_member_of_the_pair_reports_its_own_fit(nucleus):
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    from lunaNMR.core.parallel_voigt_processor import ParallelVoigtProcessor
    from lunaNMR.utils.file_manager import NMRFileManager
    it = EnhancedVoigtIntegrator()
    with contextlib.redirect_stdout(io.StringIO()):
        it.load_nmr_file(str(DNAJA1 / "2D_DNAJA1_0o0.ft"))
    it.noise_level = it._estimate_noise_level_advanced()
    it.peak_list = NMRFileManager().load_peak_list(str(DNAJA1 / "DNAJA1.txt"))
    it.gui_params = {'use_ps2d_multi_peak': True, 'use_adaptive_optimization': False}
    context = ParallelVoigtProcessor(it.enhanced_fitter, in_process=True)._peaks_context(it.peak_list)
    listed = it.peak_list.set_index('Assignment')
    for assignment, fitted in (('49SerH', SER49), ('52TyrH', TYR52)):
        x, y = float(listed.loc[assignment, 'Position_X']), float(listed.loc[assignment, 'Position_Y'])
        with contextlib.redirect_stdout(io.StringIO()):
            fit = it.fit_peak_voigt_2d(x, y, assignment, all_peaks_context=context)
        assert fit['peak_position'] == pytest.approx(fitted, abs=1e-5)
