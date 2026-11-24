# ABOUTME: The fits that teach a spectrum its linewidths are those at R-squared 0.8 or more, in the pipeline and the overlay refit.
# ABOUTME: Drives the parallel processor's statistics on fit results at chosen R-squared values.

import pytest


def _isolated_fits(r_squared, count=4):
    fit = {'success': True, 'r_squared': r_squared,
           'x_fit': {'sigma': 0.010, 'gamma': 0.012},
           'y_fit': {'sigma': 0.10, 'gamma': 0.12}}
    return [{'success': True, 'peak_results': [dict(fit)]} for _ in range(count)]


@pytest.fixture
def processor():
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.core.parallel_voigt_processor import ParallelVoigtProcessor
    set_ps2d_config('15N')
    return ParallelVoigtProcessor(object())


def test_fits_at_the_gate_teach_the_linewidths(processor):
    from lunaNMR.utils.quality_categories import GOOD_FIT_MIN_R2
    assert GOOD_FIT_MIN_R2 == 0.8
    statistics = processor._collect_linewidth_statistics(_isolated_fits(0.82))
    assert statistics['n_samples'] == 4
    assert statistics['lw_f2_median'] == pytest.approx(0.022)


def test_fits_below_the_gate_do_not(processor):
    assert processor._collect_linewidth_statistics(_isolated_fits(0.79)) is None
