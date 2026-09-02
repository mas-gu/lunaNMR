# ABOUTME: The 2D fitter must honour a position bound given per peak, not one for all.
# ABOUTME: Real synthetic data and the real optimiser - a peak is driven at its neighbour.

import numpy as np
import pytest

from lunaNMR.core.ps2d_config import set_ps2d_config
from lunaNMR.core.ps2d_2d_fitter import Ps2dMultiPeakFitter2D
from lunaNMR.core.ps2d_style_fitter import voigt_profile_2d


LW_LOR_F1, LW_GAU_F1 = 0.18, 0.06
LW_LOR_F2, LW_GAU_F2 = 0.018, 0.006


@pytest.fixture(autouse=True)
def nucleus():
    set_ps2d_config('15N')


def _plane(peaks):
    """A 15N-HSQC-like patch holding exact 2D Voigts at (f1, f2, intensity)."""
    ppm_f2 = np.linspace(7.80, 7.60, 128)
    ppm_f1 = np.linspace(122.4, 120.6, 96)
    f2_grid, f1_grid = np.meshgrid(ppm_f2, ppm_f1)
    intensity = np.zeros_like(f1_grid)
    for f1, f2, amplitude in peaks:
        intensity += voigt_profile_2d(f1_grid, f2_grid, f1, LW_GAU_F1, f2,
                                      LW_GAU_F2, LW_LOR_F1, LW_LOR_F2, amplitude)
    return f1_grid, f2_grid, intensity


def _seed(f1, f2, intensity):
    return {'pos_f1': f1, 'lw_lor_f1': LW_LOR_F1, 'lw_gau_f1': LW_GAU_F1,
            'pos_f2': f2, 'lw_lor_f2': LW_LOR_F2, 'lw_gau_f2': LW_GAU_F2,
            'intensity': intensity}


# One real peak, and a second seeded 0.02 ppm away with nothing of its own left -
# the state a titration reaches when a peak decays. The optimiser has one maximum
# to explain and two peaks to explain it with, so the empty one is drawn onto it.
TRUTH = [(121.50, 7.700, 1.0)]
SEEDS = [_seed(121.50, 7.700, 1.0), _seed(121.50, 7.720, 0.05)]


def test_a_wide_window_lets_a_peak_reach_its_neighbour():
    """The premise. Without this the bound below proves nothing."""
    fitter = Ps2dMultiPeakFitter2D(verbose=False)
    result = fitter.fit_multi_peak_2d(*_plane(TRUTH), [dict(s) for s in SEEDS],
                                      pos_margin_f1=0.3, pos_margin_f2=0.03)

    assert result['success']
    moved = abs(result['peaks'][1]['pos_f2'] - 7.720)
    assert moved > 0.010, ('the empty peak stayed put, so this data does not '
                           f'reproduce the hop (moved {moved:.4f} ppm)')


def test_a_per_peak_bound_stops_it_halfway():
    """Given its own bound, the peak cannot reach the maximum that is not its own -
    so it cannot be reported at another peak's position."""
    fitter = Ps2dMultiPeakFitter2D(verbose=False)
    result = fitter.fit_multi_peak_2d(*_plane(TRUTH), [dict(s) for s in SEEDS],
                                      pos_margin_f1=0.3, pos_margin_f2=0.03,
                                      pos_margins=[(0.15, 0.0098), (0.15, 0.0098)])

    assert result['success']
    assert abs(result['peaks'][1]['pos_f2'] - 7.720) <= 0.0098 + 1e-9
    separation = abs(result['peaks'][0]['pos_f2'] - result['peaks'][1]['pos_f2'])
    assert separation > 0, 'the two peaks are coincident'


def test_a_per_peak_bound_is_exact_not_a_floor():
    """The scalar path widens a margin to at least 0.04 x FWHM. A caller that has
    computed a territory has already decided how far the peak may go, and widening
    it would let the territories overlap again."""
    fitter = Ps2dMultiPeakFitter2D(verbose=False)
    tight = 1e-5
    result = fitter.fit_multi_peak_2d(*_plane(TRUTH), [dict(s) for s in SEEDS],
                                      pos_margins=[(tight, tight), (tight, tight)])

    assert result['success']
    for fitted, seeded in zip(result['peaks'], SEEDS):
        assert abs(fitted['pos_f1'] - seeded['pos_f1']) <= tight + 1e-9
        assert abs(fitted['pos_f2'] - seeded['pos_f2']) <= tight + 1e-9
