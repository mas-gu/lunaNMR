# ABOUTME: The volume measured by summing the DATA over the peak's region, not by
# ABOUTME: integrating the fitted Voigt to infinity. Reported beside it, never instead.

"""
Region-sum volume.

`volume` is the analytic integral of the fitted Voigt over all space, so it is
`height x width_F1 x width_F2` to within a shape factor. When a peak weakens the width
stops being determined by the data - the chi-squared surface goes flat in that direction
and the Lorentzian/Gaussian split even swaps roles between adjacent points of one series -
but it still multiplies straight into the answer, and nothing downstream detects it.
Measured on the DNAJA1/HSPA8 titration, the last point of two peaks came back 2.4x and
2.9x what their heights imply, on EVERY peak.

CcpNmr's own `Peak.estimateVolume` is the same KIND of number - a Gaussian of the fitted
lineWidth, integrated over +-1 FWHM - so it shares the failure exactly. Switching to it
would not help.

This sums the spectrum instead:

    baseline    a plane fitted to the region OUTSIDE the peak mask, then subtracted.
                The 100 uM point of that titration is a tilted sheet with no peak in it;
                remove the plane and there is nothing left to sum, which is the truth.
    total       sum of the corrected data inside the mask, x the ppm area of one cell
    apportion   split between cluster members by each fitted model's integral over the
                SAME mask - the one job only a model can do

It is reported beside the analytic volume, never instead of it. For an isolated peak with
good S/N the analytic integral is the better estimator: it weights every point through the
model and recovers the real tails, where a sum truncates them and accumulates noise. Their
RATIO is the diagnostic - near 1 the model describes what is there, far from 1 it does not.
"""

import numpy as np
import pytest

from lunaNMR.core.ps2d_2d_fitter import region_sum_volumes, Ps2dMultiPeakFitter2D


def _grids(n=40, f1c=120.0, f2c=8.0, f1r=0.30, f2r=0.045):
    f1 = np.linspace(f1c - f1r, f1c + f1r, n)
    f2 = np.linspace(f2c - f2r, f2c + f2r, n)
    F2, F1 = np.meshgrid(f2, f1)
    return F1, F2


def _gaussian(F1, F2, height, f1c=120.0, f2c=8.0, fwhm1=0.16, fwhm2=0.024):
    s1 = fwhm1 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    s2 = fwhm2 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return height * np.exp(-0.5 * ((F1 - f1c) / s1) ** 2 - 0.5 * ((F2 - f2c) / s2) ** 2)


def _params(peaks):
    out = []
    for pos1, pos2, lw1, lw2, intensity in peaks:
        out.extend([pos1, 0.3 * lw1, 0.7 * lw1, pos2, 0.3 * lw2, 0.7 * lw2,
                    intensity, 0.0])
    return np.array(out, dtype=float)


class TestItMeasuresRatherThanExtrapolates:

    def test_a_peak_of_half_the_height_gives_half_the_volume(self):
        """The whole point: the number follows the DATA. Whatever the fitted width
        does, twice the signal is twice the volume."""
        F1, F2 = _grids()
        mask = np.ones(F1.shape, dtype=bool)
        params = _params([(120.0, 8.0, 0.16, 0.024, 1.0)])

        tall = region_sum_volumes(F1, F2, _gaussian(F1, F2, 1.0), mask, params, 1)
        half = region_sum_volumes(F1, F2, _gaussian(F1, F2, 0.5), mask, params, 1)
        assert half[0] == pytest.approx(0.5 * tall[0], rel=1e-6)

    def test_the_fitted_width_does_not_enter_the_answer(self):
        """The failure being fixed. Doubling every fitted linewidth doubles the
        ANALYTIC volume and must leave this one untouched."""
        F1, F2 = _grids()
        mask = np.ones(F1.shape, dtype=bool)
        data = _gaussian(F1, F2, 1.0)

        narrow = region_sum_volumes(F1, F2, data, mask,
                                    _params([(120.0, 8.0, 0.16, 0.024, 1.0)]), 1)
        wide = region_sum_volumes(F1, F2, data, mask,
                                  _params([(120.0, 8.0, 0.32, 0.048, 1.0)]), 1)
        assert wide[0] == pytest.approx(narrow[0], rel=1e-9)

    def test_it_recovers_the_true_integral_of_a_peak_inside_its_window(self):
        """Sanity against the analytic answer for a Gaussian: 2*pi*s1*s2*height."""
        F1, F2 = _grids(n=120, f1r=0.6, f2r=0.09)
        mask = np.ones(F1.shape, dtype=bool)
        s1 = 0.16 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        s2 = 0.024 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        truth = 2.0 * np.pi * s1 * s2 * 1.0

        got = region_sum_volumes(F1, F2, _gaussian(F1, F2, 1.0), mask,
                                 _params([(120.0, 8.0, 0.16, 0.024, 1.0)]), 1)
        assert got[0] == pytest.approx(truth, rel=0.02)


class TestTheBaselineIsRemoved:

    def test_a_tilted_sheet_with_no_peak_sums_to_nothing(self):
        """The 100 uM point of the real titration. The fitted Voigt reports 79% of the
        reference volume; the spectrum contains a ramp and no peak."""
        F1, F2 = _grids()
        ramp = 0.4 + 3.0 * (F1 - 120.0) + 20.0 * (F2 - 8.0)
        mask = (((F1 - 120.0) / 0.15) ** 2 + ((F2 - 8.0) / 0.0225) ** 2) <= 1.0
        params = _params([(120.0, 8.0, 0.16, 0.024, 1.0)])

        got = region_sum_volumes(F1, F2, ramp, mask, params, 1)
        scale = abs(ramp).max() * F1.size
        assert abs(got[0]) < 0.01 * scale, f'a pure ramp summed to {got[0]:.4g}'

    def test_a_peak_on_a_ramp_still_measures_the_peak(self):
        F1, F2 = _grids()
        mask = (((F1 - 120.0) / 0.15) ** 2 + ((F2 - 8.0) / 0.0225) ** 2) <= 1.0
        params = _params([(120.0, 8.0, 0.16, 0.024, 1.0)])

        peak = _gaussian(F1, F2, 1.0)
        ramp = 0.4 + 3.0 * (F1 - 120.0) + 20.0 * (F2 - 8.0)
        clean = region_sum_volumes(F1, F2, peak, mask, params, 1)
        tilted = region_sum_volumes(F1, F2, peak + ramp, mask, params, 1)
        assert tilted[0] == pytest.approx(clean[0], rel=0.05)


class TestClustersAreApportionedByTheModel:

    def test_two_equal_peaks_split_the_sum_evenly(self):
        F1, F2 = _grids(n=60)
        mask = np.ones(F1.shape, dtype=bool)
        data = (_gaussian(F1, F2, 1.0, f2c=7.99)
                + _gaussian(F1, F2, 1.0, f2c=8.01))
        params = _params([(120.0, 7.99, 0.16, 0.024, 1.0),
                          (120.0, 8.01, 0.16, 0.024, 1.0)])

        got = region_sum_volumes(F1, F2, data, mask, params, 2)
        assert got[0] == pytest.approx(got[1], rel=1e-6)

    def test_a_three_to_one_model_splits_three_to_one(self):
        F1, F2 = _grids(n=60)
        mask = np.ones(F1.shape, dtype=bool)
        data = (_gaussian(F1, F2, 3.0, f2c=7.99)
                + _gaussian(F1, F2, 1.0, f2c=8.01))
        params = _params([(120.0, 7.99, 0.16, 0.024, 3.0),
                          (120.0, 8.01, 0.16, 0.024, 1.0)])

        got = region_sum_volumes(F1, F2, data, mask, params, 2)
        assert got[0] / got[1] == pytest.approx(3.0, rel=0.05)

    def test_the_shares_add_up_to_the_measured_total(self):
        """Whatever the split, the parts must sum to what was actually in the region -
        that is the property the model can never break."""
        F1, F2 = _grids(n=60)
        mask = np.ones(F1.shape, dtype=bool)
        data = (_gaussian(F1, F2, 3.0, f2c=7.99)
                + _gaussian(F1, F2, 1.0, f2c=8.01))
        params = _params([(120.0, 7.99, 0.16, 0.024, 3.0),
                          (120.0, 8.01, 0.16, 0.024, 1.0)])

        parts = region_sum_volumes(F1, F2, data, mask, params, 2)
        whole = region_sum_volumes(F1, F2, data, mask,
                                   _params([(120.0, 8.0, 0.16, 0.048, 4.0)]), 1)
        assert sum(parts) == pytest.approx(whole[0], rel=1e-9)


class TestItRefusesRatherThanGuesses:

    def test_an_empty_mask_returns_nothing(self):
        F1, F2 = _grids()
        params = _params([(120.0, 8.0, 0.16, 0.024, 1.0)])
        got = region_sum_volumes(F1, F2, _gaussian(F1, F2, 1.0),
                                 np.zeros(F1.shape, dtype=bool), params, 1)
        assert got == [None]

    def test_a_model_with_no_intensity_anywhere_returns_nothing(self):
        """Nothing to apportion by. A silent even split would invent a measurement."""
        F1, F2 = _grids()
        mask = np.ones(F1.shape, dtype=bool)
        params = _params([(120.0, 8.0, 0.16, 0.024, 0.0),
                          (120.0, 8.0, 0.16, 0.024, 0.0)])
        assert region_sum_volumes(F1, F2, _gaussian(F1, F2, 1.0),
                                  mask, params, 2) == [None, None]

    def test_nan_data_does_not_propagate(self):
        F1, F2 = _grids()
        mask = np.ones(F1.shape, dtype=bool)
        data = _gaussian(F1, F2, 1.0)
        data[0, 0] = np.nan
        got = region_sum_volumes(F1, F2, data, mask,
                                 _params([(120.0, 8.0, 0.16, 0.024, 1.0)]), 1)
        assert got[0] is not None and np.isfinite(got[0])


class TestItReachesTheFitResult:

    def test_every_fitted_peak_carries_a_region_volume(self):
        F1, F2 = _grids(n=32)
        data = _gaussian(F1, F2, 1.0)
        peak = [{'pos_f1': 120.0, 'lw_lor_f1': 0.03, 'lw_gau_f1': 0.13,
                 'pos_f2': 8.0, 'lw_lor_f2': 0.005, 'lw_gau_f2': 0.019,
                 'intensity': 1.0}]
        result = Ps2dMultiPeakFitter2D(verbose=False).fit_multi_peak_2d(
            F1, F2, data, peak)
        assert 'region_volume' in result['peaks'][0]
        assert result['peaks'][0]['region_volume'] > 0

    def test_the_analytic_volume_is_still_there_unchanged(self):
        """Beside, never instead. For a good isolated fit the analytic integral is the
        better estimator and must not be lost."""
        F1, F2 = _grids(n=32)
        data = _gaussian(F1, F2, 1.0)
        peak = [{'pos_f1': 120.0, 'lw_lor_f1': 0.03, 'lw_gau_f1': 0.13,
                 'pos_f2': 8.0, 'lw_lor_f2': 0.005, 'lw_gau_f2': 0.019,
                 'intensity': 1.0}]
        got = Ps2dMultiPeakFitter2D(verbose=False).fit_multi_peak_2d(
            F1, F2, data, peak)['peaks'][0]
        assert got['volume'] > 0 and got['region_volume'] > 0
