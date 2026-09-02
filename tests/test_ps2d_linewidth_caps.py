# ABOUTME: The fitted linewidth TOTAL is what the volume depends on, and what must be bounded.
# ABOUTME: Covers the sum cap and the per-peak growth cap passed through fit_multi_peak_2d.

"""
Linewidth caps in the 2D fitter.

The volume PS2D reports is the ANALYTIC integral of the fitted 2D Voigt, so a profile
wider than the region it was fitted in contributes a volume made mostly of tails
nothing measured - and neither R-squared (a squared Pearson correlation, shape only)
nor the residual over the elliptical mask can see it.

Measured on a real titration (DNAJA1/HSPA8, peak 73LysH, 7 points): the height fell to
27% of the reference while both linewidth totals grew ~1.7x, and the reported volume
came back to 94% of the reference. Two defects let that happen:

  * a bound computed on the TOTAL was applied to each COMPONENT. Its median/MAD come
    from `lw_lor + lw_gau`, so a bound meant for the total permitted twice it.
  * nothing bounded the growth from one point of a series to the next. The failure was
    one jump - per-step growth ran 1.00, 1.02, 1.03, 1.15, 1.37, 1.22 - so an absolute
    cap alone does not catch it.
"""

import numpy as np
import pytest

from lunaNMR.core.ps2d_2d_fitter import Ps2dMultiPeakFitter2D, cap_linewidth_totals


NPAR = 8


def _params(pairs):
    """A parameter vector from [(lor_f1, gau_f1, lor_f2, gau_f2), ...]."""
    params = []
    for lor1, gau1, lor2, gau2 in pairs:
        params.extend([120.0, lor1, gau1, 8.0, lor2, gau2, 1.0, 0.0])
    return np.array(params, dtype=float)


class TestTheTotalIsWhatIsCapped:

    def test_a_total_inside_its_cap_is_left_exactly_alone(self):
        """The blast radius has to be the pathological fits and nothing else: a fit
        under the cap must come back bit-identical, or every existing result moves."""
        params = _params([(0.03, 0.13, 0.005, 0.019)])
        capped, _ = cap_linewidth_totals(params.copy(), 1, [0.24], [0.036])
        assert np.array_equal(capped, params)

    def test_a_total_over_its_cap_is_scaled_back_to_it(self):
        params = _params([(0.15, 0.18, 0.005, 0.019)])       # F1 total 0.33
        capped, _ = cap_linewidth_totals(params, 1, [0.24], [0.036])
        assert capped[1] + capped[2] == pytest.approx(0.24)

    def test_the_split_survives_the_scaling(self):
        """The lor/gau split is not determined by the data in a narrow window - the
        total is stiff and the split is soft - so capping must not decide it. Scaling
        the pair proportionally leaves the ratio exactly where the fit put it."""
        params = _params([(0.11, 0.22, 0.005, 0.019)])       # F1 total 0.33, 1:2
        capped, _ = cap_linewidth_totals(params, 1, [0.24], [0.036])
        assert capped[1] / capped[2] == pytest.approx(0.5)

    def test_each_dimension_is_capped_independently(self):
        params = _params([(0.15, 0.18, 0.030, 0.030)])       # F1 0.33, F2 0.06
        capped, _ = cap_linewidth_totals(params, 1, [0.24], [0.036])
        assert capped[1] + capped[2] == pytest.approx(0.24)
        assert capped[4] + capped[5] == pytest.approx(0.036)

    def test_each_peak_gets_its_own_cap(self):
        params = _params([(0.15, 0.18, 0.005, 0.019),
                          (0.02, 0.04, 0.005, 0.019)])
        capped, _ = cap_linewidth_totals(params, 2, [0.24, 0.10], [0.036, 0.036])
        assert capped[1] + capped[2] == pytest.approx(0.24)
        assert capped[NPAR + 1] + capped[NPAR + 2] == pytest.approx(0.06)

    def test_a_missing_cap_leaves_that_peak_untouched(self):
        """None means "no information", which must never be read as a cap of zero."""
        params = _params([(0.15, 0.18, 0.005, 0.019)])
        capped, _ = cap_linewidth_totals(params.copy(), 1, [None], [None])
        assert np.array_equal(capped, params)

    def test_a_zero_total_is_not_divided_by(self):
        params = _params([(0.0, 0.0, 0.005, 0.019)])
        capped, _ = cap_linewidth_totals(params, 1, [0.24], [0.036])
        assert np.isfinite(capped).all()


class TestTheCapReachesTheFitter:

    @staticmethod
    def _region(centre_f1=120.0, centre_f2=8.0, lw_f1=0.15, lw_f2=0.02, n=24):
        f1 = np.linspace(centre_f1 - 0.30, centre_f1 + 0.30, n)
        f2 = np.linspace(centre_f2 - 0.045, centre_f2 + 0.045, n)
        F2, F1 = np.meshgrid(f2, f1)
        data = np.exp(-((F1 - centre_f1) / lw_f1) ** 2
                      - ((F2 - centre_f2) / lw_f2) ** 2)
        return F1, F2, data

    @staticmethod
    def _peak():
        return [{'pos_f1': 120.0, 'lw_lor_f1': 0.03, 'lw_gau_f1': 0.13,
                 'pos_f2': 8.0, 'lw_lor_f2': 0.005, 'lw_gau_f2': 0.019,
                 'intensity': 1.0}]

    def test_lw_margins_bound_the_total_not_the_component(self):
        """The whole point: a cap of 0.10 must mean a TOTAL of 0.10, not 0.20."""
        F1, F2, data = self._region(lw_f1=0.28)          # broader than the cap allows
        fitter = Ps2dMultiPeakFitter2D(verbose=False)
        result = fitter.fit_multi_peak_2d(F1, F2, data, self._peak(),
                                          lw_margins=[(0.10, 0.030)])
        peak = result['peaks'][0]
        assert peak['lw_lor_f1'] + peak['lw_gau_f1'] <= 0.10 + 1e-9
        assert peak['lw_lor_f2'] + peak['lw_gau_f2'] <= 0.030 + 1e-9

    def test_without_lw_margins_nothing_changes(self):
        """An argument nobody passes must not alter a single existing fit."""
        F1, F2, data = self._region()
        fitter = Ps2dMultiPeakFitter2D(verbose=False)
        withNone = fitter.fit_multi_peak_2d(F1, F2, data, self._peak())
        withOut = fitter.fit_multi_peak_2d(F1, F2, data, self._peak(),
                                            lw_margins=None)
        a, b = withNone['peaks'][0], withOut['peaks'][0]
        for key in ('lw_lor_f1', 'lw_gau_f1', 'lw_lor_f2', 'lw_gau_f2', 'volume'):
            assert a[key] == pytest.approx(b[key])

    def test_a_cap_wider_than_the_peak_does_not_move_it(self):
        """A cap only ever tightens; a generous one must leave the fit where it was."""
        F1, F2, data = self._region()
        fitter = Ps2dMultiPeakFitter2D(verbose=False)
        free = fitter.fit_multi_peak_2d(F1, F2, data, self._peak())
        capped = fitter.fit_multi_peak_2d(F1, F2, data, self._peak(),
                                          lw_margins=[(5.0, 5.0)])
        a, b = free['peaks'][0], capped['peaks'][0]
        assert a['volume'] == pytest.approx(b['volume'], rel=1e-6)

    def test_a_peak_with_no_previous_width_keeps_the_learned_bound(self):
        """A width cap is a statement about ONE peak, unlike a position territory,
        which is a guarantee about the peaks as a set. So capping the members we can
        inform is strictly better than capping none, and a None entry is allowed."""
        F1, F2, data = self._region(lw_f1=0.28)
        fitter = Ps2dMultiPeakFitter2D(verbose=False)
        result = fitter.fit_multi_peak_2d(F1, F2, data, self._peak(),
                                          lw_margins=[None])
        assert result['peaks'], 'a None margin was not tolerated'


class TestCappingReDerivesWhatDependsOnTheWidth:
    """Narrowing a peak after the last optimisation and leaving `intensity` alone
    reports a model nothing fitted - intensity IS the volume, so the volume comes back
    unchanged while the height rises to keep the same area in a narrower profile.

    Measured on DNAJA1/HSPA8 73LysH when this was wrong: F1 total 0.2670 -> 0.2438,
    volume bit-identical at 1.539e9, height +32%, residual 3.0% -> 11.6%.
    """

    @staticmethod
    def _fitted(lw_margins):
        n = 28
        f1 = np.linspace(119.7, 120.3, n)
        f2 = np.linspace(7.955, 8.045, n)
        F2, F1 = np.meshgrid(f2, f1)
        # Broader than any cap below, so the cap has to bite.
        data = np.exp(-((F1 - 120.0) / 0.30) ** 2 - ((F2 - 8.0) / 0.040) ** 2)
        peak = [{'pos_f1': 120.0, 'lw_lor_f1': 0.03, 'lw_gau_f1': 0.20,
                 'pos_f2': 8.0, 'lw_lor_f2': 0.005, 'lw_gau_f2': 0.030,
                 'intensity': 1.0}]
        fitter = Ps2dMultiPeakFitter2D(verbose=False)
        return fitter.fit_multi_peak_2d(F1, F2, data, peak,
                                        lw_margins=lw_margins)['peaks'][0]

    def test_a_capped_fit_reports_a_smaller_volume_not_the_same_one(self):
        free = self._fitted(None)
        tight = self._fitted([(0.12, 0.018)])
        assert tight['lw_lor_f1'] + tight['lw_gau_f1'] < free['lw_lor_f1'] + free['lw_gau_f1']
        assert tight['volume'] < free['volume'], (
            'the volume did not follow the width down - the intensity was frozen')

    def test_the_height_does_not_rise_by_the_whole_area_ratio(self):
        """The exact signature of a frozen intensity, and the one thing that separates
        it from an honest constrained fit.

        `volume` IS the intensity parameter, so with it frozen the height rises by
        EXACTLY the ratio of the areas - the same volume packed into a narrower profile.
        A refitted intensity may still overshoot, because a model held far narrower than
        the data compromises between the apex and the wings, but it cannot track the
        area ratio.
        """
        free = self._fitted(None)
        tight = self._fitted([(0.12, 0.018)])

        areaRatio = (((free['lw_lor_f1'] + free['lw_gau_f1'])
                      * (free['lw_lor_f2'] + free['lw_gau_f2']))
                     / ((tight['lw_lor_f1'] + tight['lw_gau_f1'])
                        * (tight['lw_lor_f2'] + tight['lw_gau_f2'])))
        heightRatio = tight['height'] / free['height']
        assert areaRatio > 4, 'the cap did not bite hard enough to tell them apart'
        assert heightRatio < areaRatio / 2, (
            f'height rose {heightRatio:.1f}x against an area ratio of {areaRatio:.1f}x '
            f'- the intensity was not re-derived')


class TestHoldingLinewidthsKeepsTheIntensityFit:
    """fix_linewidths used to skip Stage 1 outright. That stage fits linewidths AND
    INTENSITIES, so gating all of it on one of the two threw the intensity fit away:
    a weak peak carried its initial intensity GUESS into Stages 2-4 and failed.

    Measured on a real T1 series - the 2000 ms point, the weakest, dropped from 3
    fitted peaks to 2 the moment linewidths were held.
    """

    # exp(-(x/w)^2) has FWHM = 2*w*sqrt(ln 2). The widths held below are the peak's
    # TRUE ones, or the held fit is being asked to describe data of another shape and
    # the test measures the fixture rather than the code.
    FWHM_F1, FWHM_F2 = 0.16, 0.024

    @classmethod
    def _weakPeak(cls, scale):
        n = 26
        f1 = np.linspace(119.7, 120.3, n)
        f2 = np.linspace(7.955, 8.045, n)
        F2, F1 = np.meshgrid(f2, f1)
        w1 = cls.FWHM_F1 / (2.0 * np.sqrt(np.log(2.0)))
        w2 = cls.FWHM_F2 / (2.0 * np.sqrt(np.log(2.0)))
        data = scale * np.exp(-((F1 - 120.0) / w1) ** 2 - ((F2 - 8.0) / w2) ** 2)
        # An initial intensity 1000x the truth, as a decayed peak's seed is: if Stage 1
        # is skipped it stays there and Stages 2-4 start from nonsense.
        peak = [{'pos_f1': 120.0, 'lw_lor_f1': 0.03, 'lw_gau_f1': 0.13,
                 'pos_f2': 8.0, 'lw_lor_f2': 0.005, 'lw_gau_f2': 0.019,
                 'intensity': 1.0}]
        return F1, F2, data, peak

    def test_a_weak_peak_still_fits_with_its_linewidths_held(self):
        F1, F2, data, peak = self._weakPeak(0.02)
        fitter = Ps2dMultiPeakFitter2D(verbose=False)
        result = fitter.fit_multi_peak_2d(F1, F2, data, peak, fix_linewidths=True)
        assert result['success'], 'holding the widths cost the fit'
        assert result['peaks'], 'no peak came back'

    def test_the_held_widths_are_actually_held(self):
        F1, F2, data, peak = self._weakPeak(1.0)
        fitter = Ps2dMultiPeakFitter2D(verbose=False)
        result = fitter.fit_multi_peak_2d(F1, F2, data, peak, fix_linewidths=True)
        got = result['peaks'][0]
        for key in ('lw_lor_f1', 'lw_gau_f1', 'lw_lor_f2', 'lw_gau_f2'):
            assert got[key] == pytest.approx(peak[0][key], rel=1e-6), (
                f'{key} moved despite fix_linewidths=True')

    def test_the_intensity_is_still_fitted(self):
        """The seed intensity is 1000x the truth. If Stage 1 is skipped it stays
        there - so the test is whether the fitted peak resembles the DATA or the seed."""
        F1, F2, data, peak = self._weakPeak(0.02)
        fitter = Ps2dMultiPeakFitter2D(verbose=False)
        held = fitter.fit_multi_peak_2d(F1, F2, data, peak,
                                        fix_linewidths=True)['peaks'][0]
        # NOT against data.max(): fit_multi_peak_2d is called directly here, so its
        # intensities never go through the denormalisation fit_overlap_group_2d
        # applies. The invariant is that the intensity moved OFF its seed at all.
        assert held['volume'] < 0.01 * peak[0]['intensity'], (
            f"volume {held['volume']:.4g} against a seed intensity of "
            f"{peak[0]['intensity']} - Stage 1 never fitted it")
