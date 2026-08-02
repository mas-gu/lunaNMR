# ABOUTME: Tests for the pure 1D Voigt math used by the 1D integration module.
# ABOUTME: Checks area normalisation, analytic height, FWHM, and the Gaussian/Lorentzian limits.

import sys
from pathlib import Path

import numpy as np
import pytest

_MODULE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_MODULE_DIR))

GAUSS_FWHM_FACTOR = 2.0 * np.sqrt(2.0 * np.log(2.0))   # 2.3548


def _numeric_fwhm(x, y):
    """FWHM measured off a dense sampled profile by linear interpolation."""
    baseline = min(y[0], y[-1])
    half = baseline + (y.max() - baseline) / 2.0
    above = np.where(y >= half)[0]
    lo, hi = above[0], above[-1]
    # interpolate the two half-max crossings for sub-grid accuracy
    left = np.interp(half, [y[lo - 1], y[lo]], [x[lo - 1], x[lo]])
    right = np.interp(half, [y[hi + 1], y[hi]], [x[hi + 1], x[hi]])
    return abs(right - left)


class TestVoigtProfile:
    # The integration windows below are enormous relative to the linewidths because a
    # Voigt has Lorentzian tails: the area missed beyond +/-L is ~(2/pi)*(gamma/L)*area,
    # which swamps a 1e-4 assertion at any sane spectral width. Verified: the shortfall
    # tracks that expression to 3 significant figures from L=60 out to L=20000.
    def test_integrates_to_area_parameter(self):
        """The `area` parameter IS the integrated intensity — the reported volume."""
        from oned_voigt import voigt_1d
        x = np.linspace(-2000.0, 2000.0, 2000001)
        y = voigt_1d(x, area=3.7, center=0.0, sigma=0.35, gamma=0.22, baseline=0.0)
        assert np.trapezoid(y, x) == pytest.approx(3.7, rel=1e-4)

    def test_area_normalisation_holds_off_centre(self):
        from oned_voigt import voigt_1d
        x = np.linspace(-2000.0, 2000.0, 2000001)
        y = voigt_1d(x, area=2.0, center=8.5, sigma=0.4, gamma=0.1, baseline=0.0)
        assert np.trapezoid(y, x) == pytest.approx(2.0, rel=1e-4)

    def test_truncated_area_matches_analytic_lorentzian_tail(self):
        """Guards the normalisation constant: a wrong prefactor would shift the
        integral by a constant factor rather than by the tail term."""
        from oned_voigt import voigt_1d
        area, gamma, L = 2.0, 0.22, 200.0
        x = np.linspace(-L, L, 2000001)
        shortfall = area - np.trapezoid(voigt_1d(x, area, 0.0, 0.35, gamma, 0.0), x)
        assert shortfall == pytest.approx((2.0 / np.pi) * (gamma / L) * area, rel=0.02)

    def test_baseline_is_an_additive_offset(self):
        from oned_voigt import voigt_1d
        x = np.linspace(-5.0, 5.0, 2001)
        a = voigt_1d(x, 1.0, 0.0, 0.3, 0.2, baseline=0.0)
        b = voigt_1d(x, 1.0, 0.0, 0.3, 0.2, baseline=1.75)
        assert np.allclose(b - a, 1.75)

    def test_peaks_at_the_centre(self):
        from oned_voigt import voigt_1d
        x = np.linspace(-5.0, 5.0, 20001)
        y = voigt_1d(x, 1.0, 1.234, 0.3, 0.2, 0.0)
        assert x[int(np.argmax(y))] == pytest.approx(1.234, abs=1e-3)

    def test_zero_sigma_is_guarded(self):
        """sigma=0 must not divide by zero — the guard the core copy lacks."""
        from oned_voigt import voigt_1d
        x = np.linspace(-1.0, 1.0, 501)
        y = voigt_1d(x, 1.0, 0.0, 0.0, 0.1, 0.0)
        assert np.all(np.isfinite(y))


class TestVoigtHeight:
    def test_matches_sampled_maximum(self):
        from oned_voigt import voigt_1d, voigt_height
        x = np.linspace(-20.0, 20.0, 800001)
        y = voigt_1d(x, area=2.5, center=0.0, sigma=0.3, gamma=0.15, baseline=0.0)
        assert voigt_height(2.5, 0.3, 0.15) == pytest.approx(y.max(), rel=1e-5)

    def test_scales_linearly_with_area(self):
        from oned_voigt import voigt_height
        assert voigt_height(6.0, 0.3, 0.15) == pytest.approx(3.0 * voigt_height(2.0, 0.3, 0.15))

    def test_gaussian_limit_height(self):
        """gamma -> 0 gives the pure Gaussian height area/(sigma*sqrt(2pi))."""
        from oned_voigt import voigt_height
        sigma = 0.4
        expected = 5.0 / (sigma * np.sqrt(2.0 * np.pi))
        assert voigt_height(5.0, sigma, 1e-12) == pytest.approx(expected, rel=1e-6)

    def test_is_below_gaussian_height_when_lorentzian_present(self):
        """Real Voigt height is strictly below the Gaussian relation area/(sigma*sqrt(2pi)).

        This is the factor the naive seed `area = H*sigma*sqrt(2pi)` ignores.
        """
        from oned_voigt import voigt_height
        sigma, gamma, area = 0.3, 0.25, 1.0
        gaussian_height = area / (sigma * np.sqrt(2.0 * np.pi))
        assert voigt_height(area, sigma, gamma) < 0.75 * gaussian_height


class TestVoigtFwhm:
    def test_gaussian_limit(self):
        from oned_voigt import voigt_fwhm
        sigma = 0.37
        assert voigt_fwhm(sigma, 1e-12) == pytest.approx(GAUSS_FWHM_FACTOR * sigma, rel=1e-3)

    def test_lorentzian_limit(self):
        from oned_voigt import voigt_fwhm
        gamma = 0.29
        assert voigt_fwhm(1e-12, gamma) == pytest.approx(2.0 * gamma, rel=1e-3)

    @pytest.mark.parametrize("sigma,gamma", [(0.30, 0.20), (0.10, 0.45), (0.50, 0.05)])
    def test_matches_numerically_measured_width(self, sigma, gamma):
        from oned_voigt import voigt_1d, voigt_fwhm
        x = np.linspace(-30.0, 30.0, 600001)
        y = voigt_1d(x, 1.0, 0.0, sigma, gamma, 0.0)
        assert voigt_fwhm(sigma, gamma) == pytest.approx(_numeric_fwhm(x, y), rel=5e-3)

    def test_monotonic_in_both_widths(self):
        from oned_voigt import voigt_fwhm
        base = voigt_fwhm(0.3, 0.2)
        assert voigt_fwhm(0.4, 0.2) > base
        assert voigt_fwhm(0.3, 0.3) > base


class TestRSquared:
    def test_perfect_fit_is_one(self):
        from oned_voigt import r_squared
        y = np.array([1.0, 4.0, 9.0, 2.0])
        assert r_squared(y, y.copy()) == pytest.approx(1.0)

    def test_mean_prediction_is_zero(self):
        from oned_voigt import r_squared
        y = np.array([1.0, 2.0, 3.0, 4.0])
        assert r_squared(y, np.full_like(y, y.mean())) == pytest.approx(0.0)

    def test_worse_than_mean_is_negative(self):
        from oned_voigt import r_squared
        y = np.array([1.0, 2.0, 3.0, 4.0])
        assert r_squared(y, np.array([10.0, -5.0, 20.0, -8.0])) < 0.0

    def test_flat_data_does_not_divide_by_zero(self):
        from oned_voigt import r_squared
        flat = np.full(5, 2.0)
        assert np.isfinite(r_squared(flat, flat.copy()))
