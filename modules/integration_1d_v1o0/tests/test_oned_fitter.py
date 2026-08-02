# ABOUTME: Tests for the 1D peak fitter - Voigt fitting and numerical region integration.
# ABOUTME: Recovers known synthetic peaks and checks the seed, bounds and failure reporting.

import sys
from pathlib import Path

import numpy as np
import pytest

_MODULE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_MODULE_DIR))


def _synthetic_peak(area=1.0, center=5.0, sigma=0.004, gamma=0.003,
                    baseline=0.0, span=0.25, n=2001, noise=0.0, seed=0):
    """A single Voigt sampled on a descending ppm axis, as a real spectrum would be."""
    from oned_voigt import voigt_1d
    ppm = np.linspace(center + span, center - span, n)
    y = voigt_1d(ppm, area, center, sigma, gamma, baseline)
    if noise:
        y = y + np.random.default_rng(seed).normal(0.0, noise, size=y.shape)
    return ppm, y


class TestVoigtFitRecovery:
    def test_recovers_area_center_and_width(self):
        from oned_fitter import fit_peak_voigt
        ppm, y = _synthetic_peak(area=2.5, center=5.0, sigma=0.004, gamma=0.003)
        r = fit_peak_voigt(ppm, y, target=5.0)
        assert r['success']
        assert r['area'] == pytest.approx(2.5, rel=1e-3)
        assert r['center'] == pytest.approx(5.0, abs=1e-4)
        assert r['r_squared'] > 0.999

    def test_recovers_peak_offset_from_target(self):
        """The reference position is a starting point, not a constraint."""
        from oned_fitter import fit_peak_voigt
        ppm, y = _synthetic_peak(area=1.0, center=5.03)
        r = fit_peak_voigt(ppm, y, target=5.0, window=0.10)
        assert r['success']
        assert r['center'] == pytest.approx(5.03, abs=1e-3)

    def test_reported_height_matches_the_data(self):
        from oned_fitter import fit_peak_voigt
        ppm, y = _synthetic_peak(area=2.0, center=5.0, baseline=0.0)
        r = fit_peak_voigt(ppm, y, target=5.0)
        assert r['height'] == pytest.approx(y.max(), rel=5e-3)

    def test_recovers_area_on_a_raised_baseline(self):
        from oned_fitter import fit_peak_voigt
        ppm, y = _synthetic_peak(area=1.5, center=5.0, baseline=12.0)
        r = fit_peak_voigt(ppm, y, target=5.0)
        assert r['area'] == pytest.approx(1.5, rel=5e-3)
        assert r['baseline'] == pytest.approx(12.0, rel=1e-2)

    def test_survives_realistic_noise(self):
        from oned_fitter import fit_peak_voigt
        ppm, y = _synthetic_peak(area=1.0, center=5.0, noise=2.0, seed=3)
        r = fit_peak_voigt(ppm, y, target=5.0)
        assert r['success']
        assert r['area'] == pytest.approx(1.0, rel=0.05)

    def test_fwhm_is_consistent_with_fitted_widths(self):
        from oned_fitter import fit_peak_voigt
        from oned_voigt import voigt_fwhm
        ppm, y = _synthetic_peak(area=1.0, center=5.0, sigma=0.004, gamma=0.003)
        r = fit_peak_voigt(ppm, y, target=5.0)
        assert r['fwhm'] == pytest.approx(voigt_fwhm(r['sigma'], r['gamma']), rel=1e-6)

    def test_reports_parameter_errors(self):
        from oned_fitter import fit_peak_voigt
        ppm, y = _synthetic_peak(area=1.0, center=5.0, noise=1.0, seed=7)
        r = fit_peak_voigt(ppm, y, target=5.0)
        assert r['param_errors'] is not None
        assert r['param_errors']['area'] > 0.0

    def test_method_is_labelled(self):
        from oned_fitter import fit_peak_voigt
        ppm, y = _synthetic_peak()
        assert fit_peak_voigt(ppm, y, target=5.0)['method'] == 'voigt'


class TestSeedQuality:
    def test_seed_area_is_close_to_truth(self):
        """The seed must land near the true area, not ~2x low.

        Guards against reintroducing the Gaussian shortcut area = H*sigma*sqrt(2pi),
        which ignores the Voigt height factor.
        """
        from oned_fitter import initial_guess
        ppm, y = _synthetic_peak(area=1.0, center=5.0, sigma=0.004, gamma=0.003)
        guess = initial_guess(ppm, y, target=5.0)
        assert 0.5 < guess['area'] / 1.0 < 2.0

    def test_seed_centre_tracks_the_maximum(self):
        from oned_fitter import initial_guess
        ppm, y = _synthetic_peak(area=1.0, center=5.02)
        assert initial_guess(ppm, y, target=5.02)['center'] == pytest.approx(5.02, abs=5e-3)

    def test_seed_baseline_tracks_an_offset(self):
        from oned_fitter import initial_guess
        ppm, y = _synthetic_peak(area=1.0, center=5.0, baseline=7.0)
        assert initial_guess(ppm, y, target=5.0)['baseline'] == pytest.approx(7.0, rel=0.2)


class TestFitFailure:
    def test_flat_data_does_not_raise(self):
        from oned_fitter import fit_peak_voigt
        ppm = np.linspace(5.25, 4.75, 501)
        r = fit_peak_voigt(ppm, np.zeros_like(ppm), target=5.0)
        assert r['success'] is False
        assert r['area'] is None

    def test_too_few_points_fails_cleanly(self):
        from oned_fitter import fit_peak_voigt
        ppm = np.linspace(5.01, 4.99, 3)
        r = fit_peak_voigt(ppm, np.array([1.0, 2.0, 1.0]), target=5.0)
        assert r['success'] is False


class TestRegionIntegration:
    def test_recovers_area_of_a_known_peak(self):
        """Region sum over a wide window approaches the true area."""
        from oned_fitter import integrate_region
        ppm, y = _synthetic_peak(area=3.0, center=5.0, sigma=0.004, gamma=0.003, span=0.5)
        r = integrate_region(ppm, y, target=5.0, window=0.5)
        assert r['success']
        assert r['area'] == pytest.approx(3.0, rel=0.05)

    def test_subtracts_a_raised_baseline(self):
        from oned_fitter import integrate_region
        ppm, y = _synthetic_peak(area=2.0, center=5.0, baseline=30.0, span=0.5)
        r = integrate_region(ppm, y, target=5.0, window=0.5)
        assert r['area'] == pytest.approx(2.0, rel=0.10)

    def test_integrates_a_multiplet_a_single_voigt_would_miss(self):
        """The reason region sum exists: a J-coupled doublet is not one Voigt."""
        from oned_fitter import integrate_region, fit_peak_voigt
        from oned_voigt import voigt_1d
        ppm = np.linspace(5.4, 4.6, 4001)
        j = 0.012                                   # ~7 Hz at 600 MHz
        y = (voigt_1d(ppm, 1.0, 5.0 + j / 2, 0.002, 0.001, 0.0)
             + voigt_1d(ppm, 1.0, 5.0 - j / 2, 0.002, 0.001, 0.0))
        region = integrate_region(ppm, y, target=5.0, window=0.4)
        assert region['area'] == pytest.approx(2.0, rel=0.05)
        # the single-Voigt fit is the inferior model here - that is the point
        assert fit_peak_voigt(ppm, y, target=5.0)['r_squared'] < 0.99

    def test_reports_no_model_parameters(self):
        from oned_fitter import integrate_region
        ppm, y = _synthetic_peak(span=0.5)
        r = integrate_region(ppm, y, target=5.0, window=0.5)
        assert r['method'] == 'region'
        for key in ('sigma', 'gamma', 'fwhm', 'r_squared', 'param_errors'):
            assert r[key] is None, f"{key} must be None - there is no model to report it from"

    def test_empty_window_fails_cleanly(self):
        from oned_fitter import integrate_region
        ppm, y = _synthetic_peak(span=0.5)
        r = integrate_region(ppm, y, target=99.0, window=0.01)
        assert r['success'] is False
