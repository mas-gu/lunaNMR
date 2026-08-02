# ABOUTME: Tests for full-spectrum 1D peak detection, reference matching and intensity extraction.
# ABOUTME: Covers the neighbour-stealing guard that otherwise makes two peaks silently report one.

import sys
from pathlib import Path

import numpy as np
import pytest

_MODULE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_MODULE_DIR))

REAL_SPECTRUM = _MODULE_DIR.parent.parent / "data_example" / "1D" / "1D_KB_GTP_xxx.ft1"
needs_real_data = pytest.mark.skipif(not REAL_SPECTRUM.exists(),
                                     reason="real 1D spectrum not present")

# The two resolved peaks in the reference spectrum, 0.0095 ppm apart.
PEAK_A, PEAK_B = 8.1847, 8.1752


def _spectrum_with_peaks(centres, areas=None, sigma=0.0012, noise=0.0,
                         lo=0.0, hi=10.0, n=40001, seed=0):
    from oned_loader import from_arrays
    from oned_voigt import voigt_1d
    ppm = np.linspace(hi, lo, n)
    areas = areas or [1.0] * len(centres)
    y = np.zeros_like(ppm)
    for c, a in zip(centres, areas):
        y = y + voigt_1d(ppm, a, c, sigma, sigma / 2.0, 0.0)
    if noise:
        y = y + np.random.default_rng(seed).normal(0.0, noise, size=y.shape)
    return from_arrays(y, ppm)


class TestDetection:
    def test_finds_all_well_separated_peaks(self):
        from oned_detector import detect_peaks
        centres = [8.5, 7.2, 4.1, 2.0]
        found = detect_peaks(_spectrum_with_peaks(centres))
        assert len(found) == len(centres)
        for c, f in zip(sorted(centres, reverse=True), found):
            assert f['ppm'] == pytest.approx(c, abs=1e-3)

    def test_rejects_noise(self):
        from oned_detector import detect_peaks
        from oned_loader import from_arrays
        rng = np.random.default_rng(5)
        noise_only = from_arrays(rng.normal(0.0, 1.0, 20000), np.linspace(10.0, 0.0, 20000))
        assert detect_peaks(noise_only, min_snr=8.0) == []

    def test_min_separation_is_expressed_in_ppm_not_points(self):
        """The old cross-section code set min-distance as npoints/50, which on a
        120k-point spectrum is 0.31 ppm and merges whole multiplets. Separation must
        be scale-free: the same spectrum sampled twice as finely finds the same peaks.
        """
        from oned_detector import detect_peaks
        centres = [6.0, 5.98]
        coarse = detect_peaks(_spectrum_with_peaks(centres, n=20001), min_separation_ppm=0.005)
        fine = detect_peaks(_spectrum_with_peaks(centres, n=80001), min_separation_ppm=0.005)
        assert len(coarse) == len(fine) == 2

    def test_min_separation_merges_closer_peaks(self):
        from oned_detector import detect_peaks
        found = detect_peaks(_spectrum_with_peaks([6.0, 5.98]), min_separation_ppm=0.05)
        assert len(found) == 1

    def test_snr_threshold_filters_weak_peaks(self):
        from oned_detector import detect_peaks
        spec = _spectrum_with_peaks([8.0, 5.0], areas=[1.0, 0.0005], noise=0.02, seed=11)
        strong = detect_peaks(spec, min_snr=20.0)
        assert len(strong) == 1
        assert strong[0]['ppm'] == pytest.approx(8.0, abs=1e-2)

    def test_reports_snr_and_height(self):
        from oned_detector import detect_peaks
        found = detect_peaks(_spectrum_with_peaks([5.0], noise=0.01, seed=2))
        assert found[0]['height'] > 0
        assert found[0]['snr'] > 10


class TestReferenceMatching:
    def test_snaps_target_to_the_nearby_maximum(self):
        from oned_detector import match_reference_peaks
        spec = _spectrum_with_peaks([5.03])
        matched = match_reference_peaks(spec, [{'assignment': 'A', 'position': 5.0}],
                                        window=0.05)
        assert matched[0]['matched'] is True
        assert matched[0]['ppm'] == pytest.approx(5.03, abs=1e-3)

    def test_reports_unmatched_when_nothing_is_in_range(self):
        from oned_detector import match_reference_peaks
        spec = _spectrum_with_peaks([5.0])
        matched = match_reference_peaks(spec, [{'assignment': 'X', 'position': 2.0}],
                                        window=0.01)
        assert matched[0]['matched'] is False
        assert matched[0]['ppm'] is None

    def test_two_targets_never_resolve_to_the_same_maximum(self):
        """The silent failure: with a window wider than the peak separation both
        targets snap to the taller peak and the intensity ratio reads 1.0.
        """
        from oned_detector import match_reference_peaks
        spec = _spectrum_with_peaks([PEAK_A, PEAK_B], areas=[1.3, 1.0])
        targets = [{'assignment': 'A', 'position': PEAK_A},
                   {'assignment': 'B', 'position': PEAK_B}]
        matched = match_reference_peaks(spec, targets, window=0.05)
        positions = [m['ppm'] for m in matched if m['matched']]
        assert len(set(positions)) == len(positions), "two targets share one maximum"

    def test_flags_the_target_that_lost_the_contest(self):
        from oned_detector import match_reference_peaks
        spec = _spectrum_with_peaks([PEAK_A, PEAK_B], areas=[1.3, 1.0])
        targets = [{'assignment': 'A', 'position': PEAK_A},
                   {'assignment': 'B', 'position': PEAK_B}]
        matched = match_reference_peaks(spec, targets, window=0.05)
        by_name = {m['assignment']: m for m in matched}
        assert by_name['A']['ppm'] == pytest.approx(PEAK_A, abs=2e-3)
        assert by_name['B']['ppm'] == pytest.approx(PEAK_B, abs=2e-3)

    def test_narrow_window_matches_each_peak_independently(self):
        from oned_detector import match_reference_peaks
        spec = _spectrum_with_peaks([PEAK_A, PEAK_B], areas=[1.3, 1.0])
        targets = [{'assignment': 'A', 'position': PEAK_A},
                   {'assignment': 'B', 'position': PEAK_B}]
        matched = match_reference_peaks(spec, targets, window=0.003)
        assert all(m['matched'] for m in matched)


class TestWindowEstimation:
    def test_window_scales_with_linewidth(self):
        from oned_detector import estimate_window
        narrow = _spectrum_with_peaks([5.0], sigma=0.0005)
        broad = _spectrum_with_peaks([5.0], sigma=0.005)
        assert estimate_window(broad, 5.0) > estimate_window(narrow, 5.0)

    def test_window_is_capped_by_the_nearest_neighbour(self):
        """A fixed 0.05 ppm window swallows a neighbour 0.0095 ppm away and the fit
        collapses (R2 0.60 on real data). The window must respect the neighbour."""
        from oned_detector import estimate_window
        spec = _spectrum_with_peaks([PEAK_A, PEAK_B])
        w = estimate_window(spec, PEAK_A, neighbours=[PEAK_B])
        assert w < abs(PEAK_A - PEAK_B) / 2.0


class TestIntensityExtraction:
    def test_recovers_known_height(self):
        from oned_fitter import measure_intensity
        from oned_voigt import voigt_height
        spec = _spectrum_with_peaks([5.0], areas=[2.0], sigma=0.0012)
        expected = voigt_height(2.0, 0.0012, 0.0006)
        r = measure_intensity(spec.ppm_axis, spec.data, target=5.0, window=0.004)
        assert r['height'] == pytest.approx(expected, rel=0.02)
        assert r['method'] == 'intensity'

    def test_subtracts_a_baseline_offset(self):
        from oned_fitter import measure_intensity
        spec = _spectrum_with_peaks([5.0], areas=[1.0])
        r_flat = measure_intensity(spec.ppm_axis, spec.data, 5.0, window=0.004)
        r_raised = measure_intensity(spec.ppm_axis, spec.data + 500.0, 5.0, window=0.004)
        assert r_raised['height'] == pytest.approx(r_flat['height'], rel=0.02)

    def test_reports_no_model_parameters(self):
        from oned_fitter import measure_intensity
        spec = _spectrum_with_peaks([5.0])
        r = measure_intensity(spec.ppm_axis, spec.data, 5.0, window=0.004)
        for key in ('area', 'sigma', 'gamma', 'fwhm', 'r_squared', 'param_errors'):
            assert r[key] is None, f"{key} must be None for a bare intensity read"

    def test_empty_window_fails_cleanly(self):
        from oned_fitter import measure_intensity
        spec = _spectrum_with_peaks([5.0])
        assert measure_intensity(spec.ppm_axis, spec.data, 99.0, window=0.001)['success'] is False


@needs_real_data
class TestRealSpectrumIntensities:
    def _spectrum(self):
        from oned_loader import load_spectrum
        return load_spectrum(REAL_SPECTRUM)

    def test_finds_the_two_peaks_of_interest(self):
        from oned_detector import detect_peaks
        spec = self._spectrum()
        found = [f for f in detect_peaks(spec, min_snr=50.0) if 8.10 <= f['ppm'] <= 8.26]
        assert len(found) == 2
        assert {round(f['ppm'], 3) for f in found} == {round(PEAK_A, 3), round(PEAK_B, 3)}

    def test_intensity_ratio_is_stable_across_sane_windows(self):
        from oned_fitter import measure_intensity
        spec = self._spectrum()
        ratios = []
        for w in (0.003, 0.004):
            a = measure_intensity(spec.ppm_axis, spec.data, PEAK_A, window=w)['height']
            b = measure_intensity(spec.ppm_axis, spec.data, PEAK_B, window=w)['height']
            ratios.append(a / b)
        assert max(ratios) / min(ratios) < 1.05
        assert 1.2 < ratios[0] < 1.5

    def test_matching_does_not_collapse_the_ratio_at_a_wide_window(self):
        """Regression on the measured failure: at +/-0.02 ppm the naive read gave a
        ratio of 1.02 because both targets found the same maximum."""
        from oned_detector import match_reference_peaks
        spec = self._spectrum()
        targets = [{'assignment': 'A', 'position': PEAK_A},
                   {'assignment': 'B', 'position': PEAK_B}]
        matched = match_reference_peaks(spec, targets, window=0.02)
        heights = [m['height'] for m in matched if m['matched']]
        assert len(heights) == 2
        assert max(heights) / min(heights) > 1.2
