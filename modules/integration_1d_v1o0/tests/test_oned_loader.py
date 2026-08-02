# ABOUTME: Tests for 1D spectrum loading and robust noise estimation.
# ABOUTME: Synthetic arrays cover the contract; the real .ft1 covers the nmrglue read path.

import sys
from pathlib import Path

import numpy as np
import pytest

_MODULE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_MODULE_DIR))

# data_example/ is gitignored, so the real-data tests skip on a fresh clone.
REAL_SPECTRUM = _MODULE_DIR.parent.parent / "data_example" / "1D" / "1D_KB_GTP_051.ft1"
needs_real_data = pytest.mark.skipif(not REAL_SPECTRUM.exists(),
                                     reason="real 1D spectrum not present")


class TestFromArrays:
    def test_round_trips_data_and_axis(self):
        from oned_loader import from_arrays
        ppm = np.linspace(10.0, 0.0, 101)
        data = np.arange(101, dtype=float)
        s = from_arrays(data, ppm)
        assert np.array_equal(s.data, data)
        assert np.array_equal(s.ppm_axis, ppm)

    def test_orients_axis_descending(self):
        """Spectra are conventionally plotted high-ppm-first; the loader guarantees it
        so downstream code never has to test the direction."""
        from oned_loader import from_arrays
        ppm = np.linspace(0.0, 10.0, 51)          # ascending input
        data = np.arange(51, dtype=float)
        s = from_arrays(data, ppm)
        assert s.ppm_axis[0] > s.ppm_axis[-1]
        assert s.data[0] == 50.0                   # data reversed with the axis

    def test_estimates_noise_when_not_given(self):
        from oned_loader import from_arrays
        rng = np.random.default_rng(0)
        s = from_arrays(rng.normal(0.0, 3.0, 8192), np.linspace(10.0, 0.0, 8192))
        assert s.noise == pytest.approx(3.0, rel=0.15)

    def test_mismatched_lengths_raise(self):
        from oned_loader import from_arrays
        with pytest.raises(ValueError):
            from_arrays(np.zeros(10), np.zeros(11))

    def test_exposes_ppm_step(self):
        from oned_loader import from_arrays
        s = from_arrays(np.zeros(101), np.linspace(10.0, 0.0, 101))
        assert s.ppm_step == pytest.approx(0.1, rel=1e-6)


class TestNoiseEstimation:
    def test_recovers_gaussian_sigma(self):
        from oned_loader import estimate_noise
        rng = np.random.default_rng(1)
        assert estimate_noise(rng.normal(0.0, 5.0, 16384)) == pytest.approx(5.0, rel=0.1)

    def test_is_not_inflated_by_strong_signal(self):
        """A global MAD over a spectrum containing peaks reports the signal, not the
        noise. On the real spectrum that error is 7.6x; this is the guard against it."""
        from oned_loader import estimate_noise
        rng = np.random.default_rng(2)
        y = rng.normal(0.0, 1.0, 32768)
        centre = np.arange(32768)
        for pos in (4000, 12000, 20000):                       # tall, wide peaks
            y = y + 500.0 * np.exp(-0.5 * ((centre - pos) / 300.0) ** 2)
        assert estimate_noise(y) == pytest.approx(1.0, rel=0.2)

    def test_is_not_inflated_by_a_negative_artefact(self):
        """Residual water suppression leaves an inverted spike; it must not count."""
        from oned_loader import estimate_noise
        rng = np.random.default_rng(3)
        y = rng.normal(0.0, 1.0, 32768)
        y[16000:16200] -= 800.0
        assert estimate_noise(y) == pytest.approx(1.0, rel=0.2)

    def test_flat_data_gives_zero(self):
        from oned_loader import estimate_noise
        assert estimate_noise(np.ones(1000)) == 0.0

    def test_short_array_does_not_raise(self):
        from oned_loader import estimate_noise
        assert np.isfinite(estimate_noise(np.array([1.0, 2.0, 3.0])))


@needs_real_data
class TestRealSpectrum:
    def test_loads_nmrpipe_1d(self):
        from oned_loader import load_spectrum
        s = load_spectrum(REAL_SPECTRUM)
        assert s.data.ndim == 1
        assert len(s.data) == 120536
        assert len(s.ppm_axis) == len(s.data)

    def test_ppm_axis_matches_acquisition(self):
        from oned_loader import load_spectrum
        s = load_spectrum(REAL_SPECTRUM)
        assert s.ppm_axis[0] == pytest.approx(12.5188, abs=1e-3)
        assert s.ppm_axis[-1] == pytest.approx(-3.1027, abs=1e-3)
        assert s.ppm_axis[0] > s.ppm_axis[-1]

    def test_noise_matches_a_hand_measured_quiet_region(self):
        """Quiet region 11.0-12.4 ppm measures 2.506e9 by MAD; a global MAD gives
        1.9e10, so this pins the estimator to the real value, not the inflated one."""
        from oned_loader import load_spectrum
        s = load_spectrum(REAL_SPECTRUM)
        assert s.noise == pytest.approx(2.506e9, rel=0.15)

    def test_metadata_carries_spectrometer_frequency(self):
        from oned_loader import load_spectrum
        s = load_spectrum(REAL_SPECTRUM)
        assert s.metadata['obs_mhz'] == pytest.approx(800.174, abs=0.01)
        assert s.metadata['nucleus'] == 'HN'

    def test_ppm_step_is_sub_hertz(self):
        from oned_loader import load_spectrum
        s = load_spectrum(REAL_SPECTRUM)
        assert s.ppm_step * s.metadata['obs_mhz'] == pytest.approx(0.1037, rel=1e-2)

    def test_intensity_scale_defaults_to_one_without_acquisition_metadata(self):
        """NMRPipe headers carry no NS/RG, so series normalisation must fall back to
        1.0 rather than silently inventing a scale."""
        from oned_loader import load_spectrum
        assert load_spectrum(REAL_SPECTRUM).intensity_scale == 1.0


class TestUnsupportedInput:
    def test_missing_file_raises(self):
        from oned_loader import load_spectrum
        with pytest.raises(FileNotFoundError):
            load_spectrum("/nonexistent/spectrum.ft1")

    def test_unreadable_format_raises(self, tmp_path):
        from oned_loader import load_spectrum
        bogus = tmp_path / "not_a_spectrum.ft1"
        bogus.write_bytes(b"this is not NMR data")
        with pytest.raises(ValueError):
            load_spectrum(bogus)
