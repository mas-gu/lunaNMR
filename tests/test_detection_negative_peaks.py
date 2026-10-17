# ABOUTME: Peak detection finds negative peaks as it finds positive ones, carrying their sign.
# ABOUTME: Real DNAJA1 spectrum: detecting on the data and on its negation must find the same peaks.

"""
A negative peak is real data - a flexible residue in a saturated hetNOE plane, a CH2 in an
edited 13C-HSQC, a folded peak. Detection looked for maxima above a positive threshold, so
a negative peak was never detected: in a peak-list run it fell back to its reference
position, marked undetected; with no peak list it did not exist. Detection is on the
magnitude, so the sign of the spectrum cannot change which peaks are found.
"""

import contextlib
import io
from pathlib import Path

import numpy as np
import pytest

SPECTRUM = (Path(__file__).resolve().parents[2] / "data_test" / "data_binding"
            / "DNAJA1_HSPA8" / "2D_DNAJA1_0o0.ft")

pytestmark = pytest.mark.skipif(not SPECTRUM.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def integrator():
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    it = EnhancedVoigtIntegrator()
    with contextlib.redirect_stdout(io.StringIO()):
        it.load_nmr_file(str(SPECTRUM))
    it.noise_level = it._estimate_noise_level_advanced()
    return it


def _detect(it, sign):
    data = it.nmr_data
    it.nmr_data = sign * data
    try:
        return it._detect_peaks_by_threshold(10 * it.noise_level)
    finally:
        it.nmr_data = data


def _by_position(peaks):
    return {(round(p['x_ppm'], 6), round(p['y_ppm'], 6)): p for p in peaks}


def test_the_sign_of_the_spectrum_does_not_change_which_peaks_are_found(integrator):
    plus, minus = _by_position(_detect(integrator, 1)), _by_position(_detect(integrator, -1))
    assert set(plus) == set(minus)
    for pos, peak in plus.items():
        assert minus[pos]['intensity'] == -peak['intensity']
        assert minus[pos]['snr'] == peak['snr']


def test_a_negative_peak_is_detected_with_its_sign(integrator):
    negative = [p for p in _detect(integrator, 1) if p['intensity'] < 0]
    assert negative, "DNAJA1's reference spectrum has a negative peak above 10 sigma"
    for p in negative:
        assert p['snr'] >= 10


def test_a_negative_peak_in_the_peak_list_is_matched(integrator):
    import pandas as pd
    peak = next(p for p in _detect(integrator, 1) if p['intensity'] < 0)
    integrator.peak_list = pd.DataFrame({'Assignment': ['neg'], 'Position_X': [peak['x_ppm']],
                                         'Position_Y': [peak['y_ppm']]})
    integrator.search_window_x = integrator.search_window_y = 0.070
    integrator.peak_list_contour_threshold = None
    integrator.peak_list_sn_threshold = 10.0
    integrator.auto_add_dummy_peaks = False
    matched = integrator._detect_peaks_reference_based()[0]
    assert matched['detected']
    assert matched['intensity'] == peak['intensity'] < 0


def test_a_detected_peak_records_its_pixel_on_the_right_axes(integrator):
    import pandas as pd
    peak = _detect(integrator, 1)[0]
    integrator.peak_list = pd.DataFrame({'Assignment': ['p'], 'Position_X': [peak['x_ppm']],
                                         'Position_Y': [peak['y_ppm']]})
    integrator.auto_add_dummy_peaks = False
    matched = integrator._detect_peaks_reference_based()[0]
    assert integrator.nmr_data[matched['y_point'], matched['x_point']] == matched['intensity']


def test_a_count_cutoff_keeps_the_strongest_peaks_whatever_their_sign(integrator):
    peaks = _detect(integrator, 1)
    negative = next(p for p in peaks if p['intensity'] < 0)
    keep = sum(1 for p in peaks if abs(p['intensity']) >= abs(negative['intensity']))
    integrator.sn_from_gui_mode = False
    integrator.sn_threshold = 10.0
    integrator.expected_peak_count = keep
    found = integrator.detect_peaks_sn_native()
    assert any(p['intensity'] == negative['intensity'] for p in found)


def test_a_negative_peak_beside_a_listed_one_is_added_to_its_fit(integrator):
    """An unmatched detected peak close to a listed one is fitted with it; its sign must
    not decide whether it is."""
    peaks = _detect(integrator, 1)
    negative_idx = next(i for i, p in enumerate(peaks) if p['intensity'] < 0)
    negative = peaks[negative_idx]
    beside = {'assignment': 'listed', 'ppm_x': negative['x_ppm'] + 0.02,
              'ppm_y': negative['y_ppm'], 'intensity': abs(negative['intensity']),
              'detected': True}
    matched = [beside]
    used = set(range(len(peaks))) - {negative_idx}
    assert integrator._add_dummy_peaks(peaks, used, matched) == 1
    assert matched[1]['intensity'] == negative['intensity']


def test_diagnose_counts_a_negative_peak_as_captured(integrator):
    from lunaNMR.validation.spectra_check import _capture
    negative = next(p for p in _detect(integrator, 1) if p['intensity'] < 0)
    cap = _capture(str(SPECTRUM), [("neg", negative['x_ppm'], negative['y_ppm'])])
    assert cap['captured'] == 1
    assert cap['weak'] == 0
