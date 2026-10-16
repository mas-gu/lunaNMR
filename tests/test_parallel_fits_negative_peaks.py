# ABOUTME: Parallel mode must fit a significant negative peak, as sequential mode does.
# ABOUTME: Skipping it wrote height 0.0 with no bound flags, which the analyses kept as data.

"""
Measured on KRASB 600 hetNOE: residue 33 is -47.7 sigma in the saturated plane. The
parallel pre-check skipped it as 'negative_intensity' and reported height 0.0 with blank
flags; blank reads as a legacy row, so the matrix carried 0.0 and `dynamixs hetnoe`
published hetNOE 0.0 +/- 0.0. Sequential mode fitted the same peak, reported it pinned on
its intensity floor, and the point was blanked. The two modes must agree.
"""

from lunaNMR.core.parallel_voigt_processor import _should_skip_peak

SAT_SIGMA = 7.58e4          # KRASB 600 saturated plane, MAD sigma
RESIDUE_33_SAT = -3614098.0  # signed sample at residue 33, -47.7 sigma


def test_a_significant_negative_peak_is_fitted():
    skip, _ = _should_skip_peak(RESIDUE_33_SAT, SAT_SIGMA, 3.0)
    assert not skip


def test_a_negative_peak_inside_the_noise_is_still_below_threshold():
    skip, reason = _should_skip_peak(-1.0 * SAT_SIGMA, SAT_SIGMA, 3.0)
    assert skip and reason == 'below_snr_threshold'


def test_a_positive_peak_is_judged_as_before():
    assert _should_skip_peak(10.0 * SAT_SIGMA, SAT_SIGMA, 3.0) == (False, None)
    assert _should_skip_peak(1.0 * SAT_SIGMA, SAT_SIGMA, 3.0) == (True, 'below_snr_threshold')


def test_zero_with_no_noise_estimate_is_not_skipped_as_negative():
    """With no noise level there is no threshold to judge by, so the fit decides."""
    assert _should_skip_peak(0.0, None, 3.0) == (False, None)
