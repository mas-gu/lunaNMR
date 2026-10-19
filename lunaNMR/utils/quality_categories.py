#!/usr/bin/env python3
# ABOUTME: fit_quality_label, the one rule a fitted peak's quality label is given by,
# ABOUTME: and usable_points: which points of a series carry a measurement of an observable.
import math

# ---------------------------------------------------------------------------
# The quality label of a fit
# ---------------------------------------------------------------------------
# A fit this far above the noise is a measurement whatever its R-squared says.
MEASURED_SNR = 5.0


def fit_quality_label(r_squared, success=True, near_noise=False, snr=None):
    """Quality label of one fitted peak: Excellent, Good, Fair, Poor or Failed.

    The one rule every fitter, the series output and the overlay refit label a fit
    with. The Kd loaders and the intensity matrices drop a Failed point, so where this
    rule puts Failed decides what is measured.

    - A fit that did not succeed, or has no R-squared, is Failed.
    - R-squared above 0.9 is Excellent, above 0.8 Good, above 0.6 Fair.
    - Below that, a peak more than MEASURED_SNR above the noise is Poor: it stands
      clear of the noise, so it is a measurement, if not a well-shaped one.
    - Otherwise R-squared of 0.5 or more is Fair and 0.2 or more Poor.
    - Below 0.2, a near-noise peak is Poor, never Failed: its R-squared is ~0 by
      construction - a held shape against noise - and its near-zero intensity is the
      measurement. Anything else below 0.2 is Failed.

    So between 0.5 and 0.6 a peak above MEASURED_SNR is Poor and one below it Fair:
    the label is not monotonic in S/N there.

    :param r_squared: the fit's R-squared (a cluster member's is its cluster's)
    :param success: whether the fit succeeded
    :param near_noise: whether the fit held this peak's shape as near-noise
    :param snr: |height| / noise of the fitted peak, or None when not known
    """
    if not success or r_squared is None or not math.isfinite(r_squared):
        return 'Failed'
    if r_squared > 0.9:
        return 'Excellent'
    if r_squared > 0.8:
        return 'Good'
    if r_squared > 0.6:
        return 'Fair'
    if snr is not None and snr > MEASURED_SNR:
        return 'Poor'
    if r_squared >= 0.5:
        return 'Fair'
    if r_squared >= 0.2 or near_noise:
        return 'Poor'
    return 'Failed'


# ---------------------------------------------------------------------------
# Which points of a series carry a measurement
# ---------------------------------------------------------------------------
# Below this many sigma a peak is not detected. The same floor holdPositionsAtNoise
# uses; 3 sigma is the usual bar for calling a signal present at all.
DETECTION_SIGMA_FLOOR = 3.0

# Observables that need a POSITION. A peak that has decayed has none, so it cannot
# contribute a chemical-shift perturbation - but it still contributes an intensity.
POSITION_OBSERVABLES = ('csp', 'shift', 'position')


def usable_points(frame, observable):
    """The rows of a series that carry a measurement of `observable`.

    Two states look identical in a series and need opposite handling, which is the
    whole reason this exists.

    A residue whose peak has DECAYED reports a near-zero intensity, and that zero IS
    the measurement: dropping it from an I/I0 series invents a plateau that is not
    there. It has no position, though, so it cannot contribute a CSP.

    A residue whose fit stopped on its INTENSITY floor reports a number that came from
    the bound rather than from the data. The floor is cluster-relative (0.1% of the
    brightest member's seed), so that number is a product of its NEIGHBOURS' seeds and
    its own seeded linewidths, carrying nothing from the data at its own position. It
    is not a measurement of anything and is dropped from every observable.

    A residue whose POSITION stopped on its margin has no measured position, so it is
    dropped from the position observables. Its intensity is kept: on a 600 MHz hetNOE pair
    most such peaks sit on the 0.020 ppm 15N margin, a quarter of one 15N point, and
    widening it to 0.05 ppm moved the freed peaks' heights by a median 0.0%.

    Neither state is visible in `success` or in R-squared. `success` is an OR whose
    first arm is the optimiser's "chi2 stopped changing", true immediately on flat
    noise; R-squared is per CLUSTER, so a member with nothing at its position rides the
    one its healthy neighbours earn - measured, 0.58 sigma under a cluster R-squared of
    0.9884.

    A row from before these flags existed carries NaN. It is kept: refusing every
    legacy series would be worse than the risk it represents, and the flags are absent
    rather than negative. Re-run the series to get the filter.

    :param frame: a DataFrame of series rows, one per (spectrum, peak)
    :param observable: what is being measured; see POSITION_OBSERVABLES
    :return: the subset of `frame` that can be used, index preserved
    """
    keep = frame.index == frame.index                 # all True, right length
    positional = str(observable).lower() in POSITION_OBSERVABLES
    columns = ('position_pinned', 'intensity_pinned') if positional else ('intensity_pinned',)
    for column in columns:
        if column in frame:
            # `== True` and not truthiness: NaN must not read as pinned, and None
            # must not read as "not pinned" by accident either - both are unknown,
            # and unknown is kept.
            keep &= ~(frame[column].fillna(False).astype(bool))

    if positional and 'detection_sigma' in frame:
        detected = frame['detection_sigma'] >= DETECTION_SIGMA_FLOOR
        keep &= detected | frame['detection_sigma'].isna()

    return frame[keep]


__all__ = [
    'fit_quality_label',
    'usable_points',
    'DETECTION_SIGMA_FLOOR',
]