# ABOUTME: Per-peak 1D measurement - Voigt fitting, region summation and bare intensity.
# ABOUTME: Depends only on numpy/scipy plus oned_voigt, so it runs without any GUI or lunaNMR import.

import numpy as np
from scipy.optimize import curve_fit

from oned_voigt import voigt_1d, voigt_height, voigt_fwhm, r_squared

# Half-width in ppm of the data window taken around a target position.
DEFAULT_WINDOW = 0.05

# Narrowest width the optimiser may return, in ppm. At 600 MHz this is 0.06 Hz -
# below any real resonance, so it only stops the optimiser collapsing a width to zero.
MIN_WIDTH = 1e-4

PARAM_NAMES = ('area', 'center', 'sigma', 'gamma', 'baseline')

# The baseline under a peak is read from a region this many times the measurement
# window, at this percentile - low enough that peaks, which are positive excursions,
# do not pull it up.
BASELINE_REGION_FACTOR = 5.0
BASELINE_PERCENTILE = 10.0

# Same thresholds as the rest of LunaNMR so the marker colours agree across viewers.
QUALITY_THRESHOLDS = ((0.9, 'Excellent'), (0.8, 'Good'), (0.5, 'Fair'))


def quality_from_r_squared(value):
    """Bucket a fit R² into the shared quality categories."""
    if value is None or not np.isfinite(value):
        return 'Unknown'
    for threshold, label in QUALITY_THRESHOLDS:
        if value >= threshold:
            return label
    return 'Poor'


def _empty_result(method, assignment, target):
    return {
        'assignment': assignment, 'position': target, 'method': method,
        'center': None, 'area': None, 'height': None, 'sigma': None,
        'gamma': None, 'fwhm': None, 'baseline': None, 'r_squared': None,
        'snr': None, 'quality': 'Unknown', 'success': False,
        'param_errors': None, 'window_bounds': None,
    }


def _select_window(ppm, y, target, window):
    """Points within +/- `window` ppm of `target`. The ppm axis is descending."""
    ppm = np.asarray(ppm, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.abs(ppm - target) <= window
    return ppm[mask], y[mask]


def _edge_baseline(ppm, y):
    """Baseline interpolated between the medians of the two window edges.

    A straight line through the edges rather than a flat offset, so a sloping
    baseline under the peak does not bias the integral.
    """
    edge = max(3, len(y) // 10)
    left_ppm, left_val = np.median(ppm[:edge]), np.median(y[:edge])
    right_ppm, right_val = np.median(ppm[-edge:]), np.median(y[-edge:])

    if left_ppm == right_ppm:
        return np.full_like(y, (left_val + right_val) / 2.0)

    slope = (right_val - left_val) / (right_ppm - left_ppm)
    return left_val + slope * (ppm - left_ppm)


def _local_baseline(ppm, y, target, window,
                    region_factor=BASELINE_REGION_FACTOR,
                    percentile=BASELINE_PERCENTILE):
    """Baseline under a peak, taken as a low percentile of a wider surrounding region.

    Taking it from the edges of the measurement window instead subtracts part of the
    peak's own tails: the height then comes out 13.7% low at a window of one FWHM and
    2.2% low at two, so the measurement depends on the window it was made in. A low
    percentile over a wider region is flat to within a few tenths of a percent at any
    window width, and a neighbouring peak does not disturb it because peaks are
    positive excursions that the percentile discards.
    """
    ppm = np.asarray(ppm, dtype=float)
    y = np.asarray(y, dtype=float)

    region = np.abs(ppm - target) <= window * region_factor
    if not region.any():
        return 0.0

    return float(np.percentile(y[region], percentile))


def _integrate(ppm, y):
    """Trapezoidal integral that is positive for a positive peak on a descending axis."""
    integral = float(np.trapezoid(y, ppm))
    return -integral if ppm[0] > ppm[-1] else integral


def _estimate_fwhm(ppm, y, baseline):
    """FWHM from the half-maximum crossings, falling back to a fraction of the window."""
    height = y.max() - baseline
    if height <= 0:
        return None

    above = np.where(y - baseline >= height / 2.0)[0]
    if len(above) < 2:
        return None

    width = abs(ppm[above[-1]] - ppm[above[0]])
    return width if width > 0 else None


def initial_guess(ppm, y, target, window=DEFAULT_WINDOW):
    """Physically consistent Voigt seed for the data around `target`.

    The area is obtained by inverting voigt_height() at the seeded widths rather
    than through the Gaussian relation area = H*sigma*sqrt(2pi); that shortcut
    omits the Voigt height factor and seeds the area roughly 2x low.
    """
    ppm_w, y_w = _select_window(ppm, y, target, window)
    if len(y_w) == 0:
        return None

    baseline = float(np.median(np.concatenate([y_w[:max(3, len(y_w) // 10)],
                                               y_w[-max(3, len(y_w) // 10):]])))
    height = float(y_w.max() - baseline)
    if height <= 0:
        return None

    center = float(ppm_w[int(np.argmax(y_w))])

    fwhm = _estimate_fwhm(ppm_w, y_w, baseline) or (window / 5.0)
    sigma = fwhm / (2.0 * 2.3548)
    gamma = fwhm / 4.0

    unit_height = voigt_height(1.0, sigma, gamma)
    area = height / unit_height if unit_height > 0 else height

    return {'area': area, 'center': center, 'sigma': sigma,
            'gamma': gamma, 'baseline': baseline}


def _bounds(guess, target, window):
    span = abs(guess['baseline']) + abs(guess['area'] / max(guess['sigma'], MIN_WIDTH))
    lower = [0.0, target - window, MIN_WIDTH, MIN_WIDTH, guess['baseline'] - span]
    upper = [np.inf, target + window, window, window, guess['baseline'] + span]
    return lower, upper


def fit_peak_voigt(ppm, y, target, window=DEFAULT_WINDOW, assignment=None):
    """Fit a single Voigt to the data around `target`.

    The returned `area` is the integrated intensity of the fitted profile - the
    full analytic area, not the area inside the window.
    """
    result = _empty_result('voigt', assignment, target)

    ppm_w, y_w = _select_window(ppm, y, target, window)
    if len(y_w) <= len(PARAM_NAMES):
        return result

    guess = initial_guess(ppm, y, target, window)
    if guess is None:
        return result

    p0 = [guess[name] for name in PARAM_NAMES]
    try:
        popt, pcov = curve_fit(voigt_1d, ppm_w, y_w, p0=p0,
                               bounds=_bounds(guess, target, window), maxfev=10000)
    except (RuntimeError, ValueError):
        return result

    fitted = voigt_1d(ppm_w, *popt)
    residual = y_w - fitted
    values = dict(zip(PARAM_NAMES, (float(v) for v in popt)))

    errors = np.sqrt(np.diag(pcov))
    noise = float(np.std(residual))

    result.update(values)
    result['height'] = voigt_height(values['area'], values['sigma'], values['gamma'])
    result['fwhm'] = voigt_fwhm(values['sigma'], values['gamma'])
    result['r_squared'] = r_squared(y_w, fitted)
    result['snr'] = result['height'] / noise if noise > 0 else None
    result['quality'] = quality_from_r_squared(result['r_squared'])
    result['window_bounds'] = (float(ppm_w.min()), float(ppm_w.max()))
    result['param_errors'] = (dict(zip(PARAM_NAMES, (float(e) for e in errors)))
                              if np.all(np.isfinite(errors)) else None)
    result['success'] = True

    return result


def measure_intensity(ppm, y, target, window=DEFAULT_WINDOW, assignment=None):
    """Baseline-corrected peak height in the window around `target`.

    The cheapest and most stable observable: no model, no integration limits, so
    none of the window sensitivity that area carries. Valid for comparing a peak
    across a series only while its linewidth stays constant - a broadening peak
    loses height at constant area.

    The window must stay narrower than half the distance to any neighbouring peak,
    or the maximum found belongs to the neighbour. See oned_detector.estimate_window.
    """
    result = _empty_result('intensity', assignment, target)

    ppm_w, y_w = _select_window(ppm, y, target, window)
    if len(y_w) < 2:
        return result

    baseline = _local_baseline(ppm, y, target, window)
    peak_idx = int(np.argmax(y_w))

    result['height'] = float(y_w[peak_idx]) - baseline
    result['center'] = float(ppm_w[peak_idx])
    result['baseline'] = baseline
    result['window_bounds'] = (float(ppm_w.min()), float(ppm_w.max()))
    result['success'] = True

    return result


def integrate_region(ppm, y, target, window=DEFAULT_WINDOW, assignment=None):
    """Baseline-corrected numerical integral over the window around `target`.

    Model-free, so it integrates multiplets and overlapping envelopes that a single
    Voigt cannot describe. Nothing model-derived (widths, R², parameter errors) is
    reported - those stay None rather than being invented.
    """
    result = _empty_result('region', assignment, target)

    ppm_w, y_w = _select_window(ppm, y, target, window)
    if len(y_w) < 2:
        return result

    baseline = _edge_baseline(ppm_w, y_w)
    corrected = y_w - baseline
    peak_idx = int(np.argmax(corrected))

    result['area'] = _integrate(ppm_w, corrected)
    result['height'] = float(corrected[peak_idx])
    result['center'] = float(ppm_w[peak_idx])
    result['baseline'] = float(np.mean(baseline))
    result['window_bounds'] = (float(ppm_w.min()), float(ppm_w.max()))
    result['success'] = True

    return result
