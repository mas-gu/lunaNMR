# ABOUTME: Full-spectrum 1D peak detection, reference-list matching and window sizing.
# ABOUTME: Thresholds are in ppm and multiples of the noise, so they do not depend on spectrum size.

import numpy as np
from scipy.signal import find_peaks

# Detection defaults. Separation and window are ppm so they are independent of how
# finely the spectrum was sampled - the old cross-section code derived min-distance
# from the point count, which on a 120k-point spectrum came to 0.31 ppm.
DEFAULT_MIN_SNR = 10.0
DEFAULT_MIN_SEPARATION_PPM = 0.002
DEFAULT_MATCH_WINDOW_PPM = 0.05

# A fit/measurement window this many times the peak FWHM, unless a neighbour is closer.
WINDOW_FWHM_MULTIPLE = 2.0

# Keep a window clear of the midpoint to a neighbouring peak by this factor.
NEIGHBOUR_SAFETY = 0.8


def _ppm_to_points(spectrum, ppm):
    step = spectrum.ppm_step
    return max(1, int(round(ppm / step))) if step > 0 else 1


def detect_peaks(spectrum, min_snr=DEFAULT_MIN_SNR,
                 min_separation_ppm=DEFAULT_MIN_SEPARATION_PPM):
    """Find maxima across the whole spectrum, ordered high ppm first.

    Thresholds are relative to the estimated noise, so the same settings work on
    spectra of different absolute intensity.
    """
    data = spectrum.data
    noise = spectrum.noise

    if noise <= 0 or data.size == 0:
        return []

    indices, _ = find_peaks(data,
                            height=min_snr * noise,
                            prominence=min_snr * noise / 2.0,
                            distance=_ppm_to_points(spectrum, min_separation_ppm))

    return [{'ppm': float(spectrum.ppm_axis[i]),
             'index': int(i),
             'height': float(data[i]),
             'snr': float(data[i] / noise)}
            for i in indices]


def estimate_fwhm(spectrum, target, search=0.05):
    """FWHM of the maximum nearest `target`, measured by half-height crossing."""
    ppm, data = spectrum.ppm_axis, spectrum.data

    mask = np.abs(ppm - target) <= search
    if not mask.any():
        return None

    offset = int(np.argmax(mask))
    segment = data[mask]
    peak = int(np.argmax(segment))
    half = segment[peak] / 2.0

    left = peak
    while left > 0 and segment[left] > half:
        left -= 1
    right = peak
    while right < len(segment) - 1 and segment[right] > half:
        right += 1

    width = abs(float(ppm[offset + right] - ppm[offset + left]))

    return width or None


def estimate_window(spectrum, target, neighbours=(), search=0.05):
    """Half-width in ppm of the data window to use around `target`.

    Driven by the peak's own linewidth rather than a fixed constant: one setting
    cannot serve a 2 Hz ligand singlet and a broad protein resonance. Capped short
    of any neighbouring peak, because a window that reaches the neighbour makes the
    fit collapse and lets a bare maximum read the wrong peak entirely.
    """
    fwhm = estimate_fwhm(spectrum, target, search=search)
    window = fwhm * WINDOW_FWHM_MULTIPLE if fwhm else search

    for neighbour in neighbours:
        separation = abs(float(neighbour) - float(target))
        if separation > 0:
            window = min(window, separation / 2.0 * NEIGHBOUR_SAFETY)

    return window


def match_reference_peaks(spectrum, targets, window=DEFAULT_MATCH_WINDOW_PPM,
                          min_snr=DEFAULT_MIN_SNR):
    """Snap each reference position to a nearby maximum, one maximum per target.

    Targets compete for maxima: the closest target wins and the others look again
    among what is left. Without this, a window wider than the peak separation makes
    every target settle on the same tallest peak, and the reported intensities are
    silently those of one peak - a ratio of 1.0 that looks entirely plausible.
    """
    detected = detect_peaks(spectrum, min_snr=min_snr,
                            min_separation_ppm=DEFAULT_MIN_SEPARATION_PPM)

    results = [{'assignment': t.get('assignment'), 'position': float(t['position']),
                'ppm': None, 'index': None, 'height': None, 'snr': None,
                'matched': False}
               for t in targets]

    # Every (target, candidate) pair inside the window, nearest pairing settled first.
    pairs = []
    for r_idx, result in enumerate(results):
        for peak in detected:
            distance = abs(peak['ppm'] - result['position'])
            if distance <= window:
                pairs.append((distance, r_idx, peak))
    pairs.sort(key=lambda p: p[0])

    claimed_targets = set()
    claimed_peaks = set()
    for _, r_idx, peak in pairs:
        if r_idx in claimed_targets or peak['index'] in claimed_peaks:
            continue
        claimed_targets.add(r_idx)
        claimed_peaks.add(peak['index'])
        results[r_idx].update(ppm=peak['ppm'], index=peak['index'],
                              height=peak['height'], snr=peak['snr'], matched=True)

    return results
