# ABOUTME: Extracts selected peaks across an ordered series of 1D spectra.
# ABOUTME: Backs the interactive picker - box select, click select, and one-shot series integration.

import numpy as np

from oned_detector import (DEFAULT_MIN_SNR, detect_peaks, estimate_fwhm,
                           estimate_window, match_reference_peaks)
from oned_fitter import fit_peak_voigt, integrate_region, measure_intensity
from oned_loader import load_spectrum

# How far a peak may move between series points and still be the same peak. Wide enough
# for real drift, far below the spacing of resolvable peaks.
DEFAULT_TRACK_WINDOW_PPM = 0.01

# A peak may move at most this many of its own linewidths between consecutive points.
DRIFT_FWHM_MULTIPLE = 2.0


def load_series(paths):
    """Load spectra in the given order."""
    return [load_spectrum(p) for p in paths]


def peaks_in_range(spectrum, ppm_from, ppm_to, min_snr=DEFAULT_MIN_SNR):
    """Peaks detected inside a ppm range, ordered high ppm first.

    Backs box selection: the user drags a rectangle over the region of interest and
    everything detectable inside it becomes a peak.
    """
    lo, hi = sorted((float(ppm_from), float(ppm_to)))

    return [p for p in detect_peaks(spectrum, min_snr=min_snr) if lo <= p['ppm'] <= hi]


def peak_at_position(spectrum, ppm, search=0.01, min_snr=DEFAULT_MIN_SNR):
    """The peak nearest a clicked position, snapped to the local maximum.

    Backs click selection. When there is no detectable peak within `search` the click
    is taken at face value and flagged undetected, so a deliberate pick in a weak
    region still works.
    """
    candidates = [p for p in detect_peaks(spectrum, min_snr=min_snr)
                  if abs(p['ppm'] - ppm) <= search]

    if candidates:
        nearest = min(candidates, key=lambda p: abs(p['ppm'] - ppm))
        return {**nearest, 'detected': True}

    index = int(np.argmin(np.abs(spectrum.ppm_axis - ppm)))

    return {'ppm': float(spectrum.ppm_axis[index]), 'index': index,
            'height': float(spectrum.data[index]),
            'snr': float(spectrum.data[index] / spectrum.noise) if spectrum.noise else None,
            'detected': False}


def _drift_limits(spectrum, peaks, track_window):
    """How far each peak may move between consecutive points.

    Bounded by the peak's own linewidth: a resonance drifting smoothly through a
    series moves a small fraction of its width per point, so a candidate several
    linewidths away is a different resonance, not the same one relocated.

    Without this a peak that decays to nothing hands its identity to whatever is
    nearest - a substrate tracked through a conversion jumps onto the product line
    the moment it disappears and reports the product's intensity, so the decay reads
    as a full recovery. The neighbour's position cannot be used for the bound because
    at the first point the product does not exist yet.
    """
    limits = []
    for peak in peaks:
        fwhm = estimate_fwhm(spectrum, float(peak['position']), search=track_window)
        limit = track_window
        if fwhm:
            limit = min(limit, DRIFT_FWHM_MULTIPLE * fwhm)
        limits.append(limit)

    return limits


def locate_peaks(spectrum, peaks, track_window=DEFAULT_TRACK_WINDOW_PPM):
    """Where each reference peak sits in this one spectrum.

    Same matching a series run uses, so what the display shows for a spectrum is what
    the integration will measure there.
    """
    if not peaks:
        return []

    return match_reference_peaks(spectrum, peaks,
                                 window=_drift_limits(spectrum, peaks, track_window))


def _windows_for(spectrum, peaks):
    """Measurement half-width per peak, linewidth-driven and clear of its neighbours.

    Neighbours are every resonance detected in the spectrum, not only the ones the
    user happened to pick. A window sized against the peak list alone makes a
    measurement depend on what else is in that list: picking one line of a close pair
    gave a window half again as wide, R2 falling to 0.64, and a fitted area 19% off
    the value obtained when both were picked.
    """
    detected = [p['ppm'] for p in detect_peaks(spectrum)]
    positions = [float(p['position']) for p in peaks]

    windows = []
    for i, position in enumerate(positions):
        # the nearest detected maximum is this peak itself, so it is not a neighbour
        nearby = sorted(detected, key=lambda found: abs(found - position))[1:]
        neighbours = positions[:i] + positions[i + 1:] + nearby
        windows.append(estimate_window(spectrum, position, neighbours=neighbours))

    return windows


def integrate_series(spectra, peaks, names=None, track_window=DEFAULT_TRACK_WINDOW_PPM,
                     progress=None):
    """Measure every selected peak in every spectrum of the series.

    Peaks are re-matched at each point within `track_window`, so a resonance that
    drifts through the series is followed rather than measured at a stale position.
    Matching is competitive, so two peaks can never collapse onto one maximum - the
    failure that silently turns an intensity ratio into 1.0.

    Both intensity and area are reported: intensity is the direct observable, area is
    the one that stays proportional to concentration when the linewidth changes.
    """
    names = names or [s.metadata.get('source', f'point_{k}') for k, s in enumerate(spectra)]

    table = []
    current = [dict(p) for p in peaks]

    limits = _drift_limits(spectra[0], current, track_window) if spectra else []

    for point, (spectrum, name) in enumerate(zip(spectra, names)):
        rows_before = len(table)
        matched = match_reference_peaks(spectrum, current, window=limits)
        windows = _windows_for(spectrum, current)
        scale = spectrum.intensity_scale

        for peak, match, window in zip(current, matched, windows):
            target = match['ppm'] if match['matched'] else peak['position']

            # an unmatched peak is read where it should be, not where the tallest
            # thing in range happens to sit
            intensity = measure_intensity(spectrum.ppm_axis, spectrum.data, target,
                                          window=window, assignment=peak.get('assignment'),
                                          locate=match['matched'])
            area = integrate_region(spectrum.ppm_axis, spectrum.data, target,
                                    window=window, assignment=peak.get('assignment'))

            # The fit is the only measurement that carries a goodness-of-fit and a
            # linewidth. Both matter here: when a resonance broadens through a series
            # its height stops tracking concentration, and nothing else reveals that.
            fit = fit_peak_voigt(spectrum.ppm_axis, spectrum.data, target,
                                 window=window, assignment=peak.get('assignment'))

            table.append({
                'point': point,
                'spectrum': name,
                'assignment': peak.get('assignment'),
                'position': peak['position'],
                'ppm': match['ppm'],
                'matched': match['matched'],
                'height': intensity['height'] * scale if intensity['success'] else None,
                'area': area['area'] * scale if area['success'] else None,
                # r_squared and fwhm describe the shape, not the amount, so the
                # acquisition scale must not touch them
                'fit_area': fit['area'] * scale if fit['success'] else None,
                'r_squared': fit['r_squared'] if fit['success'] else None,
                'fwhm': fit['fwhm'] if fit['success'] else None,
                'snr': match['snr'],
                'window': window,
                'intensity_scale': scale,
            })

            # follow the peak: the next point looks around where it was found here
            if match['matched']:
                peak['position'] = match['ppm']

        if progress is not None:
            progress(point, table[rows_before:])

    return table


def write_series_csv(table, path, value='height'):
    """Write the series as one row per spectrum and one column per peak.

    spectrum,PeakA,PeakB
    1D_KB_GTP_001,5.77e+10,1.00e+12

    A peak that could not be measured at a point leaves its cell empty rather than
    writing a zero, which would read as a real measurement of nothing.
    """
    assignments = []
    for row in table:
        if row['assignment'] not in assignments:
            assignments.append(row['assignment'])

    by_point = {}
    for row in table:
        entry = by_point.setdefault(row['point'], {'spectrum': row['spectrum']})
        entry[row['assignment']] = row.get(value)

    lines = [','.join(['spectrum'] + [str(a) for a in assignments])]
    for point in sorted(by_point):
        entry = by_point[point]
        cells = [entry['spectrum']]
        for assignment in assignments:
            measurement = entry.get(assignment)
            cells.append('' if measurement is None else repr(float(measurement)))
        lines.append(','.join(cells))

    path = str(path)
    with open(path, 'w') as handle:
        handle.write('\n'.join(lines) + '\n')

    return path


def series_matrix(table, value='height'):
    """Pivot a series table into (assignments, matrix[peak, point])."""
    assignments = sorted({r['assignment'] for r in table},
                         key=lambda a: (a is None, a))
    points = sorted({r['point'] for r in table})

    matrix = np.full((len(assignments), len(points)), np.nan)
    for row in table:
        matrix[assignments.index(row['assignment']), points.index(row['point'])] = (
            np.nan if row[value] is None else row[value])

    return assignments, matrix
