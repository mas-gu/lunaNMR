# ABOUTME: Reads a finished run's files for what its exit code hides: pinned and undetected points, Kd pool.
# ABOUTME: Counts and names only, never raw data, because everything returned reaches the model.

import csv
import json
from collections import Counter
from pathlib import Path

DETECTION_SIGMA_FLOOR = 3.0
KD_R2_FLOOR = 0.8
MAX_NAMED_RESIDUES = 30
NAME_CHARS = 64
REASON_CHARS = 200


def _cut(text, limit):
    return str(text)[:limit]


def _is_true(cell):
    return cell == 'True'


def _below_detection(cell):
    try:
        return float(cell) < DETECTION_SIGMA_FLOOR
    except ValueError:
        return False


def series_quality(out_folder):
    """Point counts a series' exit code does not show, from its tidy CSV."""
    tidy = Path(out_folder) / 'series_analysis_tidy.csv'
    if not tidy.is_file():
        return None
    with open(tidy, newline='') as handle:
        rows = list(csv.DictReader(handle))
    flagged = Counter()
    position_pinned = intensity_pinned = below_3_sigma = 0
    for row in rows:
        pinned_position = _is_true(row.get('position_pinned'))
        pinned_intensity = _is_true(row.get('intensity_pinned'))
        undetected = _below_detection(row.get('detection_sigma'))
        position_pinned += pinned_position
        intensity_pinned += pinned_intensity
        below_3_sigma += undetected
        if pinned_position or pinned_intensity or undetected:
            flagged[_cut(row['assignment'], NAME_CHARS)] += 1
    return {'points': len(rows),
            'quality': dict(Counter(row['quality'] for row in rows)),
            'position_pinned': position_pinned,
            'intensity_pinned': intensity_pinned,
            'below_3_sigma': below_3_sigma,
            'n_residues_with_flagged_points': len(flagged),
            'residues_with_flagged_points': [
                name for name, _ in flagged.most_common(MAX_NAMED_RESIDUES)]}


def _poorly_fitted(fit):
    r_squared = fit.get('r_squared')
    return r_squared is None or not r_squared >= KD_R2_FLOOR


MEANING = {
    'csp': 'Kd of a 1:1 binding isotherm, in the units of p0',
    'intensity': 'apparent decay constant of intensity against ligand, not a thermodynamic Kd '
                 '(it ignores [P]0); do not compare it directly with the CSP Kd'}


def kd_quality(out_folder):
    """The shared Kd per observable with the size of the pool behind it, from the fit JSON.

    Poor fits are counted among the residues that fitted; a residue with no fit is named under
    `unfitted` with the command's reason."""
    fit_files = sorted((Path(out_folder) / 'data').glob('*_kd_fit_data.json'))
    if not fit_files:
        return None
    fit = json.loads(fit_files[0].read_text())
    summary = {'n_total': len(fit['fits']),
               'unfitted': {_cut(name, NAME_CHARS): [_cut(reason, REASON_CHARS) for reason in reasons]
                            for name, reasons in fit['metadata'].get('unfitted', {}).items()},
               'quality_warnings': [_cut(warning, 2000)
                                    for warning in fit['metadata'].get('quality_warnings', [])]}
    for observable, shared in fit['global'].items():
        fitted = [residue[observable] for residue in fit['fits']
                  if residue.get(observable, {}).get('success')]
        summary[observable] = {
            'kd': shared.get('Kd'),
            'kd_err': shared.get('Kd_err'),
            'reliable': shared.get('reliable'),
            'meaning': MEANING.get(observable),
            'n_in_shared_fit': shared.get('n_residues'),
            f'n_fitted_below_r2_{KD_R2_FLOOR}': sum(_poorly_fitted(f) for f in fitted)}
    return summary
