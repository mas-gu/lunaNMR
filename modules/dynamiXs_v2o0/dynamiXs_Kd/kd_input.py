# ABOUTME: Reads the LunaNMR series tidy CSV into per-residue titration observables.
# ABOUTME: Positions (for CSP) and intensities (for ratio), aligned to sorted points.

import numpy as np
import pandas as pd

from kd_models import compute_csp


def _to_point(label):
    """Parse a spectrum_name/point label to a float, or None if non-numeric."""
    try:
        return float(str(label).replace('o', '.'))
    except (ValueError, TypeError):
        return None


def load_titration_tidy(csv_path):
    """Load series_analysis_tidy.csv into per-residue titration arrays.

    Returns
    -------
    points : list[float]
        Sorted titration point values (ascending).
    residues : dict[str, dict]
        assignment -> {'points', 'ppm_x', 'ppm_y', 'height', 'volume'} lists,
        each aligned to `points` (missing entries are NaN).
    """
    df = pd.read_csv(csv_path)
    df = df.copy()
    df['_pt'] = df['spectrum_name'].map(_to_point)
    df = df[df['_pt'].notna()]
    points = sorted(df['_pt'].unique())

    cols = ['ppm_x', 'ppm_y', 'height', 'volume']
    residues = {}
    for assignment, g in df.groupby('assignment'):
        by_pt = g.set_index('_pt')
        row = {'points': list(points)}
        for c in cols:
            if c in by_pt.columns:
                row[c] = [float(by_pt[c].get(p, np.nan)) if p in by_pt.index else np.nan
                          for p in points]
            else:
                row[c] = [np.nan] * len(points)
        residues[str(assignment)] = row
    return points, residues


def csp_series(residue, alpha=0.14):
    """Per-point CSP for one residue, relative to the first (reference) point."""
    px = np.asarray(residue['ppm_x'], dtype=float)
    py = np.asarray(residue['ppm_y'], dtype=float)
    dH = px - px[0]
    dN = py - py[0]
    return list(compute_csp(dH, dN, alpha=alpha))


def intensity_ratio_series(residue, value='height'):
    """Per-point intensity normalised to the first (reference) point."""
    v = np.asarray(residue[value], dtype=float)
    ref = v[0]
    if not ref:
        return list(v)
    return list(v / ref)
