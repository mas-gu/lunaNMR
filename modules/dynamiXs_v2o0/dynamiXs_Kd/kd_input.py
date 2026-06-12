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


def load_titration_tracking(csv_path):
    """Load comprehensive_peak_tracking.csv (wide, per-point columns) into the same
    per-residue format as load_titration_tidy.

    Columns look like: Assignment, {label}_Position_X, {label}_Position_Y,
    {label}_Height, {label}_Volume, ... for each titration point label.
    """
    df = pd.read_csv(csv_path)
    labels = [c[:-len('_Position_X')] for c in df.columns if c.endswith('_Position_X')]
    labels = [lbl for lbl in labels if _to_point(lbl) is not None]
    labels.sort(key=_to_point)
    points = [_to_point(lbl) for lbl in labels]

    col_map = (('Position_X', 'ppm_x'), ('Position_Y', 'ppm_y'),
               ('Height', 'height'), ('Volume', 'volume'))
    residues = {}
    for _, row in df.iterrows():
        a = str(row.get('Assignment', '')).strip()
        if not a or a.lower() in ('nan', 'dummy', ''):
            continue
        rec = {'points': list(points)}
        for suffix, key in col_map:
            vals = []
            for lbl in labels:
                col = f'{lbl}_{suffix}'
                vals.append(float(row[col]) if col in df.columns and pd.notna(row[col]) else np.nan)
            rec[key] = vals
        residues[a] = rec
    return points, residues


def load_intensity_matrix(csv_path):
    """Load a peak_intensity/volume matrix: rows = residues, columns = point labels.

    These files (peak_intensity_matrix.csv, peak_volume_matrix.csv,
    peak_intensity_detected_matrix.csv) hold only intensities — no per-point
    positions — so they support intensity Kd but not CSP.
    """
    df = pd.read_csv(csv_path)
    point_cols = [c for c in df.columns if _to_point(c) is not None]
    point_cols.sort(key=_to_point)
    points = [_to_point(c) for c in point_cols]

    residues = {}
    for _, row in df.iterrows():
        a = str(row.get('Assignment', '')).strip()
        if not a or a.lower() in ('nan', 'dummy', ''):
            continue
        vals = [float(row[c]) if pd.notna(row[c]) else np.nan for c in point_cols]
        nan = [np.nan] * len(points)
        residues[a] = {'points': list(points), 'ppm_x': list(nan), 'ppm_y': list(nan),
                       'height': vals, 'volume': list(vals)}
    return points, residues


def load_titration(csv_path):
    """Auto-detect input format (tidy / tracking / intensity matrix) and load it."""
    cols = list(pd.read_csv(csv_path, nrows=0).columns)
    if 'spectrum_name' in cols:
        return load_titration_tidy(csv_path)
    if any(c.endswith('_Position_X') for c in cols):
        return load_titration_tracking(csv_path)
    return load_intensity_matrix(csv_path)


def csp_series(residue, alpha=0.14):
    """Per-point CSP for one residue, relative to the first (reference) point.

    Undetected points are recorded with position 0.0 (a sentinel, not a real ppm);
    treat them as NaN so they are excluded from the fit instead of producing a
    huge spurious CSP (0.0 vs the ~7/120 ppm reference).
    """
    px = np.asarray(residue['ppm_x'], dtype=float)
    py = np.asarray(residue['ppm_y'], dtype=float)
    px = np.where(px == 0.0, np.nan, px)
    py = np.where(py == 0.0, np.nan, py)
    dH = px - px[0]
    dN = py - py[0]
    return list(compute_csp(dH, dN, alpha=alpha))


def intensity_ratio_series(residue, value='height'):
    """Per-point intensity normalised to the first (reference) point.

    If the reference point is missing/non-positive (NaN, 0, or negative — peak not
    detected or below noise at L=0) the ratio is undefined, so return all-NaN
    (excluded from the fit) rather than raw unnormalised values.
    """
    v = np.asarray(residue[value], dtype=float)
    ref = v[0]
    if not np.isfinite(ref) or ref <= 0.0:
        return [float('nan')] * len(v)
    return list(v / ref)
