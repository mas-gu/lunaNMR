# ABOUTME: Reads the LunaNMR series tidy CSV into per-residue titration observables.
# ABOUTME: Positions (for CSP) and intensities (for ratio), aligned to sorted points.

import re

import numpy as np
import pandas as pd

from kd_models import compute_csp, REF_MAX_RATIO, reference_usable




# Below this many sigma a peak is not detected. The same floor the series pipeline
# uses for holding a position at the noise; 3 sigma is the usual bar for calling a
# signal present at all.
DETECTION_SIGMA_FLOOR = 3.0


def _to_point(label):
    """Parse a spectrum_name/point label to a float order key, or None.

    Accepts bare numbers ('0.0', '10o5' -> 0.0, 10.5) and real spectrum names,
    from which the trailing numeric token is taken as the point index
    (e.g. '03_2D_series_ref_001' -> 1.0). Actual concentrations come from the
    caller's --conc; the point value only needs to order and count the points.
    """
    s = str(label).strip()
    try:
        return float(s.replace('o', '.'))
    except (ValueError, TypeError):
        pass
    match = re.search(r'(\d+(?:[o.]\d+)?)(?:_\d+)?$', s)
    if match:
        try:
            return float(match.group(1).replace('o', '.'))
        except ValueError:
            return None
    return None


def load_titration_tidy(csv_path):
    """Load series_analysis_tidy.csv into per-residue titration arrays.

    Returns
    -------
    points : list[float]
        Sorted titration point values (ascending).
    residues : dict[str, dict]
        assignment -> {'points', 'ppm_x', 'ppm_y', 'height', 'volume'} lists,
        each aligned to `points` (missing entries are NaN), plus 'n_failed_points'
        counting the points this residue lost to a failed fit.

    Points the series marked quality='Failed' are blanked to NaN so they drop out of
    the fit. The series writes a sentinel height (0.1, or a negative amplitude) on a
    failed fit; read as a measurement those become a confident intensity of zero and
    an unphysical ratio. Blanking is POINT-WISE: a residue that failed at one point
    keeps every other point.
    """
    df = pd.read_csv(csv_path)
    df = df.copy()
    df['_pt'] = df['spectrum_name'].map(_to_point)
    df = df[df['_pt'].notna()]
    points = sorted(df['_pt'].unique())

    cols = ['ppm_x', 'ppm_y', 'height', 'volume']
    positional = ['ppm_x', 'ppm_y']
    failed = (df['quality'].astype(str).str.strip().str.lower() == 'failed'
              if 'quality' in df.columns else pd.Series(False, index=df.index))

    # A value taken from a BOUND is not a measurement of anything. The intensity floor
    # is cluster-relative - 0.1% of the brightest member's seed - so a weak member of a
    # bright cluster pinned there reports a product of its NEIGHBOURS' seeds and its
    # own seeded linewidths, carrying nothing from the data at its own position. A
    # pinned position means the fit could not reach the peak, and the linewidth
    # stretched to cover the gap, so the intensity is wrong too. Both are blanked
    # everywhere, exactly like a failed fit.
    pinned = pd.Series(False, index=df.index)
    for column in ('position_pinned', 'intensity_pinned'):
        if column in df.columns:
            # fillna(False), not truthiness: a legacy CSV has no such column at all
            # and a legacy ROW carries NaN. Absent is not the same as pinned, and
            # refusing every series written before these existed would be worse than
            # the risk. Re-run the series to get the filter.
            pinned |= df[column].fillna(False).astype(bool)
    df.loc[failed | pinned, [c for c in cols if c in df.columns]] = np.nan

    # A peak that has DECAYED is the opposite case and needs the opposite treatment.
    # Its intensity is the measurement - blanking it invents a plateau that is not
    # there, which is the whole point of an I/I0 series - but it has no position, so
    # it cannot contribute a chemical-shift perturbation. Only the ppm columns go.
    if 'detection_sigma' in df.columns:
        sigma = df['detection_sigma']
        undetected = sigma.notna() & (sigma < DETECTION_SIGMA_FLOOR)
        df.loc[undetected, [c for c in positional if c in df.columns]] = np.nan

    residues = {}
    for assignment, g in df.groupby('assignment'):
        row = {'points': list(points),
               'n_failed_points': int(failed.reindex(g.index, fill_value=False).sum())}
        for c in cols:
            if c in g.columns:
                # Collapse any duplicate rows at the same point to one value (mean,
                # skipping NaN). A residue listed twice at a point would otherwise make
                # the per-point lookup return a Series and crash float().
                per_pt = g.groupby('_pt')[c].mean()
                row[c] = [float(per_pt.get(p, np.nan)) for p in points]
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


def intensity_ratio_series(residue, value='height', ref_max_ratio=None):
    """Per-point intensity normalised to the first (reference) point.

    If the reference point is missing/non-positive (NaN, 0, or negative — peak not
    detected or below noise at L=0) the ratio is undefined, so return all-NaN
    (excluded from the fit) rather than raw unnormalised values.

    A reference orders of magnitude below the residue's own series maximum is also
    unusable: binding lowers intensity, so a ratio far above 1 means the reference is
    wrong, not that the peak grew. Such a reference rescales the WHOLE residue rather
    than degrading one point, and it survives a quality label of 'Good'.

    Non-positive intensities at other points are undefined too, and are dropped
    point-wise. kd_models.peak_present already defines presence as finite and > 0;
    this keeps the fit path honouring the same rule.
    """
    v = np.asarray(residue[value], dtype=float)
    ref = v[0]
    if not np.isfinite(ref) or ref <= 0.0:
        return [float('nan')] * len(v)
    if not reference_usable(v, ref, ref_max_ratio):
        return [float('nan')] * len(v)
    v = np.where(np.isfinite(v) & (v > 0.0), v, np.nan)
    return list(v / ref)
