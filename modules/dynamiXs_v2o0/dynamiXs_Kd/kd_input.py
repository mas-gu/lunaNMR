# ABOUTME: Reads the LunaNMR series tidy CSV into per-residue titration observables.
# ABOUTME: Positions (for CSP) and intensities (for ratio), aligned to sorted points.

import json
import math
import os
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
    A '_<n>' repeat suffix ('1_2', the second acquisition at point 1) names the same
    point, so the bare-number read skips any label with an underscore: float()
    accepts one as a digit separator and reads '1_2' as 12.
    """
    s = str(label).strip()
    if '_' not in s:
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


def _refuse_unread_points(labels):
    """Raise if any point label has no titration point to read.

    Concentrations are matched to the points by position, so a point dropped here
    would move every concentration after it onto the wrong spectrum.
    """
    unread = sorted({str(lbl) for lbl in labels if _to_point(lbl) is None})
    if unread:
        raise ValueError(
            f"no titration point can be read from {len(unread)} spectrum label(s): "
            f"{', '.join(unread)}. Each spectrum must be named by its point "
            f"(e.g. '0.5', '0o5' or '..._001').")


# Written by `series` beside its CSVs: every column with the spectrum it came from, and the
# point it parsed from the filename (`value: null` when it could not).
SERIES_METADATA_NAME = 'series_metadata.json'


def point_labels(csv_path):
    """The titration point labels of a titration CSV, as written, in file order."""
    cols = list(pd.read_csv(csv_path, nrows=0).columns)
    if 'spectrum_name' in cols:
        names = pd.read_csv(csv_path, usecols=['spectrum_name'], dtype=str)['spectrum_name']
        return list(dict.fromkeys(names.str.strip()))
    if any(c.endswith('_Position_X') for c in cols):
        return [c[:-len('_Position_X')] for c in cols if c.endswith('_Position_X')]
    return [c for c in cols if _to_point(c) is not None]


def _unparsed_series_columns(csv_path):
    """The columns `series` named after a file stem because it read no point from the
    filename, or None when no readable series_metadata.json sits beside the CSV."""
    path = os.path.join(os.path.dirname(os.path.abspath(csv_path)), SERIES_METADATA_NAME)
    try:
        with open(path) as fh:
            columns = json.load(fh)['columns']
        return {str(c['column']) for c in columns if c.get('value') is None}
    except (OSError, ValueError, KeyError, TypeError):
        return None


def _is_number_label(label):
    """Whether a label is a number, as `series` writes titration points ('0.125', '3'),
    a '_<n>' repeat suffix allowed ('1_2'). 'o' reads as the decimal point ('0o5')."""
    s = re.sub(r'_\d+$', '', str(label).strip())
    if '_' in s:
        return False
    try:
        return math.isfinite(float(s.replace('o', '.')))
    except ValueError:
        return False


def labels_without_points(csv_path):
    """The point labels of `csv_path` that do not name a titration point.

    Those `series` marked unparsed in the series_metadata.json beside the CSV, or, with
    no such file, every label that is not a number. Such a label still orders the points
    (_to_point reads its trailing digits), but its digits are an acquisition index, not a
    concentration.
    """
    labels = point_labels(csv_path)
    unparsed = _unparsed_series_columns(csv_path)
    if unparsed is not None:
        return [lbl for lbl in labels if lbl in unparsed]
    return [lbl for lbl in labels if not _is_number_label(lbl)]


def unnamed_points(csv_path):
    """{point: [label, ...]} for the points of `csv_path` (as load_titration keys them)
    whose labels name no titration point (labels_without_points)."""
    points = {}
    for lbl in labels_without_points(csv_path):
        points.setdefault(_to_point(lbl), []).append(lbl)
    return points


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
    # Residue names are text: an all-numeric column with a blank would read as float.
    df = pd.read_csv(csv_path, dtype={'assignment': str})
    df = df.copy()
    _refuse_unread_points(df['spectrum_name'].unique())
    df['_pt'] = df['spectrum_name'].map(_to_point)
    points = sorted(df['_pt'].unique())

    cols = ['ppm_x', 'ppm_y', 'height', 'volume']
    positional = ['ppm_x', 'ppm_y']
    failed = (df['quality'].astype(str).str.strip().str.lower() == 'failed'
              if 'quality' in df.columns else pd.Series(False, index=df.index))

    # A value taken from the INTENSITY floor is not a measurement of anything. The floor
    # is cluster-relative - 0.1% of the brightest member's seed - so a weak member of a
    # bright cluster pinned there reports a product of its NEIGHBOURS' seeds and its
    # own seeded linewidths, carrying nothing from the data at its own position. It is
    # blanked everywhere, exactly like a failed fit. A POSITION pinned on its margin is
    # not a measured shift, so only the ppm columns go; the intensity does not depend
    # on it (lunaNMR.utils.quality_categories.usable_points has the measurement).
    def _flag(column):
        # fillna(False), not truthiness: a legacy CSV has no such column at all and a
        # legacy ROW carries NaN. Absent is not the same as pinned, and refusing every
        # series written before these existed would be worse than the risk. Re-run
        # the series to get the filter.
        if column not in df.columns:
            return pd.Series(False, index=df.index)
        return df[column].fillna(False).astype(bool)

    df.loc[failed | _flag('intensity_pinned'), [c for c in cols if c in df.columns]] = np.nan
    df.loc[_flag('position_pinned'), [c for c in positional if c in df.columns]] = np.nan

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


def _join_repeats(points, values):
    """Join the repeat acquisitions of a point ('1_2' beside '1'): one value per distinct
    point, the mean of its finite values (NaN when it has none), as load_titration_tidy
    joins the rows of a repeat. `points` are the per-label points, in label order."""
    distinct = sorted(set(points))
    joined = []
    for p in distinct:
        finite = [v for q, v in zip(points, values) if q == p and np.isfinite(v)]
        joined.append(float(np.mean(finite)) if finite else np.nan)
    return distinct, joined


def load_titration_tracking(csv_path):
    """Load comprehensive_peak_tracking.csv (wide, per-point columns) into the same
    per-residue format as load_titration_tidy.

    Columns look like: Assignment, {label}_Position_X, {label}_Position_Y,
    {label}_Height, {label}_Volume, ... for each titration point label. A point whose
    {label}_Quality is 'Failed' is blanked in every column, as the tidy loader does.
    """
    df = pd.read_csv(csv_path, dtype={'Assignment': str})
    labels = [c[:-len('_Position_X')] for c in df.columns if c.endswith('_Position_X')]
    _refuse_unread_points(labels)
    labels.sort(key=_to_point)
    label_points = [_to_point(lbl) for lbl in labels]
    points = sorted(set(label_points))

    col_map = (('Position_X', 'ppm_x'), ('Position_Y', 'ppm_y'),
               ('Height', 'height'), ('Volume', 'volume'))
    residues = {}
    for _, row in df.iterrows():
        a = str(row.get('Assignment', '')).strip()
        if not a or a.lower() in ('nan', 'dummy', ''):
            continue
        failed = {lbl for lbl in labels
                  if str(row.get(f'{lbl}_Quality', '')).strip().lower() == 'failed'}
        rec = {'points': list(points), 'n_failed_points': len(failed)}
        # As load_titration_tidy: a point on its intensity floor is blanked in every
        # column, a position on its margin only in the position columns. The flags are
        # written as '1.0'/'0.0'; a file from before they existed has no such columns.
        def _pinned(lbl, flag):
            col = f'{lbl}_{flag}_Pinned'
            return col in df.columns and pd.notna(row[col]) and float(row[col]) != 0.0

        for suffix, key in col_map:
            vals = []
            for lbl in labels:
                col = f'{lbl}_{suffix}'
                is_pinned = (lbl in failed or _pinned(lbl, 'Intensity')
                             or (suffix.startswith('Position') and _pinned(lbl, 'Position')))
                usable = col in df.columns and pd.notna(row[col]) and not is_pinned
                vals.append(float(row[col]) if usable else np.nan)
            rec[key] = _join_repeats(label_points, vals)[1]
        residues[a] = rec
    return points, residues


def load_intensity_matrix(csv_path):
    """Load a peak_intensity/volume matrix: rows = residues, columns = point labels.

    These files (peak_intensity_matrix.csv, peak_volume_matrix.csv,
    peak_intensity_detected_matrix.csv) hold only intensities — no per-point
    positions — so they support intensity Kd but not CSP.

    Every CSV that is neither tidy nor tracking arrives here, so one that cannot be an
    intensity matrix is refused: with no point columns or no residue rows there is
    nothing to fit, and a matrix of 0/1 detection flags fits every residue to one Kd.
    """
    df = pd.read_csv(csv_path, dtype={'Assignment': str})
    point_cols = [c for c in df.columns if _to_point(c) is not None]
    if not point_cols:
        raise ValueError(f"no titration point columns in {csv_path}: it is not a series "
                         f"tidy, peak-tracking or intensity-matrix CSV")
    point_cols.sort(key=_to_point)
    col_points = [_to_point(c) for c in point_cols]
    points = sorted(set(col_points))

    residues = {}
    cells = []
    for _, row in df.iterrows():
        a = str(row.get('Assignment', '')).strip()
        if not a or a.lower() in ('nan', 'dummy', ''):
            continue
        raw = [float(row[c]) if pd.notna(row[c]) else np.nan for c in point_cols]
        cells.extend(raw)
        vals = _join_repeats(col_points, raw)[1]
        nan = [np.nan] * len(points)
        residues[a] = {'points': list(points), 'ppm_x': list(nan), 'ppm_y': list(nan),
                       'height': vals, 'volume': list(vals)}
    if not residues:
        raise ValueError(f"no residue rows (Assignment) in {csv_path}")
    # The cells as written, before repeats are joined: two flags average to 0.5
    values = np.array(cells, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size and np.isin(finite, (0.0, 1.0)).all():
        raise ValueError(f"every value in {csv_path} is 0 or 1: these are detection flags, "
                         f"not intensities. Use the intensity or volume matrix.")
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
