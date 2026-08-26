# ABOUTME: 1:1 binding isotherm models for Kd titration fitting (CSP and intensity).
# ABOUTME: Full quadratic fraction-bound (accounts for protein concentration P0).

import re

import numpy as np


# A reference intensity this many times below the residue's own series maximum is a
# broken denominator, not a measurement: binding lowers intensity, so I/I0 far above 1
# means the L=0 point is wrong. Measured across two datasets (n=129): legitimate ratios
# top out at 1.30, the next value up is 49.9, and nothing lies between.
# Defined here so the fit path and the figure path cannot disagree about it.
REF_MAX_RATIO = 10.0


def reference_usable(v, ref, ref_max_ratio=None):
    """Whether a residue's reference intensity can serve as a denominator."""
    if not np.isfinite(ref) or ref <= 0.0:
        return False
    finite = v[np.isfinite(v)]
    cap = REF_MAX_RATIO if ref_max_ratio is None else float(ref_max_ratio)
    return not (finite.size and float(np.max(finite)) / ref > cap)


def fraction_bound(L, P, Kd):
    """Fraction of protein bound for a 1:1 interaction (full quadratic isotherm).

    Accounts for ligand depletion (free ligand != total) via the exact solution
    of P0*L0 = (Kd + free)·bound. Valid for any [P], [L]; reduces to the
    hyperbola L/(L+Kd) only when L >> P.

    Parameters
    ----------
    L : float or array
        Total ligand concentration.
    P : float
        Total protein concentration (P0).
    Kd : float
        Dissociation constant (same units as L and P).
    """
    L = np.asarray(L, dtype=float)
    s = P + L + Kd
    disc = np.maximum(s * s - 4.0 * P * L, 0.0)
    fb = (s - np.sqrt(disc)) / (2.0 * P)
    return fb if fb.ndim else float(fb)


def csp_model(L, dd_max, Kd, P):
    """Observed CSP for a 1:1 titration: Δδ_obs = Δδ_max · fraction_bound."""
    return dd_max * fraction_bound(L, P, Kd)


def intensity_decay(L, I0, I_inf, Kd):
    """Intensity decay for a titration with a plateau:

        I = I_inf + (I0 - I_inf) · exp(-L / Kd)

    Extends the lab's KD model (A·exp(-x/t2)) with a fitted floor I_inf so the
    curve levels off at the residual intensity peaks keep at saturation instead of
    being forced to zero. I0 = intensity at L=0, I_inf = plateau at high [L],
    Kd = decay constant (the [L] where the decaying part falls to 1/e).
    """
    L = np.asarray(L, dtype=float)
    out = I_inf + (I0 - I_inf) * np.exp(-L / Kd)
    return out if out.ndim else float(out)


def compute_csp(dH, dN, alpha=0.14):
    """Combined chemical-shift perturbation magnitude.

    CSP = sqrt(Δδ_H^2 + (alpha · Δδ_N)^2), where alpha scales the 15N shift to
    1H-equivalent ppm (default 0.14).
    """
    dH = np.asarray(dH, dtype=float)
    dN = np.asarray(dN, dtype=float)
    csp = np.sqrt(dH * dH + (alpha * dN) ** 2)
    return csp if csp.ndim else float(csp)


def peak_present(series, k, value='height'):
    """True if a residue's peak is genuinely present at point k: a detected position
    (finite, non-zero sentinel) AND a real intensity (finite, > 0). Used so CSP and
    intensity treat the same residues as missing."""
    px = np.asarray(series.get('ppm_x', []), float)
    py = np.asarray(series.get('ppm_y', []), float)
    v = np.asarray(series.get(value, []), float)
    if k < 0 or k >= min(len(px), len(py), len(v)):
        return False
    return (np.isfinite(px[k]) and px[k] != 0.0
            and np.isfinite(py[k]) and py[k] != 0.0
            and np.isfinite(v[k]) and v[k] > 0.0)


def pair_observable(series, i, j, obs, alpha=0.14, value='height'):
    """Observable for one residue between reference point i and point j.

    obs='csp'  -> CSP (ppm) from the position shift: sqrt(ΔδH² + (alpha·ΔδN)²).
    obs='intensity' -> ratio v[j]/v[i] of the chosen value ('height'/'volume').
    Returns NaN if an endpoint is missing — a 0.0 position is the undetected
    sentinel, a non-positive reference intensity makes the ratio undefined.
    """
    if obs == 'csp':
        px = np.asarray(series.get('ppm_x', []), float)
        py = np.asarray(series.get('ppm_y', []), float)
        if max(i, j) >= len(px) or max(i, j) >= len(py):
            return float('nan')
        xi, xj, yi, yj = px[i], px[j], py[i], py[j]
        if 0.0 in (xi, xj, yi, yj) or not np.all(np.isfinite([xi, xj, yi, yj])):
            return float('nan')
        return float(compute_csp(xj - xi, yj - yi, alpha=alpha))
    v = np.asarray(series.get(value, []), float)
    if max(i, j) >= len(v):
        return float('nan')
    ref = v[i]
    # Same rule the fit path applies. A residue the fitter dropped for a broken reference
    # must not reappear in a figure at a ratio of 1e9 — one such bar flattens every real
    # one, and the figure would be describing a residue the fit already rejected.
    if not reference_usable(v, ref):
        return float('nan')
    if not np.isfinite(v[j]) or v[j] <= 0.0:
        return float('nan')
    return float(v[j] / ref)


def ref_point_values(fits, i, j, obs, alpha=0.14, value='height'):
    """(names, vals) of the observable between reference point i and point j for each
    fit, NaN where the peak is absent (per peak_present) at either endpoint. Shared by
    the GUI viewer's ref→point bars and the CLI export."""
    names, vals = [], []
    for f in fits:
        series = f.get('series')
        if series and peak_present(series, i, value) and peak_present(series, j, value):
            v = pair_observable(series, i, j, obs, alpha=alpha, value=value)
        else:
            v = float('nan')
        names.append(f.get('residue'))
        vals.append(v)
    return names, vals


def residue_sort_key(name):
    """Order residues by sequence number, not amino-acid letter: 'K14' before
    'A17'. Names without a number sort last, alphabetically among themselves."""
    m = re.search(r'\d+', str(name))
    return (int(m.group()) if m else float('inf'), str(name))


def global_fit_table(fits, global_fit, obs, P0):
    """Per-residue table behind the global shared-Kd fit figure. Returns (header, rows).

    For each residue pooled into the global fit (present in the shared-Kd amplitude map
    and carrying >=2 observed points), reports its amplitude param(s), the single shared
    Kd, and the R^2 of that one-Kd model against the residue's observed points — the same
    per-panel R^2 shown in <obs>_global_fit.pdf. Rows are ordered by residue sequence
    number. Rows are empty if the global fit for `obs` is absent or failed.

    obs='intensity' -> header ['Residue', 'I0', 'I_inf', 'global_Kd', 'R2']
    obs='csp'       -> header ['Residue', 'dd_max', 'global_Kd', 'R2']
    """
    header = (['Residue', 'dd_max', 'global_Kd', 'R2'] if obs == 'csp'
              else ['Residue', 'I0', 'I_inf', 'global_Kd', 'R2'])
    g = global_fit.get(obs) or {}
    amp = g.get('dd_max') if obs == 'csp' else g.get('I0')
    kd = g.get('Kd')
    rows = []
    if not g.get('success') or not amp or kd is None:
        return header, rows
    for f in fits:
        res = f.get('residue')
        fit = f.get(obs) or {}
        if res not in amp or not fit.get('L'):
            continue
        L = np.asarray(fit['L'], dtype=float)
        y = np.asarray(fit['obs'], dtype=float)
        good = np.isfinite(L) & np.isfinite(y)
        if good.sum() < 2:
            continue
        if obs == 'csp':
            yhat = csp_model(L[good], amp[res], kd, P0)
            params = [float(amp[res])]
        else:
            yhat = intensity_decay(L[good], g['I0'][res], g['I_inf'][res], kd)
            params = [float(g['I0'][res]), float(g['I_inf'][res])]
        ss_res = float(np.sum((y[good] - yhat) ** 2))
        ss_tot = float(np.sum((y[good] - np.mean(y[good])) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float('nan')
        rows.append([res] + params + [float(kd), r2])
    rows.sort(key=lambda r: residue_sort_key(r[0]))
    return header, rows


def ref_vs_point_table(fits, labels, obs, ref=0, alpha=0.14, value='height'):
    """Wide table of the observable per residue across titration points, relative to
    the reference point `ref`. Returns (header, rows).

    header = ['Residue', <label0>, <label1>, ...]; each row = [residue, v0, v1, ...]
    with '' where the observable is undefined (peak absent at that point or the
    reference). Rows are ordered by residue sequence number. Built from the embedded
    per-point `series`, so it works from a saved fit JSON with no re-fit.
    """
    header = ['Residue'] + [str(lbl) for lbl in labels]
    names = [f.get('residue') for f in fits]
    columns = [ref_point_values(fits, ref, j, obs, alpha=alpha, value=value)[1]
               for j in range(len(labels))]
    rows = []
    for idx, name in enumerate(names):
        cells = [('' if columns[j][idx] is None or not np.isfinite(columns[j][idx])
                  else columns[j][idx])
                 for j in range(len(labels))]
        rows.append([name] + cells)
    rows.sort(key=lambda r: residue_sort_key(r[0]))
    return header, rows
