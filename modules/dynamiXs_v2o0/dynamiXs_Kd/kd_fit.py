# ABOUTME: Kd titration fitter - per-residue and global shared-Kd, CSP and intensity.
# ABOUTME: run_kd_analysis_with_params orchestrates load -> fit -> JSON, dynamiXs-style.

import json
import os

import numpy as np
from scipy.optimize import curve_fit, least_squares

from kd_models import csp_model, intensity_decay
from kd_input import load_titration, csp_series, intensity_ratio_series

_MIN_POINTS = 3


def json_safe(obj):
    """Replace non-finite floats (NaN/inf) with None so output is valid JSON
    (bare NaN tokens are rejected by JS/jq/strict parsers)."""
    import math
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, dict):
        return {k: json_safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [json_safe(v) for v in obj]
    return obj


def _r_squared(y, yfit):
    ss_res = np.sum((y - yfit) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return float(1.0 - ss_res / ss_tot) if ss_tot > 0 else 0.0


def _clean(L, y):
    L = np.asarray(L, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = ~(np.isnan(L) | np.isnan(y))
    return L[mask], y[mask]


def fit_residue_csp(L, csp, P0, kd_init=None, n_bootstrap=0):
    """Fit one residue's CSP isotherm: Δδ_obs = Δδ_max · fraction_bound(L,P0,Kd)."""
    L, y = _clean(L, csp)
    if len(L) < _MIN_POINTS or np.all(y == 0):
        return {'success': False, 'reason': 'too few points or no shift'}
    kd0 = kd_init if kd_init else max(np.max(L) / 2.0, 1.0)
    dd0 = max(float(np.max(y)), 1e-6)
    try:
        popt, pcov = curve_fit(
            lambda L, dd, kd: csp_model(L, dd, kd, P0), L, y,
            p0=[dd0, kd0], bounds=([0, 0], [np.inf, np.inf]), maxfev=20000)
    except Exception as e:
        return {'success': False, 'reason': str(e)}
    dd, kd = popt
    perr = np.sqrt(np.diag(pcov))
    if n_bootstrap:
        bdd, bkd = _bootstrap_csp(L, y, P0, dd, kd, n_bootstrap)
        dd_e = bdd if np.isfinite(bdd) else float(perr[0])   # fall back to covariance
        kd_e = bkd if np.isfinite(bkd) else float(perr[1])
    else:
        dd_e, kd_e = float(perr[0]), float(perr[1])
    return {'success': True, 'Kd': float(kd), 'Kd_err': kd_e,
            'dd_max': float(dd), 'dd_max_err': dd_e,
            'r_squared': _r_squared(y, csp_model(L, dd, kd, P0)),
            'L': L.tolist(), 'obs': y.tolist()}


def fit_residue_intensity(L, obs, P0=None, kd_init=None, n_bootstrap=0):
    """Fit one residue's intensity decay: I = I0·exp(-L/Kd) (lab KD-script model).

    I0 is fitted (anchors the curve to the first point); Kd is the decay constant.
    P0 is unused (kept for signature parity with the CSP fitter).
    """
    L, y = _clean(L, obs)
    if len(L) < _MIN_POINTS:
        return {'success': False, 'reason': 'too few points'}
    kd0 = kd_init if kd_init else max(np.max(L) / 2.0, 1.0)
    i0 = max(float(y[0]), 1e-12)
    inf0 = max(float(np.min(y)), 0.0)            # plateau ~ lowest observed point
    hi = max(float(np.max(y)) * 2.0, 1.0)
    try:
        popt, pcov = curve_fit(
            intensity_decay, L, y, p0=[i0, inf0, kd0],
            bounds=([0, 0, 1e-9], [hi, hi, np.inf]), maxfev=20000)
    except Exception as e:
        return {'success': False, 'reason': str(e)}
    I0, I_inf, kd = popt
    perr = np.sqrt(np.diag(pcov))
    if n_bootstrap:
        bi0, binf, bkd = _bootstrap_intensity(L, y, I0, I_inf, kd, hi, n_bootstrap)
        i0_e = bi0 if np.isfinite(bi0) else float(perr[0])    # fall back to covariance
        inf_e = binf if np.isfinite(binf) else float(perr[1])
        kd_e = bkd if np.isfinite(bkd) else float(perr[2])
    else:
        i0_e, inf_e, kd_e = float(perr[0]), float(perr[1]), float(perr[2])
    return {'success': True, 'Kd': float(kd), 'Kd_err': kd_e,
            'I0': float(I0), 'I0_err': i0_e,
            'I_inf': float(I_inf), 'I_inf_err': inf_e,
            'r_squared': _r_squared(y, intensity_decay(L, I0, I_inf, kd)),
            'L': L.tolist(), 'obs': y.tolist()}


def _bootstrap_intensity(L, y, I0, I_inf, kd, hi, n):
    try:
        yfit = intensity_decay(L, I0, I_inf, kd)
        resid = y - yfit
        i0s, infs, kds = [], [], []
        for _ in range(n):
            yb = yfit + np.random.choice(resid, size=len(resid), replace=True)
            try:
                popt, _ = curve_fit(intensity_decay, L, yb, p0=[I0, I_inf, kd],
                                    bounds=([0, 0, 1e-9], [hi, hi, np.inf]), maxfev=20000)
                i0s.append(popt[0]); infs.append(popt[1]); kds.append(popt[2])
            except Exception:
                continue
        sd = lambda a: float(np.std(a)) if a else np.nan
        return sd(i0s), sd(infs), sd(kds)
    except Exception:
        return np.nan, np.nan, np.nan


def _bootstrap_csp(L, y, P0, dd, kd, n):
    try:
        yfit = csp_model(L, dd, kd, P0)
        resid = y - yfit
        dds, kds = [], []
        for _ in range(n):
            yb = yfit + np.random.choice(resid, size=len(resid), replace=True)
            try:
                popt, _ = curve_fit(lambda L, d, k: csp_model(L, d, k, P0), L, yb,
                                    p0=[dd, kd], bounds=([0, 0], [np.inf, np.inf]), maxfev=20000)
                dds.append(popt[0]); kds.append(popt[1])
            except Exception:
                continue
        return (float(np.std(dds)) if dds else np.nan,
                float(np.std(kds)) if kds else np.nan)
    except Exception:
        return np.nan, np.nan


def fit_global_kd_csp(residues, L, P0):
    """Joint CSP fit with one shared Kd and a per-residue Δδ_max."""
    L = np.asarray(L, dtype=float)
    # drop residues with no usable points (all-NaN) so np.nanmax can't yield NaN
    names = [n for n in residues if not np.all(np.isnan(np.asarray(residues[n], dtype=float)))]
    if len(names) < 2:
        return {'success': False, 'reason': 'fewer than 2 usable residues'}
    series = {n: np.asarray(residues[n], dtype=float) for n in names}
    dd0 = [max(float(np.nanmax(series[n])), 1e-6) for n in names]
    kd0 = max(float(np.nanmax(L)) / 2.0, 1.0)

    def resid(p):
        kd = p[0]
        out = []
        for i, n in enumerate(names):
            y = series[n]
            mask = ~np.isnan(y)
            out.append(csp_model(L[mask], p[1 + i], kd, P0) - y[mask])
        return np.concatenate(out)

    lo = [0.0] * (1 + len(names))
    hi = [np.inf] * (1 + len(names))
    try:
        res = least_squares(resid, [kd0] + dd0, bounds=(lo, hi), max_nfev=20000)
    except Exception as e:
        return {'success': False, 'reason': str(e)}
    return {'success': True, 'Kd': float(res.x[0]),
            'dd_max': {n: float(res.x[1 + i]) for i, n in enumerate(names)},
            'n_residues': len(names)}


def run_kd_analysis_with_params(params, progress_callback=None):
    """Load a titration tidy CSV, fit Kd per-residue (CSP/intensity) + global, write JSON."""
    points, residues = load_titration(params['input_csv_file'])
    concs = list(params.get('concentrations') or points)
    if len(concs) != len(points):
        raise ValueError(f"{len(concs)} concentrations given but {len(points)} titration points found")
    L = np.asarray(concs, dtype=float)
    P0 = float(params['protein_conc'])
    alpha = float(params.get('alpha', 0.14))
    observables = params.get('observables', ['csp', 'intensity'])
    n_boot = int(params.get('n_bootstrap', 0))

    # Per-point intensity/volume scaling (e.g. correct a point acquired with a
    # different number of scans). Applies ONLY to height/volume, not positions.
    scales = params.get('intensity_scales')
    if scales:
        if len(scales) != len(points):
            raise ValueError(f"{len(scales)} intensity scales but {len(points)} titration points")
        for res in residues.values():
            for key in ('height', 'volume'):
                v = res.get(key)
                if v is not None:
                    res[key] = [vi * s for vi, s in zip(v, scales)]
    value = params.get('intensity_value', 'height')

    fits = []
    csp_for_global = {}
    total = len(residues)
    n_fitted = 0
    for i, (name, res) in enumerate(sorted(residues.items())):
        entry = {'residue': name}
        # Raw per-point series (positions + intensities, post-scaling), aligned to
        # metadata 'points'. Lets the viewer compute CSP / I-ratio against ANY
        # reference point and makes the JSON self-contained for exact reopen.
        entry['series'] = {k: json_safe(res.get(k))
                           for k in ('ppm_x', 'ppm_y', 'height', 'volume')}
        if 'csp' in observables:
            csp = csp_series(res, alpha=alpha)
            entry['csp'] = fit_residue_csp(L, csp, P0, n_bootstrap=n_boot)
            if entry['csp'].get('success'):
                csp_for_global[name] = csp
        if 'intensity' in observables:
            ratio = intensity_ratio_series(res, value=value)
            entry['intensity'] = fit_residue_intensity(L, ratio, P0, n_bootstrap=n_boot)
        if any(entry.get(o, {}).get('success') for o in observables):
            n_fitted += 1
        fits.append(entry)
        if progress_callback:
            progress_callback(i + 1, total, name, f"Fitted {name}")

    global_fit = {}
    if 'csp' in observables and len(csp_for_global) >= 2:
        global_fit['csp'] = fit_global_kd_csp(csp_for_global, L, P0)

    data = {
        'metadata': {'analysis': 'Kd_titration', 'protein_conc': P0, 'alpha': alpha,
                     'concentrations': concs, 'points': points,
                     'intensity_scales': list(scales) if scales else None,
                     'intensity_value': value,
                     'observables': observables, 'n_bootstrap': n_boot},
        'fits': fits, 'global': global_fit,
    }

    os.makedirs(params['output_dir'], exist_ok=True)
    prefix = params.get('output_prefix', 'kd')
    json_file = os.path.join(params['output_dir'], f"{prefix}_kd_fit_data.json")
    with open(json_file, 'w') as f:
        json.dump(json_safe(data), f, indent=2)
    results_file = os.path.join(params['output_dir'], f"{prefix}_kd_results.txt")
    _write_results_txt(results_file, fits, global_fit, observables)

    return {'n_fitted': n_fitted, 'n_total': total, 'output_dir': params['output_dir'],
            'json_file': json_file, 'results_file': results_file}


def _write_results_txt(path, fits, global_fit, observables):
    with open(path, 'w') as f:
        if global_fit.get('csp', {}).get('success'):
            f.write(f"# Global shared Kd (CSP): {global_fit['csp']['Kd']:.4g}\n")
        f.write("Residue\tCSP_Kd\tCSP_Kd_err\tCSP_ddmax\tCSP_R2\tInt_Kd\tInt_Kd_err\n")
        for e in fits:
            c = e.get('csp', {})
            it = e.get('intensity', {})
            def g(d, k):
                v = d.get(k)
                return f"{v:.4g}" if isinstance(v, (int, float)) and v is not None else "NaN"
            f.write(f"{e['residue']}\t{g(c,'Kd')}\t{g(c,'Kd_err')}\t{g(c,'dd_max')}\t"
                    f"{g(c,'r_squared')}\t{g(it,'Kd')}\t{g(it,'Kd_err')}\n")
