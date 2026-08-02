# ABOUTME: Kd titration fitter - per-residue and global shared-Kd, CSP and intensity.
# ABOUTME: run_kd_analysis_with_params orchestrates load -> fit -> JSON, dynamiXs-style.

import json
import os

import numpy as np
from scipy.optimize import curve_fit, least_squares

from kd_models import csp_model, intensity_decay
from kd_input import load_titration, csp_series, intensity_ratio_series

_MIN_POINTS = 3
_GLOBAL_R2_MIN = 0.8    # a residue must fit its own decay this well to pool into a global Kd


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


def _good_for_global(fit):
    """A per-residue fit is eligible to pool into a global shared-Kd fit only if it
    succeeded AND fits its own decay well (R² ≥ _GLOBAL_R2_MIN). A poorly-fit residue
    carries no reliable Kd and would otherwise hijack the shared fit."""
    r2 = fit.get('r_squared')
    return bool(fit.get('success')) and isinstance(r2, (int, float)) and r2 >= _GLOBAL_R2_MIN


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
    """Fit one residue's intensity decay: I = I_inf + (I0-I_inf)·exp(-L/Kd).

    I0 (value at L=0), I_inf (plateau at saturation), and Kd (decay constant) are
    all fitted. P0 is unused (kept for signature parity with the CSP fitter).
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


def _shared_kd_err(res, n_params):
    """Standard error on the shared Kd (parameter 0) from a least_squares result, via
    the covariance s²·(JᵀJ)⁻¹ with s² the residual variance. Returns None if singular."""
    J = res.jac
    m = J.shape[0]
    dof = max(m - n_params, 1)
    s2 = 2.0 * res.cost / dof
    try:
        cov = s2 * np.linalg.inv(J.T @ J)
        err = float(np.sqrt(cov[0, 0]))
        return err if np.isfinite(err) else None
    except np.linalg.LinAlgError:
        return None


def fit_global_kd_csp(residues, L, P0):
    """Joint CSP fit with one shared Kd and a per-residue Δδ_max.

    The shared Kd is bounded to one decade below the lowest nonzero ligand
    concentration through one decade above the highest. A CSP titration in the
    tens-hundreds µM range has essentially no power to distinguish Kd=1e-8 from
    Kd=1e-3, or Kd=1e5 from Kd=1e8 — left unbounded, a weak/flat-signal residue set
    can send the optimizer to such a value with a misleadingly small formal error.
    `reliable` is False when the fit is pinned at (or within 1% of) either bound, or its
    relative error (Kd_err/Kd) exceeds 30% — a weak/flat-signal residue set can also
    settle on an interior Kd without hitting either bound, but the dd_max/Kd degeneracy
    at large Kd (or the saturated-at-every-point degeneracy at small Kd) still leaves
    the value itself essentially unconstrained, which shows up as a large relative error.
    """
    L = np.asarray(L, dtype=float)
    # drop residues with no usable points (all-NaN) so np.nanmax can't yield NaN
    names = [n for n in residues if not np.all(np.isnan(np.asarray(residues[n], dtype=float)))]
    if len(names) < 2:
        return {'success': False, 'reason': 'fewer than 2 usable residues'}
    nonzero_L = L[L > 0]
    if nonzero_L.size == 0:
        return {'success': False, 'reason': 'no nonzero ligand concentrations'}
    kd_lo = float(np.min(nonzero_L)) / 10.0
    kd_hi = float(np.max(nonzero_L)) * 10.0
    series = {n: np.asarray(residues[n], dtype=float) for n in names}
    dd0 = [max(float(np.nanmax(series[n])), 1e-6) for n in names]
    kd0 = float(np.clip(max(float(np.nanmax(L)) / 2.0, 1.0), kd_lo, kd_hi))

    def resid(p):
        kd = p[0]
        out = []
        for i, n in enumerate(names):
            y = series[n]
            mask = ~np.isnan(y)
            out.append(csp_model(L[mask], p[1 + i], kd, P0) - y[mask])
        return np.concatenate(out)

    lo = [kd_lo] + [0.0] * len(names)
    hi = [kd_hi] + [np.inf] * len(names)
    try:
        res = least_squares(resid, [kd0] + dd0, bounds=(lo, hi), max_nfev=20000)
    except Exception as e:
        return {'success': False, 'reason': str(e)}
    kd = float(res.x[0])
    kd_err = _shared_kd_err(res, 1 + len(names))
    not_pinned = kd_lo * 1.01 < kd < kd_hi * 0.99
    well_determined = kd_err is not None and np.isfinite(kd_err) and kd_err <= 0.3 * kd
    reliable = bool(not_pinned and well_determined)
    return {'success': True, 'Kd': kd,
            'Kd_err': kd_err,
            'dd_max': {n: float(res.x[1 + i]) for i, n in enumerate(names)},
            'n_residues': len(names),
            'reliable': reliable}


def fit_global_kd_intensity(residues, L, n_bootstrap=0):
    """Joint intensity-decay fit with one shared apparent Kd and per-residue I0/I_inf.

    The intensity model I = I_inf + (I0-I_inf)·exp(-L/Kd) is phenomenological, so the
    shared Kd is an apparent decay constant (concentration units), not a thermodynamic
    dissociation constant like the CSP global Kd.

    Kd_err defaults to the analytic covariance error (s²·(JᵀJ)⁻¹), which assumes the
    loss surface is locally quadratic at the solution. If n_bootstrap is given, that
    assumption is cross-checked by a residual-resampling bootstrap (same convention as
    fit_residue_intensity's own n_bootstrap, scaled to the joint fit: pool every
    residue's residuals, resample with replacement, refit the whole joint problem,
    repeat) — used in place of the analytic error when it succeeds, per the same
    fall-back-on-failure pattern as the per-residue fitters.
    """
    L = np.asarray(L, dtype=float)
    names = [n for n in residues if not np.all(np.isnan(np.asarray(residues[n], dtype=float)))]
    if len(names) < 2:
        return {'success': False, 'reason': 'fewer than 2 usable residues'}
    series = {n: np.asarray(residues[n], dtype=float) for n in names}
    masks = {n: ~np.isnan(series[n]) for n in names}
    kd0 = max(float(np.nanmax(L)) / 2.0, 1.0)
    p0 = [kd0]
    for n in names:                              # per-residue [I0, I_inf] seeds
        y = series[n]
        p0 += [max(float(np.nanmax(y)), 1e-12), max(float(np.nanmin(y)), 0.0)]

    def resid(p):
        kd = p[0]
        out = []
        for i, n in enumerate(names):
            y = series[n]
            mask = masks[n]
            out.append(intensity_decay(L[mask], p[1 + 2 * i], p[2 + 2 * i], kd) - y[mask])
        return np.concatenate(out)

    lo = [1e-9] + [0.0] * (2 * len(names))
    hi = [np.inf] * (1 + 2 * len(names))
    try:
        res = least_squares(resid, p0, bounds=(lo, hi), max_nfev=20000)
    except Exception as e:
        return {'success': False, 'reason': str(e)}

    kd_err = _shared_kd_err(res, 1 + 2 * len(names))
    if n_bootstrap:
        boot_kd_err = _bootstrap_global_intensity(L, series, masks, names, res.x,
                                                   lo, hi, n_bootstrap)
        if np.isfinite(boot_kd_err):
            kd_err = boot_kd_err          # fall back to the analytic value otherwise
    return {'success': True, 'Kd': float(res.x[0]),
            'Kd_err': kd_err,
            'I0': {n: float(res.x[1 + 2 * i]) for i, n in enumerate(names)},
            'I_inf': {n: float(res.x[2 + 2 * i]) for i, n in enumerate(names)},
            'n_residues': len(names)}


def _bootstrap_global_intensity(L, series, masks, names, popt, lo, hi, n):
    """Residual-resampling bootstrap for the joint intensity fit's shared Kd: pool
    every residue's residuals into one set (matching _shared_kd_err's single pooled-s²
    assumption), resample with replacement, refit the whole joint problem, repeat.
    Returns the std of the recovered Kd across iterations, or nan if none converge."""
    try:
        yfit_parts, y_parts, sizes = [], [], []
        for i, nm in enumerate(names):
            mask = masks[nm]
            yfit_parts.append(intensity_decay(L[mask], popt[1 + 2 * i], popt[2 + 2 * i], popt[0]))
            y_parts.append(series[nm][mask])
            sizes.append(int(mask.sum()))
        yfit_flat = np.concatenate(yfit_parts)
        resid_flat = np.concatenate(y_parts) - yfit_flat
        splits = np.cumsum(sizes)[:-1]

        kds = []
        for _ in range(n):
            yb_flat = yfit_flat + np.random.choice(resid_flat, size=len(resid_flat), replace=True)
            yb_parts = np.split(yb_flat, splits)

            def resid_b(p, yb_parts=yb_parts):
                kd = p[0]
                out = []
                for i, nm in enumerate(names):
                    mask = masks[nm]
                    out.append(intensity_decay(L[mask], p[1 + 2 * i], p[2 + 2 * i], kd)
                              - yb_parts[i])
                return np.concatenate(out)

            try:
                resb = least_squares(resid_b, popt, bounds=(lo, hi), max_nfev=20000)
                kds.append(resb.x[0])
            except Exception:
                continue
        return float(np.std(kds)) if kds else float('nan')
    except Exception:
        return float('nan')


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
    intensity_for_global = {}
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
            if _good_for_global(entry['csp']):
                csp_for_global[name] = csp
        if 'intensity' in observables:
            ratio = intensity_ratio_series(res, value=value)
            entry['intensity'] = fit_residue_intensity(L, ratio, P0, n_bootstrap=n_boot)
            if _good_for_global(entry['intensity']):
                intensity_for_global[name] = ratio
        if any(entry.get(o, {}).get('success') for o in observables):
            n_fitted += 1
        fits.append(entry)
        if progress_callback:
            progress_callback(i + 1, total, name, f"Fitted {name}")

    global_fit = {}
    if 'csp' in observables and len(csp_for_global) >= 2:
        global_fit['csp'] = fit_global_kd_csp(csp_for_global, L, P0)
    if 'intensity' in observables and len(intensity_for_global) >= 2:
        global_fit['intensity'] = fit_global_kd_intensity(intensity_for_global, L)

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
    # Record the series name so a bundled reopen (generic 'fit_data.json') can recover it.
    data['metadata']['name'] = prefix
    json_file = os.path.join(params['output_dir'], f"{prefix}_kd_fit_data.json")
    with open(json_file, 'w') as f:
        json.dump(json_safe(data), f, indent=2)
    results_file = os.path.join(params['output_dir'], f"{prefix}_kd_results.txt")
    _write_results_txt(results_file, fits, global_fit, observables)

    # Importable binding-parameters JSON, so re-running needs no manual re-entry.
    from kd_params import dump_params_json, PARAMS_SUFFIX
    params_file = os.path.join(params['output_dir'], f"{prefix}{PARAMS_SUFFIX}")
    dump_params_json(params_file, {
        'points': points, 'concentrations': concs,
        'intensity_scales': list(scales) if scales else None,
        'protein_conc': P0, 'alpha': alpha, 'observables': observables,
        'intensity_value': value, 'n_bootstrap': n_boot})

    return {'n_fitted': n_fitted, 'n_total': total, 'output_dir': params['output_dir'],
            'json_file': json_file, 'results_file': results_file,
            'params_file': params_file}


def _write_results_txt(path, fits, global_fit, observables):
    with open(path, 'w') as f:
        def _gerr(g):
            e = g.get('Kd_err')
            return f" ± {e:.4g}" if isinstance(e, (int, float)) else ""
        if global_fit.get('csp', {}).get('success'):
            g = global_fit['csp']
            f.write(f"# Global shared Kd (CSP): {g['Kd']:.4g}{_gerr(g)}\n")
        if global_fit.get('intensity', {}).get('success'):
            g = global_fit['intensity']
            f.write("# Global shared apparent Kd (intensity decay): "
                    f"{g['Kd']:.4g}{_gerr(g)}\n")
        f.write("Residue\tCSP_Kd\tCSP_Kd_err\tCSP_ddmax\tCSP_R2\tInt_Kd\tInt_Kd_err\n")
        for e in fits:
            c = e.get('csp', {})
            it = e.get('intensity', {})
            def g(d, k):
                v = d.get(k)
                return f"{v:.4g}" if isinstance(v, (int, float)) and v is not None else "NaN"
            f.write(f"{e['residue']}\t{g(c,'Kd')}\t{g(c,'Kd_err')}\t{g(c,'dd_max')}\t"
                    f"{g(c,'r_squared')}\t{g(it,'Kd')}\t{g(it,'Kd_err')}\n")
