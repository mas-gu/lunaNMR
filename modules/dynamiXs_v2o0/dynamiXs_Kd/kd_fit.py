# ABOUTME: Kd titration fitter - per-residue and global shared-Kd, CSP and intensity.
# ABOUTME: run_kd_analysis_with_params orchestrates load -> fit -> JSON, dynamiXs-style.

import json
import os

import warnings

import numpy as np
from scipy.optimize import curve_fit, least_squares
from scipy.optimize import OptimizeWarning

from kd_models import csp_model, intensity_decay, residue_sort_key
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


def _singular_covariance(caught):
    """Whether curve_fit reported it could not estimate the covariance.

    Promoted from a log line to a per-fit field: a null Kd_err has several possible
    causes and only this distinguishes a singular covariance from the others. Any
    OptimizeWarning with a different message is re-raised into the log rather than
    swallowed, so a new failure mode still surfaces.
    """
    singular = False
    for w in caught:
        if issubclass(w.category, OptimizeWarning) and 'ovariance' in str(w.message):
            singular = True
        else:
            warnings.warn_explicit(w.message, w.category, w.filename, w.lineno)
    return singular


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


_CSP_TRIM_FRACTION = 0.1    # share of the largest CSPs held out of the spread estimate


def csp_significance(csp_by_name, multiple=1.0):
    """Which residues shift significantly, from the spread of CSP at the last point.

    σ is taken over a TRIMMED population: the largest CSPs are the binders, and letting
    them set the threshold meant to identify them inflates it. On one real dataset that
    is the difference between 19 and 26 residues clearing the bar.

    Judged at the literal last titration point, not each residue's own last measured
    point — a spread is only meaningful over one population measured under one
    condition. A residue with no CSP there cannot be judged at all, so it is reported
    as `unmeasured` rather than counted as insignificant.

    The population is the MEASURED residues only; unmeasured ones are not folded in as
    zeros, which would drag the spread down and admit residues that never shifted.
    """
    values = {}
    for name, csp in csp_by_name.items():
        arr = np.asarray(csp, dtype=float)
        if arr.size and np.isfinite(arr[-1]):
            values[name] = float(arr[-1])
    unmeasured = sorted(set(csp_by_name) - set(values), key=residue_sort_key)

    if len(values) < _MIN_POINTS:
        return {'sigma': None, 'threshold': None, 'multiple': multiple,
                'trim_fraction': _CSP_TRIM_FRACTION, 'n_measured': len(values),
                'significant': sorted(values, key=residue_sort_key),
                'not_significant': [],
                'unmeasured': unmeasured,
                'reason': 'too few residues measured at the last point to judge'}

    ordered = np.sort(np.fromiter(values.values(), dtype=float))
    keep = max(_MIN_POINTS, int(round(len(ordered) * (1.0 - _CSP_TRIM_FRACTION))))
    sigma = float(np.std(ordered[:keep], ddof=1))
    threshold = float(multiple) * sigma
    significant = sorted((n for n, v in values.items() if v > threshold),
                         key=residue_sort_key)
    return {'sigma': sigma, 'threshold': threshold, 'multiple': float(multiple),
            'trim_fraction': _CSP_TRIM_FRACTION, 'n_measured': len(values),
            'significant': significant,
            'not_significant': sorted(set(values) - set(significant),
                                      key=residue_sort_key),
            'unmeasured': unmeasured, 'reason': ''}


def kd_outliers(kd_by_name, z_max=3.0, conc_range=None):
    """Residues whose Kd sits absurdly far from the rest, as name -> robust z-score.

    Median and MAD on log10(Kd), NOT mean and standard deviation. Kd is a ratio spanning
    decades — 25 of them on one real dataset — so "far from typical" is only meaningful
    in log space, and a mean is dragged by the very values it should catch: the largest
    outlier there was 1.9e+05 against a median of 29.8, which inflated σ so much that it
    fell inside its own 2σ band. A robust centre cannot be captured that way.
    """
    usable = {n: float(v) for n, v in kd_by_name.items()
              if isinstance(v, (int, float)) and np.isfinite(v) and v > 0}
    if len(usable) < _MIN_POINTS:
        return {}, {}
    logs = {n: np.log10(v) for n, v in usable.items()}
    centre = float(np.median(list(logs.values())))
    median_kd = float(10 ** centre)
    stats = {'median_log10': centre, 'median_Kd': median_kd, 'z_max': float(z_max),
             'median_credible': True, 'verdict': ''}

    # A robust centre is only a centre if it is itself believable. When the TYPICAL
    # residue fits a Kd the titration could never resolve, the population has no usable
    # middle, and excluding residues for deviating from it throws out the plausible ones
    # — the ones nearest a real Kd sit furthest from a nonsense median. Refuse to gate
    # and say so loudly instead.
    lo, hi = _credible_kd_window(conc_range)
    if lo is not None and not (lo <= median_kd <= hi):
        stats['median_credible'] = False
        stats['credible_window'] = [lo, hi]
        stats['verdict'] = (
            f"THIS FIT MAKES NO SENSE: the typical residue fits Kd = {median_kd:.3g}, "
            f"outside the {lo:.3g}-{hi:.3g} range this titration can resolve. Most "
            f"per-residue fits are meaningless, not just a few outliers, so no Kd from "
            f"this dataset should be reported. Outlier exclusion was skipped — there is "
            f"no trustworthy centre to measure deviation from.")
        return {}, stats

    # Physics before statistics. A residue whose own Kd lies outside the window this
    # titration can resolve was not measured, whatever the population looks like — and a
    # spread-based test cannot catch it, because the same bad fits that produce these
    # values also widen the MAD that is supposed to exclude them. On one real dataset the
    # MAD admitted everything from 0.11 to 71000 while the titration topped out at 150.
    outside = {n: v for n, v in usable.items() if not (lo <= v <= hi)} if lo else {}
    stats['outside_window'] = {n: float(v) for n, v in outside.items()}
    if lo:
        stats['credible_window'] = [lo, hi]

    # z_max governs only the statistical test. The credibility verdict above and the
    # physics window below describe the data, not the gate, so disabling the z-test must
    # never silence them — the user most likely to disable it is the one who most needs
    # to be told their titration cannot answer the question.
    z_enabled = bool(np.isfinite(z_max) and z_max > 0)
    inside = {n: v for n, v in logs.items() if n not in outside}
    if not z_enabled or len(inside) < _MIN_POINTS:
        return dict.fromkeys(outside, float('inf')), stats
    centre = float(np.median(list(inside.values())))
    mad = float(np.median([abs(v - centre) for v in inside.values()]))
    stats['median_log10'] = centre
    stats['median_Kd'] = float(10 ** centre)
    stats['mad_log10'] = mad
    excluded = dict.fromkeys(outside, float('inf'))
    if mad <= 0:
        return excluded, stats
    # 0.6745 puts the MAD on the same scale as a standard deviation for normal data.
    z = {n: 0.6745 * (v - centre) / mad for n, v in inside.items()}
    excluded.update({n: v for n, v in z.items() if abs(v) > z_max})
    return excluded, stats


def _credible_kd_window(conc_range):
    """The Kd range a titration can resolve: one decade beyond its own concentrations.

    Same window fit_global_kd_csp bounds its shared Kd to. A titration spanning tens of
    µM has no power to distinguish Kd = 1e-8 from 1e-3, or 1e5 from 1e8.
    """
    if not conc_range:
        return None, None
    vals = [float(c) for c in conc_range
            if isinstance(c, (int, float)) and np.isfinite(c) and c > 0]
    if not vals:
        return None, None
    return min(vals) / 10.0, max(vals) * 10.0


def _is_dummy(name):
    """Placeholder rows carried in a peak list, excluded from every fitter."""
    return str(name).lower().startswith('dummy')


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
        with warnings.catch_warnings(record=True) as caught:
            # Captured, not silenced. A singular covariance is expected on a degenerate
            # fit and already surfaces as a null Kd_err — but null has more than one
            # cause, so the condition is recorded as a field rather than deleted.
            warnings.simplefilter('always', OptimizeWarning)
            popt, pcov = curve_fit(
                lambda L, dd, kd: csp_model(L, dd, kd, P0), L, y,
                p0=[dd0, kd0], bounds=([0, 0], [np.inf, np.inf]), maxfev=20000)
        singular = _singular_covariance(caught)
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
            'covariance_singular': singular,
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
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always', OptimizeWarning)
            popt, pcov = curve_fit(
                intensity_decay, L, y, p0=[i0, inf0, kd0],
                bounds=([0, 0, 1e-9], [hi, hi, np.inf]), maxfev=20000)
        singular = _singular_covariance(caught)
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
            'covariance_singular': singular,
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


def _resolvable_window(L):
    """The Kd range a titration has power to resolve: one decade below its lowest
    nonzero ligand concentration through one decade above its highest.

    Outside it the model is saturated (or untouched) at every measured point, so
    decades of Kd fit equally well and the formal error stays misleadingly small.
    Returns (None, None) when no nonzero concentration was measured.
    """
    L = np.asarray(L, dtype=float)
    nonzero = L[L > 0]
    if nonzero.size == 0:
        return None, None
    return float(np.min(nonzero)) / 10.0, float(np.max(nonzero)) * 10.0


def _well_determined(kd, kd_err, max_rel_err=0.3):
    """A shared Kd inside the resolvable window can still be unconstrained: the
    dd_max/Kd degeneracy at large Kd, or saturation at small Kd, leaves the value free
    while the fit itself succeeds. That shows up as a large relative error."""
    return (kd_err is not None and np.isfinite(kd_err) and kd > 0
            and kd_err <= max_rel_err * kd)


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
    kd_lo, kd_hi = _resolvable_window(L)
    if kd_lo is None:
        return {'success': False, 'reason': 'no nonzero ligand concentrations'}
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
    reliable = bool(not_pinned and _well_determined(kd, kd_err))
    return {'success': True, 'Kd': kd,
            'Kd_err': kd_err,
            'dd_max': {n: float(res.x[1 + i]) for i, n in enumerate(names)},
            'residues': list(names),
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
    kd = float(res.x[0])
    # Reported as a verdict, not imposed as a bound. The CSP isotherm is a physical
    # model, so clamping its shared Kd to the resolvable decade keeps the optimizer in
    # the region the data speaks to. The intensity decay is phenomenological and its
    # Kd an apparent decay constant, so clipping it would invent a number; saying the
    # titration could not resolve it is the honest report. Same window either way.
    kd_lo, kd_hi = _resolvable_window(L)
    reliable = bool(kd_lo is not None and kd_lo <= kd <= kd_hi
                    and _well_determined(kd, kd_err))
    return {'success': True, 'Kd': kd,
            'Kd_err': kd_err,
            'I0': {n: float(res.x[1 + 2 * i]) for i, n in enumerate(names)},
            'I_inf': {n: float(res.x[2 + 2 * i]) for i, n in enumerate(names)},
            'residues': list(names),
            'n_residues': len(names),
            'reliable': reliable}


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


def load_inputs(params):
    """Load a titration CSV and resolve concentrations, scaling and observable choice.

    Shared by the fit and the survey so both read a dataset the same way.
    """
    points, residues = load_titration(params['input_csv_file'])
    concs = list(params.get('concentrations') or points)
    if len(concs) != len(points):
        raise ValueError(f"{len(concs)} concentrations given but {len(points)} titration points found")
    P0 = float(params['protein_conc'])
    # Equivalents are a ratio to [P]0, so they only become concentrations once scaled.
    # Applied to an explicit --conc too: the flag describes the numbers, not their source.
    if str(params.get('conc_units', 'absolute')).lower().startswith('eq'):
        concs = [c * P0 for c in concs]

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

    return {'points': points, 'residues': residues, 'concs': concs,
            'L': np.asarray(concs, dtype=float),
            'P0': P0,
            'alpha': float(params.get('alpha', 0.14)),
            'observables': params.get('observables', ['csp', 'intensity']),
            'n_bootstrap': int(params.get('n_bootstrap', 0)),
            'value': params.get('intensity_value', 'height'),
            'ref_max_ratio': params.get('ref_max_ratio'),
            'scales': scales}


def run_kd_analysis_with_params(params, progress_callback=None):
    """Load a titration tidy CSV, fit Kd per-residue (CSP/intensity) + global, write JSON."""
    data_in = load_inputs(params)
    points, residues = data_in['points'], data_in['residues']
    concs, L, P0 = data_in['concs'], data_in['L'], data_in['P0']
    alpha, observables = data_in['alpha'], data_in['observables']
    n_boot, value, scales = data_in['n_bootstrap'], data_in['value'], data_in['scales']
    ref_max_ratio = data_in['ref_max_ratio']

    fits = []
    csp_for_global = {}
    intensity_for_global = {}
    # dummy_* rows are placeholders, not residues, and are excluded here to match every
    # other fitter in the package. The count is reported, never silently dropped.
    all_names = set(residues)
    n_excluded_dummy = sum(1 for name in residues if _is_dummy(name))
    residues = {k: v for k, v in residues.items() if not _is_dummy(k)}
    # An explicit selection restricts the fit AND the global pool. Mechanical exclusions
    # still win: selecting a dummy_* row does not resurrect it.
    selection = params.get('residues')
    if selection is not None:
        missing = sorted(set(selection) - all_names, key=residue_sort_key)
        if missing:
            raise ValueError("residue(s) not present in the input: " + ', '.join(missing))
        residues = {k: v for k, v in residues.items() if k in set(selection)}
        if not residues:
            raise ValueError("residue selection is empty — every residue is excluded")
    total = len(residues)
    n_fitted = 0
    # Significance is a property of the whole residue population at one titration point,
    # so it has to be established before any residue is judged against it.
    csp_by_name = ({n: csp_series(r, alpha=alpha) for n, r in residues.items()}
                   if 'csp' in observables else {})
    csp_sig = (csp_significance(csp_by_name,
                                multiple=float(params.get('csp_sigma_multiple', 1.0)))
               if csp_by_name else None)
    csp_significant = set(csp_sig['significant']) if csp_sig else set()

    # Sequence order, not amino-acid letter: a plain sort puts K10 before K3.
    # This sets `fits` order, which drives results.txt, every figure's x-axis
    # and the order of csp_pool_excluded.
    ordered = sorted(residues.items(), key=lambda kv: residue_sort_key(kv[0]))
    for i, (name, res) in enumerate(ordered):
        entry = {'residue': name}
        # Raw per-point series (positions + intensities, post-scaling), aligned to
        # metadata 'points'. Lets the viewer compute CSP / I-ratio against ANY
        # reference point and makes the JSON self-contained for exact reopen.
        entry['series'] = {k: json_safe(res.get(k))
                           for k in ('ppm_x', 'ppm_y', 'height', 'volume')}
        if 'csp' in observables:
            csp = csp_by_name[name]
            entry['csp'] = fit_residue_csp(L, csp, P0, n_bootstrap=n_boot)
            # Reported either way; only a significant shift earns a place in the shared
            # Kd. None means the residue had no CSP at the last point, which is not the
            # same claim as "did not shift".
            entry['csp']['significant'] = (None if name in csp_sig['unmeasured']
                                           else name in csp_significant)
            if _good_for_global(entry['csp']) and name in csp_significant:
                csp_for_global[name] = csp
        if 'intensity' in observables:
            ratio = intensity_ratio_series(res, value=value,
                                           ref_max_ratio=ref_max_ratio)
            entry['intensity'] = fit_residue_intensity(L, ratio, P0, n_bootstrap=n_boot)
            if _good_for_global(entry['intensity']):
                intensity_for_global[name] = ratio
        if any(entry.get(o, {}).get('success') for o in observables):
            n_fitted += 1
        fits.append(entry)
        if progress_callback:
            progress_callback(i + 1, total, name, f"Fitted {name}")

    # Kd outliers leave the pool last: the robust centre must be measured on residues
    # that already passed significance and R², or the values it is meant to catch are
    # part of what defines "typical".
    fit_by_name = {f['residue']: f for f in fits}
    csp_outliers, csp_outlier_stats = kd_outliers(
        {n: fit_by_name[n].get('csp', {}).get('Kd') for n in csp_for_global},
        z_max=float(params.get('kd_outlier_z', 3.0)), conc_range=concs)
    for name in csp_outliers:
        csp_for_global.pop(name, None)

    quality_warnings = []
    if csp_outlier_stats.get('verdict'):
        quality_warnings.append(csp_outlier_stats['verdict'])

    # A residue can now miss the shared CSP fit for several separate reasons. Record which
    # one, per residue: otherwise a surprisingly small pool has its explanation spread
    # across three places and nothing says which gate removed what.
    csp_pool_excluded = {}
    if 'csp' in observables:
        for f in fits:
            name = f['residue']
            if name in csp_for_global:
                continue
            if name in (csp_outlier_stats.get('outside_window') or {}):
                kd_v = csp_outlier_stats['outside_window'][name]
                w = csp_outlier_stats.get('credible_window') or [0, 0]
                why = (f"Kd {kd_v:.3g} outside the {w[0]:.3g}-{w[1]:.3g} range this "
                       f"titration can resolve")
            elif name in csp_outliers:
                why = (f"Kd outlier: robust z {csp_outliers[name]:+.1f} from the "
                       f"population median")
            elif not f.get('csp', {}).get('success'):
                why = 'csp fit failed: ' + f.get('csp', {}).get('reason', '')
            elif f['csp'].get('significant') is None:
                why = 'no CSP at the last titration point'
            elif not f['csp'].get('significant'):
                why = 'CSP below the significance threshold'
            else:
                why = f"R² {f['csp'].get('r_squared')} below {_GLOBAL_R2_MIN}"
            csp_pool_excluded[name] = why

    global_fit = {}
    if 'csp' in observables and len(csp_for_global) >= 2:
        global_fit['csp'] = fit_global_kd_csp(csp_for_global, L, P0)
    if 'intensity' in observables and len(intensity_for_global) >= 2:
        global_fit['intensity'] = fit_global_kd_intensity(intensity_for_global, L)

    data = {
        # `concentrations` here is the L the fit actually used, i.e. post-conversion.
        # load_params accepts this metadata as a params source, so conc_units must say
        # 'absolute' explicitly — inheriting the schema default would be the same value
        # by luck rather than by statement, and re-converting it is the P0-factor bug.
        'metadata': {'analysis': 'Kd_titration', 'protein_conc': P0, 'alpha': alpha,
                     'concentrations': concs, 'conc_units': 'absolute', 'points': points,
                     'intensity_scales': list(scales) if scales else None,
                     'intensity_value': value,
                     'observables': observables, 'n_bootstrap': n_boot,
                     'n_excluded_dummy': n_excluded_dummy,
                     'csp_significance': csp_sig,
                     'csp_pool_excluded': csp_pool_excluded,
                     'quality_warnings': quality_warnings,
                     'csp_kd_outliers': {'excluded': csp_outliers,
                                         **(csp_outlier_stats or {})},
                     # Residues that were selected but yielded no usable fit, named per
                     # observable. A residue whose CSP fits while its intensity does not
                     # is why the two global pools differ in size, so a partial failure
                     # is recorded here too — otherwise that difference has no
                     # explanation anywhere in the output.
                     'unfitted': {f['residue']: [f"{o}: {f.get(o, {}).get('reason', '')}"
                                                 for o in observables
                                                 if not f.get(o, {}).get('success')]
                                  for f in fits
                                  if any(not f.get(o, {}).get('success')
                                         for o in observables)},
                     'n_failed_points': {name: int(res['n_failed_points'])
                                         for name, res in residues.items()
                                         if res.get('n_failed_points')}},
        'fits': fits, 'global': global_fit,
    }

    os.makedirs(params['output_dir'], exist_ok=True)
    prefix = params.get('output_prefix', 'kd')
    # Record the series name so a bundled reopen (generic 'fit_data.json') can recover it.
    data['metadata']['name'] = prefix
    # Machine-readable output one level down; the top level stays figures + report.
    data_dir = os.path.join(params['output_dir'], 'data')
    os.makedirs(data_dir, exist_ok=True)
    json_file = os.path.join(data_dir, f"{prefix}_kd_fit_data.json")
    with open(json_file, 'w') as f:
        json.dump(json_safe(data), f, indent=2)
    results_file = os.path.join(params['output_dir'], f"{prefix}_kd_results.txt")
    _write_results_txt(results_file, fits, global_fit, observables)

    # Importable binding-parameters JSON, so re-running needs no manual re-entry.
    # The run's own params, in the units the caller gave: `concs` has already had
    # equivalents multiplied by P0, so persisting it under the caller's conc_units would
    # convert twice on read-back. normalize_params keeps the schema keys and drops the
    # rest, so a setting added later round-trips without being enumerated here.
    from kd_params import dump_params_json, PARAMS_SUFFIX
    params_file = os.path.join(data_dir, f"{prefix}{PARAMS_SUFFIX}")
    dump_params_json(params_file, {**params, 'points': points})

    return {'n_fitted': n_fitted, 'n_total': total, 'output_dir': params['output_dir'],
            'json_file': json_file, 'results_file': results_file,
            'params_file': params_file, 'n_excluded_dummy': n_excluded_dummy,
            'quality_warnings': quality_warnings}


def _write_results_txt(path, fits, global_fit, observables):
    with open(path, 'w') as f:
        def _gerr(g):
            e = g.get('Kd_err')
            return f" ± {e:.4g}" if isinstance(e, (int, float)) else ""

        def _verdict(g):
            """A shared Kd is not quotable without saying whether the titration could
            resolve it, so the flag travels with the number rather than living only in
            the JSON."""
            return "  [reliable]" if g.get('reliable') else (
                "  [NOT reliable: outside the resolvable range, or relative error > 30%]")
        if global_fit.get('csp', {}).get('success'):
            g = global_fit['csp']
            f.write(f"# Global shared Kd (CSP): {g['Kd']:.4g}{_gerr(g)}{_verdict(g)}\n")
        if global_fit.get('intensity', {}).get('success'):
            g = global_fit['intensity']
            f.write("# Global shared apparent Kd (intensity decay): "
                    f"{g['Kd']:.4g}{_gerr(g)}{_verdict(g)}\n")
        f.write("Residue\tCSP_Kd\tCSP_Kd_err\tCSP_ddmax\tCSP_R2\tInt_Kd\tInt_Kd_err\n")
        for e in fits:
            c = e.get('csp', {})
            it = e.get('intensity', {})
            def g(d, k):
                v = d.get(k)
                return f"{v:.4g}" if isinstance(v, (int, float)) and v is not None else "NaN"
            f.write(f"{e['residue']}\t{g(c,'Kd')}\t{g(c,'Kd_err')}\t{g(c,'dd_max')}\t"
                    f"{g(c,'r_squared')}\t{g(it,'Kd')}\t{g(it,'Kd_err')}\n")
