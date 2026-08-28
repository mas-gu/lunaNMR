#!/usr/bin/env python3
"""
Enhanced T1/T2 NMR Relaxation Fitting Script with Configurable File Handling

This script performs exponential decay fitting for T1/T2 relaxation data with bootstrap error estimation.
RE VERSION: Enhanced with configurable input/output file management in main function.

Features:
- Exponential decay fitting using lmfit
- Bootstrap error estimation
- Configurable file paths in main function
- Multi-figure PDF output with customizable plots per figure
- Comprehensive results output
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for thread safety
import matplotlib.pyplot as plt
from lmfit import Model
import pandas as pd
import os
import json
import re
from pathlib import Path

try:
    from delay_parser import parse_delay_column, parse_delay_columns, require_delay_start
except ImportError:  # imported as dynamiXs_T1_T2.fit_Tx_NMRRE (parent dir on path)
    from dynamiXs_T1_T2.delay_parser import parse_delay_column, parse_delay_columns, require_delay_start


# A fitted time constant beyond this multiple of the longest delay means the decay
# is not measurable in the sampled window (flat / no-signal / bad peak) — the value
# is meaningless, so the residue is flagged unreliable and left out of the summary.
DEGENERATE_T2_OVER_TMAX = 100.0


def _is_reliable_t2(t2, t2_err, x, amplitude=None, y=None):
    """True when a fitted time constant is measurable within the delay window.

    Two ways to fail. The time constant can exceed the window so far that nothing
    decays within it. Or there can be no signal to decay: a flat or empty trace fits
    an amplitude of zero, at which the time constant is unconstrained and whatever
    the optimiser last held is meaningless. The amplitude check matters because a
    fixed baseline no longer forces such a fit to a visibly absurd time constant.
    """
    if not np.isfinite(t2):
        return False
    scale = 0.0
    if y is not None and np.size(y):
        finite = np.asarray(y, float)[np.isfinite(np.asarray(y, float))]
        scale = float(np.max(np.abs(finite))) if finite.size else 0.0
        if scale <= 0:
            return False
    if amplitude is not None:
        # Compared against the data's own scale: an all-zero trace fits an amplitude of
        # ~1e-12, which is positive but carries no decay.
        if not np.isfinite(amplitude) or amplitude <= 1e-6 * max(scale, 0.0):
            return False
        if scale == 0.0 and amplitude <= 0:
            return False
    x_max = float(np.max(x)) if np.size(x) else 0.0
    return x_max <= 0 or t2 <= DEGENERATE_T2_OVER_TMAX * x_max


# A free per-residue baseline is not identifiable at realistic sampling. At t_max/T ~ 1.7
# it correlates with the time constant at -0.92 to -0.97 and inflates sigma(T) 2-4x, and on
# data whose baseline is genuinely zero it lands negative for 40-95% of residues. It is
# therefore fixed at zero by default; whether a real offset exists is a question about the
# whole series, answered by shared_baseline_test() over all residues at once.
FIT_BASELINE_BY_DEFAULT = False

# Below this ratio of longest delay to fitted time constant the decay is too incomplete
# for the window to pin T down; the fit still returns, but the value is weakly determined.
MARGINAL_TMAX_OVER_T = 2.0


def exp_decay(x, A, t2, C=0.0):
    """Mono-exponential decay with baseline offset: I(t) = A * exp(-t/t2) + C."""
    return A * np.exp(-x / t2) + C


def bootstrap_errors(x, y, model, params, n_bootstrap=1000, fit_baseline=None):
    """
    Classic residual bootstrap error estimation.

    Uses simple residual resampling with replacement - the standard,
    well-understood bootstrap approach for regression. Assumes homoscedastic
    noise (constant variance), which is appropriate for NMR thermal noise.

    Parameters
    ----------
    x : array
        Time points
    y : array
        Signal intensities
    model : lmfit.Model
        Fitting model
    params : lmfit.Parameters
        Initial parameters
    n_bootstrap : int
        Number of bootstrap iterations (default: 1000)

    Returns
    -------
    tuple : (A_error, t2_error)

    Notes
    -----
    Algorithm:
    1. Fit data to get y_fit and residuals
    2. For each bootstrap iteration:
       - Resample residuals WITH REPLACEMENT
       - Create synthetic data: y_synthetic = y_fit + resampled_residuals
       - Refit to get parameter estimates
    3. Return std of bootstrap parameter distributions
    """
    # Initial fit (unweighted)
    result = model.fit(y, params, x=x)
    y_fit = result.best_fit
    residuals = y - y_fit
    n = len(residuals)

    a_values = []
    t2_values = []
    c_values = []

    for _ in range(n_bootstrap):
        # Resample residuals WITH REPLACEMENT
        resampled_indices = np.random.randint(0, n, size=n)
        resampled_residuals = residuals[resampled_indices]

        # Create synthetic data
        y_synthetic = y_fit + resampled_residuals

        # Ensure positive values
        y_synthetic = np.maximum(y_synthetic, 1e-10)

        # Refit (unweighted)
        params_boot = model.make_params(
            A=params['A'].value, t2=params['t2'].value, C=params['C'].value
        )
        params_boot['A'].min = 0
        params_boot['t2'].min = 0
        # make_params always returns a free C; without this the resampled fits would carry
        # a baseline the original fit did not have.
        vary_c = params['C'].vary if fit_baseline is None else bool(fit_baseline)
        params_boot['C'].set(vary=vary_c)

        try:
            res = model.fit(y_synthetic, params_boot, x=x)
            if res.success:
                a_values.append(res.params['A'].value)
                t2_values.append(res.params['t2'].value)
                c_values.append(res.params['C'].value)
        except Exception:
            continue

    # Return standard deviation of bootstrap distribution
    if len(a_values) < 10:
        # Not enough successful fits, fall back to covariance estimate
        a_err = result.params['A'].stderr if result.params['A'].stderr else np.nan
        t2_err = result.params['t2'].stderr if result.params['t2'].stderr else np.nan
        c_err = result.params['C'].stderr if result.params['C'].stderr else np.nan
        return a_err, t2_err, c_err

    return np.std(a_values), np.std(t2_values), np.std(c_values)


def _median_window_ratio(data, t_values):
    """Median t_max/T over residues, for whichever fitted T is passed in."""
    ratios = [float(np.max(x)) / t for (x, _), t in zip(data, t_values)
              if np.isfinite(t) and t > 0]
    return float(np.nanmedian(ratios)) if ratios else float('nan')


def shared_baseline_test(series, f_max=0.45, f_min=-0.15, n_grid=61):
    """Is there a baseline COMMON to every residue in this series?

    A per-residue baseline cannot be measured at realistic sampling, but a shared one can:
    it is a single parameter against hundreds of points. The model is

        I_i(t) = A_i * [ (1 - f) * exp(-t / T_i) + f ]

    with A_i and T_i free per residue and one global plateau fraction f, profiled over f so
    each inner fit stays a well-conditioned two-parameter problem.

    `series` is either {residue: (x, y)} or a sequence of (residue, x, y). Prefer the
    sequence: a dict cannot hold two rows with the same name, and dropping one silently
    changes the answer.

    Two things make this robust, and both were learned the hard way on real data:

    * Each residue contributes its SSR **relative to its own scale**, and the residues are
      combined by **median**, not sum. A summed raw SSR is dominated by whichever rows fit
      worst: on one real series seven `dummy_*` rows held 2.6% of the signal but 44-54% of
      the SSR, and removing them moved the estimate from 0.00 to 0.32.
    * Significance is a **per-residue sign test**, not an F-test on the pooled SSR. With
      hundreds of residues an F-test turns a 1% SSR improvement into p ~ 0, which says
      nothing about whether the effect is real in each residue.

    `dummy_*` rows are excluded, matching every other fitter in the package.

    Read `significant` and, when it is False, `reason` — which names the gate that stopped
    it. Do NOT read `well_determined` as evidence a baseline exists: it asks only whether f
    is resolved above the grid step, and a precisely located NEGATIVE f is precisely no
    baseline. Real 800 T1rho is the case to remember: f = -0.08, eight grid steps from
    zero, so well_determined is True and significant is False.

    `profile_depth` orders candidates WITHIN one series but is NOT comparable ACROSS series:
    it is not scale-free in t_max/T, and truncating a window moves it by as much as the whole
    significant-versus-rejected range, in either direction (6.25% -> 15.71% on one real
    series, 25.81% -> 15.08% on another). `window_ratio_at_f` and
    `window_ratio_at_zero` are both reported so that confound is visible; check the windows
    match before comparing two depths. They differ by 37-43% wherever a baseline is claimed,
    and can straddle the t_max/T < 2 line, so always read the one matching the model you are
    quoting — the fitter's own per-residue `window_ratio` is the C = 0 one. `f` is quantised to
    `f_grid_step`; quote it to two digits.
    """
    from scipy.optimize import curve_fit

    items = (list(series) if not isinstance(series, dict)
             else [(k, x, y) for k, (x, y) in series.items()])
    n_dummy = sum(1 for name, _, _ in items if str(name).lower().startswith('dummy'))
    data = [(np.asarray(x, float), np.asarray(y, float))
            for name, x, y in items
            if not str(name).lower().startswith('dummy') and np.size(x) >= 3]

    base = {'f': 0.0, 'significant': False, 'well_determined': False,
            'n_improved': 0, 'frac_improved': 0.0, 'p_value': 1.0,
            'n_residues': len(data), 'n_excluded_dummy': n_dummy, 'reason': ''}
    if len(data) < 5:
        return dict(base, reason='too few residues')

    def scaled_ssr(x, y, f, want_t=False):
        """This residue's residual sum of squares at trial f, relative to its own scale."""
        shape = lambda t, A, T: A * ((1.0 - f) * np.exp(-t / T) + f)
        try:
            p_, _ = curve_fit(shape, x, y,
                              p0=[float(np.max(y)), max(float(np.max(x)) / 2, 1e-9)],
                              bounds=([0, 1e-9], [np.inf, np.inf]), maxfev=20000)
        except Exception:
            return (np.nan, np.nan) if want_t else np.nan
        scale = float(np.sum(y ** 2))
        ssr = float(np.sum((y - shape(x, *p_)) ** 2)) / scale if scale > 0 else np.nan
        return (ssr, float(p_[1])) if want_t else ssr

    # The grid reaches below zero on purpose: when the true f is zero the estimate must be
    # free to land either side of it, or the null is truncated by construction.
    grid = np.linspace(f_min, f_max, n_grid)
    per_f = np.array([[scaled_ssr(x, y, f) for x, y in data] for f in grid])
    profile = np.nanmedian(per_f, axis=1)
    if not np.any(np.isfinite(profile)):
        return dict(base, reason='no residue could be fitted')
    # f is quantised to the grid, and the floor is reported as `f_grid_step`. The default
    # grid gives +/-0.005, small against a real baseline of ~0.15. Parabolic refinement was
    # tried and rejected: the median profile is not locally parabolic, and on real series it
    # returned no improvement over the grid point. Pass a larger n_grid if you need more.
    # Rounded to the grid's own precision: comparing an unrounded 0.010000000000000009
    # against a rounded step of 0.01 made the size gate miss by one ulp, and the row was
    # then rejected for a different reason than the one that actually applied.
    f_hat = round(float(grid[int(np.nanargmin(profile))]), 12)

    # Per-residue sign test: does this residue's own fit actually improve at f_hat?
    zero = [scaled_ssr(x, y, 0.0, want_t=True) for x, y in data]
    at_zero = np.array([v[0] for v in zero])
    t_zero = np.array([v[1] for v in zero])
    hat = [scaled_ssr(x, y, f_hat, want_t=True) for x, y in data]
    at_hat = np.array([v[0] for v in hat])
    t_hat = np.array([v[1] for v in hat])
    ok = np.isfinite(at_zero) & np.isfinite(at_hat)
    n_ok = int(ok.sum())
    n_improved = int(np.sum(at_hat[ok] < at_zero[ok]))
    frac = n_improved / n_ok if n_ok else 0.0
    try:
        from scipy.stats import binomtest
        p_value = float(binomtest(n_improved, n_ok, 0.5, alternative='greater').pvalue) if n_ok else 1.0
    except Exception:
        p_value = float('nan')

    grid_step = round(float(grid[1] - grid[0]), 12)
    # `well_determined` asks ONE question: is f resolved above the measurement floor?
    # It deliberately excludes frac_improved, which asks whether residues AGREE — a
    # different question, reported separately. Folding the two together described a real
    # series whose f sat 8 grid steps from zero as "not well determined" because 59% of
    # residues improved rather than 60%, which is not what the name says.
    well_determined = bool(abs(f_hat) > grid_step)

    # Every rejection states its own reason, in the order the gates apply. A caller should
    # not have to read the fields in a particular sequence to avoid the wrong conclusion.
    if f_hat < 0:
        reason = ('f is negative: a plateau cannot sit below zero, so this is scatter '
                  'about no baseline, not a baseline')
    elif abs(f_hat) <= grid_step:
        reason = f'f is within one grid step ({grid_step:g}) of zero — too small to call'
    elif frac < 0.6:
        reason = (f'only {100 * frac:.0f}% of residues individually improve; a shared '
                  f'baseline should help most of them')
    elif not (p_value < 0.05):
        reason = f'per-residue sign test not significant (p = {p_value:.3g})'
    else:
        reason = ''
    return dict(base, f=f_hat, n_residues=n_ok, n_improved=n_improved,
                frac_improved=frac, p_value=p_value, well_determined=well_determined,
                significant=(reason == ''), reason=reason,
                f_grid_step=grid_step,
                # Median t_max/T, reported under BOTH models because they disagree by
                # 37-43% wherever a baseline is claimed: a plateau absorbs part of the tail,
                # so T shortens and the ratio rises. The names say which model each is
                # conditioned on, because the two can straddle the t_max/T < 2 line and give
                # opposite readings of whether the window was adequate. The fitter's own
                # per-residue `window_ratio` is the C = 0 one.
                window_ratio_at_f=_median_window_ratio(data, t_hat),
                window_ratio_at_zero=_median_window_ratio(data, t_zero),
                # How much better the profile is at f_hat than at zero. Orders WITHIN a
                # series; not comparable ACROSS series with different windows — it is not
                # scale-free in t_max/T. Read alongside `window_ratio`.
                profile_depth=(float((np.nanmedian(at_zero[ok]) - np.nanmin(profile))
                                     / np.nanmedian(at_zero[ok]))
                               if n_ok and np.nanmedian(at_zero[ok]) > 0 else float('nan')),
                profile_min=float(np.nanmin(profile)),
                profile_at_zero=float(np.nanmedian(at_zero[ok])) if n_ok else float('nan'))


def fit_single_residue(x, y, residue_name, initial_A=None, initial_t2=None,
                       initial_C=None, n_bootstrap=1000, error_method='analytical',
                       fit_baseline=FIT_BASELINE_BY_DEFAULT):
    """
    Fit exponential decay with baseline offset to single residue data.

    Model: I(t) = A * exp(-t/t2) + C

    Parameters
    ----------
    x : array
        Time points
    y : array
        Signal intensities
    residue_name : str
        Residue identifier
    initial_A : float
        Initial amplitude estimate
    initial_t2 : float
        Initial time constant estimate
    initial_C : float or None
        Initial baseline offset estimate. If None, uses min(y) as a
        data-driven default.
    n_bootstrap : int
        Number of bootstrap iterations (only used if error_method='bootstrap')
    error_method : str
        Error estimation method: 'analytical' (fast, from covariance matrix)
        or 'bootstrap' (robust, resampling-based)

    Returns
    -------
    dict : Fitting results
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    # Data-driven initial guesses (fixed defaults like A=5/t2=100 fail to converge
    # when the amplitude or time-scale differs by orders of magnitude, e.g. real
    # intensities ~1e6 over a 0-2 s delay range).
    if initial_C is None:
        initial_C = float(np.min(y))
    if initial_A is None:
        initial_A = float(np.max(y) - np.min(y)) or 1.0
    if initial_t2 is None:
        x_max = float(np.max(x))
        initial_t2 = x_max / 2.0 if x_max > 0 else 1.0

    # Create model and parameters
    model = Model(exp_decay)
    if not fit_baseline:
        initial_C = 0.0
    params = model.make_params(A=initial_A, t2=initial_t2, C=initial_C)
    params['A'].min = 0
    params['t2'].min = 0
    params['C'].set(vary=bool(fit_baseline))

    # Unweighted least squares fit
    result = model.fit(y, params, x=x)

    a = result.params['A'].value
    t2 = result.params['t2'].value
    c = result.params['C'].value

    # Error estimation based on selected method
    if error_method == 'bootstrap':
        a_err, t2_err, c_err = bootstrap_errors(x, y, model, params, n_bootstrap,
                                               fit_baseline=fit_baseline)
    else:
        # Analytical: use covariance matrix from lmfit
        a_err = result.params['A'].stderr if result.params['A'].stderr else np.nan
        t2_err = result.params['t2'].stderr if result.params['t2'].stderr else np.nan
        c_err = result.params['C'].stderr if result.params['C'].stderr else np.nan
    # lmfit reports no stderr both when the covariance is singular and when the residuals
    # vanish, but those mean opposite things: an exact fit has ~zero uncertainty, a
    # singular one has unknown uncertainty. Only the first is safe to call zero, and
    # leaving it NaN drops the residue from every downstream table.
    _resid = np.asarray(getattr(result, 'residual', []), float).ravel()
    _scale = float(np.max(np.abs(np.asarray(y, float)))) if np.size(y) else 0.0
    if _resid.size and _scale > 0 and float(np.sqrt(np.mean(_resid ** 2))) < 1e-9 * _scale:
        a_err = 0.0 if not np.isfinite(a_err) else a_err
        t2_err = 0.0 if not np.isfinite(t2_err) else t2_err
    if not fit_baseline:
        # Fixed, not unmeasured: NaN here would be indistinguishable from a singular
        # covariance, which downstream code treats as a failed error estimate.
        c_err = 0.0

    x_max = float(np.max(x)) if np.size(x) else 0.0
    window_ratio = (x_max / t2) if (np.isfinite(t2) and t2 > 0) else np.nan
    window_marginal = bool(np.isfinite(window_ratio) and window_ratio < MARGINAL_TMAX_OVER_T)

    success = _is_reliable_t2(t2, t2_err, x, amplitude=a, y=y)
    if success:
        print(f"Fitting residue: {residue_name}  ->  {t2:.4g}")
    else:
        print(f"Fitting residue: {residue_name}  ->  unreliable (no decay in window), excluded")

    return {
        'residue': residue_name,
        'A': a,
        't2': t2,
        'C': c,
        'A_err': a_err,
        't2_err': t2_err,
        'C_err': c_err,
        'baseline_fixed': not bool(fit_baseline),
        'window_ratio': window_ratio,
        'window_marginal': window_marginal,
        'x': x,
        'y': y,
        'result': result,
        'success': success
    }


def create_plots(results_list, output_prefix, n_plots_per_figure=20, 
                 experiment_type="T1", time_units="ms", signal_units="Intensity"):
    """
    Create multi-figure PDF plots of fitting results
    
    Parameters:
    -----------
    results_list : list
        List of fitting results dictionaries
    output_prefix : str
        Prefix for output PDF files
    n_plots_per_figure : int
        Number of plots per figure
    experiment_type : str
        Type of experiment (T1, T2, etc.)
    time_units : str
        Units for time axis
    signal_units : str
        Units for signal axis
    """
    n_datasets = len(results_list)
    n_figures = int(np.ceil(n_datasets / n_plots_per_figure))
    
    for fig_idx in range(n_figures):
        fig, axes = plt.subplots(5, 4, figsize=(20, 25), sharex=True)
        axes = axes.flatten()
        start_idx = fig_idx * n_plots_per_figure
        end_idx = min((fig_idx + 1) * n_plots_per_figure, n_datasets)
        
        for idx, result_idx in enumerate(range(start_idx, end_idx)):
            if result_idx >= len(results_list):
                break
                
            result = results_list[result_idx]
            x = result['x']
            y = result['y']
            residue = result['residue']
            a = result['A']
            t2 = result['t2']
            t2_err = result['t2_err']
            
            c = result.get('C', 0.0)
            # Generate fit curve
            x_fit = np.linspace(0, max(x)*1.2, 50)
            y_fit = a * np.exp(-x_fit / t2) + c
            
            ax = axes[idx]
            ax.plot(x, y, 'ko', lw=2, ms=8, label="Data")
            ax.plot(x_fit, y_fit, 'r-', linewidth=2, label="Fit")
            
            # Add fitting results text
            textstr = f"{experiment_type} = {t2:.2f} ± {t2_err:.2f} {time_units}"
            ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
                    fontsize=10, verticalalignment='top', horizontalalignment='left',
                    bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.8))
            
            ax.set_title(f"Residue: {residue}")
            ax.set_xlabel(f"Time ({time_units})")
            ax.set_ylabel(f"Signal ({signal_units})")
            ax.legend()
            ax.set_xlim(0, max(x)*1.1)
            ax.set_ylim(min(y)*0.9, max(y)*1.1)
        
        # Hide unused subplots
        for idx in range(end_idx - start_idx, len(axes)):
            axes[idx].set_visible(False)
        
        plt.tight_layout()
        fig.savefig(f"{output_prefix}_fit_results_fig{fig_idx + 1}.pdf", format='pdf', dpi=300)
        # The figure is already saved; showing it only makes sense on an interactive
        # backend. Under Agg (CLI, tests) plt.show() warns and does nothing.
        if matplotlib.get_backend().lower() not in ('agg', 'pdf', 'ps', 'svg', 'template'):
            plt.show()
        plt.close(fig)


def save_results(results_list, output_file, experiment_type="T1"):
    """
    Save fitting results to text file

    Parameters:
    -----------
    results_list : list
        List of fitting results dictionaries
    output_file : str
        Output filename
    experiment_type : str
        Type of experiment for headers
    """
    with open(output_file, "w") as f:
        f.write(
            f"Residue\tA\t{experiment_type}\tC\tA_err\t{experiment_type}_err"
            f"\tC_err\tSuccess\tWindowRatio\n"
        )
        for result in results_list:
            # fitting_wrapper.py reads row['Success']; the multicore writer has always
            # emitted it and this one did not, so the single-core path raised KeyError.
            ok = 'Yes' if result.get('reliable', result.get('success', True)) else 'No'
            wr = result.get('window_ratio', float('nan'))
            f.write(
                f"{result['residue']}\t{result['A']:.6e}\t{result['t2']:.6e}\t"
                f"{result['C']:.6e}\t{result['A_err']:.6e}\t"
                f"{result['t2_err']:.6e}\t{result['C_err']:.6e}\t{ok}\t{wr:.3f}\n"
            )


def save_fit_data_json(results_list, output_file, experiment_type, time_units,
                       signal_units, n_bootstrap, field_freq):
    """
    Save complete fitting data as JSON for interactive visualization

    Parameters:
    -----------
    results_list : list
        List of fitting results dictionaries
    output_file : str
        Output JSON filename
    experiment_type : str
        Type of experiment (T1, T2)
    time_units : str
        Units for time axis
    signal_units : str
        Units for signal axis
    n_bootstrap : int
        Number of bootstrap iterations used
    field_freq : float
        Magnetic field frequency in MHz
    """
    # Get time points from first result (same for all residues)
    if not results_list:
        raise ValueError("No results to save")

    time_points = results_list[0]['x'].tolist()

    # Build fit curve points (dense sampling for smooth plotting)
    max_time = max(time_points)
    fit_time_dense = np.linspace(0, max_time * 1.2, 100)

    fits_data = []
    for result in results_list:
        c = result.get('C', 0.0)
        c_err = result.get('C_err', np.nan)
        # Calculate fit curve using fitted parameters
        fit_intensity = result['A'] * np.exp(-fit_time_dense / result['t2']) + c

        fits_data.append({
            'residue': str(result['residue']),
            'A': float(result['A']),
            't2': float(result['t2']),
            'C': float(c),
            'A_err': float(result['A_err']),
            't2_err': float(result['t2_err']),
            'C_err': float(c_err),
            # Fit provenance: whether the baseline was fixed, and whether the delay window
            # was long enough to pin the time constant down. Computed by the fitter and
            # previously dropped here, which made them unreachable from the CLI.
            'baseline_fixed': bool(result.get('baseline_fixed', True)),
            'window_ratio': float(result.get('window_ratio', float('nan'))),
            'window_marginal': bool(result.get('window_marginal', False)),
            'intensities': [float(val) for val in result['y']],
            'fit_curve': {
                'time': [float(t) for t in fit_time_dense],
                'intensity': [float(i) for i in fit_intensity]
            }
        })

    output_data = {
        'metadata': {
            'experiment_type': experiment_type,
            'field_freq': float(field_freq),
            'time_units': time_units,
            'signal_units': signal_units,
            'n_bootstrap': n_bootstrap,
            'n_residues': len(results_list),
            'time_points': [float(t) for t in time_points]
        },
        'fits': fits_data
    }

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)

    print(f"Saved JSON fit data to: {output_file}")


def main():
    """
    Main function with configurable file handling
    """
    
    # ========================
    # CONFIGURATION SECTION
    # ========================
    
    # Input file configuration
    input_csv_file = "T1_data_600.csv"  # Input CSV file path
    
    # Output file configuration
    output_prefix = "600_T1"  # Prefix for output files
    results_txt_file = "600_T1_fit_results.txt"  # Results text file
    
    # Experiment configuration
    experiment_type = "T1"  # T1, T2, etc. (for labels and headers)
    time_units = "ms"  # Time axis units
    signal_units = "Intensity"  # Signal axis units
    
    # Fitting parameters
    initial_A = None  # Data-driven (max-min) unless overridden
    initial_t2 = None  # Data-driven (from delay range) unless overridden
    n_bootstrap = 1000  # Number of bootstrap iterations for error estimation
    error_method = 'analytical'  # 'analytical' (fast) or 'bootstrap' (robust)

    # Plot configuration
    n_plots_per_figure = 20  # Number of plots per PDF page

    # Optional: Time scaling factor (uncomment if needed)
    # time_scaling_factor = 1000  # Multiply time values by this factor

    # ========================
    # END CONFIGURATION
    # ========================

    # Validate input file exists
    if not os.path.exists(input_csv_file):
        raise FileNotFoundError(f"Input file not found: {input_csv_file}")

    print(f"Starting {experiment_type} relaxation fitting analysis...")
    print(f"Input file: {input_csv_file}")
    print(f"Output prefix: {output_prefix}")
    print(f"Error method: {error_method}")
    if error_method == 'bootstrap':
        print(f"Bootstrap iterations: {n_bootstrap}")

    # Load CSV using pandas to handle mixed data types
    print("Loading data...")
    raw_df = pd.read_csv(input_csv_file, header=None)

    # Detect format: LunaNMR Fit Series vs simple DynamiXs format
    header_row = raw_df.iloc[0].astype(str).tolist()
    lunaNMR_columns = {'Peak_Number', 'Assignment', 'Reference_X', 'Reference_Y'}
    # Columns whose header carries no delay: reported so a spectrum cannot go
    # missing from a series without anyone noticing.
    dropped_columns = []
    detected_lunaNMR_cols = [col for col in header_row if col in lunaNMR_columns]

    if detected_lunaNMR_cols:
        print(f"Detected LunaNMR Fit Series format (columns: {detected_lunaNMR_cols})")
        delay_start_idx = require_delay_start(header_row, lunaNMR_columns, input_csv_file)
        assignment_idx = header_row.index('Assignment') if 'Assignment' in header_row else 1
        print(f"  Using 'Assignment' column (index {assignment_idx}) for residue names")
        print(f"  Delay columns start at index {delay_start_idx}")
        residue_names = raw_df.iloc[1:, assignment_idx].to_numpy()
        # Parse delay columns - handles duplicates like "300_2" -> 300.0
        delay_headers = raw_df.iloc[0, delay_start_idx:].tolist()
        x, _keep, _dropped = parse_delay_columns(delay_headers)
        dropped_columns.extend(_dropped)
        y_data = raw_df.iloc[1:, delay_start_idx:].astype(float).to_numpy()[:, _keep]
    else:
        print("Detected simple DynamiXs format")
        residue_names = raw_df.iloc[1:, 0].to_numpy()  # Column 1 (residue names), skip header
        # Parse delay columns - handles duplicates like "300_2" -> 300.0
        delay_headers = raw_df.iloc[0, 1:].tolist()
        x, _keep, _dropped = parse_delay_columns(delay_headers)
        dropped_columns.extend(_dropped)
        y_data = raw_df.iloc[1:, 1:].astype(float).to_numpy()[:, _keep]  # Values to fit (rows 2+, columns 2+)

    # Filter out dummy residues (case-insensitive)
    valid_mask = np.array([not str(name).lower().startswith('dummy') for name in residue_names])
    n_dummies = np.sum(~valid_mask)
    if n_dummies > 0:
        print(f"Excluding {n_dummies} dummy residue(s) from fitting")
        residue_names = residue_names[valid_mask]
        y_data = y_data[valid_mask]

    print(f"Loaded {len(residue_names)} residues with {len(x)} time points")
    print(f"Time range: {x.min():.3f} to {x.max():.3f} {time_units}")

    # Fit all residues
    print("Starting fitting process...")
    results_list = []

    for col_idx in range(len(residue_names)):
        y = y_data[col_idx, :]
        residue = residue_names[col_idx]

        result = fit_single_residue(x, y, residue, initial_A, initial_t2,
                                    n_bootstrap=n_bootstrap, error_method=error_method)
        results_list.append(result)

    print(f"Completed fitting for {len(results_list)} residues")
    
    # Create plots
    print("Generating plots...")
    create_plots(results_list, output_prefix, n_plots_per_figure, 
                 experiment_type, time_units, signal_units)
    
    # Save results
    print(f"Saving results to {results_txt_file}...")
    save_results(results_list, results_txt_file, experiment_type)
    
    # Summary statistics
    t2_values = [r['t2'] for r in results_list]
    t2_errors = [r['t2_err'] for r in results_list]
    
    print(f"\n{experiment_type} Analysis Summary:")
    print(f"  Number of residues fitted: {len(results_list)}")
    print(f"  {experiment_type} range: {min(t2_values):.2f} to {max(t2_values):.2f} {time_units}")
    print(f"  Mean {experiment_type}: {np.mean(t2_values):.2f} ± {np.std(t2_values):.2f} {time_units}")
    print(f"  Mean fitting error: {np.mean(t2_errors):.2f} {time_units}")
    print(f"  Results saved to: {results_txt_file}")
    print(f"  Plots saved with prefix: {output_prefix}")
    
    print("Analysis completed successfully!")


def run_analysis_with_params(params, progress_callback=None):
    """
    Run analysis with parameters provided by GUI

    Parameters:
    -----------
    params : dict
        Dictionary containing all analysis parameters
    progress_callback : callable, optional
        Function(completed, total, residue_name, message) called after each residue

    Returns:
    --------
    dict : Results summary
    """

    # Callers may choose; the CLI never does, because C = 0 is the convention.
    fit_baseline = bool(params.get('fit_baseline', FIT_BASELINE_BY_DEFAULT))
    
    # Input file configuration
    input_csv_file = params['input_csv_file']
    
    # Output file configuration
    output_prefix = params['output_prefix']
    results_txt_file = params['results_txt_file']
    
    # Experiment configuration
    experiment_type = params['experiment_type']
    time_units = params.get('time_units', 'ms')
    signal_units = params.get('signal_units', 'Intensity')
    
    # Fitting parameters
    initial_A = params.get('initial_A')   # None -> data-driven per residue
    initial_t2 = params.get('initial_t2')  # None -> data-driven per residue
    n_bootstrap = params.get('n_bootstrap', 1000)
    error_method = params.get('error_method', 'analytical')

    # Plot configuration
    n_plots_per_figure = params.get('n_plots_per_figure', 20)

    # Validate input file exists
    if not os.path.exists(input_csv_file):
        raise FileNotFoundError(f"Input file not found: {input_csv_file}")

    print(f"Starting {experiment_type} relaxation fitting analysis...")
    print(f"Input file: {input_csv_file}")
    print(f"Output prefix: {output_prefix}")
    print(f"Error method: {error_method}")
    if error_method == 'bootstrap':
        print(f"Bootstrap iterations: {n_bootstrap}")

    # Load CSV using pandas to handle mixed data types
    print("Loading data...")
    raw_df = pd.read_csv(input_csv_file, header=None)

    # Detect format: LunaNMR Fit Series vs simple DynamiXs format
    header_row = raw_df.iloc[0].astype(str).tolist()
    lunaNMR_columns = {'Peak_Number', 'Assignment', 'Reference_X', 'Reference_Y'}
    # Columns whose header carries no delay: reported so a spectrum cannot go
    # missing from a series without anyone noticing.
    dropped_columns = []
    detected_lunaNMR_cols = [col for col in header_row if col in lunaNMR_columns]

    if detected_lunaNMR_cols:
        print(f"Detected LunaNMR Fit Series format (columns: {detected_lunaNMR_cols})")
        delay_start_idx = require_delay_start(header_row, lunaNMR_columns, input_csv_file)
        assignment_idx = header_row.index('Assignment') if 'Assignment' in header_row else 1
        print(f"  Using 'Assignment' column (index {assignment_idx}) for residue names")
        print(f"  Delay columns start at index {delay_start_idx}")
        residue_names = raw_df.iloc[1:, assignment_idx].to_numpy()
        # Parse delay columns - handles duplicates like "300_2" -> 300.0
        delay_headers = raw_df.iloc[0, delay_start_idx:].tolist()
        x, _keep, _dropped = parse_delay_columns(delay_headers)
        dropped_columns.extend(_dropped)
        y_data = raw_df.iloc[1:, delay_start_idx:].astype(float).to_numpy()[:, _keep]
    else:
        print("Detected simple DynamiXs format")
        residue_names = raw_df.iloc[1:, 0].to_numpy()  # Column 1 (residue names), skip header
        # Parse delay columns - handles duplicates like "300_2" -> 300.0
        delay_headers = raw_df.iloc[0, 1:].tolist()
        x, _keep, _dropped = parse_delay_columns(delay_headers)
        dropped_columns.extend(_dropped)
        y_data = raw_df.iloc[1:, 1:].astype(float).to_numpy()[:, _keep]  # Values to fit (rows 2+, columns 2+)

    # Filter out dummy residues (case-insensitive)
    valid_mask = np.array([not str(name).lower().startswith('dummy') for name in residue_names])
    n_dummies = np.sum(~valid_mask)
    if n_dummies > 0:
        print(f"Excluding {n_dummies} dummy residue(s) from fitting")
        residue_names = residue_names[valid_mask]
        y_data = y_data[valid_mask]

    print(f"Loaded {len(residue_names)} residues with {len(x)} time points")
    print(f"Time range: {x.min():.3f} to {x.max():.3f} {time_units}")

    # Fit all residues
    print("Starting fitting process...")
    results_list = []
    total_residues = len(residue_names)

    for col_idx in range(total_residues):
        y = y_data[col_idx, :]
        residue = residue_names[col_idx]

        result = fit_single_residue(x, y, residue, initial_A, initial_t2,
                                    n_bootstrap=n_bootstrap, error_method=error_method)
        results_list.append(result)

        # Report progress
        completed = col_idx + 1
        if progress_callback:
            progress_callback(completed, total_residues, residue,
                              f"Fitted {residue} ({completed}/{total_residues})")

    print(f"Completed fitting for {len(results_list)} residues")
    
    # Create plots
    print("Generating plots...")
    create_plots(results_list, output_prefix, n_plots_per_figure, 
                 experiment_type, time_units, signal_units)
    
    # Save results
    print(f"Saving results to {results_txt_file}...")
    save_results(results_list, results_txt_file, experiment_type)

    # Save JSON data for interactive visualization
    json_folder = params.get('json_folder')
    json_file = None
    if json_folder:
        # Create json folder if it doesn't exist
        os.makedirs(json_folder, exist_ok=True)

        # Construct JSON filename: field{N}_{T1|T2}_fit_data.json
        field_name = params.get('field_name', 'field1')
        json_filename = f"{field_name}_{experiment_type}_fit_data.json"
        json_file = os.path.join(json_folder, json_filename)

        print(f"Saving JSON fit data to {json_file}...")
        save_fit_data_json(
            results_list=results_list,
            output_file=json_file,
            experiment_type=experiment_type,
            time_units=time_units,
            signal_units=signal_units,
            n_bootstrap=n_bootstrap,
            field_freq=params.get('field_freq', 600.0)
        )

    # Summary statistics over reliable residues only (degenerate/no-decay excluded)
    reliable = [r for r in results_list if r.get('success', True)]
    n_excluded = len(results_list) - len(reliable)
    t2_values = [r['t2'] for r in reliable]
    t2_errors = [r['t2_err'] for r in reliable]

    print(f"\n{experiment_type} Analysis Summary:")
    print(f"  Number of residues fitted: {len(reliable)}")
    if n_excluded:
        print(f"  Excluded (no measurable decay): {n_excluded}")
    if t2_values:
        print(f"  {experiment_type} range: {min(t2_values):.2f} to {max(t2_values):.2f} {time_units}")
        print(f"  Mean {experiment_type}: {np.mean(t2_values):.2f} ± {np.std(t2_values):.2f} {time_units}")
        print(f"  Mean fitting error: {np.nanmean(t2_errors):.2f} {time_units}")
    else:
        print("  No residues with a measurable decay.")
    print(f"  Results saved to: {results_txt_file}")
    print(f"  Plots saved with prefix: {output_prefix}")
    if json_file:
        print(f"  JSON data saved to: {json_file}")

    print("Analysis completed successfully!")

    # Return results summary
    return {
        'n_fitted': len(reliable),
        'dropped_columns': dropped_columns,
        'n_excluded': n_excluded,
        'results_file': results_txt_file,
        'plots_prefix': output_prefix,
        'json_file': json_file,
        't2_range': (min(t2_values), max(t2_values)) if t2_values else (np.nan, np.nan),
        'mean_t2': np.mean(t2_values) if t2_values else np.nan,
        'std_t2': np.std(t2_values) if t2_values else np.nan,
        'mean_error': np.nanmean(t2_errors) if t2_errors else np.nan,
        # Fit provenance, so a caller sees the convention and the sampling quality
        # without opening the per-residue files.
        'baseline_fixed': not FIT_BASELINE_BY_DEFAULT,
        'n_window_marginal': sum(1 for r in results_list
                                 if isinstance(r, dict) and r.get('window_marginal')),
    }


if __name__ == "__main__":
    main()