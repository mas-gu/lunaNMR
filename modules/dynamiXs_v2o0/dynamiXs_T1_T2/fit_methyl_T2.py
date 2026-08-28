# ABOUTME: Bi-exponential T2 fitting for methyl groups: I(t) = 0.5*A*[exp(-t/T2a) + exp(-t/T2b)].
# ABOUTME: Fitted in (A, T2_avg, dT2) basis to keep T2_avg well-determined when bi-exp is unidentifiable.

import json
import os
import re

import numpy as np
import pandas as pd
from lmfit import Model

# Initial-guess defaults; the GUI exposes T2 starting values as user-overridable.
DEFAULT_INITIAL_T2A = 100.0  # ms (slow component)
DEFAULT_INITIAL_T2B = 20.0   # ms (fast component)

# A bi-exp fit is flagged unidentifiable when |dT2|/sigma(dT2) is below this many
# standard deviations — i.e. dT2 cannot be distinguished from zero, so the data
# is effectively mono-exponential and the (T2a, T2b) decomposition is meaningless.
UNIDENTIFIABLE_DT2_SIGMA = 3.0

T2_LOWER_BOUND_MS = 0.1


try:
    from delay_parser import parse_delay_column, parse_delay_columns, require_delay_start
except ImportError:  # imported as dynamiXs_T1_T2.fit_methyl_T2 (parent dir on path)
    from dynamiXs_T1_T2.delay_parser import parse_delay_column, parse_delay_columns, require_delay_start


def biexp_decay(x, A, t2_a, t2_b):
    """Methyl bi-exponential T2 decay (shared amplitude, equal-weight 0.5).

    I(t) = 0.5 * A * (exp(-t/t2_a) + exp(-t/t2_b))
    """
    return 0.5 * A * (np.exp(-x / t2_a) + np.exp(-x / t2_b))


def biexp_decay_avg(x, A, T2_avg, dT2):
    """Equivalent bi-exp model in (T2_avg, dT2) coordinates.

    t2_a = T2_avg + dT2/2  (slow component, dT2 >= 0 enforces ordering)
    t2_b = T2_avg - dT2/2  (fast component)
    I(t) = 0.5 * A * (exp(-t/t2_a) + exp(-t/t2_b))

    The reparameterization diagonalizes the Hessian: T2_avg is the well-determined
    direction, dT2 is the soft direction that becomes unidentifiable as dT2 -> 0.
    """
    t2_a = T2_avg + dT2 / 2.0
    t2_b = T2_avg - dT2 / 2.0
    if t2_b < T2_LOWER_BOUND_MS:
        t2_b = T2_LOWER_BOUND_MS
    if t2_a < T2_LOWER_BOUND_MS:
        t2_a = T2_LOWER_BOUND_MS
    return 0.5 * A * (np.exp(-x / t2_a) + np.exp(-x / t2_b))


def _derive_secondary_errors(T2_avg, dT2, covar, var_names):
    """Compute (t2_a_err, t2_b_err, eta_err) from the (A, T2_avg, dT2) covariance.

    Change of variables:
        t2_a = T2_avg + dT2/2,  t2_b = T2_avg - dT2/2
        eta  = t2_a / t2_b

    Applies the standard sigma^2 = J * Cov * J^T propagation. Returns NaN entries
    when the covariance is missing (e.g. lmfit could not estimate it).
    """
    if covar is None or var_names is None:
        return float("nan"), float("nan"), float("nan")

    try:
        i_avg = list(var_names).index("T2_avg")
        i_dT2 = list(var_names).index("dT2")
    except ValueError:
        return float("nan"), float("nan"), float("nan")

    var_avg = covar[i_avg, i_avg]
    var_dT2 = covar[i_dT2, i_dT2]
    cov_avg_dT2 = covar[i_avg, i_dT2]

    var_t2_a = var_avg + var_dT2 / 4.0 + cov_avg_dT2
    var_t2_b = var_avg + var_dT2 / 4.0 - cov_avg_dT2
    t2_a_err = float(np.sqrt(max(var_t2_a, 0.0)))
    t2_b_err = float(np.sqrt(max(var_t2_b, 0.0)))

    t2_b = T2_avg - dT2 / 2.0
    if abs(t2_b) < 1e-12:
        return t2_a_err, t2_b_err, float("nan")

    inv_t2b_sq = 1.0 / (t2_b * t2_b)
    deta_dT2avg = -dT2 * inv_t2b_sq
    deta_ddT2 = T2_avg * inv_t2b_sq
    var_eta = (
        deta_dT2avg ** 2 * var_avg
        + deta_ddT2 ** 2 * var_dT2
        + 2.0 * deta_dT2avg * deta_ddT2 * cov_avg_dT2
    )
    eta_err = float(np.sqrt(max(var_eta, 0.0)))

    return t2_a_err, t2_b_err, eta_err


def _mono_exp_refit(x, y, A_init, T2_init, A_max, t2_max):
    """Fit a single-component model `A * exp(-t/T2)` to (x, y).

    Used when the bi-exp's t2_b hits the lower bound: in that regime the data
    only constrains the slow component, and the bi-exp's reported A is twice
    the visible amplitude (the 'fast component' carries half but decays before
    the first measurement — a fitting artifact). Refitting as mono-exp returns
    physically-meaningful (A, T2) and their stderrs.

    Returns (A, T2, A_err, T2_err).
    """
    def _mono(x, A, T2):
        return A * np.exp(-x / T2)

    model = Model(_mono)
    params = model.make_params(A=A_init, T2=T2_init)
    params["A"].min, params["A"].max = 0.0, A_max
    params["T2"].min, params["T2"].max = T2_LOWER_BOUND_MS, t2_max
    result = model.fit(y, params, x=x)

    def _se(name):
        v = result.params[name].stderr
        return float(v) if v is not None else float("nan")

    return (
        float(result.params["A"].value),
        float(result.params["T2"].value),
        _se("A"),
        _se("T2"),
    )


def _is_unidentifiable(dT2, dT2_err):
    """True when the data does not significantly distinguish two components.

    Flag is set when dT2 cannot be told apart from zero at ~UNIDENTIFIABLE_DT2_SIGMA
    standard deviations (the bound-active or vanishing-curvature regime).
    """
    if not np.isfinite(dT2_err):
        return True
    if dT2_err <= 0:
        return False
    return bool(dT2 < UNIDENTIFIABLE_DT2_SIGMA * dT2_err)


def bootstrap_errors_methyl(x, y, model, params, n_bootstrap=1000):
    """Residual-resampling bootstrap for the (A, T2_avg, dT2) fit.

    Returns a dict with stderrs and a sample-derived covariance matrix. The
    covariance is reused by `_derive_secondary_errors` to give consistent
    t2_a/t2_b/eta errors. Falls back to lmfit analytical results when fewer
    than 10 bootstrap fits succeed.
    """
    result = model.fit(y, params, x=x)
    y_fit = result.best_fit
    residuals = y - y_fit
    n = len(residuals)

    A_vals, T2_avg_vals, dT2_vals = [], [], []

    for _ in range(n_bootstrap):
        idx = np.random.randint(0, n, size=n)
        y_synth = y_fit + residuals[idx]
        y_synth = np.maximum(y_synth, 1e-10)

        params_boot = params.copy()
        try:
            res = model.fit(y_synth, params_boot, x=x)
            if res.success:
                A_vals.append(res.params["A"].value)
                T2_avg_vals.append(res.params["T2_avg"].value)
                dT2_vals.append(res.params["dT2"].value)
        except Exception:
            continue

    if len(A_vals) < 10:
        def _se(name):
            v = result.params[name].stderr
            return float(v) if v is not None else float("nan")
        return {
            "A_err": _se("A"),
            "T2_avg_err": _se("T2_avg"),
            "dT2_err": _se("dT2"),
            "covar": result.covar,
            "var_names": result.var_names,
        }

    A_arr = np.asarray(A_vals)
    T2_avg_arr = np.asarray(T2_avg_vals)
    dT2_arr = np.asarray(dT2_vals)
    sample_covar = np.cov(np.vstack([A_arr, T2_avg_arr, dT2_arr]))
    return {
        "A_err": float(np.std(A_arr)),
        "T2_avg_err": float(np.std(T2_avg_arr)),
        "dT2_err": float(np.std(dT2_arr)),
        "covar": sample_covar,
        "var_names": ["A", "T2_avg", "dT2"],
    }


def fit_single_residue_methyl(x, y, residue_name,
                              initial_A=None,
                              initial_t2_a=DEFAULT_INITIAL_T2A,
                              initial_t2_b=DEFAULT_INITIAL_T2B,
                              n_bootstrap=1000, error_method="analytical"):
    """Fit shared-amplitude bi-exponential T2 decay in (A, T2_avg, dT2) basis.

    Model: I(t) = 0.5 * A * (exp(-t/t2_a) + exp(-t/t2_b)) with
    t2_a = T2_avg + dT2/2, t2_b = T2_avg - dT2/2. dT2 >= 0 enforces slow-first
    ordering automatically (no post-fit swap).

    Returns t2_a, t2_b, eta as derived quantities alongside the natively-fit
    T2_avg and dT2. Errors on the derived quantities use the full covariance
    of (A, T2_avg, dT2) via change of variables.

    `bi_exp_unidentifiable` is True when dT2 cannot be distinguished from zero
    (mono-exp regime) — in that case T2_avg is well-determined but T2a, T2b
    individually are not.
    """
    print(f"Fitting methyl residue: {residue_name}")

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    span = float(np.max(y) - np.min(y))
    if initial_A is None:
        initial_A = span if span > 0 else float(np.max(y))

    y_max = float(np.max(y))
    x_max = float(np.max(x)) if len(x) else 1.0
    t2_max = max(5.0 * x_max, 1.0)
    A_max = max(5.0 * y_max, 1.0)

    initial_T2_avg = 0.5 * (float(initial_t2_a) + float(initial_t2_b))
    initial_dT2 = abs(float(initial_t2_a) - float(initial_t2_b))

    initial_T2_avg = float(np.clip(initial_T2_avg, T2_LOWER_BOUND_MS, t2_max))
    initial_dT2 = float(np.clip(initial_dT2, 0.0, t2_max))

    model = Model(biexp_decay_avg)
    params = model.make_params(A=initial_A, T2_avg=initial_T2_avg, dT2=initial_dT2)
    params["A"].min, params["A"].max = 0.0, A_max
    params["T2_avg"].min, params["T2_avg"].max = T2_LOWER_BOUND_MS, t2_max
    params["dT2"].min, params["dT2"].max = 0.0, t2_max

    result = model.fit(y, params, x=x)

    A = result.params["A"].value
    T2_avg = result.params["T2_avg"].value
    dT2 = result.params["dT2"].value

    if error_method == "bootstrap":
        # Warm-start each bootstrap iteration from the converged solution, not
        # the initial guess — every resample is fitting data that's only a
        # residual-shuffle away from the converged y_fit, so the global minimum
        # is right next to result.params. Starting from the initial guess wastes
        # iterations and inflates the failure rate.
        boot = bootstrap_errors_methyl(x, y, model, result.params, n_bootstrap)
        A_err = boot["A_err"]
        T2_avg_err = boot["T2_avg_err"]
        dT2_err = boot["dT2_err"]
        covar = boot["covar"]
        var_names = boot["var_names"]
    else:
        def _se(name):
            v = result.params[name].stderr
            return float(v) if v is not None else float("nan")
        A_err = _se("A")
        T2_avg_err = _se("T2_avg")
        dT2_err = _se("dT2")
        covar = result.covar
        var_names = result.var_names

    t2_a = T2_avg + dT2 / 2.0
    t2_b = T2_avg - dT2 / 2.0

    t2_b_bound_active = t2_b < T2_LOWER_BOUND_MS
    if t2_b_bound_active:
        # The optimizer pushed ΔT2 past 2·T2_avg; biexp_decay_avg internally
        # clipped t2_b at T2_LOWER_BOUND_MS, but lmfit's reported parameters
        # carry the unclipped values. The bi-exp model is over-parameterized
        # for this data — its A is twice the visible amplitude (the "fast
        # component" decayed before the first measurement). Refit as a clean
        # mono-exp to recover physically-meaningful (A, T2) and stderrs.
        slow_T2 = T2_avg + dT2 / 2.0
        A_mono, T2_mono, A_mono_err, T2_mono_err = _mono_exp_refit(
            x, y, A_init=0.5 * A, T2_init=slow_T2,
            A_max=A_max, t2_max=t2_max,
        )
        A = A_mono
        A_err = A_mono_err
        T2_avg = T2_mono
        T2_avg_err = T2_mono_err
        t2_a = T2_mono
        t2_b = T2_mono
        dT2 = 0.0
        # Bi-exp covariance no longer applies — secondary errs derived below
        # will fall back to NaN, which we override appropriately.
        covar, var_names = None, None

    t2_a_err, t2_b_err, eta_err = _derive_secondary_errors(
        T2_avg, dT2, covar, var_names
    )

    # When dT2 sits at its 0 lower bound (true mono-exp regime, OR forced
    # there above), the Hessian is exactly singular along the dT2 direction;
    # lmfit's pseudo-inverse spits out a finite-but-meaningless huge stderr
    # (e.g. 1e6 ms) from the unconstrained Jacobian. Suppress dT2_err / eta_err
    # to NaN — a parameter pinned at a constraint has no curvature-based
    # uncertainty. T2_avg_err and A_err remain meaningful (orthogonal
    # well-determined directions, or replaced by the mono-exp refit above).
    # Optimizer convergence usually leaves dT2 at ~1e-6 instead of exactly 0
    # when it lands on the lower bound. 1e-3 ms is three orders of magnitude
    # smaller than any plausible methyl T2 separation, so it's a safe "at bound"
    # threshold that won't catch real-but-small bi-exp character.
    DT2_AT_BOUND_EPS = 1e-3
    dT2_at_lower_bound = dT2 < DT2_AT_BOUND_EPS
    bound_active = t2_b_bound_active or dT2_at_lower_bound
    if bound_active:
        dT2_err = float("nan")
        eta_err = float("nan")
        # In the bound-active case t2_a == t2_b == T2_avg (collapsed to mono-exp),
        # so their stderrs follow T2_avg_err — finite when the mono-exp refit
        # ran above; NaN otherwise (true mono-exp case where the bi-exp covar
        # in dT2 was singular).
        if t2_b_bound_active:
            t2_a_err = T2_avg_err
            t2_b_err = T2_avg_err
        else:
            t2_a_err = float("nan")
            t2_b_err = float("nan")

    if t2_b != 0:
        eta = t2_a / t2_b
    else:
        eta = float("nan")

    bi_exp_unidentifiable = bound_active or _is_unidentifiable(dT2, dT2_err)

    return {
        "residue": residue_name,
        "A": float(A),
        "T2_avg": float(T2_avg),
        "dT2": float(dT2),
        "t2_a": float(t2_a),
        "t2_b": float(t2_b),
        "eta": float(eta),
        "A_err": float(A_err),
        "T2_avg_err": float(T2_avg_err),
        "dT2_err": float(dT2_err),
        "t2_a_err": float(t2_a_err),
        "t2_b_err": float(t2_b_err),
        "eta_err": float(eta_err),
        "bi_exp_unidentifiable": bi_exp_unidentifiable,
        "x": x, "y": y,
        "result": result,
    }


def save_results_methyl(results_list, output_file):
    """Write methyl T2 bi-exp results to a tab-delimited text file.

    Columns:
        Residue \t A \t T2a \t T2b \t T2_avg \t dT2 \t eta \t
        A_err \t T2a_err \t T2b_err \t T2_avg_err \t dT2_err \t eta_err \t
        bi_exp_unidentifiable \t Success
    """
    with open(output_file, "w") as f:
        f.write(
            "Residue\tA\tT2a\tT2b\tT2_avg\tdT2\teta\t"
            "A_err\tT2a_err\tT2b_err\tT2_avg_err\tdT2_err\teta_err\t"
            "bi_exp_unidentifiable\tSuccess\n"
        )
        for result in results_list:
            success = "Yes" if result.get("success", False) else "No"
            unident = "Yes" if result.get("bi_exp_unidentifiable", False) else "No"
            f.write(
                f"{result['residue']}\t"
                f"{result['A']:.6e}\t{result['t2_a']:.6e}\t{result['t2_b']:.6e}\t"
                f"{result.get('T2_avg', float('nan')):.6e}\t"
                f"{result.get('dT2', float('nan')):.6e}\t"
                f"{result['eta']:.6e}\t"
                f"{result['A_err']:.6e}\t{result['t2_a_err']:.6e}\t{result['t2_b_err']:.6e}\t"
                f"{result.get('T2_avg_err', float('nan')):.6e}\t"
                f"{result.get('dT2_err', float('nan')):.6e}\t"
                f"{result['eta_err']:.6e}\t"
                f"{unident}\t{success}\n"
            )


def save_fit_data_json_methyl(results_list, output_file, time_units, signal_units,
                              n_bootstrap, field_freq):
    """Save methyl shared-amplitude bi-exp fit data as JSON for the methyl viewer.

    Each `fits[i]` entry carries A, T2_avg, dT2, t2_a, t2_b, eta, their errors,
    and the `bi_exp_unidentifiable` flag.
    """
    if not results_list:
        raise ValueError("No results to save")

    first_x = results_list[0]["x"]
    time_points = first_x.tolist() if hasattr(first_x, "tolist") else list(first_x)
    max_time = max(time_points)
    fit_time_dense = np.linspace(0.0, max_time * 1.2, 100)

    fits_data = []
    for result in results_list:
        fit_intensity = 0.5 * result["A"] * (
            np.exp(-fit_time_dense / result["t2_a"])
            + np.exp(-fit_time_dense / result["t2_b"])
        )
        y = result["y"]
        intensities = y.tolist() if hasattr(y, "tolist") else list(y)

        fits_data.append({
            "residue": str(result["residue"]),
            "A": float(result["A"]),
            "T2_avg": float(result.get("T2_avg", float("nan"))),
            "dT2": float(result.get("dT2", float("nan"))),
            "t2_a": float(result["t2_a"]),
            "t2_b": float(result["t2_b"]),
            "eta": float(result["eta"]),
            "A_err": float(result["A_err"]),
            "T2_avg_err": float(result.get("T2_avg_err", float("nan"))),
            "dT2_err": float(result.get("dT2_err", float("nan"))),
            "t2_a_err": float(result["t2_a_err"]),
            "t2_b_err": float(result["t2_b_err"]),
            "eta_err": float(result["eta_err"]),
            "bi_exp_unidentifiable": bool(result.get("bi_exp_unidentifiable", False)),
            "intensities": [float(v) for v in intensities],
            "fit_curve": {
                "time": [float(t) for t in fit_time_dense],
                "intensity": [float(i) for i in fit_intensity],
            },
        })

    output_data = {
        "metadata": {
            "experiment_type": "methylT2",
            "field_freq": float(field_freq),
            "time_units": time_units,
            "signal_units": signal_units,
            "n_bootstrap": n_bootstrap,
            "n_residues": len(results_list),
            "time_points": [float(t) for t in time_points],
        },
        "fits": fits_data,
    }

    out_dir = os.path.dirname(output_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(output_file, "w") as f:
        json.dump(output_data, f, indent=2)


def run_methyl_t2_analysis_with_params(params, progress_callback=None):
    """Run methyl shared-amplitude bi-exp T2 analysis on a CSV input.

    Required `params` keys:
        input_csv_file, output_prefix, results_txt_file, json_folder,
        field_name, field_freq
    Optional:
        time_units (default 'ms'), signal_units (default 'Intensity'),
        initial_A, initial_t2_a, initial_t2_b,
        n_bootstrap (default 1000), error_method (default 'analytical')
    """
    input_csv_file = params["input_csv_file"]
    output_prefix = params["output_prefix"]
    results_txt_file = params["results_txt_file"]
    json_folder = params["json_folder"]
    field_name = params["field_name"]
    field_freq = float(params["field_freq"])

    time_units = params.get("time_units", "ms")
    signal_units = params.get("signal_units", "Intensity")

    initial_A = params.get("initial_A")
    initial_t2_a = params.get("initial_t2_a", DEFAULT_INITIAL_T2A)
    initial_t2_b = params.get("initial_t2_b", DEFAULT_INITIAL_T2B)
    n_bootstrap = int(params.get("n_bootstrap", 1000))
    error_method = params.get("error_method", "analytical")

    if not os.path.exists(input_csv_file):
        raise FileNotFoundError(f"Input file not found: {input_csv_file}")

    print(f"Starting methyl T2 (shared-amp bi-exp) fitting...")
    print(f"Input file: {input_csv_file}")
    print(f"Output prefix: {output_prefix}")
    print(f"Error method: {error_method}")

    raw_df = pd.read_csv(input_csv_file, header=None)

    header_row = raw_df.iloc[0].astype(str).tolist()
    lunaNMR_columns = {"Peak_Number", "Assignment", "Reference_X", "Reference_Y"}
    detected = [col for col in header_row if col in lunaNMR_columns]
    # Columns whose header carries no delay: reported so a spectrum cannot go
    # missing from a series without anyone noticing.
    dropped_columns = []

    if detected:
        delay_start_idx = require_delay_start(header_row, lunaNMR_columns, input_csv_file)
        assignment_idx = header_row.index("Assignment") if "Assignment" in header_row else 1
        residue_names = raw_df.iloc[1:, assignment_idx].to_numpy()
        delay_headers = raw_df.iloc[0, delay_start_idx:].tolist()
        x, _keep, _dropped = parse_delay_columns(delay_headers)
        dropped_columns.extend(_dropped)
        y_data = raw_df.iloc[1:, delay_start_idx:].astype(float).to_numpy()[:, _keep]
    else:
        residue_names = raw_df.iloc[1:, 0].to_numpy()
        delay_headers = raw_df.iloc[0, 1:].tolist()
        x, _keep, _dropped = parse_delay_columns(delay_headers)
        dropped_columns.extend(_dropped)
        y_data = raw_df.iloc[1:, 1:].astype(float).to_numpy()[:, _keep]

    valid_mask = np.array([not str(name).lower().startswith("dummy") for name in residue_names])
    n_dummies = int(np.sum(~valid_mask))
    if n_dummies:
        residue_names = residue_names[valid_mask]
        y_data = y_data[valid_mask]

    total = len(residue_names)
    print(f"Loaded {total} residues with {len(x)} time points")

    results_list = []
    for idx in range(total):
        y = y_data[idx, :]
        residue = residue_names[idx]
        try:
            res = fit_single_residue_methyl(
                x, y, residue,
                initial_A=initial_A,
                initial_t2_a=initial_t2_a, initial_t2_b=initial_t2_b,
                n_bootstrap=n_bootstrap, error_method=error_method,
            )
            res["success"] = True
        except Exception as e:
            print(f"Warning: fit failed for {residue}: {e}")
            res = {
                "residue": residue, "success": False,
                "A": np.nan,
                "T2_avg": np.nan, "dT2": np.nan,
                "t2_a": np.nan, "t2_b": np.nan, "eta": np.nan,
                "A_err": np.nan,
                "T2_avg_err": np.nan, "dT2_err": np.nan,
                "t2_a_err": np.nan, "t2_b_err": np.nan, "eta_err": np.nan,
                "bi_exp_unidentifiable": False,
                "x": x, "y": y,
            }
        results_list.append(res)
        if progress_callback:
            progress_callback(idx + 1, total, str(residue),
                              f"Fitted {residue} ({idx + 1}/{total})")

    save_results_methyl(results_list, results_txt_file)

    os.makedirs(json_folder, exist_ok=True)
    json_filename = f"{field_name}_methylT2_fit_data.json"
    json_path = os.path.join(json_folder, json_filename)
    successful = [r for r in results_list if r.get("success", False)]
    if successful:
        save_fit_data_json_methyl(
            successful, json_path,
            time_units=time_units, signal_units=signal_units,
            n_bootstrap=n_bootstrap, field_freq=field_freq,
        )

    return {
        "n_fitted": len(successful),
        "dropped_columns": dropped_columns,
        "n_total": total,
        "output_dir": os.path.dirname(output_prefix),
        "results_file": results_txt_file,
        "json_file": json_path if successful else None,
    }
