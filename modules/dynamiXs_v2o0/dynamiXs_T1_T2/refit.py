# ABOUTME: Headless helpers for interactive outlier rejection in dynamiXs T1/T2 fits.
# ABOUTME: Sidecar (.outliers.json) IO, single-residue refit, and surgical JSON/TSV row updates.

import json
from pathlib import Path
from typing import Iterable

import numpy as np

from fit_Tx_NMRRE import fit_single_residue
from fit_methyl_T2 import fit_single_residue_methyl


def sidecar_path_for(json_path) -> Path:
    """Return the .outliers.json path that pairs with a *_fit_data.json file."""
    p = Path(json_path)
    return p.with_suffix(".outliers.json") if p.suffix == ".json" else p.parent / (p.name + ".outliers.json")


def load_outliers(json_path) -> dict:
    """Load per-residue exclusion indices from the sidecar next to json_path.

    Returns an empty dict if the sidecar is absent or malformed.
    """
    sc = sidecar_path_for(json_path)
    if not sc.exists():
        return {}
    try:
        payload = json.loads(sc.read_text())
    except (json.JSONDecodeError, OSError):
        return {}
    exclusions = payload.get("exclusions", {})
    return {str(k): list(v) for k, v in exclusions.items()}


def save_outliers(json_path, exclusions: dict) -> None:
    """Write exclusions to the sidecar. Overwrites any previous content."""
    sc = sidecar_path_for(json_path)
    payload = {
        "version": 1,
        "exclusions": {str(k): list(v) for k, v in exclusions.items()},
    }
    sc.write_text(json.dumps(payload, indent=2))


def refit_residue(fit_entry: dict, metadata: dict,
                  excluded_indices: Iterable[int],
                  error_method: str = "analytical") -> dict:
    """Re-fit one residue with the given time-point indices excluded.

    The original `intensities` array is preserved verbatim in the returned entry —
    exclusions live in the sidecar, never in the data record.

    Initial guesses are warm-started from the prior fit's A/t2/C.
    The dense `fit_curve` is regenerated over the original full time range so plots
    remain visually consistent after points are dropped.
    """
    excluded = sorted(set(int(i) for i in excluded_indices))

    time_points = np.asarray(metadata["time_points"], dtype=float)
    intensities = np.asarray(fit_entry["intensities"], dtype=float)

    keep_mask = np.ones(len(time_points), dtype=bool)
    for i in excluded:
        if 0 <= i < len(keep_mask):
            keep_mask[i] = False

    x_fit = time_points[keep_mask]
    y_fit = intensities[keep_mask]

    n_bootstrap = int(metadata.get("n_bootstrap", 1000))
    result = fit_single_residue(
        x_fit, y_fit, fit_entry.get("residue", ""),
        initial_A=float(fit_entry.get("A", 1.0)),
        initial_t2=float(fit_entry.get("t2", 100.0)),
        initial_C=float(fit_entry.get("C", 0.0)),
        n_bootstrap=n_bootstrap,
        error_method=error_method,
    )

    # Dense curve over the original full time range, so plots stay consistent.
    max_t = float(time_points.max()) if len(time_points) else 0.0
    fit_time_dense = np.linspace(0.0, max_t * 1.2, 100)
    fit_intensity_dense = (
        result["A"] * np.exp(-fit_time_dense / result["t2"]) + result["C"]
    )

    return {
        "residue": fit_entry.get("residue", ""),
        "A": float(result["A"]),
        "t2": float(result["t2"]),
        "C": float(result["C"]),
        "A_err": float(result["A_err"]),
        "t2_err": float(result["t2_err"]),
        "C_err": float(result["C_err"]),
        "intensities": list(map(float, intensities)),
        "fit_curve": {
            "time": [float(t) for t in fit_time_dense],
            "intensity": [float(i) for i in fit_intensity_dense],
        },
    }


def update_json_fit_entry(json_path, residue: str, new_fit_entry: dict) -> None:
    """Replace the entry for `residue` in a *_fit_data.json file in place.

    Preserves metadata and all other fits. Raises KeyError if the residue is absent.
    """
    p = Path(json_path)
    payload = json.loads(p.read_text())
    fits = payload.get("fits", [])
    for i, entry in enumerate(fits):
        if str(entry.get("residue")) == str(residue):
            fits[i] = new_fit_entry
            payload["fits"] = fits
            p.write_text(json.dumps(payload, indent=2))
            return
    raise KeyError(f"Residue {residue!r} not found in {p}")


def update_tsv_row(tsv_path, residue: str, new_fit: dict, experiment_type: str) -> None:
    """Replace one Residue row in a fit_results.txt TSV.

    The TSV is written by `fit_Tx_NMRRE.save_results` with header:
        Residue \t A \t {experiment_type} \t C \t A_err \t {experiment_type}_err \t C_err
    (optionally followed by `\tSuccess` from the multicore variant).

    Raises KeyError if the residue is absent.
    """
    p = Path(tsv_path)
    lines = p.read_text().splitlines()
    if not lines:
        raise KeyError(f"TSV at {p} is empty")

    header = lines[0]
    cols = header.split("\t")
    has_success = cols and cols[-1] == "Success"

    out_lines = [header]
    found = False
    for line in lines[1:]:
        if not line.strip():
            out_lines.append(line)
            continue
        first = line.split("\t", 1)[0]
        if first == str(residue):
            cells = [
                str(residue),
                f"{new_fit['A']:.6e}",
                f"{new_fit['t2']:.6e}",
                f"{new_fit['C']:.6e}",
                f"{new_fit['A_err']:.6e}",
                f"{new_fit['t2_err']:.6e}",
                f"{new_fit['C_err']:.6e}",
            ]
            if has_success:
                cells.append("Yes")
            out_lines.append("\t".join(cells))
            found = True
        else:
            out_lines.append(line)

    if not found:
        raise KeyError(f"Residue {residue!r} not found in {p}")

    p.write_text("\n".join(out_lines) + "\n")


# ---------- methyl bi-exp variants ----------

def refit_residue_methyl(fit_entry: dict, metadata: dict,
                         excluded_indices: Iterable[int],
                         error_method: str = "analytical") -> dict:
    """Re-fit one residue's methyl bi-exp data (shared-amplitude form).

    Model: I(t) = 0.5 * A * (exp(-t/t2_a) + exp(-t/t2_b)).
    Original `intensities` are preserved verbatim; the dense `fit_curve` is
    regenerated over the full original time range so plots stay consistent.

    Uses fresh data-driven initial guesses, NOT warm-start from the prior fit:
    a refit is normally triggered because the prior fit was unsatisfactory, and
    warm-starting from those parameters can leave the L-M optimizer trapped in
    the same poor local minimum.
    """
    excluded = sorted(set(int(i) for i in excluded_indices))

    time_points = np.asarray(metadata["time_points"], dtype=float)
    intensities = np.asarray(fit_entry["intensities"], dtype=float)

    keep_mask = np.ones(len(time_points), dtype=bool)
    for i in excluded:
        if 0 <= i < len(keep_mask):
            keep_mask[i] = False

    x_fit = time_points[keep_mask]
    y_fit = intensities[keep_mask]

    n_bootstrap = int(metadata.get("n_bootstrap", 1000))
    result = fit_single_residue_methyl(
        x_fit, y_fit, fit_entry.get("residue", ""),
        n_bootstrap=n_bootstrap,
        error_method=error_method,
    )

    max_t = float(time_points.max()) if len(time_points) else 0.0
    fit_time_dense = np.linspace(0.0, max_t * 1.2, 100)
    fit_intensity_dense = 0.5 * result["A"] * (
        np.exp(-fit_time_dense / result["t2_a"])
        + np.exp(-fit_time_dense / result["t2_b"])
    )

    return {
        "residue": fit_entry.get("residue", ""),
        "A": float(result["A"]),
        "T2_avg": float(result["T2_avg"]),
        "dT2": float(result["dT2"]),
        "t2_a": float(result["t2_a"]), "t2_b": float(result["t2_b"]),
        "eta": float(result["eta"]) if not np.isnan(result["eta"]) else result["eta"],
        "A_err": float(result["A_err"]),
        "T2_avg_err": float(result["T2_avg_err"]),
        "dT2_err": float(result["dT2_err"]),
        "t2_a_err": float(result["t2_a_err"]), "t2_b_err": float(result["t2_b_err"]),
        "eta_err": float(result["eta_err"]) if not np.isnan(result["eta_err"]) else result["eta_err"],
        "bi_exp_unidentifiable": bool(result["bi_exp_unidentifiable"]),
        "intensities": list(map(float, intensities)),
        "fit_curve": {
            "time": [float(t) for t in fit_time_dense],
            "intensity": [float(i) for i in fit_intensity_dense],
        },
    }


def update_tsv_row_methyl(tsv_path, residue: str, new_fit: dict) -> None:
    """Replace one Residue row in a methyl shared-amp bi-exp fit_results.txt TSV.

    Header layout written by `fit_methyl_T2.save_results_methyl`:
        Residue \t A \t T2a \t T2b \t T2_avg \t dT2 \t eta
                \t A_err \t T2a_err \t T2b_err \t T2_avg_err \t dT2_err \t eta_err
                \t bi_exp_unidentifiable \t Success
    """
    p = Path(tsv_path)
    lines = p.read_text().splitlines()
    if not lines:
        raise KeyError(f"TSV at {p} is empty")

    header = lines[0]
    cols = header.split("\t")
    has_success = bool(cols) and cols[-1] == "Success"
    has_unident = "bi_exp_unidentifiable" in cols

    out_lines = [header]
    found = False
    nan = float("nan")
    for line in lines[1:]:
        if not line.strip():
            out_lines.append(line)
            continue
        first = line.split("\t", 1)[0]
        if first == str(residue):
            cells = [
                str(residue),
                f"{new_fit['A']:.6e}",
                f"{new_fit['t2_a']:.6e}",
                f"{new_fit['t2_b']:.6e}",
            ]
            if has_unident:
                cells += [
                    f"{new_fit.get('T2_avg', nan):.6e}",
                    f"{new_fit.get('dT2', nan):.6e}",
                ]
            cells += [
                f"{new_fit['eta']:.6e}",
                f"{new_fit['A_err']:.6e}",
                f"{new_fit['t2_a_err']:.6e}",
                f"{new_fit['t2_b_err']:.6e}",
            ]
            if has_unident:
                cells += [
                    f"{new_fit.get('T2_avg_err', nan):.6e}",
                    f"{new_fit.get('dT2_err', nan):.6e}",
                ]
            cells.append(f"{new_fit['eta_err']:.6e}")
            if has_unident:
                cells.append(
                    "Yes" if new_fit.get("bi_exp_unidentifiable", False) else "No"
                )
            if has_success:
                cells.append("Yes")
            out_lines.append("\t".join(cells))
            found = True
        else:
            out_lines.append(line)

    if not found:
        raise KeyError(f"Residue {residue!r} not found in {p}")

    p.write_text("\n".join(out_lines) + "\n")
