# ABOUTME: The methyl bi-exponential T2 fitter gives the same decay whatever unit its delays
# ABOUTME: are stated in: its bounds and starting values follow that unit.

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_FITTER_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_T1_T2"
sys.path.insert(0, str(_FITTER_DIR))

DELAYS_MS = np.array([2.0, 4.0, 8.0, 12.0, 18.0, 25.0, 35.0, 50.0, 75.0, 100.0,
                      150.0, 200.0, 250.0, 300.0])
PER_MS = {"s": 1e-3, "us": 1e3}


def _decay(seed=0, t2_a=80.0, t2_b=15.0):
    """A noisy 80 / 15 ms shared-amplitude decay sampled at DELAYS_MS."""
    rng = np.random.default_rng(seed)
    y = 0.5 * 2.0e6 * (np.exp(-DELAYS_MS / t2_a) + np.exp(-DELAYS_MS / t2_b))
    return y + rng.normal(0.0, 4.0e3, size=y.size)


@pytest.mark.parametrize("units", ["s", "us"])
def test_a_decay_fits_the_same_in_any_delay_unit(units):
    from fit_methyl_T2 import fit_single_residue_methyl
    y = _decay()
    in_ms = fit_single_residue_methyl(DELAYS_MS, y, "M")
    in_units = fit_single_residue_methyl(DELAYS_MS * PER_MS[units], y, "M", time_units=units)

    assert not in_ms["bi_exp_unidentifiable"]
    assert in_ms["t2_a"] == pytest.approx(80.0, rel=0.1)
    assert in_ms["t2_b"] == pytest.approx(15.0, rel=0.1)
    assert in_units["bi_exp_unidentifiable"] == in_ms["bi_exp_unidentifiable"]
    for key in ("T2_avg", "dT2", "t2_a", "t2_b", "T2_avg_err", "dT2_err", "t2_a_err", "t2_b_err"):
        assert in_units[key] == pytest.approx(in_ms[key] * PER_MS[units], rel=1e-4), key
    for key in ("A", "eta"):
        assert in_units[key] == pytest.approx(in_ms[key], rel=1e-4), key


def test_a_single_component_decay_in_seconds_keeps_its_time_constant():
    """One component only: the fit collapses to mono-exponential, at the decay's own
    time constant and not at the smallest one the fitter allows."""
    from fit_methyl_T2 import fit_single_residue_methyl
    rng = np.random.default_rng(1)
    y = 1.0e6 * np.exp(-DELAYS_MS / 60.0) + rng.normal(0.0, 2.0e3, size=DELAYS_MS.size)
    in_ms = fit_single_residue_methyl(DELAYS_MS, y, "MONO")
    in_s = fit_single_residue_methyl(DELAYS_MS * 1e-3, y, "MONO", time_units="s")

    assert in_ms["bi_exp_unidentifiable"] and in_s["bi_exp_unidentifiable"]
    assert in_s["T2_avg"] == pytest.approx(0.060, rel=0.05)
    assert in_s["T2_avg"] == pytest.approx(in_ms["T2_avg"] * 1e-3, rel=1e-3)


def _matrix(path, per_ms):
    rows = [{"Assignment": f"M{seed}", **{f"{d * per_ms:g}": v
                                          for d, v in zip(DELAYS_MS, _decay(seed))}}
            for seed in range(3)]
    pd.DataFrame(rows).to_csv(path, index=False)


def _run(tmp_path, units, per_ms):
    from fit_methyl_T2 import run_methyl_t2_analysis_with_params
    out = tmp_path / units
    out.mkdir()
    _matrix(out / "matrix.csv", per_ms)
    result = run_methyl_t2_analysis_with_params({
        "input_csv_file": str(out / "matrix.csv"), "output_prefix": str(out / "p"),
        "results_txt_file": str(out / "p_fit_results.txt"), "json_folder": str(out),
        "field_name": "field1", "field_freq": 600.0, "time_units": units})
    assert result["n_fitted"] == 3
    return result, pd.read_csv(result["results_file"], sep="\t")


def test_a_matrix_in_seconds_is_fitted_as_the_same_matrix_in_ms(tmp_path):
    result_ms, in_ms = _run(tmp_path, "ms", 1.0)
    result_s, in_s = _run(tmp_path, "s", 1e-3)

    assert result_s["n_bi_exp_unidentifiable"] == result_ms["n_bi_exp_unidentifiable"] == 0
    assert list(in_s["bi_exp_unidentifiable"]) == ["No", "No", "No"]
    for column in ("T2a", "T2b"):
        assert in_s[column].to_numpy() == pytest.approx(in_ms[column].to_numpy() * 1e-3, rel=1e-4)
    assert in_s["T2a"].between(0.07, 0.09).all() and in_s["T2b"].between(0.012, 0.018).all()


def test_a_refit_reads_the_units_the_fit_json_records(tmp_path):
    from refit import refit_residue_methyl
    result, _ = _run(tmp_path, "s", 1e-3)
    payload = json.loads(Path(result["json_file"]).read_text())
    entry = payload["fits"][0]

    refitted = refit_residue_methyl(entry, payload["metadata"], excluded_indices=[])

    assert payload["metadata"]["time_units"] == "s"
    assert not refitted["bi_exp_unidentifiable"]
    assert refitted["t2_a"] == pytest.approx(entry["t2_a"], rel=1e-3)
    assert refitted["t2_b"] == pytest.approx(entry["t2_b"], rel=1e-3)
