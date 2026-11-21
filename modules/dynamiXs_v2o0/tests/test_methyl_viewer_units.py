# ABOUTME: The methyl T2 viewer labels its across-residue bars in the delay units the fit
# ABOUTME: JSON records. A real fit of a synthetic matrix in seconds, opened in the viewer.

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DIR))
sys.path.insert(0, str(_DIR / "dynamiXs_T1_T2"))

DELAYS_MS = np.array([2.0, 4.0, 8.0, 12.0, 18.0, 25.0, 35.0, 50.0, 75.0, 100.0,
                      150.0, 200.0, 250.0, 300.0])


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


def _fit_json(folder, units, per_ms):
    from fit_methyl_T2 import run_methyl_t2_analysis_with_params
    folder.mkdir()
    rows = []
    for seed in range(3):
        rng = np.random.default_rng(seed)
        y = 0.5 * 2.0e6 * (np.exp(-DELAYS_MS / 80.0) + np.exp(-DELAYS_MS / 15.0))
        y = y + rng.normal(0.0, 4.0e3, size=y.size)
        rows.append({"Assignment": f"M{seed}",
                     **{f"{d * per_ms:g}": v for d, v in zip(DELAYS_MS, y)}})
    pd.DataFrame(rows).to_csv(folder / "matrix.csv", index=False)
    result = run_methyl_t2_analysis_with_params({
        "input_csv_file": str(folder / "matrix.csv"), "output_prefix": str(folder / "p"),
        "results_txt_file": str(folder / "p_fit_results.txt"), "json_folder": str(folder / "json"),
        "field_name": "p_field1", "field_freq": 600.0, "time_units": units})
    return result["json_file"]


@pytest.mark.parametrize("units, per_ms", [("s", 1e-3), ("ms", 1.0)])
def test_the_summary_bars_are_labelled_in_the_fits_units(app, tmp_path, monkeypatch, units, per_ms):
    import visualization.methyl_t2_fit_viewer as module
    for box in ("show_info", "show_warning", "show_error"):
        monkeypatch.setattr(module, box, lambda *a, **k: None)
    json_file = _fit_json(tmp_path / units, units, per_ms)
    viewer = module.MethylT2FitViewer(summary_only=True, json_file=json_file)
    try:
        labels = {}
        for key, _, _ in viewer.SUMMARY_QUANTITIES:
            viewer.summary_quantity = key
            viewer._update_plot()
            labels[key] = viewer.figure.axes[0].get_ylabel()
        assert labels["T2_avg"] == f"T2_avg ({units})"
        assert labels["dT2"] == f"ΔT2 ({units})"
        assert labels["t2_a"] == f"t2_a ({units})"
        assert labels["t2_b"] == f"t2_b ({units})"
        assert labels["eta"] == "η = T2a / T2b"
    finally:
        viewer.close()
