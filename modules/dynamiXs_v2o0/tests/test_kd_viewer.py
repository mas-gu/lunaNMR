# ABOUTME: Smoke test for the Kd titration viewer - loads a fit JSON and renders.
# ABOUTME: Headless (offscreen Qt); exercises load + view switching.

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_MOD = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_MOD))
sys.path.insert(0, str(_MOD / "dynamiXs_Kd"))

from PySide6.QtWidgets import QApplication


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication(sys.argv)


def _make_fit_json(tmp_path):
    from kd_models import csp_model
    from kd_fit import run_kd_analysis_with_params
    P0 = 50.0
    pts = [0.0, 10.0, 25.0, 60.0, 150.0, 300.0]
    rows = []
    for name, dd_max, kd in (("A", 0.2, 15.0), ("B", 0.3, 15.0)):
        d = csp_model(np.array(pts), dd_max, kd, P0)
        for p, dv in zip(pts, d):
            rows.append((str(p), name, 8.0 + dv, 120.0, 1000.0 - 2 * p, 2000.0))
    df = pd.DataFrame(rows, columns=["spectrum_name", "assignment",
                                     "ppm_x", "ppm_y", "height", "volume"])
    csv = tmp_path / "series_analysis_tidy.csv"
    df.to_csv(csv, index=False)
    r = run_kd_analysis_with_params({
        "input_csv_file": str(csv), "output_dir": str(tmp_path), "output_prefix": "kd",
        "concentrations": pts, "protein_conc": P0, "alpha": 0.14,
        "observables": ["csp", "intensity"], "n_bootstrap": 0})
    return r["json_file"]


def test_viewer_loads_and_switches_views(app, tmp_path):
    from visualization.kd_titration_fit_viewer import open_kd_titration_viewer
    jf = _make_fit_json(tmp_path)
    v = open_kd_titration_viewer(parent=None, json_file=jf)
    try:
        assert v.residue_list.count() == 2
        assert "Global shared Kd" in v.global_label.text()   # 2 residues -> global fit
        for i in range(3):           # per-residue / Kd-bar / ddmax-bar
            v.view_combo.setCurrentIndex(i)
            v._refresh()
    finally:
        v.deleteLater()


def test_exclude_point_and_refit_persists(app, tmp_path):
    import json
    from visualization.kd_titration_fit_viewer import open_kd_titration_viewer
    jf = _make_fit_json(tmp_path)
    v = open_kd_titration_viewer(parent=None, json_file=jf)
    try:
        v.obs_combo.setCurrentIndex(0)           # CSP
        v.residue_list.setCurrentRow(0)
        f = v._current_residue()
        kd_before = f["csp"]["Kd"]
        # exclude the last point and refit
        v.excluded["csp"][f["residue"]] = {len(f["csp"]["L"]) - 1}
        v._on_refit()
        # full display points are preserved; exclusion recorded; params recomputed
        assert len(f["csp"]["L"]) == 6
        assert f["csp"]["excluded"] == [5]
        # persisted to disk
        saved = json.loads(open(jf).read())
        assert saved["fits"][0]["csp"]["excluded"] == [5]
        # reset clears the exclusion
        v._on_reset()
        assert f["csp"].get("excluded") == []
        assert f["csp"]["Kd"] == pytest.approx(kd_before, rel=1e-6)
    finally:
        v.deleteLater()
