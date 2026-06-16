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


def test_residue_sort_key_orders_by_number_not_letter():
    from visualization.kd_titration_fit_viewer import _residue_sort_key
    names = ["A17", "K14", "C76", "E20", "E21", "A25"]
    assert sorted(names, key=_residue_sort_key) == [
        "K14", "A17", "E20", "E21", "A25", "C76"]


def test_residue_list_sorted_by_number(app, tmp_path):
    # Build a JSON whose fits are alphabetical (A.. before K..) and confirm the
    # viewer reorders the residue list numerically (K14 before A17).
    import json
    from visualization.kd_titration_fit_viewer import open_kd_titration_viewer
    jf = _make_fit_json(tmp_path)
    data = json.loads(Path(jf).read_text())
    base = data["fits"][0]
    data["fits"] = [dict(base, residue=r) for r in ("A17", "K14", "C76")]
    Path(jf).write_text(json.dumps(data))
    v = open_kd_titration_viewer(parent=None, json_file=jf)
    try:
        order = [v.residue_list.item(i).text() for i in range(v.residue_list.count())]
        assert order == ["K14", "A17", "C76"]
    finally:
        v.deleteLater()


def test_pair_observable_csp_and_ratio():
    import math
    from visualization.kd_titration_fit_viewer import _pair_observable
    series = {"ppm_x": [8.0, 8.1, 8.3], "ppm_y": [120.0, 120.5, 121.0],
              "height": [1000.0, 500.0, 250.0], "volume": [2000.0, 1000.0, 500.0]}
    # CSP point0 -> point2: sqrt(0.3^2 + (0.14*1.0)^2)
    assert _pair_observable(series, 0, 2, "csp", alpha=0.14) == pytest.approx(
        math.hypot(0.3, 0.14 * 1.0))
    # intensity ratio point0 -> point1
    assert _pair_observable(series, 0, 1, "intensity", value="height") == pytest.approx(0.5)
    # undefined reference intensity -> NaN
    assert math.isnan(_pair_observable({"height": [0.0, 5.0]}, 0, 1, "intensity"))
    # undetected position sentinel (0.0) -> NaN
    assert math.isnan(_pair_observable(
        {"ppm_x": [0.0, 8.1], "ppm_y": [0.0, 120.5]}, 0, 1, "csp"))


def test_comparison_view_renders_with_point_selectors(app, tmp_path):
    from visualization.kd_titration_fit_viewer import open_kd_titration_viewer
    jf = _make_fit_json(tmp_path)
    v = open_kd_titration_viewer(parent=None, json_file=jf)
    try:
        # selectors hidden until the comparison view; default ref=0, compare=1
        assert v.cmp_box.isVisible() is False
        assert v.ref_point_combo.count() == 6
        assert v.ref_point_combo.currentIndex() == 0
        assert v.cmp_point_combo.currentIndex() == 1
        v.view_combo.setCurrentIndex(3)          # "Ref vs point (bars)"
        assert v.cmp_box.isVisible() is True
        # renders bars for both observables without error
        v.obs_combo.setCurrentIndex(0)           # CSP
        assert len(v.figure.axes[0].patches) >= 1
        v.obs_combo.setCurrentIndex(1)           # Intensity
        assert len(v.figure.axes[0].patches) >= 1
    finally:
        v.deleteLater()
