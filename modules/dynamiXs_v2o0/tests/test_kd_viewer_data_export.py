# ABOUTME: Tests the Kd fit viewer writes CSV + JSON data alongside the ref->point bars.
# ABOUTME: The vs-sequence figures must be reproducible from a table, not PDF-only.

import os
import sys
import json
import csv
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_MOD = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_MOD))


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _fit(res, xs, ys, hs):
    L = [0.0, 6.25, 12.5]
    return {"residue": res,
            "series": {"ppm_x": xs, "ppm_y": ys, "height": hs, "volume": hs},
            "csp": {"success": True, "Kd": 10.0, "Kd_err": 1.0, "dd_max": 0.1,
                    "r_squared": 0.9, "L": L, "obs": [0.0, 0.02, 0.03]},
            "intensity": {"success": True, "Kd": 10.0, "Kd_err": 1.0, "I0": 1.0,
                          "I_inf": 0.2, "r_squared": 0.9, "L": L, "obs": [1.0, 0.8, 0.5]}}


def _write_json(tmp_path):
    data = {"metadata": {"analysis": "Kd_titration", "protein_conc": 50.0, "alpha": 0.14,
                        "concentrations": [0.0, 6.25, 12.5], "points": [0.0, 1.0, 2.0],
                        "intensity_value": "height", "observables": ["csp", "intensity"],
                        "n_bootstrap": 0, "name": "MINI"},
            "fits": [_fit("A17", [8.0, 8.02, 8.05], [120.0, 120.1, 120.3], [100.0, 80.0, 50.0]),
                     _fit("K14", [7.5, 7.51, 0.0], [110.0, 110.05, 0.0], [200.0, 150.0, 0.0])],
            "global": {}}
    jf = tmp_path / "MINI_kd_fit_data.json"
    jf.write_text(json.dumps(data))
    return jf


def test_viewer_writes_csv_and_json(app, tmp_path):
    from visualization.kd_titration_fit_viewer import KdTitrationFitViewer
    jf = _write_json(tmp_path)
    viewer = KdTitrationFitViewer(json_file=str(jf))

    paths = viewer._write_ref_vs_point_data(str(tmp_path), "MINI", "csp")
    assert sorted(Path(p).suffix for p in paths) == ['.csv', '.json']

    rows = list(csv.reader(open(str(tmp_path / "MINI_csp_ref_vs_point.csv"))))
    assert rows[0][0] == 'Residue'
    assert [r[0] for r in rows[1:]] == ['K14', 'A17']    # sequence-number order


def test_viewer_exports_kd_vs_residue_pdf(app, tmp_path):
    from visualization.kd_titration_fit_viewer import KdTitrationFitViewer
    jf = _write_json(tmp_path)
    viewer = KdTitrationFitViewer(json_file=str(jf))
    for obs in ("csp", "intensity"):
        out = tmp_path / f"MINI_{obs}_kd_vs_residue.pdf"
        assert viewer._export_kd_vs_residue(str(out), obs) == 1
        assert out.exists() and out.stat().st_size > 0


def _write_json_with_global(tmp_path):
    data = json.loads(_write_json(tmp_path).read_text())
    data["global"] = {"intensity": {"success": True, "Kd": 9.0,
                                     "I0": {"A17": 1.0, "K14": 1.0},
                                     "I_inf": {"A17": 0.2, "K14": 0.3},
                                     "n_residues": 2},
                      "csp": {"success": True, "Kd": 11.0,
                              "dd_max": {"A17": 0.1, "K14": 0.08},
                              "n_residues": 2}}
    jf = tmp_path / "GLOB_kd_fit_data.json"
    jf.write_text(json.dumps(data))
    return jf


def test_viewer_exports_global_fit_pdf(app, tmp_path):
    from visualization.kd_titration_fit_viewer import KdTitrationFitViewer
    viewer = KdTitrationFitViewer(json_file=str(_write_json_with_global(tmp_path)))
    for obs in ("csp", "intensity"):
        out = tmp_path / f"GLOB_{obs}_global_fit.pdf"
        assert viewer._export_global_fit(str(out), obs) >= 1
        assert out.exists() and out.stat().st_size > 0


def test_global_label_tracks_observable(app, tmp_path):
    from visualization.kd_titration_fit_viewer import KdTitrationFitViewer
    viewer = KdTitrationFitViewer(json_file=str(_write_json_with_global(tmp_path)))
    viewer.obs_combo.setCurrentIndex(0)          # CSP
    viewer._update_global_label()
    assert "CSP" in viewer.global_label.text()
    assert "11" in viewer.global_label.text()    # csp global Kd = 11.0
    viewer.obs_combo.setCurrentIndex(1)          # Intensity
    viewer._update_global_label()
    txt = viewer.global_label.text()
    assert "intensity" in txt.lower() and "CSP" not in txt
    assert "9" in txt                            # intensity global Kd = 9.0
    assert "avg" in txt.lower()                  # per-residue average also shown


def test_mean_kd_helper(app, tmp_path):
    from visualization.kd_titration_fit_viewer import KdTitrationFitViewer
    viewer = KdTitrationFitViewer(json_file=str(_write_json_with_global(tmp_path)))
    mean, n = viewer._mean_kd("intensity")       # both residues have per-residue Kd=10
    assert mean == pytest.approx(10.0) and n == 2
    viewer.excluded_residues["intensity"].add("A17")   # exclude one -> average over the rest
    mean2, n2 = viewer._mean_kd("intensity")
    assert n2 == 1


def test_mean_kd_skips_low_r2_outliers(app, tmp_path):
    from visualization.kd_titration_fit_viewer import KdTitrationFitViewer
    data = json.loads(_write_json_with_global(tmp_path).read_text())
    # inject a low-R² residue with a wildly off Kd — must NOT enter the mean
    data["fits"].append({"residue": "BAD99",
                         "series": {"ppm_x": [8.0], "ppm_y": [120.0], "height": [1.0],
                                    "volume": [1.0]},
                         "intensity": {"success": True, "Kd": 900.0, "Kd_err": 500.0,
                                       "I0": 1.0, "I_inf": 0.5, "r_squared": 0.2,
                                       "L": [0.0, 6.25, 12.5], "obs": [1.0, 0.9, 0.8]}})
    jf = tmp_path / "OUT_kd_fit_data.json"
    jf.write_text(json.dumps(data))
    viewer = KdTitrationFitViewer(json_file=str(jf))
    mean, n = viewer._mean_kd("intensity")
    assert n == 2 and mean == pytest.approx(10.0)     # BAD99 (R²=0.2) excluded, not ~306


def test_mean_kd_hampel_rejects_high_r2_outlier(app, tmp_path):
    from visualization.kd_titration_fit_viewer import KdTitrationFitViewer
    data = json.loads(_write_json_with_global(tmp_path).read_text())
    L = [0.0, 6.25, 12.5]
    # 4 well-fit residues near Kd=10, plus one well-fit (R²=0.97) but Kd off by 30x
    for i, kd in enumerate([10.1, 9.8, 10.2, 900.0]):
        data["fits"].append({"residue": f"R{i}",
                             "series": {"ppm_x": [8.0], "ppm_y": [120.0],
                                        "height": [1.0], "volume": [1.0]},
                             "intensity": {"success": True, "Kd": kd, "Kd_err": 0.5,
                                           "I0": 1.0, "I_inf": 0.2, "r_squared": 0.97,
                                           "L": L, "obs": [1.0, 0.6, 0.4]}})
    jf = tmp_path / "HAMP_kd_fit_data.json"
    jf.write_text(json.dumps(data))
    viewer = KdTitrationFitViewer(json_file=str(jf))
    kept = viewer._clean_kds("intensity")
    assert 900.0 not in kept                          # high-R² but 30x median -> dropped
    assert all(k < 20 for k in kept)


def test_viewer_global_fit_absent_writes_nothing(app, tmp_path):
    from visualization.kd_titration_fit_viewer import KdTitrationFitViewer
    viewer = KdTitrationFitViewer(json_file=str(_write_json(tmp_path)))   # global == {}
    out = tmp_path / "GLOB_intensity_global_fit.pdf"
    assert viewer._export_global_fit(str(out), "intensity") == 0
    assert not out.exists()


def test_viewer_kd_vs_residue_no_fits_writes_nothing(app, tmp_path):
    from visualization.kd_titration_fit_viewer import KdTitrationFitViewer
    data = json.loads(_write_json(tmp_path).read_text())
    for f in data["fits"]:                       # strip all intensity fits
        f["intensity"] = {"success": False}
    jf = tmp_path / "NOINT_kd_fit_data.json"
    jf.write_text(json.dumps(data))
    viewer = KdTitrationFitViewer(json_file=str(jf))
    out = tmp_path / "NOINT_intensity_kd_vs_residue.pdf"
    assert viewer._export_kd_vs_residue(str(out), "intensity") == 0
    assert not out.exists()
