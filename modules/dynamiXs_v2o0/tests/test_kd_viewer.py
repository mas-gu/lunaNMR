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


def test_comparison_shows_missing_residue_as_grey_bar(app, tmp_path):
    import json
    from visualization.kd_titration_fit_viewer import open_kd_titration_viewer
    jf = _make_fit_json(tmp_path)
    data = json.loads(Path(jf).read_text())
    good = dict(data["fits"][0], residue="K14")     # has embedded series
    bad = {"residue": "A17"}                          # no series -> missing
    data["fits"] = [good, bad]
    Path(jf).write_text(json.dumps(data))
    v = open_kd_titration_viewer(parent=None, json_file=jf)
    try:
        v.view_combo.setCurrentIndex(3)               # "Ref vs point (bars)"
        bars = v.figure.axes[0].patches
        assert len(bars) == 2                          # both shown; missing one is grey
        # the missing residue's bar is grey, the good one is not
        colors = {tuple(round(c, 2) for c in b.get_facecolor()) for b in bars}
        grey = tuple(round(c, 2) for c in __import__("matplotlib").colors.to_rgba("#d0d0d0"))
        assert grey in colors
    finally:
        v.deleteLater()


def test_exclude_residue_via_bar_click_persists_and_reloads(app, tmp_path):
    import json
    from visualization.kd_titration_fit_viewer import open_kd_titration_viewer
    jf = _make_fit_json(tmp_path)
    v = open_kd_titration_viewer(parent=None, json_file=jf)
    try:
        v.view_combo.setCurrentIndex(3)          # "Ref vs point (bars)"
        v._toggle_edit(True)                      # edit mode -> bars pickable
        patch = v.figure.axes[0].patches[0]
        resid = v._pick_registry[id(patch)][0]    # residue behind that bar

        class _Evt:
            pass
        v.obs_combo.setCurrentIndex(0)            # CSP
        e = _Evt(); e.artist = patch; e.ind = [0]
        v._on_pick(e)                             # click the bar -> exclude residue (CSP)
        assert resid in v.excluded_residues["csp"]
        saved = json.loads(Path(jf).read_text()).get("excluded_residues", {})
        assert resid in saved["csp"]

        # reload restores the exclusion
        v2 = open_kd_titration_viewer(parent=None, json_file=jf)
        try:
            assert resid in v2.excluded_residues["csp"]
        finally:
            v2.deleteLater()

        # clicking again un-excludes
        v._toggle_edit(True)
        patch2 = next(p for p in v.figure.axes[0].patches
                      if v._pick_registry.get(id(p), (None,))[0] == resid)
        e2 = _Evt(); e2.artist = patch2; e2.ind = [0]
        v._on_pick(e2)
        assert resid not in v.excluded_residues["csp"]
    finally:
        v.deleteLater()


def test_residue_exclusion_is_independent_per_observable(app, tmp_path):
    # Excluding a residue in CSP must NOT exclude it in Intensity, and vice versa.
    from visualization.kd_titration_fit_viewer import open_kd_titration_viewer
    jf = _make_fit_json(tmp_path)
    v = open_kd_titration_viewer(parent=None, json_file=jf)
    try:
        v.view_combo.setCurrentIndex(3)
        v.obs_combo.setCurrentIndex(0)            # CSP
        v._toggle_edit(True)
        patch = v.figure.axes[0].patches[0]
        resid = v._pick_registry[id(patch)][0]

        class _Evt:
            pass
        e = _Evt(); e.artist = patch; e.ind = [0]
        v._on_pick(e)                             # exclude in CSP
        assert resid in v.excluded_residues["csp"]
        assert resid not in v.excluded_residues["intensity"]   # independent
    finally:
        v.deleteLater()


def test_missing_grey_is_consistent_across_csp_and_intensity(app, tmp_path):
    # A residue with valid positions but a vanished peak (zero intensity) must be
    # treated as missing (grey) in BOTH CSP and Intensity, not just intensity.
    import json
    import matplotlib
    from visualization.kd_titration_fit_viewer import open_kd_titration_viewer
    jf = _make_fit_json(tmp_path)
    data = json.loads(Path(jf).read_text())
    good = dict(data["fits"][0], residue="K14")
    gone = {"residue": "A17", "series": {
        "ppm_x": [8.0, 8.0, 8.0], "ppm_y": [120.0, 120.0, 120.0],
        "height": [1000.0, 0.0, 0.0], "volume": [2000.0, 0.0, 0.0]}}
    data["fits"] = [good, gone]
    Path(jf).write_text(json.dumps(data))
    grey = tuple(round(c, 2) for c in matplotlib.colors.to_rgba("#d0d0d0"))

    v = open_kd_titration_viewer(parent=None, json_file=jf)
    try:
        v.view_combo.setCurrentIndex(3)
        for obs_idx in (0, 1):                       # 0=CSP, 1=Intensity
            v.obs_combo.setCurrentIndex(obs_idx)
            cols = [tuple(round(c, 2) for c in b.get_facecolor())
                    for b in v.figure.axes[0].patches]
            assert grey in cols, f"missing residue not grey for obs index {obs_idx}"
    finally:
        v.deleteLater()


def test_intensity_comparison_yaxis_is_zero_to_one(app, tmp_path):
    # I/I₀ is a ratio → the intensity comparison y-axis stays 0..1 (not data-scaled),
    # while CSP (ppm) auto-scales.
    from visualization.kd_titration_fit_viewer import open_kd_titration_viewer
    jf = _make_fit_json(tmp_path)
    v = open_kd_titration_viewer(parent=None, json_file=jf)
    try:
        v.view_combo.setCurrentIndex(3)
        v.obs_combo.setCurrentIndex(1)            # Intensity
        top = v.figure.axes[0].get_ylim()[1]
        assert top == pytest.approx(1.0, abs=1e-9)
        v.obs_combo.setCurrentIndex(0)            # CSP → auto-scaled (well below 1 ppm)
        assert v.figure.axes[0].get_ylim()[1] < 0.5
    finally:
        v.deleteLater()


def test_export_intensity_fits_writes_vector_pdf(app, tmp_path):
    from visualization.kd_titration_fit_viewer import open_kd_titration_viewer
    jf = _make_fit_json(tmp_path)
    v = open_kd_titration_viewer(parent=None, json_file=jf)
    try:
        out = tmp_path / "int_fits.pdf"
        assert v._export_intensity_fits(str(out)) >= 1   # ≥1 page of per-residue curves
        data = out.read_bytes()
        assert data[:4] == b"%PDF"
        # pdf.fonttype=42 embeds TrueType (FontFile2) → text stays selectable/editable
        assert b"/FontFile2" in data
    finally:
        v.deleteLater()


def test_export_ref_vs_point_page_per_point(app, tmp_path):
    from visualization.kd_titration_fit_viewer import open_kd_titration_viewer
    jf = _make_fit_json(tmp_path)                        # 6 titration points
    v = open_kd_titration_viewer(parent=None, json_file=jf)
    try:
        for obs, fname in (("csp", "csp_ref.pdf"), ("intensity", "int_ref.pdf")):
            out = tmp_path / fname
            assert v._export_ref_vs_point(str(out), obs) == 5   # ref(0) → points 1..5
            assert out.read_bytes()[:4] == b"%PDF"
    finally:
        v.deleteLater()


def test_export_ref_vs_point_skips_when_too_few_points(app, tmp_path):
    from visualization.kd_titration_fit_viewer import open_kd_titration_viewer
    jf = _make_fit_json(tmp_path)
    v = open_kd_titration_viewer(parent=None, json_file=jf)
    try:
        v.data["metadata"]["concentrations"] = [0.0]     # only one point → nothing to compare
        out = tmp_path / "none.pdf"
        assert v._export_ref_vs_point(str(out), "csp") == 0
        assert not out.exists()
    finally:
        v.deleteLater()


def test_export_all_writes_three_pdfs_next_to_json(app, tmp_path, monkeypatch):
    from visualization import kd_titration_fit_viewer as mod
    jf = _make_fit_json(tmp_path)
    v = mod.open_kd_titration_viewer(parent=None, json_file=jf)
    try:
        # accept the choice dialog with all three (default-checked) selected
        monkeypatch.setattr(mod._ExportChoiceDialog, "exec", lambda self: 1)
        monkeypatch.setattr(mod.QMessageBox, "information", lambda *a, **k: None)
        v._export()
        dest = Path(jf).parent
        assert (dest / "kd_intensity_titration_fits.pdf").exists()
        assert (dest / "kd_intensity_ref_vs_point.pdf").exists()
        assert (dest / "kd_csp_ref_vs_point.pdf").exists()
    finally:
        v.deleteLater()


def test_export_base_falls_back_to_metadata_name_for_generic_filename(app, tmp_path, monkeypatch):
    # A project-bundled fit is stored as generic 'fit_data.json'; the export prefix then
    # comes from the series name recorded in metadata, not the filename.
    import json
    from visualization import kd_titration_fit_viewer as mod
    jf = _make_fit_json(tmp_path)
    data = json.loads(Path(jf).read_text())
    data["metadata"]["name"] = "HSPA8"
    generic = Path(jf).parent / "fit_data.json"
    generic.write_text(json.dumps(data))
    v = mod.open_kd_titration_viewer(parent=None, json_file=str(generic))
    try:
        monkeypatch.setattr(mod._ExportChoiceDialog, "exec", lambda self: 1)
        monkeypatch.setattr(mod.QMessageBox, "information", lambda *a, **k: None)
        v._export()
        assert (generic.parent / "HSPA8_intensity_titration_fits.pdf").exists()
        assert (generic.parent / "HSPA8_csp_ref_vs_point.pdf").exists()
    finally:
        v.deleteLater()


def test_fit_json_metadata_records_name(tmp_path):
    # run_kd_analysis_with_params stamps the output prefix into metadata['name'] so a
    # bundled reopen can still recover the series name.
    import json
    from kd_fit import run_kd_analysis_with_params
    from kd_models import csp_model
    pts = [0.0, 10.0, 25.0, 60.0]
    rows = []
    for name in ("A17", "K14"):
        d = csp_model(np.array(pts), 0.2, 15.0, 50.0)
        for p, dv in zip(pts, d):
            rows.append((str(p), name, 8.0 + dv, 120.0, 1000.0 - 2 * p, 2000.0))
    df = pd.DataFrame(rows, columns=["spectrum_name", "assignment",
                                     "ppm_x", "ppm_y", "height", "volume"])
    csv = tmp_path / "series_analysis_tidy.csv"
    df.to_csv(csv, index=False)
    r = run_kd_analysis_with_params({
        "input_csv_file": str(csv), "output_dir": str(tmp_path), "output_prefix": "HSPA8",
        "concentrations": pts, "protein_conc": 50.0, "alpha": 0.14,
        "observables": ["csp"], "n_bootstrap": 0})
    meta = json.loads(Path(r["json_file"]).read_text())["metadata"]
    assert meta["name"] == "HSPA8"


def test_export_filenames_carry_series_prefix(app, tmp_path, monkeypatch):
    # A fit saved as '<series>_kd_fit_data.json' exports PDFs prefixed with the series name.
    from visualization import kd_titration_fit_viewer as mod
    jf = _make_fit_json(tmp_path)
    named = Path(jf).parent / "HSPA8_kd_fit_data.json"
    Path(jf).rename(named)
    v = mod.open_kd_titration_viewer(parent=None, json_file=str(named))
    try:
        monkeypatch.setattr(mod._ExportChoiceDialog, "exec", lambda self: 1)
        monkeypatch.setattr(mod.QMessageBox, "information", lambda *a, **k: None)
        v._export()
        dest = named.parent
        assert (dest / "HSPA8_intensity_titration_fits.pdf").exists()
        assert (dest / "HSPA8_intensity_ref_vs_point.pdf").exists()
        assert (dest / "HSPA8_csp_ref_vs_point.pdf").exists()
    finally:
        v.deleteLater()


# The CLI's `export kd` (lunaNMR/cli.py) uses a 4-col x 5-row grid with figsize=(16, 15)
# inches (1152x1080 pt). These tests pin the GUI viewer's multi-panel PDFs to the exact
# same page size, so figures from either surface can be directly overlaid/compared.
_CLI_MEDIABOX = b"/MediaBox [ 0 0 1152 1080 ]"


def test_export_intensity_fits_page_size_matches_cli(app, tmp_path):
    from visualization.kd_titration_fit_viewer import open_kd_titration_viewer
    jf = _make_fit_json(tmp_path)
    v = open_kd_titration_viewer(parent=None, json_file=jf)
    try:
        out = tmp_path / "int_fits.pdf"
        assert v._export_intensity_fits(str(out)) >= 1
        assert _CLI_MEDIABOX in out.read_bytes()
    finally:
        v.deleteLater()


def test_export_global_fit_page_size_matches_cli(app, tmp_path):
    from visualization.kd_titration_fit_viewer import open_kd_titration_viewer
    jf = _make_fit_json(tmp_path)
    v = open_kd_titration_viewer(parent=None, json_file=jf)
    try:
        out = tmp_path / "global_fit.pdf"
        assert v._export_global_fit(str(out), "intensity") >= 1
        assert _CLI_MEDIABOX in out.read_bytes()
    finally:
        v.deleteLater()


def test_draw_intensity_fit_normalizes_raw_scale_and_shares_ylim(app, tmp_path):
    # Same fix as the CLI's _build_kd_panel/_INTENSITY_YLIM: a raw-scale residue
    # (I0 in the hundreds) must plot as I/I(0) on the shared 0-1 axis, not auto-scale
    # to its own raw magnitude.
    from visualization.kd_titration_fit_viewer import open_kd_titration_viewer, _INTENSITY_YLIM
    jf = _make_fit_json(tmp_path)
    v = open_kd_titration_viewer(parent=None, json_file=jf)
    try:
        ax = v.figure.axes[0] if v.figure.axes else v.figure.add_subplot(111)
        ax.clear()
        f = {'residue': 'RAW', 'intensity': {
            'L': np.array([0.0, 10.0, 50.0]), 'obs': np.array([500.0, 300.0, 100.0]),
            'I0': 500.0, 'I_inf': 50.0, 'Kd': 40.0, 'Kd_err': 1.0, 'r_squared': 0.9}}
        v._draw_intensity_fit(ax, f)
        assert ax.get_ylim() == _INTENSITY_YLIM
    finally:
        v.deleteLater()


def test_draw_global_fit_panel_normalizes_intensity_and_shares_ylim(app, tmp_path):
    from visualization.kd_titration_fit_viewer import open_kd_titration_viewer, _INTENSITY_YLIM
    jf = _make_fit_json(tmp_path)
    v = open_kd_titration_viewer(parent=None, json_file=jf)
    try:
        ax = v.figure.axes[0] if v.figure.axes else v.figure.add_subplot(111)
        ax.clear()
        res = v.fits[0]['residue']
        f = {'residue': res, 'intensity': {
            'L': np.array([0.0, 10.0, 50.0]), 'obs': np.array([500.0, 300.0, 100.0])}}
        g = {'Kd': 40.0, 'I0': {res: 500.0}, 'I_inf': {res: 50.0}}
        v._draw_global_fit_panel(ax, f, 'intensity', g)
        assert ax.get_ylim() == _INTENSITY_YLIM
    finally:
        v.deleteLater()


def test_plot_curve_intensity_shares_ylim_with_exports(app, tmp_path):
    # The interactive single-residue view (shown before export) must match what gets
    # exported: I/I(0) on the same shared 0-1 axis, not raw-scale auto-scaling.
    from visualization.kd_titration_fit_viewer import open_kd_titration_viewer, _INTENSITY_YLIM
    jf = _make_fit_json(tmp_path)
    v = open_kd_titration_viewer(parent=None, json_file=jf)
    try:
        v.view_combo.setCurrentIndex(0)   # per-residue curve view
        v.obs_combo.setCurrentIndex(1)    # Intensity
        v._refresh()
        assert v.figure.axes[0].get_ylim() == _INTENSITY_YLIM
    finally:
        v.deleteLater()


def test_draw_global_fit_panel_leaves_csp_autoscaled(app, tmp_path):
    from visualization.kd_titration_fit_viewer import open_kd_titration_viewer, _INTENSITY_YLIM
    jf = _make_fit_json(tmp_path)
    v = open_kd_titration_viewer(parent=None, json_file=jf)
    try:
        ax = v.figure.axes[0] if v.figure.axes else v.figure.add_subplot(111)
        ax.clear()
        res = v.fits[0]['residue']
        f = {'residue': res, 'csp': {
            'L': np.array([0.0, 10.0, 50.0]), 'obs': np.array([0.0, 0.1, 0.2])}}
        g = v.data['global']['csp']
        v._draw_global_fit_panel(ax, f, 'csp', g)
        assert ax.get_ylim() != _INTENSITY_YLIM
    finally:
        v.deleteLater()
