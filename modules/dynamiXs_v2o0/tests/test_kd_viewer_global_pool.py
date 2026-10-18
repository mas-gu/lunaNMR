# ABOUTME: The Kd viewer's Kd-vs-residue figure tells the global fit's pool from the rest:
# ABOUTME: export plots the pool only (as `export kd` does), the view marks the others. Real fit.

import contextlib
import io
import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_MOD = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_MOD))
sys.path.insert(0, str(_MOD / "dynamiXs_Kd"))

DNAJA1 = (Path(__file__).resolve().parents[4] / "data_test" / "data_binding" / "tidy"
          / "DNAJA1_HSPA8_tidy.csv")

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


@pytest.fixture(scope="module")
def fit_json(tmp_path_factory):
    """DNAJA1 fitted in equivalents by `kd`: a shared CSP Kd over a pool far smaller
    than the residues fitted."""
    from lunaNMR.cli import main
    folder = tmp_path_factory.mktemp("dnaja1")
    csv = folder / DNAJA1.name
    shutil.copy(DNAJA1, csv)
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        code = main(["kd", "--format", "json", "--p0", "50", "--conc-units", "equivalents",
                     "--input", str(csv), "--out", str(folder / "kd")])
    assert code == 0, out.getvalue()
    return next((folder / "kd" / "data").glob("*_kd_fit_data.json"))


def _viewer(fit_json):
    from visualization.kd_titration_fit_viewer import KdTitrationFitViewer
    return KdTitrationFitViewer(json_file=str(fit_json))


def _pool(fit_json, obs):
    return set(json.loads(fit_json.read_text())["global"][obs]["residues"])


def _export_axes(viewer, obs):
    from matplotlib.figure import Figure
    ax = Figure().add_subplot(111)
    viewer._plot_kd_vs_residue_export(ax, obs)
    return ax


def test_the_pool_is_smaller_than_the_fit(app, fit_json):
    data = json.loads(fit_json.read_text())
    assert len(_pool(fit_json, "csp")) < len(data["fits"])


@pytest.mark.parametrize("obs", ["csp", "intensity"])
def test_the_export_plots_only_the_pool(app, fit_json, obs):
    ax = _export_axes(_viewer(fit_json), obs)
    names = {t.get_text() for t in ax.get_xticklabels()}
    assert names == _pool(fit_json, obs)


def test_the_export_says_how_many_of_how_many(app, fit_json):
    data = json.loads(fit_json.read_text())
    ax = _export_axes(_viewer(fit_json), "csp")
    assert f"{len(_pool(fit_json, 'csp'))} of {len(data['fits'])} residues" in ax.get_title()


def test_the_export_writes_its_pdf(app, fit_json, tmp_path):
    out = tmp_path / "csp_kd_vs_residue.pdf"
    assert _viewer(fit_json)._export_kd_vs_residue(str(out), "csp") == 1
    assert out.stat().st_size > 0


def test_the_view_keeps_every_residue_and_marks_the_rest(app, fit_json):
    data = json.loads(fit_json.read_text())
    viewer = _viewer(fit_json)
    viewer.obs_combo.setCurrentIndex(0)
    viewer.view_combo.setCurrentIndex(1)                  # Kd vs residue
    ax = viewer.figure.axes[0]
    names = [t.get_text() for t in ax.get_xticklabels()]
    assert len(names) == len(data["fits"])

    pool = _pool(fit_json, "csp")
    finite = {f["residue"] for f in data["fits"]
              if f["csp"].get("success") and np.isfinite(f["csp"].get("Kd") or np.nan)}
    hollow = {n for n, bar in zip(names, ax.patches[:len(names)])
              if bar.get_facecolor()[3] == 0}
    assert hollow == finite - pool
    assert any("pool" in t.get_text() for t in ax.get_legend().get_texts())


def test_the_average_is_over_the_pool(app, fit_json):
    data = json.loads(fit_json.read_text())
    pool = _pool(fit_json, "csp")
    expected = [f["csp"]["Kd"] for f in data["fits"] if f["residue"] in pool]
    mean, n = _viewer(fit_json)._mean_kd("csp")
    assert n == len(expected)
    assert mean == pytest.approx(float(np.mean(expected)))


def test_a_user_exclusion_still_leaves_the_average(app, fit_json):
    viewer = _viewer(fit_json)
    viewer.excluded_residues["csp"].add(sorted(_pool(fit_json, "csp"))[0])
    assert viewer._mean_kd("csp")[1] == len(_pool(fit_json, "csp")) - 1
