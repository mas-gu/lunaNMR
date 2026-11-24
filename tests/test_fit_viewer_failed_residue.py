# ABOUTME: The Fit Viewer lists a residue whose fit failed, marks it as the peak navigators mark a
# ABOUTME: failed fit, draws its points with no curve or value, and lets a refit make it a fitted residue.

import json
import os
import sys
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

MATRIX = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
          / "run_20260825_1405" / "series" / "700_T1" / "peak_intensity_matrix.csv")
# These two peaks of the real 700 MHz T1 series do not decay within the delays.
FAILED = ["10ValH", "43GluH"]
KEY = "field1_T1"

pytestmark = pytest.mark.skipif(not MATRIX.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def run(tmp_path):
    """A `dynamixs t1t2` fit of the real matrix: (fit JSON, results file)."""
    from lunaNMR.cli import main
    assert main(["dynamixs", "t1t2", "--input", str(MATRIX), "--out", str(tmp_path),
                 "--exp", "T1", "--field-freq", "700", "--prefix", "P"]) == 0
    return tmp_path / "P_field1_T1_fit_data.json", tmp_path / "P_fit_results.txt"


@pytest.fixture
def viewer(app, run, monkeypatch):
    import lunaNMR.gui.dialogs.dynamixs_dialog  # noqa: F401 - puts the module on sys.path
    import visualization.fit_viewer as fit_viewer
    for box in ("show_info", "show_warning", "show_error"):
        monkeypatch.setattr(fit_viewer, box, lambda *a, **k: None)
    viewer = fit_viewer.FitViewer(json_files=[str(run[0])])
    yield viewer
    viewer.close()
    viewer.deleteLater()


def _item(viewer, residue):
    from PySide6.QtCore import Qt
    items = [viewer.residue_list.item(i) for i in range(viewer.residue_list.count())]
    (found,) = [item for item in items if item.data(Qt.UserRole) == residue]
    return found


def _t1_panel(viewer):
    """The T1 panel: the first drawn, above the T2 one this run has no data for."""
    return viewer.axes[0]


def _show(viewer, residue):
    viewer.residue_list.setCurrentItem(_item(viewer, residue))
    return _t1_panel(viewer)


def _failed_mark():
    from lunaNMR.gui.styles.quality_colors import navigator_row_display
    return navigator_row_display("Failed")


def test_a_failed_residue_is_listed_with_the_navigators_failed_mark(viewer):
    symbol, background = _failed_mark()
    for residue in FAILED:
        item = _item(viewer, residue)
        assert item.text() == f"{symbol} {residue}"
        assert item.background().color() == background
        assert "not fitted" in item.toolTip()


def test_a_fitted_residue_is_listed_plain(viewer):
    fitted = next(f["residue"] for f in viewer.data[KEY]["fits"] if f["success"])
    item = _item(viewer, fitted)
    assert item.text() == fitted
    assert item.toolTip() == ""


def test_a_failed_residue_shows_its_points_and_no_fit(viewer):
    ax = _show(viewer, FAILED[0])
    n_points = len(viewer.data[KEY]["metadata"]["time_points"])
    assert sum(len(points.get_offsets()) for points in ax.collections) == n_points
    assert len(ax.lines) == 0
    texts = [t.get_text() for t in ax.texts]
    assert texts == ["T1 (700.0 MHz): not fitted"]


def test_a_fitted_residue_shows_its_curve_and_value(viewer):
    fit = next(f for f in viewer.data[KEY]["fits"] if f["success"])
    ax = _show(viewer, fit["residue"])
    assert len(ax.lines) == 1
    assert [t.get_text() for t in ax.texts] == [
        f"T1 (700.0 MHz) = {fit['t2']:.2f} ± {fit['t2_err']:.2f} ms"]


def test_inspecting_a_series_hands_over_no_fit_for_a_failed_residue(viewer):
    handed = viewer._get_all_fit_data(KEY)
    fitted = next(f for f in viewer.data[KEY]["fits"] if f["success"])
    assert handed[fitted["residue"]]["t_value"] == fitted["t2"]
    assert handed[fitted["residue"]]["fit_curve"]["time"]
    for residue in FAILED:
        assert handed[residue]["t_value"] is None and handed[residue]["t_error"] is None
        assert handed[residue]["fit_curve"] == {}
        assert handed[residue]["intensities"]


def test_a_refit_that_succeeds_makes_the_residue_a_fitted_one(viewer, run):
    json_file, results_file = run
    residue = FAILED[0]
    time_points = viewer.data[KEY]["metadata"]["time_points"]
    two_longest = sorted(range(len(time_points)), key=time_points.__getitem__)[-2:]
    _show(viewer, residue)
    viewer.exclusions.setdefault(KEY, {})[residue] = two_longest

    viewer._on_refit_clicked()

    entry = next(f for f in json.loads(json_file.read_text())["fits"] if f["residue"] == residue)
    assert entry["success"] is True
    item = _item(viewer, residue)
    assert item.text() == residue and item.toolTip() == ""
    assert item.background().color() != _failed_mark()[1]
    ax = _t1_panel(viewer)
    assert len(ax.lines) == 1
    assert [t.get_text() for t in ax.texts] == [
        f"T1 (700.0 MHz) = {entry['t2']:.2f} ± {entry['t2_err']:.2f} ms"]
    row = next(line.split("\t") for line in results_file.read_text().splitlines()
               if line.split("\t")[0] == residue)
    assert row[7] == "Yes" and float(row[2]) == pytest.approx(entry["t2"], rel=1e-6)
    assert _item(viewer, FAILED[1]).text() == f"{_failed_mark()[0]} {FAILED[1]}"


def test_a_refit_that_fails_leaves_the_residue_marked(viewer, run):
    _, results_file = run
    fitted = next(f["residue"] for f in viewer.data[KEY]["fits"] if f["success"])
    n_points = len(viewer.data[KEY]["metadata"]["time_points"])
    _show(viewer, fitted)
    viewer.exclusions.setdefault(KEY, {})[fitted] = list(range(2, n_points))

    viewer._on_refit_clicked()

    assert _item(viewer, fitted).text() == f"{_failed_mark()[0]} {fitted}"
    ax = _t1_panel(viewer)
    assert len(ax.lines) == 0
    row = next(line.split("\t") for line in results_file.read_text().splitlines()
               if line.split("\t")[0] == fitted)
    assert row[1:] == ["NaN"] * 6 + ["No", "NaN"]


def test_the_two_field_overlay_draws_only_the_field_that_fitted(app, run, tmp_path, monkeypatch):
    """The same residue at 600 MHz, where its peak decays, as the viewer's second field."""
    from lunaNMR.cli import main
    import visualization.fit_viewer as fit_viewer
    for box in ("show_info", "show_warning", "show_error"):
        monkeypatch.setattr(fit_viewer, box, lambda *a, **k: None)
    other = tmp_path / "600"
    assert main(["dynamixs", "t1t2", "--input", str(MATRIX.parents[1] / "600_T1" / MATRIX.name),
                 "--out", str(other), "--exp", "T1", "--field-freq", "600", "--prefix", "Q",
                 "--field-name", "field2"]) == 0
    viewer = fit_viewer.FitViewer(json_files=[str(run[0]),
                                              str(other / "Q_field2_T1_fit_data.json")])
    try:
        assert viewer.field_mode == "overlay"
        residue = FAILED[0]
        at_600 = next(f for f in viewer.data["field2_T1"]["fits"] if f["residue"] == residue)
        assert at_600["success"] is True

        ax = _show(viewer, residue)

        assert len(ax.lines) == 1
        (text,) = [t.get_text() for t in ax.texts]
        assert text.split("\n") == [
            "T1 (700.0 MHz): not fitted",
            f"T1 (600.0 MHz) = {at_600['t2']:.2f} ± {at_600['t2_err']:.2f} ms"]
        assert _item(viewer, residue).toolTip() == "T1 (700.0 MHz): not fitted"
    finally:
        viewer.close()
        viewer.deleteLater()
