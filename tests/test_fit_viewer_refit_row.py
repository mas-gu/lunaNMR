# ABOUTME: A Fit Viewer refit writes a full results-file row and a full JSON entry, with its own
# ABOUTME: Success, WindowRatio and window flags over the points it kept: CLI and page T1/T2 runs.

import json
import math
import os
import sys
import types
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

SERIES = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
          / "run_20260825_1405" / "series")
MATRIX = SERIES / "600_T1" / "peak_intensity_matrix.csv"

pytestmark = pytest.mark.skipif(not MATRIX.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def quiet_viewers(monkeypatch):
    """The fit viewer reports a refit or a load problem in a message box; nothing closes one offscreen."""
    import lunaNMR.gui.dialogs.dynamixs_dialog  # noqa: F401 - puts the module on sys.path
    import visualization.fit_viewer as fit_viewer
    for box in ("show_info", "show_warning", "show_error"):
        monkeypatch.setattr(fit_viewer, box, lambda *a, **k: None)


def _row(results_file, residue):
    lines = Path(results_file).read_text().splitlines()
    header = lines[0].split("\t")
    for line in lines[1:]:
        cells = line.split("\t")
        if cells[0] == residue:
            assert len(cells) == len(header), (header, cells)
            return dict(zip(header, cells))
    raise KeyError(residue)


def _refit_without_the_last_delay(json_file, results_file):
    """Refit the first residue in the viewer with its longest delay excluded; return the
    results-file row, the JSON entry on disk, an entry left alone, and the delays kept."""
    from visualization.fit_viewer import FitViewer
    viewer = FitViewer(json_files=[str(json_file)])
    try:
        key = "field1_T1"
        fits = viewer.data[key]["fits"]
        residue = fits[0]["residue"]
        time_points = viewer.data[key]["metadata"]["time_points"]
        last = max(range(len(time_points)), key=time_points.__getitem__)
        viewer.exclusions.setdefault(key, {})[residue] = [last]
        viewer._refit_one(key, residue)
        on_disk = json.loads(Path(json_file).read_text())["fits"]
        kept = [t for i, t in enumerate(time_points) if i != last]
        return _row(results_file, residue), on_disk[0], on_disk[1], kept
    finally:
        viewer.close()


def _assert_the_refit_reports_itself(row, entry, untouched, kept):
    from fit_Tx_NMRRE import MARGINAL_TMAX_OVER_T
    assert set(entry) == set(untouched)
    ratio = max(kept) / entry["t2"]
    assert entry["window_ratio"] == pytest.approx(ratio)
    assert entry["window_marginal"] is (ratio < MARGINAL_TMAX_OVER_T)
    assert entry["baseline_fixed"] is True
    assert isinstance(entry["success"], bool)
    assert row["Success"] == ("Yes" if entry["success"] else "No")
    assert row["WindowRatio"] == f"{ratio:.3f}"
    assert float(row["T1"]) == pytest.approx(entry["t2"], rel=1e-6)


def test_a_cli_runs_refit_writes_a_full_row_and_entry(app, tmp_path, quiet_viewers):
    from lunaNMR.cli import main
    out = tmp_path / "cli"
    assert main(["dynamixs", "t1t2", "--input", str(MATRIX), "--out", str(out),
                 "--exp", "T1", "--field-freq", "600", "--prefix", "P"]) == 0
    _assert_the_refit_reports_itself(*_refit_without_the_last_delay(
        out / "P_field1_T1_fit_data.json", out / "P_fit_results.txt"))


@pytest.fixture
def page_run(app, tmp_path):
    """A T1 fit by the real T1/T2 page, which fits on every core: (json, results)."""
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    d = DynamiXsDialog(parent=None, main_window=types.SimpleNamespace())
    try:
        page = d.t1t2_page
        page.output_dir = str(tmp_path)
        page.field1_t1_file = str(MATRIX)
        page.analysis_combo.setCurrentIndex(0)          # Field 1 - T1
        page.prefix_edit.setText("G")
        page._start_analysis()
        page.worker.wait()
        app.processEvents()
        entry = page.session_results["field1"]["t1"]
        yield entry["json_file"], entry["results_file"]
    finally:
        d.deleteLater()


def test_a_page_run_reports_each_fits_window(page_run):
    json_file, results_file = page_run
    fits = json.loads(Path(json_file).read_text())["fits"]
    assert all("success" in f for f in fits)
    assert all(math.isfinite(f["window_ratio"]) for f in fits)
    row = _row(results_file, fits[0]["residue"])
    assert row["WindowRatio"] == f"{fits[0]['window_ratio']:.3f}"


def test_a_page_runs_refit_writes_a_full_row_and_entry(page_run, quiet_viewers):
    _assert_the_refit_reports_itself(*_refit_without_the_last_delay(*page_run))
