# ABOUTME: A Fit Viewer refit rewrites the results file of the run that wrote the fit JSON,
# ABOUTME: for a CLI `dynamixs t1t2` run and a T1/T2 page run alike, and inside a bundle.

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
    for line in Path(results_file).read_text().splitlines()[1:]:
        if line.split("\t", 1)[0] == residue:
            return line
    raise KeyError(residue)


def _refit_rewrites_the_results_file(json_file, results_file):
    """Exclude the first delay of the first residue, refit it in the viewer, and return
    (row before, row after, the refit's T1) from `results_file`."""
    from visualization.fit_viewer import FitViewer
    viewer = FitViewer(json_files=[str(json_file)])
    try:
        key = "field1_T1"
        residue = viewer.data[key]["fits"][0]["residue"]
        before = _row(results_file, residue)
        viewer.exclusions.setdefault(key, {})[residue] = [0]
        entry = viewer._refit_one(key, residue)
        return before, _row(results_file, residue), entry["t2"]
    finally:
        viewer.close()


def test_a_cli_t1t2_runs_results_file_takes_the_refit(app, tmp_path, quiet_viewers):
    from lunaNMR.cli import main
    out = tmp_path / "cli"
    assert main(["dynamixs", "t1t2", "--input", str(MATRIX), "--out", str(out),
                 "--exp", "T1", "--field-freq", "600", "--prefix", "P"]) == 0
    json_file = out / "P_field1_T1_fit_data.json"
    results = out / "P_fit_results.txt"
    before, after, t1 = _refit_rewrites_the_results_file(json_file, results)
    assert after != before
    assert float(after.split("\t")[2]) == pytest.approx(t1, rel=1e-6)


def test_a_t1t2_page_runs_results_file_takes_the_refit(app, tmp_path, quiet_viewers):
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
        before, after, t1 = _refit_rewrites_the_results_file(entry["json_file"],
                                                             entry["results_file"])
        assert after != before
        assert float(after.split("\t")[2]) == pytest.approx(t1, rel=1e-6)
    finally:
        d.deleteLater()


def test_a_built_bundles_results_file_takes_the_refit(app, tmp_path, quiet_viewers):
    """`project build` stages a CLI run's JSON and results file under new names."""
    from lunaNMR.cli import main
    run = tmp_path / "run"
    assert main(["dynamixs", "t1t2", "--input", str(MATRIX), "--out",
                 str(run / "fits" / "600_T1"), "--exp", "T1", "--field-freq", "600",
                 "--prefix", "600_T1"]) == 0
    bundle = tmp_path / "B.lunaNMR"
    assert main(["project", "build", str(run), "--out", str(bundle)]) == 0
    from lunaNMR.utils.project_manager import ProjectManager
    window = types.SimpleNamespace(dynamixs_dialog=None, dynamixs_analyses=None)
    assert ProjectManager(window)._load_dynamixs_state(bundle)
    (analysis,) = window.dynamixs_analyses.values()
    entry = analysis["state"]["t1t2"]["session_results"]["field1"]["t1"]
    assert Path(entry["json_file"]).is_relative_to(bundle)
    before, after, t1 = _refit_rewrites_the_results_file(entry["json_file"],
                                                         entry["results_file"])
    assert after != before
    assert float(after.split("\t")[2]) == pytest.approx(t1, rel=1e-6)
