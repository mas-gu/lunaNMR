# ABOUTME: A methyl Fit Viewer refit rewrites the results file the fit JSON records, for a CLI
# ABOUTME: `dynamixs methyl-t2` run, a methyl T2 page run and a built bundle, with a full row.

import json
import os
import sys
import types
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

METHYL = (Path(__file__).resolve().parents[2] / "data_test" / "methyl" / "sampleNBD"
          / "series_results_sampleNBD_20260818_124619" / "peak_intensity_matrix.csv")

pytestmark = pytest.mark.skipif(not METHYL.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def quiet_viewer(monkeypatch):
    """The methyl viewer reports each load in a message box; nothing closes one offscreen."""
    import lunaNMR.gui.dialogs.dynamixs_dialog  # noqa: F401 - puts the module on sys.path
    import visualization.methyl_t2_fit_viewer as methyl_viewer
    for box in ("show_info", "show_warning", "show_error"):
        monkeypatch.setattr(methyl_viewer, box, lambda *a, **k: None)


def _row(results_file, residue):
    lines = Path(results_file).read_text().splitlines()
    header = lines[0].split("\t")
    for line in lines[1:]:
        cells = line.split("\t")
        if cells[0] == residue:
            assert len(cells) == len(header), (header, cells)
            return dict(zip(header, cells))
    raise KeyError(residue)


def _refit_rewrites_the_results_file(json_file, results_file):
    """Open the run's JSON in the viewer, check it finds `results_file`, refit the first
    residue without its longest delay, and check the row and entry it writes."""
    from visualization.methyl_t2_fit_viewer import MethylT2FitViewer
    viewer = MethylT2FitViewer(json_file=str(json_file))
    try:
        (key,) = viewer.data
        assert Path(viewer._find_tsv_for(key)).resolve() == Path(results_file).resolve()
        fits = viewer.data[key]["fits"]
        residue = fits[0]["residue"]
        time_points = viewer.data[key]["metadata"]["time_points"]
        last = max(range(len(time_points)), key=time_points.__getitem__)
        before = _row(results_file, residue)
        viewer.exclusions.setdefault(key, {})[residue] = [last]
        entry = viewer._refit_one(key, residue)
        after = _row(results_file, residue)
        on_disk = json.loads(Path(json_file).read_text())["fits"]
    finally:
        viewer.close()
    assert after != before
    assert float(after["T2a"]) == pytest.approx(entry["t2_a"], rel=1e-6)
    assert after["Success"] == "Yes"
    assert after["bi_exp_unidentifiable"] == ("Yes" if entry["bi_exp_unidentifiable"] else "No")
    assert set(on_disk[0]) == set(on_disk[1])


def test_a_cli_methyl_runs_results_file_takes_the_refit(app, tmp_path, quiet_viewer):
    from lunaNMR.cli import main
    out = tmp_path / "cli"
    assert main(["dynamixs", "methyl-t2", "--input", str(METHYL), "--out", str(out),
                 "--prefix", "M"]) == 0
    json_file = out / "M_field1_methylT2_fit_data.json"
    assert json.loads(json_file.read_text())["metadata"]["results_file"] == "M_fit_results.txt"
    _refit_rewrites_the_results_file(json_file, out / "M_fit_results.txt")


def test_a_methyl_page_runs_results_file_takes_the_refit(app, tmp_path, quiet_viewer):
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    d = DynamiXsDialog(parent=None, main_window=types.SimpleNamespace())
    try:
        page = d.methyl_t2_page
        page.output_dir = str(tmp_path)
        page._set_input_file(str(METHYL))
        page.prefix_edit.setText("G")
        page._start_analysis()
        page.worker.wait()
        app.processEvents()
        _refit_rewrites_the_results_file(page.last_json_file, page.last_results_file)
    finally:
        d.deleteLater()


def test_a_built_bundles_results_file_takes_the_refit(app, tmp_path, quiet_viewer):
    """`project build` stages a CLI run's JSON and results file inside the bundle."""
    from lunaNMR.cli import main
    run = tmp_path / "run"
    assert main(["dynamixs", "methyl-t2", "--input", str(METHYL),
                 "--out", str(run / "fits" / "methyl"), "--prefix", "M"]) == 0
    bundle = tmp_path / "B.lunaNMR"
    assert main(["project", "build", str(run), "--out", str(bundle)]) == 0
    from lunaNMR.utils.project_manager import ProjectManager
    window = types.SimpleNamespace(dynamixs_dialog=None, dynamixs_analyses=None)
    assert ProjectManager(window)._load_dynamixs_state(bundle)
    (analysis,) = window.dynamixs_analyses.values()
    page = analysis["state"]["methyl_t2"]
    (json_file,) = Path(page["last_json_folder"]).glob("*_methylT2_fit_data.json")
    assert json_file.is_relative_to(bundle)
    assert Path(page["last_results_file"]).is_relative_to(bundle)
    _refit_rewrites_the_results_file(json_file, page["last_results_file"])
