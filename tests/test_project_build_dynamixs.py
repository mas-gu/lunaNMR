# ABOUTME: `project build` reopens a run's DynamiXs analyses (T1/T2 fits, model-free) in the GUI.
# ABOUTME: Real JA1 relaxation matrices -> CLI fits -> build -> real loader -> real DynamiXs dialog.

"""
The CLI fits each relaxation experiment separately, each into its own folder; the GUI's
T1/T2 page holds both fields and every experiment as ONE analysis. Each `dynamixs`
command records what it ran beside its outputs, and `project build` turns those
records into the page states the GUI saves. What matters is the far end: the dialog,
opened on the bundle, finds every result inside the bundle - not at the run folder,
which a bundle must not depend on.
"""

import json
import os
import shutil
import sys
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

RUN = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
       / "run_20260825_1405")
NOE = RUN / "results" / "hetnoe"

pytestmark = pytest.mark.skipif(not RUN.exists(), reason="data_test not present")


def _matrix(name):
    return str(RUN / "series" / name / "peak_intensity_matrix.csv")


@pytest.fixture(scope="module")
def run_dir(tmp_path_factory):
    """A relaxation run: four T1/T2 fits and a dual-field model-free."""
    from lunaNMR.cli import main
    run = tmp_path_factory.mktemp("JA1run")
    for name, exp, freq in (("600_T1", "T1", "600"), ("600_T2", "T2", "600"),
                            ("700_T1", "T1", "700"), ("700_T2", "T2", "700")):
        assert main(["dynamixs", "t1t2", "--input", _matrix(name), "--out",
                     str(run / "fits" / name), "--exp", exp, "--field-freq", freq,
                     "--prefix", name]) == 0
    assert main(["dynamixs", "modelfree", "--dual",
                 "--f1-t1", _matrix("600_T1"), "--f1-t2", _matrix("600_T2"),
                 "--f1-noe-sat", str(NOE / "f600_noe_sat.csv"),
                 "--f1-noe-unsat", str(NOE / "f600_noe_unsat.csv"), "--field1-freq", "600",
                 "--f2-t1", _matrix("700_T1"), "--f2-t2", _matrix("700_T2"),
                 "--f2-noe-sat", str(NOE / "f700_noe_sat.csv"),
                 "--f2-noe-unsat", str(NOE / "f700_noe_unsat.csv"), "--field2-freq", "700",
                 "--out", str(run / "modelfree" / "dual"), "--prefix", "dual"]) == 0
    return run


@pytest.fixture(scope="module")
def bundle(run_dir):
    """Built, then the run folder deleted: the bundle must stand on its own."""
    from lunaNMR.cli import main
    out = run_dir.parent / "JA1.lunaNMR"
    assert main(["project", "build", str(run_dir), "--out", str(out)]) == 0
    shutil.rmtree(run_dir / "fits")
    shutil.rmtree(run_dir / "modelfree")
    return out


@pytest.fixture(scope="module")
def loaded(bundle):
    import types
    from lunaNMR.utils.project_manager import ProjectManager
    window = types.SimpleNamespace(integrator=types.SimpleNamespace(peak_list=None))
    assert ProjectManager(window).load_project(bundle)
    return window


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


def _dialog(app, loaded, name):
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    dialog = DynamiXsDialog(parent=None, main_window=loaded)
    dialog.open_analysis(loaded.dynamixs_analyses[name], name)
    return dialog


def _inside(path, bundle):
    return Path(path).resolve().is_relative_to(Path(bundle).resolve())


T1T2 = "JA1run0_T1T2"
MODELFREE = "dual_modelfree"


def test_each_command_records_what_it_ran(run_dir):
    record = json.loads((run_dir / "fits" / "600_T1" / "600_T1_dynamixs.json").read_text())
    assert record["command"].startswith("dynamixs t1t2")
    assert Path(record["args"]["input"]).is_absolute()
    record = json.loads((run_dir / "modelfree" / "dual" / "dual_dynamixs.json").read_text())
    assert record["command"].startswith("dynamixs modelfree")


def test_the_four_fits_become_one_t1t2_analysis_and_modelfree_its_own(loaded):
    assert set(loaded.dynamixs_analyses) == {T1T2, MODELFREE}


def test_the_t1t2_page_reopens_both_fields_from_the_bundle(app, loaded, bundle):
    dialog = _dialog(app, loaded, T1T2)
    try:
        page = dialog.t1t2_page
        assert dialog.stack.currentWidget() is page
        assert page.field2_enabled
        assert page.field1_freq_spin.value() == 600 and page.field2_freq_spin.value() == 700
        assert page.fitted_experiments == {"field1_t1", "field1_t2", "field2_t1", "field2_t2"}
        for field in ("field1", "field2"):
            for exp in ("t1", "t2"):
                entry = page.session_results[field][exp]
                folder = Path(entry["json_folder"])
                assert _inside(folder, bundle), folder
                assert list(folder.glob("*_fit_data.json")), f"no fit JSON in {folder}"
                assert _inside(entry["results_file"], bundle)
                assert Path(entry["results_file"]).exists()
    finally:
        dialog.deleteLater()


def test_the_modelfree_page_reopens_its_results_from_the_bundle(app, loaded, bundle):
    dialog = _dialog(app, loaded, MODELFREE)
    try:
        page = dialog.integrated_page
        assert dialog.stack.currentWidget() is page
        assert page.prefix_edit.text() == "dual"
        assert _inside(page.output_dir, bundle)
        assert (Path(page.output_dir) / "dual_spectral_density_detailed.csv").exists()
    finally:
        dialog.deleteLater()
