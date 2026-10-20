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
METHYL = (Path(__file__).resolve().parents[2] / "data_test" / "methyl" / "sampleNBD"
          / "series_results_sampleNBD_20260818_124619" / "peak_intensity_matrix.csv")

pytestmark = pytest.mark.skipif(not (RUN.exists() and METHYL.exists()),
                                reason="data_test not present")


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
    assert main(["dynamixs", "methyl-t2", "--input", str(METHYL),
                 "--out", str(run / "fits" / "NBD_methyl"), "--prefix", "NBD"]) == 0
    assert main(["dynamixs", "hetnoe", "--sat", str(NOE / "f600_noe_sat.csv"),
                 "--unsat", str(NOE / "f600_noe_unsat.csv"),
                 "--out", str(run / "hetnoe" / "f600"), "--prefix", "f600"]) == 0
    # The per-field spectral-density table the model-free run built, from its own folder.
    density_in = run / "density" / "f600_input.csv"
    density_in.parent.mkdir()
    shutil.copy(run / "modelfree" / "dual" / "dual_field1_input.csv", density_in)
    assert main(["dynamixs", "density", "--input", str(density_in), "--field1-freq", "600",
                 "--out", str(run / "density" / "f600"), "--prefix", "f600",
                 "--no-plot"]) == 0
    # A second model-free record of the same prefix, elsewhere in the run: both must be
    # bundled, not the later one over the earlier.
    shutil.copytree(run / "modelfree" / "dual", run / "modelfree" / "dual_rerun")
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
MODELFREE_RERUN = "modelfree_dual_rerun_modelfree"
METHYL_T2 = "NBD_methylT2"
HETNOE = "f600_hetNOE"
DENSITY = "f600_spectralDensity"


def test_modelfree_keeps_its_intermediates_in_its_own_folder(run_dir):
    """A folder shared by every run would let two model-free runs on one prefix
    overwrite each other's density inputs mid-run."""
    for field in ("field1", "field2"):
        assert (run_dir / "modelfree" / "dual" / f"dual_{field}_input.csv").exists()


def test_each_command_records_what_it_ran(run_dir):
    record = json.loads((run_dir / "fits" / "600_T1" / "600_T1_dynamixs.json").read_text())
    assert record["command"].startswith("dynamixs t1t2")
    assert Path(record["args"]["input"]).is_absolute()
    record = json.loads((run_dir / "modelfree" / "dual" / "dual_dynamixs.json").read_text())
    assert record["command"].startswith("dynamixs modelfree")


def test_the_four_fits_become_one_t1t2_analysis_and_modelfree_its_own(loaded):
    assert set(loaded.dynamixs_analyses) == {T1T2, MODELFREE, MODELFREE_RERUN, METHYL_T2,
                                             HETNOE, DENSITY}


def test_the_hetnoe_page_reopens_its_results_from_the_bundle(app, loaded, bundle):
    dialog = _dialog(app, loaded, HETNOE)
    try:
        page = dialog.hetnoe_page
        assert dialog.stack.currentWidget() is page
        assert page.sat_file == str(NOE / "f600_noe_sat.csv")
        assert _inside(page.last_results_file, bundle)
        assert page.results_table.rowCount() > 0
    finally:
        dialog.deleteLater()


def test_the_spectral_density_page_reopens_from_the_bundle(app, loaded, bundle):
    dialog = _dialog(app, loaded, DENSITY)
    try:
        page = dialog.spectral_page
        assert dialog.stack.currentWidget() is page
        assert page.prefix_edit.text() == "f600"
        assert page.field1_spin.value() == 600
        results = page.session_results["results_file"]
        assert _inside(results, bundle) and Path(results).exists()
    finally:
        dialog.deleteLater()


def test_the_methyl_t2_page_reopens_its_fits_from_the_bundle(app, loaded, bundle):
    dialog = _dialog(app, loaded, METHYL_T2)
    try:
        page = dialog.methyl_t2_page
        assert dialog.stack.currentWidget() is page
        assert page.input_file == str(METHYL)
        folder = Path(page.last_json_folder)
        assert _inside(folder, bundle), folder
        assert list(folder.glob("*_methylT2_fit_data.json"))
        assert _inside(page.last_results_file, bundle)
        assert Path(page.last_results_file).exists()
        assert page.viewer_btn.isEnabled()
    finally:
        dialog.deleteLater()


def test_the_t1t2_results_hold_one_entry_per_field(app, loaded):
    """The page's viewers walk session_results as {field: {t1, t2, t1rho}}; a stray
    top-level value there stops both from opening."""
    dialog = _dialog(app, loaded, T1T2)
    try:
        session = dialog.t1t2_page.session_results
        assert set(session) == {"field1", "field2"}
        assert all(v is None or isinstance(v, dict) for v in session.values())
    finally:
        dialog.deleteLater()


def test_an_analysis_keeps_the_command_it_came_from(loaded):
    """Inputs the page has no field for (an R2 table, say) are still on record."""
    meta = loaded.dynamixs_analyses[MODELFREE]["state"]["analysis_metadata"]
    assert meta["cli_runs"][0]["args"]["f1_t1"] == _matrix("600_T1")


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


PAGE_OF = {T1T2: "t1t2", MODELFREE: "integrated", METHYL_T2: "methyl_t2",
           HETNOE: "hetnoe", DENSITY: "spectral"}


@pytest.mark.parametrize("name", list(PAGE_OF))
def test_a_built_analysis_names_its_page_by_key(loaded, name):
    """A stack index moves whenever a page is added or removed; the key does not."""
    assert loaded.dynamixs_analyses[name]["state"]["active_page"] == PAGE_OF[name]


@pytest.mark.parametrize("name", list(PAGE_OF))
def test_a_gui_saved_analysis_reopens_on_its_page(app, loaded, name, tmp_path):
    """Reopened from the built bundle, saved again by the dialog through a real project
    save, loaded and reopened: it lands on its own page."""
    import types
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    from lunaNMR.utils.project_manager import ProjectManager
    window = types.SimpleNamespace(integrator=types.SimpleNamespace(peak_list=None),
                                   dynamixs_analyses={})
    first = DynamiXsDialog(parent=None, main_window=window)
    try:
        first.open_analysis(loaded.dynamixs_analyses[name], name)
        assert first.ensure_current_saved() == name
        assert window.dynamixs_analyses[name]["state"]["active_page"] == PAGE_OF[name]
        assert ProjectManager(window).save_project(tmp_path / "again.lunaNMR")[0]
    finally:
        first.deleteLater()
    reloaded = types.SimpleNamespace(integrator=types.SimpleNamespace(peak_list=None))
    assert ProjectManager(reloaded).load_project(tmp_path / "again.lunaNMR")
    second = DynamiXsDialog(parent=None, main_window=reloaded)
    try:
        second.open_analysis(reloaded.dynamixs_analyses[name], name)
        assert second.stack.currentWidget() is second._pages()[PAGE_OF[name]]
        assert second.page_notice.isHidden()
    finally:
        second.deleteLater()


@pytest.mark.parametrize("active_page", [6, "plotting"])
def test_a_page_it_cannot_read_opens_no_page_and_says_so(app, loaded, active_page):
    """An integer is a stack index from another page layout; whatever page it named then
    is not the one it names now, so none is guessed."""
    import copy
    import types
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    entry = copy.deepcopy(loaded.dynamixs_analyses[HETNOE])
    entry["state"]["active_page"] = active_page
    # Its own window: closing the dialog keeps its state on the window it was given
    dialog = DynamiXsDialog(parent=None, main_window=types.SimpleNamespace())
    try:
        dialog.show()
        dialog.open_analysis(entry, HETNOE)
        assert dialog.stack.currentWidget() is dialog.main_menu_page
        assert not dialog.page_notice.isHidden()
        assert "page layout this version does not read" in dialog.page_notice.text()
        dialog.show_hetnoe()
        assert dialog.page_notice.isHidden()
    finally:
        dialog.close()
        dialog.deleteLater()


def test_the_menu_offers_no_plot_data_page(app):
    import types
    from PySide6.QtWidgets import QPushButton
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    dialog = DynamiXsDialog(parent=None, main_window=types.SimpleNamespace())
    try:
        labels = [b.text() for b in dialog.main_menu_page.findChildren(QPushButton)]
        assert "Plot Data" not in labels
        assert not hasattr(dialog, "plotting_page")
        assert dialog.stack.count() == 1 + len(dialog._pages())
        opened = []
        for button in dialog.main_menu_page.findChildren(QPushButton):
            if button.text() in ("Change Directory", "Open selected run", "Delete selected"):
                continue
            dialog.show_main_menu()
            button.click()
            opened.append(dialog._active_page_key())
        assert sorted(opened) == sorted(dialog._pages())
    finally:
        dialog.deleteLater()


@pytest.fixture
def default_prefix_run(tmp_path):
    """Two fits on the default --prefix, relative paths, from a folder that holds names
    equal to their argument values - as a run made by hand is laid out."""
    from lunaNMR.cli import main
    work = tmp_path / "work"
    for name in ("field1", "T1", "T2"):
        (work / name).mkdir(parents=True)
    here = os.getcwd()
    os.chdir(work)
    try:
        for out, exp in (("t1", "T1"), ("t2", "T2")):
            assert main(["dynamixs", "t1t2", "--input", _matrix(f"600_{exp}"), "--out", out,
                         "--exp", exp, "--field-freq", "600"]) == 0
    finally:
        os.chdir(here)
    return work


def test_a_record_keeps_values_that_are_not_paths(default_prefix_run):
    record = json.loads((default_prefix_run / "t1" / "field1_dynamixs.json").read_text())
    assert record["args"]["prefix"] == "field1"
    assert record["args"]["exp"] == "T1"
    assert record["args"]["out"] == str(default_prefix_run / "t1")


def test_experiments_on_one_prefix_keep_their_own_results(default_prefix_run, tmp_path):
    import types
    from lunaNMR.cli import main
    from lunaNMR.utils.project_manager import ProjectManager
    out = tmp_path / "x.lunaNMR"
    assert main(["project", "build", str(default_prefix_run), "--out", str(out)]) == 0
    window = types.SimpleNamespace(integrator=types.SimpleNamespace(peak_list=None))
    assert ProjectManager(window).load_project(out)
    (state,) = [a["state"] for a in window.dynamixs_analyses.values()]
    field = state["t1t2"]["session_results"]["field1"]
    t1, t2 = Path(field["t1"]["results_file"]), Path(field["t2"]["results_file"])
    assert t1 != t2 and t1.read_text() != t2.read_text()
    assert Path(field["t1"]["json_folder"]) != Path(field["t2"]["json_folder"])


def test_the_same_experiment_twice_is_refused_not_overwritten(default_prefix_run, tmp_path):
    from lunaNMR.cli import main
    shutil.copytree(default_prefix_run / "t1", default_prefix_run / "t1_again")
    out = tmp_path / "x.lunaNMR"
    assert main(["project", "build", str(default_prefix_run), "--out", str(out)]) == 1
    assert not out.exists(), "a failed build must not leave a bundle behind"


def test_a_malformed_record_is_an_error_not_a_traceback(default_prefix_run, tmp_path):
    from lunaNMR.cli import main
    (default_prefix_run / "t1" / "broken_dynamixs.json").write_text("{")
    out = tmp_path / "x.lunaNMR"
    assert main(["project", "build", str(default_prefix_run), "--out", str(out)]) == 1
    assert not out.exists()


def test_a_cli_modelfree_run_writes_its_fit_json(run_dir):
    """The page's View Fits reads the T1/T2 fit JSONs a run writes, as a GUI run does."""
    names = sorted(p.name for p in (run_dir.parent / "JA1.lunaNMR").rglob(
        "dynamixs/analyses/dual_modelfree/**/*_fit_data.json"))
    assert names == ["dual_field1_T1_fit_data.json", "dual_field1_T2_fit_data.json",
                     "dual_field2_T1_fit_data.json", "dual_field2_T2_fit_data.json"]


def test_the_modelfree_page_views_its_fits_from_the_bundle(app, loaded, bundle, monkeypatch):
    from visualization import fit_viewer
    for box in ("show_info", "show_warning", "show_error"):
        monkeypatch.setattr(fit_viewer, box, lambda *a, **k: None)
    dialog = _dialog(app, loaded, MODELFREE)
    try:
        page = dialog.integrated_page
        json_files = page.session_results.get('json_files') or []
        assert len(json_files) == 4
        assert all(_inside(f, bundle) and Path(f).exists() for f in json_files)
        page._open_fit_viewer()
        assert getattr(page, 'fit_viewer', None) is not None
        assert "nothing to view" not in page.progress_text.toPlainText()
        page.fit_viewer.close()
    finally:
        dialog.deleteLater()
