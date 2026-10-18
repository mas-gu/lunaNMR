# ABOUTME: A project save writes exactly the session, and a project load starts from a clean one.
# ABOUTME: Real DNAJA1 series + Kd fit + methyl T2 fit -> bundle -> the real ProjectManager.

"""
The GUI holds one project at a time on its main window. Opening a second project, or
saving after deleting something, must leave the bundle and the session saying the same
thing: what one project lacks must not come back from another, and what the user removed
must not come back on reopen.
"""

import json
import shutil
import types
from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parents[2] / "data_test"
METHYL = (DATA / "methyl" / "sampleNBD" / "series_results_sampleNBD_20260818_124619"
          / "peak_intensity_matrix.csv")

pytestmark = pytest.mark.skipif(not METHYL.exists(), reason="data_test not present")

SERIES, KD, METHYL_T2 = "DNAJA1", "DNAJA1", "NBD_methylT2"


@pytest.fixture(scope="module")
def run_dir(dnaja1_titration_series, tmp_path_factory):
    """A run folder: a 3-point titration series, its Kd fit and a methyl T2 fit."""
    from lunaNMR.cli import main
    run = tmp_path_factory.mktemp("session_run")
    series = run / "series" / SERIES
    shutil.copytree(dnaja1_titration_series, series)
    assert main(["kd", "--input", str(series / "series_analysis_tidy.csv"),
                 "--out", str(run / "kd" / KD), "--prefix", KD,
                 "--p0", "50", "--conc", "0,25,150"]) == 0
    assert main(["dynamixs", "methyl-t2", "--input", str(METHYL),
                 "--out", str(run / "fits" / "NBD"), "--prefix", "NBD"]) == 0
    return run


@pytest.fixture(scope="module")
def bundle_a(run_dir):
    """Project A: the run's series, Kd and DynamiXs analyses, plus the series' first
    spectrum's per-peak fits as the project's fit results."""
    from lunaNMR.cli import main
    from lunaNMR.utils.project_manager import ProjectManager
    built = run_dir.parent / "built.lunaNMR"
    assert main(["project", "build", str(run_dir), "--out", str(built)]) == 0
    window = _window()
    pm = ProjectManager(window)
    assert pm.load_project(built)[0]
    first = next(iter(window.saved_series[SERIES].results.values()))
    window.last_fitting_results = [dict(row) for row in first["integration_results"]]
    out = run_dir.parent / "A.lunaNMR"
    ok, _ = pm.save_project(out)
    assert ok
    return out


@pytest.fixture(scope="module")
def empty_bundle(tmp_path_factory):
    """Project B: saved from a session holding nothing."""
    from lunaNMR.utils.project_manager import ProjectManager
    out = tmp_path_factory.mktemp("empty") / "B.lunaNMR"
    ok, _ = ProjectManager(_window()).save_project(out)
    assert ok
    return out


def _window():
    return types.SimpleNamespace(integrator=types.SimpleNamespace(peak_list=None))


def _copy(bundle, tmp_path, name="P.lunaNMR"):
    dest = tmp_path / name
    shutil.copytree(bundle, dest)
    return dest


def _loaded(bundle):
    from lunaNMR.utils.project_manager import ProjectManager
    window = _window()
    pm = ProjectManager(window)
    ok, _, summary = pm.load_project(bundle)
    assert ok
    return window, pm, summary


def test_project_a_holds_every_kind_of_result(bundle_a):
    window, _, _ = _loaded(bundle_a)
    assert list(window.saved_series) == [SERIES]
    assert list(window.kd_analyses) == [KD]
    assert list(window.dynamixs_analyses) == [METHYL_T2]
    assert window.last_fitting_results


# ---------- a reopened project saves onto itself ----------

def test_a_reopened_project_saves_onto_itself(bundle_a, tmp_path):
    """A loaded DynamiXs analysis points at its bundled results; saving back to the same
    bundle copies them onto themselves."""
    bundle = _copy(bundle_a, tmp_path)
    window, pm, _ = _loaded(bundle)
    ok, summary = pm.save_project(bundle)
    assert ok
    assert summary["dynamixs"] == {"analyses": [METHYL_T2]}
    page = bundle / "dynamixs" / "analyses" / METHYL_T2 / "results" / "methyl_t2"
    assert json.loads((page / "NBD_field1_methylT2_fit_data.json").read_text())
    assert _loaded(bundle)[0].dynamixs_analyses[METHYL_T2]


def test_a_bundled_json_folder_survives_a_save_onto_itself(bundle_a, tmp_path):
    """The GUI's methyl and T1/T2 pages keep their fit JSON in a json* folder."""
    bundle = _copy(bundle_a, tmp_path)
    page = bundle / "dynamixs" / "analyses" / METHYL_T2 / "results" / "methyl_t2"
    (page / "json").mkdir()
    shutil.move(page / "NBD_field1_methylT2_fit_data.json", page / "json")
    window, pm, _ = _loaded(bundle)
    assert pm.save_project(bundle)[0]
    assert (page / "json" / "NBD_field1_methylT2_fit_data.json").exists()


# ---------- opening a project starts from a clean session ----------

def test_opening_a_project_keeps_nothing_of_the_one_before(bundle_a, empty_bundle):
    window, pm, _ = _loaded(bundle_a)
    window.current_voigt_result = window.last_fitting_results[0]
    ok, _, summary = pm.load_project(empty_bundle)
    assert ok
    assert not window.kd_analyses
    assert not window.dynamixs_analyses
    assert not window.saved_series
    assert window.batch_results is None
    assert window.last_fitting_results is None
    assert window.current_voigt_result is None
    assert window.integrator.peak_list is None
    assert summary["dynamixs"] == {}


def test_saving_the_second_project_writes_nothing_of_the_first(bundle_a, empty_bundle,
                                                               tmp_path):
    bundle = _copy(empty_bundle, tmp_path)
    window, pm, _ = _loaded(bundle_a)
    assert pm.load_project(bundle)[0]
    assert pm.save_project(bundle)[0]
    for held in ("kd", "dynamixs", "series_results", "fit_results", "peak_list.csv"):
        assert not (bundle / held).exists(), held


# ---------- a save prunes what the session no longer holds ----------

def test_deleting_everything_and_saving_removes_it_from_the_bundle(bundle_a, tmp_path):
    bundle = _copy(bundle_a, tmp_path)
    window, pm, _ = _loaded(bundle)
    window.kd_analyses = {}
    window.dynamixs_analyses = {}
    window.saved_series = {}
    window.batch_results = None
    window.last_fitting_results = None
    window.integrator.peak_list = None
    assert pm.save_project(bundle)[0]
    window, _, _ = _loaded(bundle)
    assert not window.kd_analyses
    assert not window.dynamixs_analyses
    assert not window.saved_series
    assert window.last_fitting_results is None
    assert window.integrator.peak_list is None


def test_a_renamed_series_leaves_no_folder_under_its_old_name(bundle_a, tmp_path):
    bundle = _copy(bundle_a, tmp_path)
    window, pm, _ = _loaded(bundle)
    window.saved_series = {"renamed": window.saved_series[SERIES]}
    assert pm.save_project(bundle)[0]
    assert sorted(p.name for p in (bundle / "series_results").iterdir()) == \
        ["renamed", "series_manifest.json"]
    assert list(_loaded(bundle)[0].saved_series) == ["renamed"]


def test_save_as_onto_another_bundle_replaces_what_it_held(bundle_a, empty_bundle, tmp_path):
    target = _copy(bundle_a, tmp_path, "target.lunaNMR")
    window, pm, _ = _loaded(empty_bundle)
    assert pm.save_project(target)[0]
    for held in ("kd", "dynamixs", "series_results", "fit_results", "peak_list.csv"):
        assert not (target / held).exists(), held


def test_a_slim_save_leaves_no_region_arrays_of_an_embedded_one(bundle_a, tmp_path):
    """fits.json no longer names them; left behind they are only weight."""
    bundle = _copy(bundle_a, tmp_path)
    window, pm, _ = _loaded(bundle)
    for row in window.last_fitting_results:
        row["region_2d"] = {"f1_ppm": [1.0], "f2_ppm": [2.0], "intensity": [[3.0]]}
    assert pm.save_project(bundle, embed_region_data=True)[0]
    assert list((bundle / "fit_results" / "arrays").glob("*.npz"))
    window, pm, _ = _loaded(bundle)
    assert pm.save_project(bundle)[0]
    assert not list((bundle / "fit_results" / "arrays").glob("*.npz"))


# ---------- the main window: what still shows the previous project ----------

@pytest.fixture(scope="module")
def app():
    import os
    import sys
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


def test_opening_a_project_closes_the_windows_showing_the_previous_one(app, bundle_a,
                                                                       empty_bundle,
                                                                       tmp_path):
    """The Kd and DynamiXs windows, and the Spectral Inspector, show the previous
    project; a save captures what they show."""
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    spectrum = DATA / "data_binding" / "DNAJA1_HSPA8" / "2D_DNAJA1_0o0.ft"
    window = LunaNMRMainWindow()
    try:
        pm = window._get_project_manager()
        assert pm.load_project(bundle_a)[0]
        window._refresh_after_project_load()
        assert window.integrator.fitted_peaks
        window.launch_kd_titration()
        window.launch_dynamixs()
        window.launch_spectral_inspector()
        window._spectral_inspector.load_state({"groups": [{"name": "A", "spectra": [
            {"display_name": "A", "file_path": str(spectrum)}]}]})
        assert window._spectral_inspector.groups

        bundle = _copy(empty_bundle, tmp_path)
        assert pm.load_project(bundle)[0]
        window._refresh_after_project_load()
        assert window.kd_titration_dialog is None
        assert window.dynamixs_dialog is None
        assert not window.integrator.fitted_peaks
        assert not window._spectral_inspector.groups
        assert pm.save_project(bundle)[0]
        for held in ("kd", "dynamixs", "series_results", "fit_results",
                     "spectral_inspector"):
            assert not (bundle / held).exists(), held
    finally:
        if window._spectral_inspector is not None:
            window._spectral_inspector.close()
            window._spectral_inspector.deleteLater()
        window.close()
        window.deleteLater()


# ---------- a series' CSVs: its bundled copy, or its own run's ----------

SERIES_CSVS = ("series_analysis_tidy.csv", "peak_intensity_matrix.csv")


def _point_series_at(bundle, folder):
    """Rewrite the bundled series' record so its run folder is `folder`."""
    batch = bundle / "series_results" / SERIES / "batch_results.json"
    data = json.loads(batch.read_text())
    data["metadata"]["output_folder"] = str(folder)
    data["metadata"]["csv_path"] = str(folder / "series_analysis_tidy.csv")
    batch.write_text(json.dumps(data))


def _rerun_folder(bundle, folder):
    """A re-run's output folder: the bundled CSVs, each marked so the two are told apart."""
    folder.mkdir()
    for name in SERIES_CSVS:
        text = (bundle / "series_results" / SERIES / name).read_text()
        (folder / name).write_text(text + "\n")
    return folder


def test_save_as_carries_the_bundled_csvs_once_the_run_has_moved(bundle_a, tmp_path):
    bundle = _copy(bundle_a, tmp_path)
    _point_series_at(bundle, tmp_path / "gone")
    window, pm, _ = _loaded(bundle)
    other = tmp_path / "P2.lunaNMR"
    assert pm.save_project(other)[0]
    for name in SERIES_CSVS:
        assert ((other / "series_results" / SERIES / name).read_text()
                == (bundle / "series_results" / SERIES / name).read_text()), name


def test_a_renamed_series_keeps_its_bundled_csvs(bundle_a, tmp_path):
    bundle = _copy(bundle_a, tmp_path)
    _point_series_at(bundle, tmp_path / "gone")
    original = {name: (bundle / "series_results" / SERIES / name).read_text()
                for name in SERIES_CSVS}
    window, pm, _ = _loaded(bundle)
    window.saved_series = {"renamed": window.saved_series[SERIES]}
    assert pm.save_project(bundle)[0]
    for name in SERIES_CSVS:
        assert (bundle / "series_results" / "renamed" / name).read_text() == original[name]


def _dynamixs_dialog(window, bundle):
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    window.current_project_path = bundle
    return DynamiXsDialog(parent=None, main_window=window)


def test_a_series_rerun_under_its_name_reads_its_own_csvs(app, bundle_a, tmp_path):
    """A re-run that overwrites a loaded series' name is a new series; the project's
    copy of the old one must not stand in for it."""
    from lunaNMR.utils.project_manager import batch_results_from_dict
    bundle = _copy(bundle_a, tmp_path)
    window, _, _ = _loaded(bundle)
    rerun = _rerun_folder(bundle, tmp_path / "rerun")
    data = json.loads((bundle / "series_results" / SERIES / "batch_results.json").read_text())
    batch = batch_results_from_dict(data)
    batch.metadata["output_folder"] = str(rerun)
    batch.metadata["csv_path"] = str(rerun / "series_analysis_tidy.csv")
    window.saved_series[SERIES] = batch
    dialog = _dynamixs_dialog(window, bundle)
    try:
        page = dialog.hetnoe_page
        [series] = page._get_available_series(include=lambda br: True)
        assert series["csv_path"] == str(rerun / "series_analysis_tidy.csv")
        assert page._series_intensity_matrix(SERIES) == str(rerun / "peak_intensity_matrix.csv")
    finally:
        dialog.deleteLater()


def test_a_loaded_series_reads_its_bundled_csvs(app, bundle_a, tmp_path):
    """The CSVs that belong to a loaded series are the ones beside the batch_results.json
    it was read from, whatever its recorded run folder holds now."""
    bundle = _copy(bundle_a, tmp_path)
    rerun = _rerun_folder(bundle, tmp_path / "rerun")
    _point_series_at(bundle, rerun)
    window, _, _ = _loaded(bundle)
    dialog = _dynamixs_dialog(window, bundle)
    try:
        page = dialog.hetnoe_page
        held = bundle / "series_results" / SERIES
        [series] = page._get_available_series(include=lambda br: True)
        assert series["csv_path"] == str(held / "series_analysis_tidy.csv")
        assert page._series_intensity_matrix(SERIES) == str(held / "peak_intensity_matrix.csv")
    finally:
        dialog.deleteLater()


# ---------- a saved DynamiXs analysis keeps the results it was saved with ----------

METHYL_JSON = "NBD_field1_methylT2_fit_data.json"


def _mark_run(folder, run):
    path = folder / METHYL_JSON
    data = json.loads(path.read_text())
    data["run"] = run
    path.write_text(json.dumps(data))


def test_a_later_run_into_the_same_folder_leaves_a_saved_analysis_alone(run_dir, tmp_path):
    """A page writes every run into its one output folder under the same names; the
    analysis saved from the first run must not take the second's fits on the next save."""
    from lunaNMR.utils.project_manager import ProjectManager
    live = tmp_path / "methyl_out"
    shutil.copytree(run_dir / "fits" / "NBD", live)
    _mark_run(live, "first")
    window = _window()
    window.dynamixs_analyses = {"first": {"state": {"methyl_t2": {"output_dir": str(live)}},
                                          "file_refs": {}}}
    pm = ProjectManager(window)
    bundle = tmp_path / "P.lunaNMR"
    assert pm.save_project(bundle)[0]

    _mark_run(live, "second")
    window.dynamixs_analyses["second"] = {"state": {"methyl_t2": {"output_dir": str(live)}},
                                          "file_refs": {}}
    assert pm.save_project(bundle)[0]
    for name in ("first", "second"):
        page = bundle / "dynamixs" / "analyses" / name / "results" / "methyl_t2"
        assert json.loads((page / METHYL_JSON).read_text())["run"] == name

    other = tmp_path / "P2.lunaNMR"
    assert pm.save_project(other)[0]
    page = other / "dynamixs" / "analyses" / "first" / "results" / "methyl_t2"
    assert json.loads((page / METHYL_JSON).read_text())["run"] == "first"


# ---------- a bundle saved with its region data says so when reopened ----------

REGION = {"f1_ppm": [1.0], "f2_ppm": [2.0], "intensity": [[3.0]]}


def test_a_slim_bundle_reopens_as_slim(bundle_a):
    window, _, _ = _loaded(bundle_a)
    assert window.region_data_embedded is False


def test_a_bundle_with_embedded_fit_regions_reopens_as_embedded(bundle_a, tmp_path):
    bundle = _copy(bundle_a, tmp_path)
    window, pm, _ = _loaded(bundle)
    for row in window.last_fitting_results:
        row["region_2d"] = dict(REGION)
    assert pm.save_project(bundle, embed_region_data=True)[0]
    assert window.region_data_embedded is True
    assert _loaded(bundle)[0].region_data_embedded is True


def test_a_bundle_with_embedded_series_regions_reopens_as_embedded(bundle_a, tmp_path):
    bundle = _copy(bundle_a, tmp_path)
    window, pm, _ = _loaded(bundle)
    window.last_fitting_results = None
    for result in window.saved_series[SERIES].results.values():
        for row in result["integration_results"]:
            row["region_2d"] = dict(REGION)
    assert pm.save_project(bundle, embed_region_data=True)[0]
    assert _loaded(bundle)[0].region_data_embedded is True
    window, pm, _ = _loaded(bundle)
    assert pm.save_project(bundle)[0]
    assert window.region_data_embedded is False
    assert _loaded(bundle)[0].region_data_embedded is False


@pytest.mark.parametrize("embedded", [False, True])
def test_the_save_box_starts_from_how_the_project_was_saved(app, monkeypatch, tmp_path,
                                                            embedded):
    """Unticked by default, a re-save of an embedded bundle opened without its spectra
    strips the regions it can no longer rebuild."""
    from PySide6.QtWidgets import QMessageBox
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    shown = []

    def answer(box):
        shown.append(box.checkBox().isChecked())
        return QMessageBox.Cancel

    monkeypatch.setattr(QMessageBox, "exec", answer)
    window = LunaNMRMainWindow()
    try:
        window.region_data_embedded = embedded
        window._do_save_project(tmp_path / "P.lunaNMR")
        assert shown == [embedded]
        assert not (tmp_path / "P.lunaNMR").exists()
    finally:
        window.close()
        window.deleteLater()


# ---------- the load summary says what loaded and what did not ----------

def test_a_bundle_without_series_does_not_report_series(empty_bundle):
    assert "Series results" not in _loaded(empty_bundle)[2]["loaded_items"]


def test_what_could_not_be_read_is_reported(bundle_a, tmp_path):
    bundle = _copy(bundle_a, tmp_path)
    (bundle / "series_results" / SERIES / "batch_results.json").unlink()
    (bundle / "kd" / "analyses" / KD / "fit_data.json").write_text("{not json")
    (bundle / "dynamixs" / "analyses" / METHYL_T2 / "state.json").write_text("{not json")
    window, _, summary = _loaded(bundle)
    assert "Series results" not in summary["loaded_items"]
    assert "Kd titration" not in summary["loaded_items"]
    assert summary["dynamixs"] == {}
    assert sorted(summary["skipped"]) == sorted([
        f"Series '{SERIES}' (no batch_results.json)",
        f"Kd analysis '{KD}' (unreadable)",
        f"DynamiXs analysis '{METHYL_T2}' (unreadable)",
    ])


def _popup_text(monkeypatch, show):
    from PySide6.QtWidgets import QMessageBox
    shown = []

    def answer(box):
        shown.append(box.text())
        return QMessageBox.Ok

    monkeypatch.setattr(QMessageBox, "exec", answer)
    show()
    [text] = shown
    return text


def test_the_popups_name_the_dynamixs_analyses_and_what_was_skipped(app, monkeypatch,
                                                                     tmp_path):
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    window = LunaNMRMainWindow()
    try:
        loaded = _popup_text(monkeypatch, lambda: window._show_load_summary(
            tmp_path / "P.lunaNMR",
            {"loaded_items": ["Project manifest"], "dynamixs": {"analyses": [METHYL_T2]},
             "skipped": [f"Kd analysis '{KD}' (unreadable)"]}))
        assert METHYL_T2 in loaded
        assert f"Kd analysis '{KD}' (unreadable)" in loaded
        saved = _popup_text(monkeypatch, lambda: window._show_save_summary(
            tmp_path / "P.lunaNMR",
            {"saved_items": ["Project manifest"], "dynamixs": {"analyses": [METHYL_T2]}}))
        assert METHYL_T2 in saved
    finally:
        window.close()
        window.deleteLater()
