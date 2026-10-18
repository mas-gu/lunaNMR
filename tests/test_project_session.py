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
