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
