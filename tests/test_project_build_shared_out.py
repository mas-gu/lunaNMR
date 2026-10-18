# ABOUTME: `project build` bundles each DynamiXs run with its own output files only,
# ABOUTME: even when several runs of different --prefix were written into one --out folder.

"""
The project save copies a page's results from its output folder by name pattern, so a
folder shared by two runs gave each named analysis both runs' files - and the methyl
viewer opens every fit JSON in its folder, so reopening one analysis showed the other's
fits as its own. Real CLI runs on real matrices, built, then the run folder deleted.
"""

import shutil
import sys
import types
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

DATA = Path(__file__).resolve().parents[2] / "data_test"
METHYL = (DATA / "methyl" / "sampleNBD" / "series_results_sampleNBD_20260818_124619"
          / "peak_intensity_matrix.csv")
NOE = DATA / "data_JA1_relaxation" / "run_20260825_1405" / "results" / "hetnoe"

pytestmark = pytest.mark.skipif(not (METHYL.exists() and NOE.exists()),
                                reason="data_test not present")

PREFIXES = ("first", "second")


@pytest.fixture(scope="module")
def loaded(tmp_path_factory):
    from lunaNMR.cli import main
    from lunaNMR.utils.project_manager import ProjectManager
    run = tmp_path_factory.mktemp("run")
    for prefix in PREFIXES:
        assert main(["dynamixs", "methyl-t2", "--input", str(METHYL),
                     "--out", str(run / "methyl"), "--prefix", prefix]) == 0
        assert main(["dynamixs", "hetnoe", "--sat", str(NOE / "f600_noe_sat.csv"),
                     "--unsat", str(NOE / "f600_noe_unsat.csv"),
                     "--out", str(run / "noe"), "--prefix", prefix]) == 0
    bundle = run.parent / "shared.lunaNMR"
    assert main(["project", "build", str(run), "--out", str(bundle)]) == 0
    shutil.rmtree(run)
    window = types.SimpleNamespace(integrator=types.SimpleNamespace(peak_list=None))
    assert ProjectManager(window).load_project(bundle)
    return bundle, window


def _bundled(bundle, name, page):
    return sorted(p.name for p in
                  (bundle / "dynamixs" / "analyses" / name / "results" / page).iterdir())


@pytest.mark.parametrize("suffix,page", [("methylT2", "methyl_t2"), ("hetNOE", "hetnoe")])
def test_each_analysis_holds_only_its_own_run_files(loaded, suffix, page):
    bundle, window = loaded
    for prefix in PREFIXES:
        other = next(p for p in PREFIXES if p != prefix)
        files = _bundled(bundle, f"{prefix}_{suffix}", page)
        assert files, f"{prefix}_{suffix}: nothing bundled"
        assert all(f.startswith(f"{prefix}_") for f in files), files
        assert not any(f.startswith(f"{other}_") for f in files), files


def test_the_methyl_viewer_folder_holds_one_fit(loaded):
    bundle, window = loaded
    for prefix in PREFIXES:
        page = window.dynamixs_analyses[f"{prefix}_methylT2"]["state"]["methyl_t2"]
        folder = Path(page["last_json_folder"])
        assert folder.resolve().is_relative_to(bundle.resolve())
        assert [p.name for p in folder.glob("*_methylT2_fit_data.json")] == [
            f"{prefix}_field1_methylT2_fit_data.json"]
        assert Path(page["last_results_file"]).name == f"{prefix}_fit_results.txt"
        assert Path(page["last_results_file"]).exists()


def test_the_hetnoe_results_file_is_its_own_inside_the_bundle(loaded):
    bundle, window = loaded
    for prefix in PREFIXES:
        page = window.dynamixs_analyses[f"{prefix}_hetNOE"]["state"]["hetnoe"]
        results = Path(page["last_results_file"])
        assert results.name == f"{prefix}_hetnoe.csv"
        assert results.exists() and results.resolve().is_relative_to(bundle.resolve())


def test_prefixes_that_extend_one_another_are_told_apart(tmp_path):
    """`a` must not claim `a_b`'s files: a file belongs to the longest prefix it starts with."""
    from lunaNMR.cli import main
    from lunaNMR.utils.project_builder import dynamixs_analyses
    run = tmp_path / "run"
    for prefix in ("a", "a_b"):
        assert main(["dynamixs", "hetnoe", "--sat", str(NOE / "f600_noe_sat.csv"),
                     "--unsat", str(NOE / "f600_noe_unsat.csv"),
                     "--out", str(run / "noe"), "--prefix", prefix]) == 0
    analyses, _ = dynamixs_analyses(run, tmp_path / "staging")
    for name, prefix in (("a_hetNOE", "a"), ("a_b_hetNOE", "a_b")):
        out = Path(analyses[name]["state"]["hetnoe"]["output_dir"])
        assert sorted(p.name for p in out.iterdir() if p.suffix == ".csv") == [
            f"{prefix}_hetnoe.csv"]
