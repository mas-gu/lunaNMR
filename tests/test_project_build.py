# ABOUTME: `project build` turns an automated run's output folders into a .lunaNMR bundle.
# ABOUTME: Round trip on real data: series + kd -> build -> load with the real ProjectManager.

"""
A full automated run writes series and Kd results to plain folders; the GUI opens
.lunaNMR bundles. `project build` assembles one through ProjectManager.save_project
itself, so the bundle is exactly what the GUI would have saved. The test that matters is
the round trip: what the loader reads back.
"""

import json
import sys
import types
from pathlib import Path

import pandas as pd
import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
P0, CONC = "50", "0,6.25,12.5,25,50,100,150"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def run_dir(tmp_path_factory):
    """A run folder laid out as the automated runs lay theirs out: series/ and kd/."""
    from lunaNMR.cli import main
    run = tmp_path_factory.mktemp("run")
    series = run / "series" / "DNAJA1"
    assert main(["series", "--spectra", str(DNAJA1), "--peaks", str(DNAJA1 / "DNAJA1.txt"),
                 "--out", str(series), "--mode", "titration",
                 "--peak-source", "cascade"]) == 0
    assert main(["kd", "--input", str(series / "series_analysis_tidy.csv"),
                 "--out", str(run / "kd" / "DNAJA1"), "--prefix", "DNAJA1",
                 "--p0", P0, "--conc", CONC]) == 0
    return run


@pytest.fixture(scope="module")
def bundle(run_dir):
    from lunaNMR.cli import main
    out = run_dir.parent / "DNAJA1.lunaNMR"
    assert main(["project", "build", str(run_dir), "--out", str(out)]) == 0
    return out


@pytest.fixture(scope="module")
def loaded(bundle):
    """What the GUI gets back: the real loader on a bare main window."""
    from lunaNMR.utils.project_manager import ProjectManager
    window = types.SimpleNamespace(integrator=types.SimpleNamespace(peak_list=None))
    assert ProjectManager(window).load_project(bundle)
    return window


def test_series_writes_its_batch_results(run_dir):
    """The loader rebuilds a series from batch_results.json; nothing else carries it."""
    data = json.loads((run_dir / "series" / "DNAJA1" / "batch_results.json").read_text())
    assert len(data["results"]) == 7
    assert data["series_mode"] == "titration"


def test_the_series_reopens_with_every_spectrum(loaded):
    series = loaded.saved_series["DNAJA1"]
    assert series.series_mode == "titration"
    assert len(series.results) == 7
    first = next(iter(series.results.values()))
    assert first["integration_results"], "the overlay viewer reads the per-peak rows"


def test_the_tidy_csv_is_inside_the_bundle(bundle):
    """The Kd page resolves a titration's input here first."""
    assert (bundle / "series_results" / "DNAJA1" / "series_analysis_tidy.csv").exists()


def test_the_spectrum_and_peak_list_reopen(loaded):
    """current_nmr_file is the one attribute that reloads a spectrum; it must be absolute."""
    spectrum = Path(loaded.current_nmr_file)
    assert spectrum.is_absolute() and spectrum.exists()
    assert spectrum.parent == DNAJA1
    assert set(loaded.integrator.peak_list["Assignment"].astype(str)) >= {"3LysH", "29LeuH"}


def test_the_kd_analysis_reopens_as_the_cli_wrote_it(loaded, run_dir):
    written = json.loads((run_dir / "kd" / "DNAJA1" / "data" / "DNAJA1_kd_fit_data.json")
                         .read_text())
    entry = loaded.kd_analyses["DNAJA1"]
    assert entry["fit_data"] == written
    assert Path(entry["fit_data_path"]).exists()


def test_the_kd_page_parameters_come_from_the_fit(loaded):
    meta = loaded.kd_analyses["DNAJA1"]["meta"]
    assert meta["p0"] == pytest.approx(50.0)
    assert meta["concentrations"] == pytest.approx([0, 6.25, 12.5, 25, 50, 100, 150])
    assert len(meta["detected_points"]) == 7
    assert meta["observable"] == 0          # csp + intensity, the combo's first entry


def test_an_empty_run_is_refused(tmp_path):
    from lunaNMR.cli import main
    (tmp_path / "run").mkdir()
    assert main(["project", "build", str(tmp_path / "run"),
                 "--out", str(tmp_path / "x.lunaNMR")]) == 1
    assert not (tmp_path / "x.lunaNMR").exists()


def test_a_dry_run_writes_nothing(run_dir, tmp_path):
    from lunaNMR.cli import main
    out = tmp_path / "plan.lunaNMR"
    assert main(["project", "build", str(run_dir), "--out", str(out), "--dry-run"]) == 0
    assert not out.exists()


def test_a_series_run_on_relative_paths_records_absolute_ones(tmp_path):
    """A bundle is built later, from anywhere; a path relative to where the series ran
    resolves to nothing there."""
    import os
    from lunaNMR.cli import main
    here = os.getcwd()
    os.chdir(tmp_path)
    try:
        assert main(["series", "--spectra", os.path.relpath(DNAJA1, tmp_path),
                     "--peaks", os.path.relpath(DNAJA1 / "DNAJA1.txt", tmp_path),
                     "--out", "rel", "--mode", "titration", "--peak-source", "cascade"]) == 0
    finally:
        os.chdir(here)
    meta = json.loads((tmp_path / "rel" / "batch_results.json").read_text())["metadata"]
    for key in ("data_folder", "output_folder", "csv_path", "peak_list_file"):
        assert Path(meta[key]).is_absolute(), key
    assert Path(meta["data_folder"]) == DNAJA1


SERIES_FILES = ("series_analysis_tidy.csv", "peak_intensity_matrix.csv",
                "comprehensive_peak_tracking.csv", "series_metadata.json")


def _copy_of_run(run_dir, dest):
    """A run folder copied elsewhere, with each series CSV marked so a copy taken from
    the original folder is told apart from one taken from this folder."""
    import shutil
    shutil.copytree(run_dir, dest)
    for name in SERIES_FILES:
        path = dest / "series" / "DNAJA1" / name
        path.write_text(path.read_text() + "\n")
    return dest


def test_a_copied_run_bundles_its_own_csvs_not_the_originals(run_dir, tmp_path):
    """batch_results.json records where the series first wrote; a copy of the run still
    names the original folder, which may since hold a different run."""
    from lunaNMR.cli import main
    moved = _copy_of_run(run_dir, tmp_path / "moved")
    out = tmp_path / "B.lunaNMR"
    assert main(["project", "build", str(moved), "--out", str(out)]) == 0
    for name in SERIES_FILES:
        assert ((out / "series_results" / "DNAJA1" / name).read_text()
                == (moved / "series" / "DNAJA1" / name).read_text()), name


def test_a_moved_run_still_bundles_its_csvs(run_dir, tmp_path):
    """A moved run's recorded folder is gone; its CSVs sit beside its batch_results.json."""
    from lunaNMR.cli import main
    moved = _copy_of_run(run_dir, tmp_path / "moved")
    batch = moved / "series" / "DNAJA1" / "batch_results.json"
    data = json.loads(batch.read_text())
    gone = tmp_path / "gone"
    data["metadata"]["output_folder"] = str(gone)
    data["metadata"]["csv_path"] = str(gone / "series_analysis_tidy.csv")
    batch.write_text(json.dumps(data))
    out = tmp_path / "B.lunaNMR"
    assert main(["project", "build", str(moved), "--out", str(out)]) == 0
    for name in SERIES_FILES:
        assert (out / "series_results" / "DNAJA1" / name).exists(), name


@pytest.mark.parametrize("folders", [
    ("600/titration", "titration"),
    ("a/x", "b/x", "x"),
    ("a/x", "b/x", "b_x"),
])
def test_series_folders_of_one_name_are_all_found(tmp_path, folders):
    """A series is named after its folder, then after its place under the run; neither
    may be one another series already holds, or that series leaves the bundle unsaid."""
    from lunaNMR.utils.project_builder import find_series
    for folder in folders:
        (tmp_path / folder).mkdir(parents=True)
        (tmp_path / folder / "batch_results.json").write_text("{}")
    found = find_series(tmp_path)
    assert sorted(p.parent.relative_to(tmp_path).as_posix() for p in found.values()) \
        == sorted(folders)


@pytest.mark.parametrize("dry_run", [False, True])
def test_a_spectrum_that_does_not_exist_is_refused(run_dir, tmp_path, dry_run):
    """The bundle reopens its spectrum by path; a path to nothing is a bundle that cannot."""
    from lunaNMR.cli import main
    out = tmp_path / "C.lunaNMR"
    assert main(["project", "build", str(run_dir), "--out", str(out),
                 "--spectrum", str(tmp_path / "nonexistent.ft"),
                 *(["--dry-run"] if dry_run else [])]) == 1
    assert not out.exists()


def _run_with_data_folder(run_dir, dest, data_folder):
    """A copy of the run whose series records `data_folder` as where its spectra are."""
    import shutil
    shutil.copytree(run_dir, dest)
    batch = dest / "series" / "DNAJA1" / "batch_results.json"
    data = json.loads(batch.read_text())
    data["metadata"]["data_folder"] = str(data_folder)
    batch.write_text(json.dumps(data))
    return dest


@pytest.mark.parametrize("dry_run", [False, True])
def test_a_default_spectrum_that_does_not_exist_is_refused(run_dir, tmp_path, dry_run):
    """With no --spectrum the bundle records the series' reference spectrum; when the
    spectra are no longer where the series read them, that is a path to nothing."""
    import subprocess
    moved = _run_with_data_folder(run_dir, tmp_path / "moved", tmp_path / "gone")
    out = tmp_path / "D.lunaNMR"
    p = subprocess.run([sys.executable, "-m", "lunaNMR", "project", "build", str(moved),
                        "--out", str(out), "--format", "json",
                        *(["--dry-run"] if dry_run else [])],
                       cwd=str(Path(__file__).resolve().parents[1]),
                       capture_output=True, text=True)
    assert p.returncode == 1, p.stdout
    payload = json.loads(p.stdout)
    assert payload["ok"] is False
    assert "--spectrum" in payload["error"]["message"]
    assert str(tmp_path / "gone") in payload["error"]["message"]
    assert not out.exists()


def test_the_same_run_builds_when_given_the_spectrum(run_dir, tmp_path):
    from lunaNMR.cli import main
    moved = _run_with_data_folder(run_dir, tmp_path / "moved", tmp_path / "gone")
    out = tmp_path / "D.lunaNMR"
    spectrum = DNAJA1 / "2D_DNAJA1_0o0.ft"
    assert main(["project", "build", str(moved), "--out", str(out),
                 "--spectrum", str(spectrum)]) == 0
    manifest = json.loads((out / "project.json").read_text())
    assert manifest["file_references"]["nmr_spectrum"] == str(spectrum.resolve())
