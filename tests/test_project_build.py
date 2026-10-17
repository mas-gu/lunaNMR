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
