# ABOUTME: The MCP tools that plan a series, survey a titration and export a Kd fit, run on real data.
# ABOUTME: Each tool's argument vector is also pinned to the real CLI parser, so a renamed flag fails here.

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("mcp")

from lunaNMR.cli import build_parser
from lunaNMR.mcp import cli_call
from mcp_helpers import call_tool

BINDING = Path(__file__).resolve().parents[2] / "data_test" / "data_binding"
needs_binding_data = pytest.mark.skipif(
    not BINDING.is_dir(), reason="data_test/data_binding is absent")


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    """A private copy of the DNAJA1 titration: the tools write beside their inputs."""
    root = tmp_path_factory.mktemp("dataset")
    series = root / "series"
    series.mkdir()
    source = BINDING / "DNAJA1_HSPA8"
    for spectrum in source.glob("*.ft"):
        shutil.copy(spectrum, series / spectrum.name)
    shutil.copy(source / "DNAJA1.txt", series / "DNAJA1.txt")
    shutil.copy(BINDING / "tidy" / "DNAJA1_HSPA8_tidy.csv", root / "tidy.csv")
    return root


@pytest.fixture(scope="module")
def fit_json(dataset):
    """A Kd fit made by the real CLI, for the export tool to read."""
    completed = subprocess.run(
        [sys.executable, "-m", "lunaNMR", "kd", "--input", str(dataset / "tidy.csv"),
         "--out", str(dataset / "fit"), "--p0", "50", "--conc-units", "equivalents",
         "--format", "json"],
        capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr[-500:]
    return next((dataset / "fit" / "data").glob("*_kd_fit_data.json"))


def test_every_argument_vector_is_accepted_by_the_cli_parser():
    parser = build_parser()
    for mode, source, deep in [("time", "reference", True), ("titration", "cascade", False)]:
        parser.parse_args(cli_call.plan_series_args(
            "spectra", "peaks.txt", "out", mode=mode, peak_source=source, deep=deep))
    for conc in (None, [0, 25, 50]):
        parser.parse_args(cli_call.kd_survey_args(
            "tidy.csv", 50, "absolute", "out", conc=conc))
    for kinds in (None, ["curves", "kd-bars"]):
        parser.parse_args(cli_call.export_kd_args("fit.json", "out", kinds=kinds))


@needs_binding_data
def test_plan_series_checks_the_inputs_and_creates_nothing(dataset):
    summary = call_tool(dataset, "plan_series", {
        "spectra": "series", "peaks": "series/DNAJA1.txt", "out": "plan",
        "mode": "titration", "peak_source": "cascade"})
    assert summary["ok"] is True
    assert summary["dry_run"] is True
    assert summary["spectra_found"] == 7
    assert not (dataset / "plan").exists()


@needs_binding_data
def test_plan_series_leaves_out_the_parameter_block_and_the_assignment_list(dataset):
    summary = call_tool(dataset, "plan_series", {
        "spectra": "series", "peaks": "series/DNAJA1.txt", "out": "plan_trimmed"})
    assert "params" not in summary
    assert summary["checks"]["peak_list"]["n_peaks"] == 68
    assert "assignments" not in summary["checks"]["peak_list"]
    assert summary["spectra_found"] == 7


@needs_binding_data
def test_plan_series_refuses_an_output_folder_that_exists(dataset):
    summary = call_tool(dataset, "plan_series", {
        "spectra": "series", "peaks": "series/DNAJA1.txt", "out": "series"})
    assert summary["ok"] is False
    assert summary["error"]["type"] == "PathRefused"
    assert "already exists" in summary["error"]["message"]


@needs_binding_data
def test_plan_series_refuses_a_peak_list_outside_the_root(dataset):
    summary = call_tool(dataset, "plan_series", {
        "spectra": "series", "peaks": "../DNAJA1.txt", "out": "plan"})
    assert summary["ok"] is False
    assert "outside the data root" in summary["error"]["message"]


@needs_binding_data
def test_kd_survey_selects_residues_and_writes_beside_the_input(dataset):
    summary = call_tool(dataset, "kd_survey", {
        "input": "tidy.csv", "p0": 50, "conc_units": "equivalents", "out": "survey"})
    assert summary["ok"] is True
    assert summary["n_total"] > 0
    assert Path(summary["residues_file"]).parent == dataset.resolve()
    assert (dataset / "survey").is_dir()


@needs_binding_data
def test_kd_survey_refuses_an_output_folder_that_exists(dataset):
    summary = call_tool(dataset, "kd_survey", {
        "input": "tidy.csv", "p0": 50, "conc_units": "equivalents", "out": "series"})
    assert summary["ok"] is False
    assert summary["error"]["type"] == "PathRefused"


@needs_binding_data
def test_export_kd_writes_figures_for_the_requested_kinds(dataset, fit_json):
    summary = call_tool(dataset, "export_kd", {
        "fit_json": str(fit_json.relative_to(dataset)), "out": "figures",
        "kinds": ["kd-bars"]})
    assert summary["ok"] is True
    assert any(figure.endswith("kd_vs_residue.pdf") for figure in summary["figures"])
    assert all("ref_vs_point" not in figure for figure in summary["figures"])


@needs_binding_data
def test_export_kd_refuses_an_output_folder_that_exists(dataset, fit_json):
    summary = call_tool(dataset, "export_kd", {
        "fit_json": str(fit_json.relative_to(dataset)), "out": "series"})
    assert summary["ok"] is False
    assert summary["error"]["type"] == "PathRefused"
