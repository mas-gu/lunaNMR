# ABOUTME: The project-bundle tools of the MCP server: build a bundle from a run, list it, export it.
# ABOUTME: Real series and Kd runs feed the bundle; remove and overwrite are deliberately not exposed.

import pytest

pytest.importorskip("mcp")

from lunaNMR.cli import build_parser
from lunaNMR.mcp import cli_call
from mcp_helpers import (BINDING, call_result, call_tool, list_tools, stage_binding_dataset,
                         wait_for_job)

needs_binding_data = pytest.mark.skipif(
    not BINDING.is_dir(), reason="data_test/data_binding is absent")


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    root = stage_binding_dataset(tmp_path_factory.mktemp("dataset"))
    series = call_tool(root, "start_series", {
        "spectra": "series", "peaks": "series/DNAJA1.txt", "out": "runs/s1",
        "mode": "titration", "peak_source": "cascade"}, wait_seconds=1)
    wait_for_job(root, series["job"], 300)
    kd = call_tool(root, "kd_fit", {
        "input": "runs/s1/series_analysis_tidy.csv", "p0": 50, "conc_units": "equivalents",
        "out": "runs/kd1"}, wait_seconds=1)
    wait_for_job(root, kd["job"], 300)
    return root


BUILD = {"run_dir": "runs", "out": "bundle.lunaNMR", "spectrum": "series/2D_DNAJA1_0o0.ft",
         "peaks": "series/DNAJA1.txt"}


@pytest.fixture(scope="module")
def built(dataset):
    """The bundle a build produced, so each test that reads it can run on its own."""
    return call_tool(dataset, "project_build", BUILD)


def test_argument_vectors_are_accepted_by_the_cli_parser():
    parser = build_parser()
    parser.parse_args(cli_call.project_build_args("run", "b.lunaNMR", "s.ft", "p.txt"))
    parser.parse_args(cli_call.project_inventory_args("b.lunaNMR"))
    parser.parse_args(cli_call.project_export_args("b.lunaNMR", "out"))


def test_remove_and_overwrite_are_not_offered(tmp_path):
    tools = {tool.name: tool for tool in list_tools(tmp_path)}
    assert "project_remove" not in tools
    assert "force" not in tools["project_export"].input_schema["properties"]


@needs_binding_data
def test_a_bundle_is_built_from_a_series_and_a_kd_run(dataset, built):
    summary = built
    assert summary["ok"] is True
    assert summary["series"] == ["s1"]
    assert summary["kd"] == ["kd"]
    assert (dataset / "bundle.lunaNMR" / "project.json").is_file()


@needs_binding_data
def test_building_into_an_existing_bundle_is_refused(dataset, built):
    summary = call_tool(dataset, "project_build", {**BUILD, "out": "bundle.lunaNMR"})
    assert summary["error"]["type"] == "PathRefused"


@needs_binding_data
@pytest.mark.parametrize("field", ["spectrum", "peaks"])
def test_the_recorded_spectrum_and_peak_list_must_be_inside_the_root(dataset, tmp_path, field):
    outside = tmp_path / "outside.txt"
    outside.write_text("x")
    summary = call_tool(dataset, "project_build", {
        **BUILD, "out": f"bundle_{field}.lunaNMR", field: str(outside)})
    assert summary["error"]["type"] == "PathRefused"
    assert not (dataset / f"bundle_{field}.lunaNMR").exists()


@needs_binding_data
def test_the_spectrum_and_peak_list_cannot_be_left_to_the_runs_own_record(dataset):
    arguments = {key: value for key, value in BUILD.items() if key != "peaks"}
    result = call_result(dataset, "project_build", {**arguments, "out": "bundle_nopeaks.lunaNMR"})
    assert result.is_error is True


@needs_binding_data
def test_the_inventory_lists_what_the_bundle_holds(dataset, built):
    summary = call_tool(dataset, "project_inventory", {"bundle": "bundle.lunaNMR"})
    assert summary["ok"] is True
    labels = [category["label"] for category in summary["categories"]]
    assert "Project manifest" in labels


@needs_binding_data
def test_an_export_copies_the_bundle_to_a_new_folder_only(dataset, built):
    summary = call_tool(dataset, "project_export", {"bundle": "bundle.lunaNMR", "out": "plain"})
    assert summary["ok"] is True
    assert summary["n_files"] > 0 and (dataset / "plain").is_dir()
    again = call_tool(dataset, "project_export", {"bundle": "bundle.lunaNMR", "out": "plain"})
    assert again["error"]["type"] == "PathRefused"
