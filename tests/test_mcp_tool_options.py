# ABOUTME: The options of the series, Kd and export tools beyond the basics, each with an observable effect.
# ABOUTME: Every argument vector is also pinned to the real parser, so a renamed flag fails here.

import json

import pytest

pytest.importorskip("mcp")

from lunaNMR.cli import build_parser
from lunaNMR.mcp import cli_call
from mcp_helpers import BINDING, call_tool, stage_binding_dataset, wait_for_job

needs_binding_data = pytest.mark.skipif(
    not BINDING.is_dir(), reason="data_test/data_binding is absent")

KD = {"input": "tidy.csv", "p0": 50, "conc_units": "equivalents"}


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    return stage_binding_dataset(tmp_path_factory.mktemp("dataset"))


def test_every_added_argument_vector_is_accepted_by_the_cli_parser():
    parser = build_parser()
    parser.parse_args(cli_call.diagnose_args("d", quick=True, mode="titration", sample=True,
                                             strict=True))
    parser.parse_args(cli_call.plan_series_args("s", "p", "o", quick=True, sample=True,
                                                params_file="params.json"))
    parser.parse_args(cli_call.start_series_args("s", "p", "o", params_file="params.json"))
    kd_options = dict(prefix="run1", alpha=0.2, intensity_from="volume", bootstrap=50,
                      csp_sigma_multiple=2.0, kd_outlier_z=0.0, noise_quantile=0.3,
                      dd_runaway_ratio=5.0, ref_max_ratio=8.0)
    parser.parse_args(cli_call.kd_survey_args(
        "t.csv", 50, "absolute", "o", conc=[0, 1], selection_file="sel.txt", observables=["csp"],
        **kd_options))
    parser.parse_args(cli_call.kd_fit_args(
        "t.csv", 50, "absolute", "o", conc=[0, 1], observables=["csp", "intensity"], **kd_options))
    parser.parse_args(cli_call.export_kd_args(
        "f.json", "o", kinds=["curves"], observables=["csp"], fig_format=["pdf", "png"],
        per_page=10, prefix="p", summary_only=True))


@needs_binding_data
def test_diagnose_sample_assesses_only_the_first_and_last_spectrum(dataset):
    summary = call_tool(dataset, "diagnose", {
        "path": "series", "quick": True, "sample": True, "mode": "titration"})
    assert summary["ok"] is True
    assert len(summary["experiments"][0]["spectra"]) == 2


@needs_binding_data
def test_diagnose_strict_fails_on_a_warning(dataset):
    arguments = {"path": "series", "quick": True, "mode": "time"}
    assert call_tool(dataset, "diagnose", arguments)["ok"] is True
    assert call_tool(dataset, "diagnose", {**arguments, "strict": True})["ok"] is False


@needs_binding_data
def test_plan_series_reads_a_parameter_file_inside_the_root(dataset):
    (dataset / "good.json").write_text(json.dumps({"gui_params": {"centroid_window_x_ppm": 0.03}}))
    (dataset / "bad.json").write_text(json.dumps({"no_such_group": {"x": 1}}))
    base = {"spectra": "series", "peaks": "series/DNAJA1.txt", "mode": "titration",
            "peak_source": "cascade", "quick": True, "sample": True}
    good = call_tool(dataset, "plan_series", {**base, "out": "plan_good", "params_file": "good.json"})
    bad = call_tool(dataset, "plan_series", {**base, "out": "plan_bad", "params_file": "bad.json"})
    assert good["ok"] is True
    assert bad["ok"] is False and "no_such_group" in json.dumps(bad)


@needs_binding_data
def test_a_parameter_file_outside_the_root_is_refused(dataset, tmp_path):
    outside = tmp_path / "outside.json"
    outside.write_text("{}")
    summary = call_tool(dataset, "plan_series", {
        "spectra": "series", "peaks": "series/DNAJA1.txt", "out": "plan_out",
        "params_file": str(outside)})
    assert summary["error"]["type"] == "PathRefused"


@needs_binding_data
def test_kd_options_reach_the_fit_and_the_prefix_names_its_files(dataset):
    started = call_tool(dataset, "kd_fit", {
        **KD, "out": "runs/kd_opts", "prefix": "mykd", "alpha": 0.2, "kd_outlier_z": 0.0,
        "csp_sigma_multiple": 2.0, "intensity_from": "volume"})
    status = wait_for_job(dataset, started["job"])
    assert status["status"] == "done"
    fit_file = dataset / "runs" / "kd_opts" / "data" / "mykd_kd_fit_data.json"
    metadata = json.loads(fit_file.read_text())["metadata"]
    assert metadata["alpha"] == 0.2
    assert metadata["intensity_value"] == "volume"
    assert metadata["csp_significance"]["multiple"] == 2.0


@needs_binding_data
def test_kd_survey_writes_its_selection_where_told_inside_the_root(dataset):
    summary = call_tool(dataset, "kd_survey", {
        **KD, "out": "survey_sel", "selection_file": "my_selection.txt", "prefix": "sv"})
    assert summary["ok"] is True
    assert (dataset / "my_selection.txt").is_file()


@needs_binding_data
def test_kd_survey_refuses_a_selection_file_outside_the_root(dataset, tmp_path):
    summary = call_tool(dataset, "kd_survey", {
        **KD, "out": "survey_out", "selection_file": str(tmp_path / "elsewhere.txt")})
    assert summary["error"]["type"] == "PathRefused"


@needs_binding_data
def test_export_summary_only_writes_no_figures_and_prefix_names_the_files(dataset):
    started = call_tool(dataset, "kd_fit", {**KD, "out": "runs/kd_for_export"})
    wait_for_job(dataset, started["job"])
    fit = next((dataset / "runs" / "kd_for_export").glob("data/*_kd_fit_data.json"))
    summary = call_tool(dataset, "export_kd", {
        "fit_json": str(fit.relative_to(dataset)), "out": "export_summary",
        "summary_only": True, "prefix": "only"})
    assert summary["ok"] is True
    names = sorted(p.name for p in (dataset / "export_summary").rglob("*") if p.is_file())
    assert names and all(name.startswith("only") for name in names)
    assert not any(name.endswith(".pdf") for name in names)
