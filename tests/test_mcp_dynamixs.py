# ABOUTME: The relaxation tools of the MCP server (T1/T2, methyl T2, hetNOE, T1rho, density, model-free).
# ABOUTME: Real relaxation test inputs from an earlier series run; each result is checked against the files it wrote.

import csv
import json

import pytest

pytest.importorskip("mcp")

from lunaNMR.cli import build_parser
from lunaNMR.mcp import cli_call
from mcp_helpers import call_tool, relaxation_data_present, stage_relaxation_dataset, wait_for_job

needs_relaxation_data = pytest.mark.skipif(
    not relaxation_data_present(), reason="the relaxation relaxation test data is absent")

WAIT = 180


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    return stage_relaxation_dataset(tmp_path_factory.mktemp("relaxation"))


def _run(dataset, tool, arguments, wait_seconds=WAIT):
    summary = call_tool(dataset, tool, arguments, wait_seconds=wait_seconds)
    if summary.get("status") == "running":
        summary = wait_for_job(dataset, summary["job"], 600)
    return summary


@pytest.fixture(scope="module")
def modelfree_single(dataset):
    return _run(dataset, "dynamixs_modelfree", {
        "f1_t1": "in/600_T1_matrix.csv", "f1_t2": "in/600_T2_matrix.csv",
        "f1_noe_sat": "in/600_noe_sat.csv", "f1_noe_unsat": "in/600_noe_unsat.csv",
        "field1_freq": 600, "out": "mf600"})


def test_every_argument_vector_is_accepted_by_the_cli_parser():
    parser = build_parser()
    parser.parse_args(cli_call.dynamixs_t1t2_args(
        "m.csv", "o", "T1", prefix="p", field_name="f", field_freq=700, time_units="s",
        allow_unparsed_delays=True, error_method="bootstrap", bootstrap=200))
    parser.parse_args(cli_call.dynamixs_methyl_t2_args(
        "m.csv", "o", prefix="p", field_name="f", field_freq=700, time_units="ms",
        error_method="analytical"))
    parser.parse_args(cli_call.dynamixs_hetnoe_args("s.csv", "u.csv", "o", prefix="p"))
    parser.parse_args(cli_call.dynamixs_t1rho_args(
        "r.csv", "t.csv", "p.txt", 2000.0, 118.0, 90.0, "o", prefix="p", field_freq=600,
        time_units="us", allow_unparsed_delays=True, error_method="analytical"))
    parser.parse_args(cli_call.dynamixs_density_args(
        "t.csv", "o", input2="t2.csv", prefix="p", dual=True, use_087=False, monte_carlo=True,
        n_samples=50, plot=False, parallel=False, field1_freq=600, field2_freq=700, rnh=1.02,
        csa=-170.0))
    parser.parse_args(cli_call.dynamixs_modelfree_args(
        "o", f1_t1="a", f1_t2="b", f1_noe_sat="c", f1_noe_unsat="d", f2_t1="e", f2_r2_table="f",
        f2_noe_sat="g", f2_noe_unsat="h", f1_t1_units="s", f1_t2_units="ms", f2_t1_units="us",
        f2_t2_units="ms", allow_unparsed_delays=True, prefix="p", dual=True, method="dual_jwh",
        init_amp=0.5, init_t1=800.0, init_t2=100.0, n_bootstrap=100, error_method="bootstrap",
        n_monte_carlo=10, field1_freq=600, field2_freq=700, rnh=1.015, csa=-172.0))


@needs_relaxation_data
def test_t1_fit_reports_the_residues_that_fail_the_usual_gate(dataset):
    summary = _run(dataset, "dynamixs_t1t2", {
        "input": "in/600_T1_matrix.csv", "out": "t1_600", "exp": "T1", "field_freq": 600})
    assert summary["status"] == "done"
    assert summary["result"]["n_fitted"] == 60
    table = dataset / "t1_600" / "field1_fit_results.txt"
    rows = list(csv.DictReader(open(table), delimiter="\t"))
    failing = [r["Residue"] for r in rows
               if float(r["T1_err"]) >= float(r["T1"]) or r["Success"] != "Yes"]
    quality = summary["relaxation_quality"]
    assert quality["n_residues"] == len(rows)
    assert quality["n_failing_gate"] == len(failing)
    assert sorted(quality["residues_failing_gate"]) == sorted(failing)[:30]


@needs_relaxation_data
def test_t1_fit_names_the_residues_whose_delay_window_is_marginal(dataset):
    summary = _run(dataset, "dynamixs_t1t2", {
        "input": "in/600_T1_matrix.csv", "out": "t1_600_window", "exp": "T1", "field_freq": 600})
    table = dataset / "t1_600_window" / "field1_fit_results.txt"
    rows = list(csv.DictReader(open(table), delimiter="\t"))
    marginal = [r["Residue"] for r in rows if float(r["WindowRatio"]) < 2]
    assert marginal, "the data must hold a marginal residue for this test to mean anything"
    quality = summary["relaxation_quality"]
    assert quality["n_window_marginal"] == summary["result"]["n_window_marginal"]
    assert sorted(quality["residues_window_marginal"]) == sorted(marginal)[:30]


@needs_relaxation_data
def test_a_tidy_csv_is_refused_with_the_reason(dataset):
    summary = _run(dataset, "dynamixs_t1t2", {
        "input": "in/tidy.csv", "out": "t1_tidy", "exp": "T1"})
    assert summary["status"] == "failed"
    assert "matrix" in json.dumps(summary["result"]["error"]).lower()


@needs_relaxation_data
@pytest.mark.parametrize("field", ["prefix", "field_name"])
def test_names_that_become_filenames_cannot_name_a_path(dataset, field):
    summary = call_tool(dataset, "dynamixs_t1t2", {
        "input": "in/600_T1_matrix.csv", "out": f"t1_bad_{field}", "exp": "T1", field: "../x"})
    assert summary["error"]["type"] == "InvalidName"
    assert not (dataset / f"t1_bad_{field}").exists()


@needs_relaxation_data
def test_hetnoe_reports_its_median_and_the_implausible_residues(dataset):
    summary = _run(dataset, "dynamixs_hetnoe", {
        "sat": "in/600_noe_sat.csv", "unsat": "in/600_noe_unsat.csv", "out": "noe_600"})
    assert summary["status"] == "done"
    rows = list(csv.DictReader(open(dataset / "noe_600" / "field1_hetnoe.csv")))
    values = sorted(float(r["hetNOE"]) for r in rows)
    quality = summary["hetnoe_quality"]
    assert quality["n_residues"] == len(values)
    assert quality["median"] == pytest.approx(values[len(values) // 2], abs=0.02)
    assert quality["n_outside_plausible"] == sum(1 for v in values if not -0.2 <= v <= 1.05)


@needs_relaxation_data
def test_single_field_model_free_runs_and_reports_its_attrition(modelfree_single, dataset):
    assert modelfree_single["status"] == "done"
    assert modelfree_single["result"]["n_successful"] == 55
    assert modelfree_single["result"]["n_at_iteration_cap"] == 0
    assert any(name.endswith("spectral_density_basic.csv") for name in modelfree_single["files"])


@needs_relaxation_data
def test_spectral_density_runs_on_a_table_a_model_free_run_wrote(modelfree_single, dataset):
    summary = _run(dataset, "dynamixs_density", {
        "input": "mf600/modelfree_field1_input.csv", "out": "density_600", "field1_freq": 600,
        "plot": False})
    assert summary["status"] == "done"
    assert summary["result"]["n_residues"] == 55


@needs_relaxation_data
def test_dual_field_model_free_outlasts_the_wait_and_finishes_as_a_job(dataset):
    started = call_tool(dataset, "dynamixs_modelfree", {
        "dual": True, "f1_t1": "in/600_T1_matrix.csv", "f1_t2": "in/600_T2_matrix.csv",
        "f1_noe_sat": "in/600_noe_sat.csv", "f1_noe_unsat": "in/600_noe_unsat.csv",
        "field1_freq": 600, "f2_t1": "in/700_T1_matrix.csv", "f2_t2": "in/700_T2_matrix.csv",
        "f2_noe_sat": "in/700_noe_sat.csv", "f2_noe_unsat": "in/700_noe_unsat.csv",
        "field2_freq": 700, "out": "mf_dual"}, wait_seconds=1)
    assert started["status"] == "running"
    status = wait_for_job(dataset, started["job"], 600)
    assert status["status"] == "done"
    assert status["result"]["method"] == "dual_087"
    assert status["result"]["n_at_iteration_cap"] == 0


@needs_relaxation_data
def test_methyl_t2_relays_how_many_fits_were_unidentifiable(dataset):
    summary = _run(dataset, "dynamixs_methyl_t2", {
        "input": "in/methyl_matrix.csv", "out": "methyl", "field_freq": 600, "time_units": "ms"})
    assert summary["status"] == "done"
    assert summary["result"]["n_total"] == 77
    assert summary["result"]["n_bi_exp_unidentifiable"] > 0


@needs_relaxation_data
def test_t1rho_runs_its_plumbing_with_a_stand_in_matrix(dataset):
    """No real T1rho series exists in the test data: a real T2 decay matrix stands in for it,
    so this checks the inputs and the output files, not the physics."""
    summary = _run(dataset, "dynamixs_t1rho", {
        "input": "in/600_T2_matrix.csv", "t1": "in/600_T1_matrix.csv", "peaks": "in/peaks.txt",
        "omega1": 2000.0, "carrier": 118.0, "theta": 90.0, "out": "t1rho"})
    assert summary["status"] == "done", summary
    assert any(name.endswith("r2_from_t1rho.csv") for name in summary["files"])


@needs_relaxation_data
@pytest.mark.parametrize("field", ["f1_t1", "f1_t2", "f1_noe_sat", "f1_noe_unsat", "f1_r2_table",
                                   "f2_t1", "f2_t2", "f2_noe_sat", "f2_noe_unsat", "f2_r2_table"])
def test_every_model_free_input_path_is_confined_to_the_root(dataset, tmp_path, field):
    outside = tmp_path / "outside.csv"
    outside.write_text("Residue\n")
    arguments = {"f1_t1": "in/600_T1_matrix.csv", "f1_noe_sat": "in/600_noe_sat.csv",
                 "f1_noe_unsat": "in/600_noe_unsat.csv", "out": f"mf_guard_{field}"}
    arguments[field] = str(outside)
    summary = call_tool(dataset, "dynamixs_modelfree", arguments)
    assert summary["error"]["type"] == "PathRefused"
    assert not (dataset / f"mf_guard_{field}").exists()
