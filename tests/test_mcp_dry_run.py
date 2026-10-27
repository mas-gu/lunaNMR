# ABOUTME: `dry_run` on the analysis tools checks the inputs and reports the plan, running and writing nothing.
# ABOUTME: Real staged inputs; a dry run answers at once with no job and no output folder.

import pytest

pytest.importorskip("mcp")

from mcp_helpers import (BINDING, call_tool, relaxation_data_present, stage_binding_dataset,
                         stage_relaxation_dataset)

pytestmark = pytest.mark.skipif(
    not (BINDING.is_dir() and relaxation_data_present()),
    reason="the binding and relaxation test data are both needed")

KD = {"input": "tidy.csv", "p0": 50, "conc_units": "equivalents"}
RELAXATION = {"f1_t1": "in/600_T1_matrix.csv", "f1_t2": "in/600_T2_matrix.csv",
              "f1_noe_sat": "in/600_noe_sat.csv", "f1_noe_unsat": "in/600_noe_unsat.csv"}

CASES = {
    "kd_fit": KD,
    "kd_survey": KD,
    "peaks_shift": {"peaks": "series/DNAJA1.txt", "dx": 0.01},
    "dynamixs_t1t2": {"input": "in/600_T1_matrix.csv", "exp": "T1"},
    "dynamixs_methyl_t2": {"input": "in/methyl_matrix.csv"},
    "dynamixs_hetnoe": {"sat": "in/600_noe_sat.csv", "unsat": "in/600_noe_unsat.csv"},
    "dynamixs_modelfree": RELAXATION,
}


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    root = tmp_path_factory.mktemp("dryrun")
    stage_binding_dataset(root)
    stage_relaxation_dataset(root)
    return root


@pytest.mark.parametrize("tool", sorted(CASES))
def test_a_dry_run_answers_at_once_and_writes_nothing(dataset, tool):
    out = f"planned_{tool}.txt" if tool == "peaks_shift" else f"planned_{tool}"
    summary = call_tool(dataset, tool, {**CASES[tool], "out": out, "dry_run": True})
    assert summary["ok"] is True, summary
    assert summary.get("dry_run") is True
    assert "job" not in summary
    assert not (dataset / out).exists()


def test_a_dry_run_still_refuses_what_a_real_run_would_refuse(dataset):
    summary = call_tool(dataset, "dynamixs_t1t2", {
        "input": "in/does_not_exist.csv", "exp": "T1", "out": "planned_missing", "dry_run": True})
    assert summary["ok"] is False
    assert summary["error"]["type"] == "PathRefused"


def test_named_residues_are_checked_but_not_sent_to_a_dry_run(dataset):
    summary = call_tool(dataset, "kd_fit", {
        **KD, "out": "planned_names", "residues": ["A17", "A25"], "dry_run": True})
    assert summary["ok"] is True and summary.get("dry_run") is True
    bad = call_tool(dataset, "kd_fit", {
        **KD, "out": "planned_bad", "residues": ["../x"], "dry_run": True})
    assert bad["error"]["type"] == "InvalidResidueName"
