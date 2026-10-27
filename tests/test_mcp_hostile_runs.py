# ABOUTME: MCP server against run folders a dataset controls: records, links and sidecars that point outside the root.
# ABOUTME: Nothing a tool reads back from a folder may open a file beyond --root, nor leak its contents.

import json

import pytest

pytest.importorskip("mcp")

from mcp_helpers import call_tool


@pytest.fixture
def root(tmp_path):
    folder = tmp_path / "dataset"
    folder.mkdir()
    return folder


@pytest.fixture
def outside(tmp_path):
    folder = tmp_path / "outside"
    folder.mkdir()
    return folder


def _finished_job(folder, kind):
    (folder / "mcp_job").mkdir(parents=True)
    (folder / "mcp_job" / "job.json").write_text(json.dumps(
        {"kind": kind, "args": [], "pid": 2 ** 22 + 7, "started": 1.0}))
    (folder / "mcp_job" / "result.json").write_text(json.dumps({"exit_code": 0, "finished": 2.0}))


def test_project_build_refuses_a_record_whose_files_lie_outside_the_root(root, outside):
    secret = outside / "secret.json"
    secret.write_text(json.dumps({"metadata": {}, "SECRET_JSON": "LEAKED"}))
    (root / "run").mkdir()
    (root / "run" / "x_dynamixs.json").write_text(json.dumps({
        "command": "dynamixs t1t2", "args": {"field_freq": 600, "exp": "T1", "input": "m.csv"},
        "summary": {"json_file": str(secret), "results_file": str(secret)}}))
    (root / "spectrum.ft2").write_text("x")
    (root / "peaks.txt").write_text("x")
    summary = call_tool(root, "project_build", {
        "run_dir": "run", "out": "b.lunaNMR", "spectrum": "spectrum.ft2", "peaks": "peaks.txt"})
    assert summary["ok"] is False
    assert summary["error"]["type"] == "PathRefused"
    assert not (root / "b.lunaNMR").exists()


def test_project_build_dry_run_refuses_the_same_record(root, outside):
    (root / "run").mkdir()
    (root / "run" / "x_dynamixs.json").write_text(json.dumps({
        "command": "dynamixs t1t2", "args": {}, "summary": {"r2_table": str(outside / "t.csv")}}))
    (root / "spectrum.ft2").write_text("x")
    (root / "peaks.txt").write_text("x")
    summary = call_tool(root, "project_build", {
        "run_dir": "run", "out": "b.lunaNMR", "spectrum": "spectrum.ft2", "peaks": "peaks.txt",
        "dry_run": True})
    assert summary["error"]["type"] == "PathRefused"


def test_job_status_does_not_read_a_results_table_that_is_a_link_out(root, outside):
    (outside / "secret.tsv").write_text("Residue\tT1\tT1_err\tSuccess\nSECRETRES\t5\t9\tYes\n")
    _finished_job(root / "evil", "t1t2")
    (root / "evil" / "x_fit_results.txt").symlink_to(outside / "secret.tsv")
    summary = call_tool(root, "job_status", {"job": "evil"})
    assert summary["ok"] is False
    assert "SECRETRES" not in json.dumps(summary)


def test_job_status_does_not_open_an_input_a_record_names_outside_the_root(root, outside):
    (outside / "secret.csv").write_text("GUESSME,1\n")
    _finished_job(root / "evil", "modelfree")
    (root / "evil" / "p_dynamixs.json").write_text(json.dumps({
        "args": {"f1_t1": str(outside / "secret.csv")},
        "summary": {"excluded_residues": {"field1": ["GUESSME", "NOPE"]}}}))
    summary = call_tool(root, "job_status", {"job": "evil"})
    assert "GUESSME" not in json.dumps(summary.get("modelfree_quality", {}))
    assert "NOPE" not in json.dumps(summary.get("modelfree_quality", {}))
    assert "summary_error" in summary


def test_a_malformed_record_does_not_hide_the_finished_job(root):
    _finished_job(root / "broken", "modelfree")
    (root / "broken" / "p_dynamixs.json").write_text("{not json")
    summary = call_tool(root, "job_status", {"job": "broken"})
    assert summary["ok"] is True
    assert summary["status"] == "done"
    assert "summary_error" in summary


def test_a_series_sidecar_that_is_a_link_out_blocks_a_relaxation_input(root, outside):
    (outside / "meta.json").write_text(json.dumps({"n_value_unparsed": 7}))
    (root / "m.csv").write_text("Assignment,0,10\nA1,1,0.5\n")
    (root / "series_metadata.json").symlink_to(outside / "meta.json")
    summary = call_tool(root, "dynamixs_t1t2", {"input": "m.csv", "out": "o", "exp": "T1"})
    assert summary["ok"] is False
    assert summary["error"]["type"] == "PathRefused"
    assert not (root / "o").exists()
