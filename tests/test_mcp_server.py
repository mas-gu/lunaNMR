# ABOUTME: Drives the lunaNMR MCP server through the SDK's in-process client against real spectra.
# ABOUTME: Also pins every CLI argument vector the server builds to the real parser, so a renamed flag fails here.

import asyncio
import json
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("mcp")

from mcp.client import Client
from mcp.client.stdio import StdioServerParameters

from lunaNMR.cli import build_parser
from lunaNMR.mcp import cli_call
from lunaNMR.mcp.data_root import DataRoot
from lunaNMR.mcp.server import build_server
from mcp_helpers import call_result, call_tool, list_tools

REPO_ROOT = Path(__file__).resolve().parents[1]
BINDING_ROOT = REPO_ROOT.parent / "data_test" / "data_binding"
needs_binding_data = pytest.mark.skipif(
    not BINDING_ROOT.is_dir(), reason="data_test/data_binding is absent")


def test_every_argument_vector_the_server_builds_is_accepted_by_the_cli_parser(tmp_path):
    parser = build_parser()
    for quick in (False, True):
        for mode in ("time", "titration"):
            argv = cli_call.diagnose_args(tmp_path, quick=quick, mode=mode)
            parser.parse_args(argv)


def test_the_server_lists_diagnose(tmp_path):
    async def run():
        async with Client(build_server(DataRoot(tmp_path))) as client:
            return [tool.name for tool in (await client.list_tools()).tools]
    assert "diagnose" in asyncio.run(run())


def test_a_path_outside_the_root_comes_back_as_an_error_envelope(tmp_path):
    root = tmp_path / "dataset"
    root.mkdir()
    summary = call_tool(root, "diagnose", {"path": "../elsewhere"})
    assert summary["ok"] is False
    assert summary["error"]["type"] == "PathRefused"
    assert "outside the data root" in summary["error"]["message"]


def test_a_missing_path_comes_back_as_an_error_envelope(tmp_path):
    summary = call_tool(tmp_path, "diagnose", {"path": "nothing_here"})
    assert summary["ok"] is False
    assert summary["error"]["type"] == "PathRefused"


def test_the_cli_failure_envelope_is_relayed(tmp_path):
    summary = call_tool(tmp_path, "diagnose", {"path": "."})
    assert summary["ok"] is False
    assert summary["error"]["message"]


@needs_binding_data
def test_diagnose_runs_the_real_cli_on_real_spectra():
    summary = call_tool(BINDING_ROOT, "diagnose", {"path": "DNAJA1_HSPA8", "quick": True})
    assert summary["ok"] is True
    assert summary["command"] == "diagnose"
    assert summary["experiments"][0]["n_spectra"] == 7


def test_the_installed_entry_point_serves_over_stdio(tmp_path):
    async def run():
        server = StdioServerParameters(
            command=sys.executable,
            args=['-m', 'lunaNMR.mcp.server', '--root', str(tmp_path)])
        async with Client(server) as client:
            return [tool.name for tool in (await client.list_tools()).tools]
    assert "diagnose" in asyncio.run(run())


def test_no_lunaNMR_module_imports_the_mcp_package():
    package = REPO_ROOT / "lunaNMR"
    importers = [
        path for path in package.rglob("*.py")
        if "mcp" not in path.relative_to(package).parts[:1]
        and "lunaNMR.mcp" in path.read_text(encoding="utf-8", errors="ignore")]
    assert importers == []


def test_a_missing_root_is_a_usage_error():
    completed = subprocess.run(
        [sys.executable, '-m', 'lunaNMR.mcp.server', '--root', '/no/such/folder'],
        capture_output=True, text=True)
    assert completed.returncode == 2
    assert "data root" in completed.stderr


@needs_binding_data
def test_diagnose_leaves_out_the_assignment_list_but_keeps_the_count():
    summary = call_tool(BINDING_ROOT, "diagnose", {"path": "DNAJA1_HSPA8", "quick": True})
    peak_list = summary["experiments"][0]["peak_list"]
    assert peak_list["n_peaks"] == 68
    assert "assignments" not in peak_list


def test_a_refused_path_sets_the_protocol_error_flag_and_keeps_the_reason(tmp_path):
    result = call_result(tmp_path, "diagnose", {"path": "../elsewhere"})
    assert result.is_error is True
    assert json.loads(result.content[0].text)["error"]["type"] == "PathRefused"


def test_a_failed_cli_run_sets_the_protocol_error_flag(tmp_path):
    result = call_result(tmp_path, "diagnose", {"path": "."})
    assert result.is_error is True
    assert json.loads(result.content[0].text)["ok"] is False


@needs_binding_data
def test_a_successful_run_is_not_an_error_and_carries_structured_content():
    result = call_result(BINDING_ROOT, "diagnose", {"path": "DNAJA1_HSPA8", "quick": True})
    assert result.is_error is False
    assert result.structured_content["command"] == "diagnose"
    assert result.structured_content == json.loads(result.content[0].text)


def test_tools_declare_whether_they_write(tmp_path):
    hints = {tool.name: tool.annotations for tool in list_tools(tmp_path)}
    assert hints["diagnose"].read_only_hint is True
    assert hints["plan_series"].read_only_hint is True
    assert hints["kd_survey"].read_only_hint is False
    assert hints["kd_survey"].destructive_hint is True
    assert hints["export_kd"].read_only_hint is False
    assert hints["export_kd"].destructive_hint is False
    assert hints["job_status"].read_only_hint is True
    for starter in ("start_series", "kd_fit"):
        assert hints[starter].read_only_hint is False
        assert hints[starter].destructive_hint is False
    assert hints["cancel_job"].destructive_hint is True
    assert all(annotation.open_world_hint is False for annotation in hints.values())


@needs_binding_data
def test_diagnose_reads_titration_filenames_as_points_when_asked():
    arguments = {"path": "DNAJA1_HSPA8", "quick": True}
    as_time = call_tool(BINDING_ROOT, "diagnose", arguments)
    as_titration = call_tool(BINDING_ROOT, "diagnose", {**arguments, "mode": "titration"})
    assert as_time["experiments"][0]["delays"]["n_parsed"] == 0
    assert as_titration["experiments"][0]["delays"]["n_parsed"] == 7
