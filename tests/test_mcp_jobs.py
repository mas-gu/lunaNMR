# ABOUTME: Background jobs of the MCP server: a series fit and a Kd fit started, polled, cancelled and read back.
# ABOUTME: Real runs on a private copy of the DNAJA1 titration; job state lives on disk, not in the server.

import asyncio
import csv
import json
import os
import time
from pathlib import Path

import pytest

pytest.importorskip("mcp")

from mcp.client import Client

from lunaNMR.cli import build_parser
from lunaNMR.mcp import cli_call
from lunaNMR.mcp.data_root import DataRoot
from lunaNMR.mcp.server import build_server
from mcp_helpers import BINDING, call_tool, stage_binding_dataset, wait_for_job

needs_binding_data = pytest.mark.skipif(
    not BINDING.is_dir(), reason="data_test/data_binding is absent")

SERIES = {"spectra": "series", "peaks": "series/DNAJA1.txt",
          "mode": "titration", "peak_source": "cascade"}
KD = {"input": "tidy.csv", "p0": 50, "conc_units": "equivalents"}


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    return stage_binding_dataset(tmp_path_factory.mktemp("dataset"))


@pytest.fixture(scope="module")
def series_job(dataset):
    started = call_tool(dataset, "start_series", {**SERIES, "out": "runs/series_1"})
    assert started["ok"] is True
    return started, wait_for_job(dataset, started["job"])


@pytest.fixture(scope="module")
def kd_job(dataset):
    started = call_tool(dataset, "kd_fit", {**KD, "out": "runs/kd_1"})
    assert started["ok"] is True
    return started, wait_for_job(dataset, started["job"])


def _pid_is_gone(pid, timeout_s=15):
    deadline = time.monotonic() + timeout_s
    while time.monotonic() < deadline:
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            return True
        time.sleep(0.5)
    return False


def test_every_job_argument_vector_is_accepted_by_the_cli_parser():
    parser = build_parser()
    for parallel in (False, True):
        parser.parse_args(cli_call.start_series_args(
            "spectra", "peaks.txt", "out", mode="time", peak_source="reference",
            parallel=parallel))
    parser.parse_args(cli_call.kd_fit_args(
        "tidy.csv", 50, "absolute", "out", conc=[0, 25, 50], observables=["csp"]))
    parser.parse_args(cli_call.kd_fit_args(
        "tidy.csv", 50, "equivalents", "out", residues_file="kd_residues.txt",
        intensity_scales=[1, 1, 0.5]))


@needs_binding_data
def test_a_started_series_is_a_running_job_with_a_record_on_disk(dataset, series_job):
    started, _ = series_job
    assert started["status"] == "running"
    record = json.loads((dataset / "runs" / "series_1" / "mcp_job" / "job.json").read_text())
    assert record["kind"] == "series"


@needs_binding_data
def test_a_finished_series_job_carries_the_cli_summary(series_job):
    _, status = series_job
    assert status["status"] == "done"
    assert status["result"]["spectra_fitted"] == 7
    assert "params" not in status["result"]


@needs_binding_data
def test_a_finished_series_job_reports_what_exit_zero_hides(dataset, series_job):
    _, status = series_job
    with open(dataset / "runs" / "series_1" / "series_analysis_tidy.csv", newline="") as handle:
        rows = list(csv.DictReader(handle))
    quality = status["series_quality"]
    assert quality["points"] == len(rows) == 483
    assert quality["position_pinned"] == sum(r["position_pinned"] == "True" for r in rows) > 0
    assert quality["intensity_pinned"] == sum(r["intensity_pinned"] == "True" for r in rows) > 0
    assert quality["below_3_sigma"] == sum(float(r["detection_sigma"]) < 3 for r in rows)
    assert quality["residues_with_flagged_points"]


@needs_binding_data
def test_a_job_can_be_read_back_by_a_server_that_did_not_start_it(dataset, series_job):
    _, status = series_job
    again = call_tool(dataset, "job_status", {"job": "runs/series_1"})
    assert again["status"] == "done"
    assert again["series_quality"] == status["series_quality"]


@needs_binding_data
def test_start_series_refuses_an_output_folder_that_exists(dataset):
    summary = call_tool(dataset, "start_series", {**SERIES, "out": "series"})
    assert summary["ok"] is False
    assert summary["error"]["type"] == "PathRefused"


@needs_binding_data
def test_job_status_of_a_folder_that_is_not_a_job_is_an_error(dataset):
    summary = call_tool(dataset, "job_status", {"job": "series"})
    assert summary["ok"] is False
    assert summary["error"]["type"] == "NotAJob"


@needs_binding_data
def test_a_finished_kd_job_reports_the_shared_kd_and_its_pool(dataset, kd_job):
    _, status = kd_job
    assert status["status"] == "done"
    fit = json.loads(next((dataset / "runs" / "kd_1").glob("data/*_kd_fit_data.json")).read_text())
    quality = status["kd_quality"]
    assert quality["n_total"] == len(fit["fits"]) == status["result"]["n_total"]
    for observable in ("csp", "intensity"):
        shared = fit["global"][observable]
        assert quality[observable]["kd"] == shared["Kd"]
        assert quality[observable]["reliable"] == shared["reliable"]
        assert quality[observable]["n_in_shared_fit"] == shared["n_residues"]
        fitted = [f[observable] for f in fit["fits"] if f[observable].get("success")]
        assert quality[observable]["n_fitted_below_r2_0.8"] == sum(
            1 for f in fitted if f.get("r_squared") is None or f["r_squared"] < 0.8)


@needs_binding_data
def test_the_intensity_kd_is_labelled_as_an_apparent_decay_constant(kd_job):
    _, status = kd_job
    assert "apparent" in status["kd_quality"]["intensity"]["meaning"]
    assert "1:1" in status["kd_quality"]["csp"]["meaning"]


@needs_binding_data
def test_an_unfittable_residue_is_named_not_counted_as_fitted(dataset):
    started = call_tool(dataset, "kd_fit", {
        **KD, "out": "runs/kd_unfittable", "residues": ["A17", "A25", "D36"]})
    status = wait_for_job(dataset, started["job"])
    assert status["status"] == "done"
    assert status["result"]["n_fitted"] == 2
    assert status["result"]["n_total"] == status["kd_quality"]["n_total"] == 3
    assert "D36" in status["kd_quality"]["unfitted"]


@needs_binding_data
def test_scan_count_scales_reach_the_fit_and_are_recorded(dataset):
    scales = [1, 1, 1, 1, 0.5, 0.5, 0.5]
    started = call_tool(dataset, "kd_fit", {
        **KD, "out": "runs/kd_scaled", "intensity_scales": scales})
    status = wait_for_job(dataset, started["job"])
    assert status["status"] == "done"
    fit = json.loads(next((dataset / "runs" / "kd_scaled").glob("data/*_kd_fit_data.json")).read_text())
    assert fit["metadata"]["intensity_scales"] == scales


@needs_binding_data
def test_kd_fit_can_be_limited_to_named_residues(dataset):
    started = call_tool(dataset, "kd_fit", {
        **KD, "out": "runs/kd_subset", "residues": ["A17", "A25", "A30"]})
    status = wait_for_job(dataset, started["job"])
    assert status["status"] == "done"
    assert status["result"]["n_fitted"] == status["kd_quality"]["n_total"] == 3


@needs_binding_data
def test_a_failed_job_keeps_the_commands_reason(dataset):
    started = call_tool(dataset, "kd_fit", {
        **KD, "out": "runs/kd_unknown", "residues": ["NoSuchResidue"]})
    status = wait_for_job(dataset, started["job"])
    assert status["status"] == "failed"
    assert "NoSuchResidue" in status["result"]["error"]["message"]
    assert "kd_quality" not in status


@needs_binding_data
def test_kd_fit_refuses_residues_and_a_residue_file_together(dataset):
    summary = call_tool(dataset, "kd_fit", {
        **KD, "out": "runs/kd_both", "residues": ["A17"], "residues_file": "tidy.csv"})
    assert summary["ok"] is False
    assert not (dataset / "runs" / "kd_both").exists()


@needs_binding_data
def test_a_running_job_can_be_cancelled_and_its_process_dies(dataset):
    async def run():
        async with Client(build_server(DataRoot(dataset), wait_seconds=0)) as client:
            started = json.loads((await client.call_tool(
                "start_series", {**SERIES, "out": "runs/series_cancel"})).content[0].text)
            cancelled = json.loads((await client.call_tool(
                "cancel_job", {"job": started["job"]})).content[0].text)
            return started, cancelled
    started, cancelled = asyncio.run(run())
    pid = json.loads((dataset / "runs" / "series_cancel" / "mcp_job" / "job.json").read_text())["pid"]
    assert cancelled["status"] == "cancelled"
    assert call_tool(dataset, "job_status", {"job": started["job"]})["status"] == "cancelled"
    assert _pid_is_gone(pid)


@needs_binding_data
def test_a_job_can_only_be_cancelled_by_the_server_that_started_it(dataset):
    async def run():
        async with Client(build_server(DataRoot(dataset), wait_seconds=0)) as client:
            return json.loads((await client.call_tool(
                "start_series", {**SERIES, "out": "runs/series_foreign"})).content[0].text)
    started = asyncio.run(run())
    refused = call_tool(dataset, "cancel_job", {"job": started["job"]})
    assert refused["ok"] is False
    assert refused["error"]["type"] == "NotCancellable"
    pid = json.loads((dataset / "runs" / "series_foreign" / "mcp_job" / "job.json").read_text())["pid"]
    os.killpg(pid, 15)


@needs_binding_data
def test_the_server_refuses_to_start_more_jobs_than_it_allows(dataset):
    async def run():
        async with Client(build_server(DataRoot(dataset), max_running_jobs=1, wait_seconds=0)) as client:
            first = await client.call_tool("start_series", {**SERIES, "out": "runs/limit_1"})
            second = await client.call_tool("start_series", {**SERIES, "out": "runs/limit_2"})
            await client.call_tool("cancel_job", {"job": "runs/limit_1"})
            return first, second
    first, second = asyncio.run(run())
    assert first.is_error is False
    assert second.is_error is True
    assert json.loads(second.content[0].text)["error"]["type"] == "TooManyJobs"
    assert not (dataset / "runs" / "limit_2").exists()
