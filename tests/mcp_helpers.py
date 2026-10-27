# ABOUTME: Shared helper for the MCP server tests: one tool call through the SDK's in-process client.
# ABOUTME: Returns the tool's JSON envelope exactly as a client would receive it.

import asyncio
import json
import shutil
import time
from pathlib import Path

from mcp.client import Client

from lunaNMR.mcp.data_root import DataRoot
from lunaNMR.mcp.server import build_server


def call_result(root, tool, arguments, **server_options):
    """The protocol-level result: content, structured content and the error flag.

    Jobs return at once unless `wait_seconds` is given, so tests that poll stay deterministic.
    """
    server_options.setdefault("wait_seconds", 0)

    async def run():
        async with Client(build_server(DataRoot(root), **server_options)) as client:
            return await client.call_tool(tool, arguments)
    return asyncio.run(run())


def call_tool(root, tool, arguments, **server_options):
    return json.loads(call_result(root, tool, arguments, **server_options).content[0].text)


def list_tools(root):
    async def run():
        async with Client(build_server(DataRoot(root))) as client:
            return (await client.list_tools()).tools
    return asyncio.run(run())


BINDING = Path(__file__).resolve().parents[2] / "data_test" / "data_binding"


def stage_binding_dataset(root):
    """A private copy of the DNAJA1 titration under `root`: the tools write beside their inputs."""
    series = root / "series"
    series.mkdir()
    source = BINDING / "DNAJA1_HSPA8"
    for spectrum in source.glob("*.ft"):
        shutil.copy(spectrum, series / spectrum.name)
    shutil.copy(source / "DNAJA1.txt", series / "DNAJA1.txt")
    shutil.copy(BINDING / "tidy" / "DNAJA1_HSPA8_tidy.csv", root / "tidy.csv")
    return root


def wait_for_job(root, job, timeout_s=600):
    deadline = time.monotonic() + timeout_s
    while time.monotonic() < deadline:
        status = call_tool(root, "job_status", {"job": job})
        if status["status"] != "running":
            return status
        time.sleep(2)
    raise TimeoutError(f"job {job} still running after {timeout_s} s")


DATA_TEST = BINDING.parent
PREVIOUS_RUN = DATA_TEST / "_JA1_previous_run" / "analysis"
METHYL_MATRIX = (DATA_TEST / "methyl" / "sampleNBD" / "series_results_sampleNBD_20260818_124619"
                 / "peak_intensity_matrix.csv")
T1_PEAK_LIST = DATA_TEST / "data_JA1_relaxation" / "600_WT" / "T1" / "600_T1_assignment.txt"


def relaxation_data_present():
    return PREVIOUS_RUN.is_dir() and METHYL_MATRIX.is_file() and T1_PEAK_LIST.is_file()


def stage_relaxation_dataset(root):
    """A private copy of the real relaxation test inputs: matrices from an earlier series run,
    the hetNOE planes, a real methyl matrix and the T1 peak list, all under `root/in`."""
    inputs = root / "in"
    inputs.mkdir()
    for field in ("600", "700"):
        for experiment in ("T1", "T2"):
            shutil.copy(PREVIOUS_RUN / "series" / f"{field}_{experiment}" / "peak_intensity_matrix.csv",
                        inputs / f"{field}_{experiment}_matrix.csv")
        for plane in ("sat", "unsat"):
            shutil.copy(PREVIOUS_RUN / "results" / f"{field}_noe_{plane}.csv",
                        inputs / f"{field}_noe_{plane}.csv")
    shutil.copy(METHYL_MATRIX, inputs / "methyl_matrix.csv")
    shutil.copy(T1_PEAK_LIST, inputs / "peaks.txt")
    shutil.copy(BINDING / "tidy" / "DNAJA1_HSPA8_tidy.csv", inputs / "tidy.csv")
    return root
