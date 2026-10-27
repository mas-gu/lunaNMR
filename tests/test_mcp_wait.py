# ABOUTME: A job tool waits a short while for the run to finish, then hands back a job to poll.
# ABOUTME: Quick analyses therefore answer in one call; slow ones never block the client.

import pytest

pytest.importorskip("mcp")

from mcp_helpers import BINDING, call_tool, stage_binding_dataset

needs_binding_data = pytest.mark.skipif(
    not BINDING.is_dir(), reason="data_test/data_binding is absent")


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    return stage_binding_dataset(tmp_path_factory.mktemp("dataset"))


@needs_binding_data
def test_a_quick_job_comes_back_finished_in_one_call(dataset):
    summary = call_tool(dataset, "kd_fit", {
        "input": "tidy.csv", "p0": 50, "conc_units": "equivalents", "out": "wait_quick"},
        wait_seconds=120)
    assert summary["status"] == "done"
    assert summary["kd_quality"]["csp"]["reliable"] in (True, False)


@needs_binding_data
def test_a_job_that_outlasts_the_wait_comes_back_running(dataset):
    summary = call_tool(dataset, "start_series", {
        "spectra": "series", "peaks": "series/DNAJA1.txt", "out": "wait_slow",
        "mode": "titration", "peak_source": "cascade"}, wait_seconds=1)
    assert summary["status"] == "running"
    call_tool(dataset, "cancel_job", {"job": summary["job"]})
