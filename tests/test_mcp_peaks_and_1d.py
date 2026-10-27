# ABOUTME: The peak-list shift and 1D integration tools of the MCP server, run on real spectra.
# ABOUTME: The 1D series is written as NMRPipe files by the same helper the CLI's own tests use.

import json

import pytest

pytest.importorskip("mcp")
pytest.importorskip("nmrglue")

from lunaNMR.cli import build_parser
from lunaNMR.mcp import cli_call
from mcp_helpers import BINDING, call_tool, stage_binding_dataset, wait_for_job
from test_cli_integrate_1d import _write_series

needs_binding_data = pytest.mark.skipif(
    not BINDING.is_dir(), reason="data_test/data_binding is absent")


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    return stage_binding_dataset(tmp_path_factory.mktemp("dataset"))


@pytest.fixture(scope="module")
def series_1d(tmp_path_factory):
    folder = tmp_path_factory.mktemp("oned")
    spectra = folder / "spectra"
    spectra.mkdir()
    _write_series(spectra, n=4)
    return folder


def test_argument_vectors_are_accepted_by_the_cli_parser():
    parser = build_parser()
    parser.parse_args(cli_call.peaks_shift_args("p.txt", "o.txt", dx=0.01, dy=-0.1))
    parser.parse_args(cli_call.peaks_shift_args("p.txt", "o.txt", auto=True, spectrum="s.ft",
                                                quick=True))
    parser.parse_args(cli_call.integrate_1d_args("s", "o", peaks=[5.0, 5.3], prefix="x",
                                                 min_snr=5.0, value="height"))
    parser.parse_args(cli_call.integrate_1d_args("s", "o", detect=True))


def _positions(path):
    rows = [[cell.strip() for cell in line.replace("\t", ",").split(",")]
            for line in path.read_text().splitlines()[1:] if line.strip()]
    return {row[0]: (float(row[1]), float(row[2])) for row in rows}


@needs_binding_data
def test_a_manual_shift_writes_the_moved_list_to_a_new_file(dataset):
    summary = call_tool(dataset, "peaks_shift", {
        "peaks": "series/DNAJA1.txt", "out": "shifted.txt", "dx": 0.01, "dy": -0.1})
    assert summary["ok"] is True
    before = _positions(dataset / "series" / "DNAJA1.txt")
    after = _positions(dataset / "shifted.txt")
    name = next(iter(before))
    assert after[name][0] == pytest.approx(before[name][0] + 0.01, abs=1e-4)
    assert after[name][1] == pytest.approx(before[name][1] - 0.1, abs=1e-3)


@needs_binding_data
def test_a_shift_will_not_overwrite_a_file(dataset):
    summary = call_tool(dataset, "peaks_shift", {
        "peaks": "series/DNAJA1.txt", "out": "series/DNAJA1.txt", "dx": 0.01})
    assert summary["ok"] is False
    assert summary["error"]["type"] == "PathRefused"


@needs_binding_data
def test_a_shift_will_not_write_outside_the_root(dataset, tmp_path):
    summary = call_tool(dataset, "peaks_shift", {
        "peaks": "series/DNAJA1.txt", "out": str(tmp_path / "out.txt"), "dx": 0.01})
    assert summary["error"]["type"] == "PathRefused"


@needs_binding_data
def test_an_automatic_shift_measures_a_registered_list_as_unmoved(dataset):
    summary = call_tool(dataset, "peaks_shift", {
        "peaks": "series/DNAJA1.txt", "out": "auto_shifted.txt", "auto": True,
        "spectrum": "series/2D_DNAJA1_0o0.ft", "quick": True})
    assert summary["ok"] is True
    assert (dataset / "auto_shifted.txt").is_file()


def test_the_1d_series_is_measured_and_its_files_are_listed(series_1d):
    started = call_tool(series_1d, "integrate_1d", {
        "spectra": "spectra", "out": "measured", "peaks": [5.0, 5.3], "prefix": "run1"})
    status = wait_for_job(series_1d, started["job"], 120)
    assert status["status"] == "done"
    assert "run1_series.csv" in status["files"]
    assert (series_1d / "measured" / "run1_series.csv").is_file()


def test_the_1d_prefix_cannot_name_a_path(series_1d):
    summary = call_tool(series_1d, "integrate_1d", {
        "spectra": "spectra", "out": "measured_bad", "detect": True, "prefix": "../escape"})
    assert summary["error"]["type"] == "InvalidName"
    assert not (series_1d / "measured_bad").exists()
