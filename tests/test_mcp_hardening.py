# ABOUTME: Hardening of the MCP server: nothing a client names may read or reveal files outside the data root.
# ABOUTME: Also bounds what dataset-controlled text can put in front of the model, and survives hostile job records.

import asyncio
import csv
import json
import subprocess
import sys
import time

import pytest

pytest.importorskip("mcp")

from lunaNMR.mcp import cli_call, jobs, server, summaries
from lunaNMR.mcp.data_root import DataRoot, PathRefused
from mcp_helpers import call_result, call_tool

TITRATION_FIELDS = ["spectrum_name", "assignment", "peak_number", "ppm_x", "ppm_y", "height",
                    "volume", "snr", "quality", "r_squared", "detection_sigma",
                    "position_pinned", "intensity_pinned"]


@pytest.fixture
def root(tmp_path):
    folder = tmp_path / "dataset"
    folder.mkdir()
    return folder


@pytest.fixture
def outside(tmp_path):
    folder = tmp_path / "outside"
    folder.mkdir()
    (folder / "secret.txt").write_text("SECRET_MARKER\n")
    return folder


def _tidy_csv(path, assignment="A17"):
    rows = [dict(zip(TITRATION_FIELDS, [point, assignment, 1, 8.0, 120.0, 100.0 / (1 + point),
                                        10.0, 50.0, "Excellent", 0.99, 50.0, "False", "False"]))
            for point in range(5)]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TITRATION_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


KD = {"input": "tidy.csv", "p0": 50, "conc_units": "equivalents"}


class TestResidueNamesAreNeverFilePaths:
    @pytest.fixture(autouse=True)
    def _input(self, root):
        _tidy_csv(root / "tidy.csv")

    def test_an_absolute_path_is_refused_and_nothing_is_started(self, root, outside):
        summary = call_tool(root, "kd_fit", {
            **KD, "out": "run", "residues": [str(outside / "secret.txt")]})
        assert summary["ok"] is False
        assert summary["error"]["type"] == "InvalidResidueName"
        assert not (root / "run").exists()

    @pytest.mark.parametrize("name", ["a/b", "..\\x", "line\nbreak", "#comment", "", "x" * 65, "-flag"])
    def test_names_that_are_not_residue_names_are_refused(self, root, name):
        summary = call_tool(root, "kd_fit", {**KD, "out": "run", "residues": [name]})
        assert summary["error"]["type"] == "InvalidResidueName"

    def test_a_name_equal_to_a_file_in_the_working_directory_is_still_only_a_name(
            self, root, tmp_path, monkeypatch):
        decoy = tmp_path / "workdir"
        decoy.mkdir()
        (decoy / "decoy.txt").write_text("DECOY_CONTENT\n")
        monkeypatch.chdir(decoy)
        started = call_tool(root, "kd_fit", {**KD, "out": "run", "residues": ["decoy.txt"]})
        assert started["ok"] is True
        deadline = 120
        status = None
        for _ in range(deadline):
            status = call_tool(root, "job_status", {"job": "run"})
            if status["status"] != "running":
                break
            time.sleep(1)
        assert status["status"] == "failed"
        assert "decoy.txt" in status["result"]["error"]["message"]
        assert "DECOY_CONTENT" not in json.dumps(status)


class TestSymlinksCannotLeadOutOfTheRoot:
    def test_an_escaping_link_blocks_a_scan_of_the_folder_holding_it(self, root, outside):
        (root / "linked").symlink_to(outside)
        summary = call_tool(root, "diagnose", {"path": "."})
        assert summary["ok"] is False
        assert summary["error"]["type"] == "PathRefused"
        assert "linked" in summary["error"]["message"]

    def test_a_link_to_a_file_outside_blocks_the_folder_too(self, root, outside):
        (root / "series").mkdir()
        (root / "series" / "peaks.txt").symlink_to(outside / "secret.txt")
        summary = call_tool(root, "plan_series", {
            "spectra": "series", "peaks": "series/peaks.txt", "out": "plan"})
        assert summary["error"]["type"] == "PathRefused"

    def test_a_link_beside_the_kd_input_blocks_it(self, root, outside):
        _tidy_csv(root / "tidy.csv")
        (root / "kd_params.json").symlink_to(outside / "secret.txt")
        summary = call_tool(root, "kd_survey", {**KD, "out": "survey"})
        assert summary["error"]["type"] == "PathRefused"

    def test_a_link_that_stays_inside_the_root_is_fine(self, root):
        (root / "real").mkdir()
        (root / "alias").symlink_to(root / "real")
        assert DataRoot(root).refuse_escaping_links(root) is None


class TestWhatReachesTheModelIsBounded:
    def test_a_long_log_line_is_cut(self, tmp_path):
        out = tmp_path / "job"
        (out / "mcp_job").mkdir(parents=True)
        (out / "mcp_job" / "job.json").write_text(json.dumps(
            {"kind": "kd", "args": [], "pid": 2 ** 22 + 7, "started": 1.0}))
        (out / "mcp_job" / "result.json").write_text(json.dumps({"exit_code": 1, "finished": 2.0}))
        (out / "mcp_job" / "stderr.log").write_text(("x" * 5000 + "\n") * 3)
        status = jobs.JobManager().status(out)
        assert status["status"] == "failed"
        assert max(len(line) for line in status["stderr_tail"]) <= jobs.LINE_CHARS

    def test_long_residue_names_are_cut(self, tmp_path):
        _tidy_csv(tmp_path / "series_analysis_tidy.csv", assignment="R" * 500)
        rows = list(csv.DictReader(open(tmp_path / "series_analysis_tidy.csv")))
        for row in rows:
            row["detection_sigma"] = "1.0"
        with open(tmp_path / "series_analysis_tidy.csv", "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=TITRATION_FIELDS)
            writer.writeheader()
            writer.writerows(rows)
        names = summaries.series_quality(tmp_path)["residues_with_flagged_points"]
        assert names and max(len(name) for name in names) <= summaries.NAME_CHARS

    def test_the_text_of_a_run_that_printed_no_json_is_cut(self):
        envelope = asyncio.run(cli_call.run_cli(["definitely-not-a-subcommand-" + "x" * 3000]))
        assert envelope["ok"] is False
        assert len(envelope["error"]["message"]) <= cli_call.MESSAGE_CHARS + 200

    def test_an_oversized_result_becomes_an_error_not_a_flood(self):
        result = server._result({"ok": True, "blob": "y" * (server.MAX_RESULT_CHARS + 1)})
        assert result.is_error is True
        assert json.loads(result.content[0].text)["error"]["type"] == "ResultTooLarge"


class TestHostileJobRecords:
    def test_a_job_record_missing_its_fields_is_not_a_job(self, root):
        (root / "a" / "mcp_job").mkdir(parents=True)
        (root / "a" / "mcp_job" / "job.json").write_text("{}")
        summary = call_tool(root, "job_status", {"job": "a"})
        assert summary["ok"] is False
        assert summary["error"]["type"] == "NotAJob"

    def test_a_result_record_without_an_exit_code_is_treated_as_no_result(self, root):
        finished = subprocess.Popen([sys.executable, "-c", "pass"])
        finished.wait()
        (root / "b" / "mcp_job").mkdir(parents=True)
        (root / "b" / "mcp_job" / "job.json").write_text(json.dumps(
            {"kind": "kd", "args": [], "pid": finished.pid, "started": 1.0}))
        (root / "b" / "mcp_job" / "result.json").write_text(json.dumps({"finished": 2.0}))
        summary = call_tool(root, "job_status", {"job": "b"})
        assert summary["ok"] is True
        assert summary["status"] == "failed"

    def test_a_second_start_on_the_same_output_folder_is_refused_not_a_crash(self, root):
        out = root / "run"
        (out / "mcp_job").mkdir(parents=True)
        with pytest.raises(PathRefused, match="already exists"):
            jobs.JobManager().start("kd", ["kd"], out)


class TestPackaging:
    def test_the_extra_pins_the_sdk_major_version_and_the_python_it_needs(self):
        from pathlib import Path
        import re
        text = (Path(__file__).resolve().parents[1] / "pyproject.toml").read_text()
        requirement = re.search(r'(?m)^mcp = \["([^"]+)"\]', text).group(1)
        assert "<3" in requirement
        assert "python_version" in requirement and "3.10" in requirement
