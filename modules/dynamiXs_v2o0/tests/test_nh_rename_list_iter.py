# ABOUTME: End-to-end test for the NH_rename_list_iter peak-list transform script,
# ABOUTME: running it as a subprocess against real files in a temporary data/ folder.
"""The script is top-level code with no entry point, so it is exercised by running it.

It reads `data/` relative to the working directory, so the test supplies a real
tab-separated peak list there and checks the transformed output rather than
mocking pandas or the filesystem.
"""

import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = (Path(__file__).resolve().parent.parent
          / "dynamiXs_format" / "NH_rename_list_iter.py")

INPUT_ROWS = "\n".join([
    "Assignment\tPosition_X\tPosition_Y\tHeight",
    "(<NA:A.12.SER.H>, <NA:A.12.SER.N>)\t8.51\t121.4\t20000",
    "(<NA:A.10.VAL.H>, <NA:A.10.VAL.N>)\t8.33\t128.3\t10000",
]) + "\n"


@pytest.fixture
def data_dir(tmp_path):
    d = tmp_path / "data"
    d.mkdir()
    (d / "T1_sample_0o0.txt").write_text(INPUT_ROWS, encoding="utf-8")
    return tmp_path


def _run(cwd):
    return subprocess.run([sys.executable, str(SCRIPT)], cwd=cwd,
                          capture_output=True, text=True)


def test_script_runs_without_error(data_dir):
    result = _run(data_dir)
    assert result.returncode == 0, result.stderr
    assert "NameError" not in result.stderr


def test_assignments_are_transformed_and_sorted_by_residue(data_dir):
    _run(data_dir)
    out = (data_dir / "data" / "T1_sample_0o0_transformed.txt").read_text(encoding="utf-8")
    lines = out.strip().splitlines()
    assert lines[0] == "Assignment, Position_X, Position_Y, Height"
    # residue 10 sorts ahead of 12 even though it is second in the input
    assert lines[1].startswith("10.VAL, ")
    assert lines[2].startswith("12.SER, ")


def test_missing_data_directory_is_reported(tmp_path):
    result = _run(tmp_path)
    assert "does not exist" in (result.stdout + result.stderr)
