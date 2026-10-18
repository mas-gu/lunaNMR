# ABOUTME: `lunaNMR batch` keeps the CLI's contract: a dry run that cannot list is a failure,
# ABOUTME: and --format json gives the one envelope on stdout. Real DNAJA1 spectra.

"""
batch delegates to its own CLI, so the rest of the contract never reached it: a dry run
that could not list the folder printed the error and exited 0, and there was no
--format json - a caller parsing stdout got batch's progress text.
"""

import contextlib
import io
import json
import os
from pathlib import Path

import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


def _batch(*argv):
    from lunaNMR.cli import main
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
        code = main(["batch", *argv])
    return code, out.getvalue()


@pytest.fixture
def folder(tmp_path):
    spectra = tmp_path / "spectra"
    spectra.mkdir()
    (spectra / "2D_DNAJA1_0o0.ft").symlink_to(DNAJA1 / "2D_DNAJA1_0o0.ft")
    return spectra


def test_a_dry_run_that_finds_nothing_to_process_fails(folder):
    """An unreadable folder lists as empty: a run would process nothing."""
    os.chmod(folder, 0)
    try:
        code, _ = _batch(str(folder), "--dry-run")
    finally:
        os.chmod(folder, 0o755)
    assert code == 1


def test_a_dry_run_reports_in_the_envelope(folder):
    code, out = _batch(str(folder), "--dry-run", "--format", "json")
    envelope = json.loads(out)
    assert code == 0
    assert envelope["ok"] is True
    assert envelope["command"] == "batch"
    assert envelope["dry_run"] is True
    assert envelope["spectra_found"] == 1


def test_a_failed_run_reports_why_in_the_envelope(tmp_path):
    """A folder with no spectra: the run's failure reaches the envelope with its reason,
    and its non-zero exit makes it ok:false. (A run with spectra detects and fits every
    peak it finds, which takes minutes per spectrum.)"""
    (tmp_path / "empty").mkdir()
    code, out = _batch(str(tmp_path / "empty"), "--format", "json",
                       "--output", str(tmp_path / "out"))
    envelope = json.loads(out)
    assert code == 1
    assert envelope["ok"] is False
    assert envelope["error"]["type"] == "BatchFailed"
    assert "No spectrum files found" in envelope["error"]["message"]


def test_a_rejected_command_line_returns_instead_of_raising(tmp_path):
    code, out = _batch("--no-such-flag", "--format", "json")
    assert code == 2
    assert json.loads(out)["error"]["type"] == "UsageError"
