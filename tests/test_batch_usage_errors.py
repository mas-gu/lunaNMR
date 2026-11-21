# ABOUTME: A command line `lunaNMR batch` refuses exits 2 as a UsageError naming the reason,
# ABOUTME: and a folder or config file that does not exist exits 1 with the path in the envelope.

"""
batch validated its flag combinations after argparse and returned 1 for all of them, so
under --format json `--quiet --verbose` was a `BatchFailed` whose message was "batch
exited 1; its report is on stderr": the exit code of a failed run, and no reason.
"""

import contextlib
import io
import json

import pytest


def _batch(*argv):
    from lunaNMR.cli import main
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
        code = main(["batch", *argv, "--format", "json"])
    return code, json.loads(out.getvalue())


@pytest.fixture
def config(tmp_path):
    path = tmp_path / "config.json"
    path.write_text(json.dumps({"processing": {"sn_thresholds": {"15N1H": 3.0}}}))
    return path


@pytest.mark.parametrize("flags, named", [
    (["--quiet", "--verbose"], "--quiet and --verbose"),
    (["--skip-errors", "--strict"], "--skip-errors and --strict"),
    (["--preset", "15N1H", "--config", "CONFIG"], "--preset and --config"),
])
def test_a_flag_combination_batch_refuses_is_a_usage_error(tmp_path, config, flags, named):
    flags = [str(config) if flag == "CONFIG" else flag for flag in flags]
    code, envelope = _batch(str(tmp_path), *flags)
    assert code == 2
    assert envelope["ok"] is False
    assert envelope["error"]["type"] == "UsageError"
    assert named in envelope["error"]["message"]


def test_a_command_line_without_a_folder_is_a_usage_error():
    code, envelope = _batch()
    assert code == 2
    assert envelope["error"]["type"] == "UsageError"
    assert "folder" in envelope["error"]["message"]


def test_a_folder_that_does_not_exist_is_bad_input_and_is_named(tmp_path):
    code, envelope = _batch(str(tmp_path / "nowhere"))
    assert code == 1
    assert envelope["error"]["type"] == "BatchFailed"
    assert "nowhere" in envelope["error"]["message"]


def test_a_config_file_that_does_not_exist_is_bad_input_and_is_named(tmp_path):
    code, envelope = _batch(str(tmp_path), "--config", str(tmp_path / "absent.json"))
    assert code == 1
    assert envelope["error"]["type"] == "BatchFailed"
    assert "absent.json" in envelope["error"]["message"]
