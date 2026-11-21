# ABOUTME: --quiet, --verbose, --log-level and a config's logging.level set what `lunaNMR batch`
# ABOUTME: prints of its processing log; the log file keeps every line. Real runs on a DNAJA1 spectrum.

"""
The processing log has its own console handler, fixed at INFO: --quiet still printed every
"Processing: ..." line, --verbose printed no DEBUG line, and a config file's logging.level
was replaced by the --log-level default before nothing read it.
"""

import contextlib
import io
import json
from pathlib import Path

import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
SPECTRUM = "2D_DNAJA1_0o0.ft"
INFO_LINE = f"INFO - Processing: {SPECTRUM}"
WARNING_LINE = f"WARNING -   No peaks detected in {SPECTRUM}"
DEBUG_LINE = "DEBUG - Parameters for 15N1H"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


def _run(tmp_path, *flags, logging_level=None):
    """One `lunaNMR batch` run with a threshold no peak reaches, so it ends at once.
    Returns (what it printed, its log file)."""
    from lunaNMR.cli import main
    spectra = tmp_path / "spectra"
    spectra.mkdir()
    (spectra / SPECTRUM).symlink_to(DNAJA1 / SPECTRUM)
    content = {"processing": {"sn_thresholds": {"15N1H": 1000.0}}}
    if logging_level:
        content["logging"] = {"level": logging_level}
    config = tmp_path / "config.json"
    config.write_text(json.dumps(content))
    out = tmp_path / "out"
    console = io.StringIO()
    with contextlib.redirect_stdout(console), contextlib.redirect_stderr(console):
        code = main(["batch", str(spectra), "--config", str(config), "--output", str(out), *flags])
    assert code == 1
    return console.getvalue(), next(out.glob("batch_processing_*.log")).read_text()


def test_without_a_flag_the_processing_log_is_printed_from_info_up(tmp_path):
    console, _ = _run(tmp_path)
    assert INFO_LINE in console
    assert WARNING_LINE in console
    assert DEBUG_LINE not in console


def test_quiet_prints_no_line_below_an_error_and_the_file_keeps_them(tmp_path):
    console, log = _run(tmp_path, "--quiet")
    assert "Processing:" not in console
    assert "No peaks detected in" not in console
    assert INFO_LINE in log
    assert WARNING_LINE in log


def test_log_level_warning_prints_the_warnings_only(tmp_path):
    console, log = _run(tmp_path, "--log-level", "WARNING")
    assert "Processing:" not in console
    assert WARNING_LINE in console
    assert INFO_LINE in log


def test_verbose_prints_and_writes_the_debug_lines(tmp_path):
    console, log = _run(tmp_path, "--verbose")
    assert DEBUG_LINE in console
    assert DEBUG_LINE in log


def test_the_level_of_a_config_file_is_the_one_printed_from(tmp_path):
    console, log = _run(tmp_path, logging_level="WARNING")
    assert "Processing:" not in console
    assert WARNING_LINE in console
    assert INFO_LINE in log


def test_log_level_on_the_command_line_replaces_the_config_files(tmp_path):
    console, _ = _run(tmp_path, "--log-level", "INFO", logging_level="WARNING")
    assert INFO_LINE in console
