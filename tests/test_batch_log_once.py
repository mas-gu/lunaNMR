# ABOUTME: A `lunaNMR batch` log line reaches the console once and the batch log file once.
# ABOUTME: Real run of the CLI in a fresh process on the DNAJA1 spectrum, as a user starts it.

"""
Each batch line was printed three times: by the batch logger's own console handler, and
again by the two handlers on the root logger it propagated to - the console handler the
batch CLI adds, and the one an import-time logging.basicConfig adds to an empty root. The
other batch modules' lines reached the two root handlers: twice.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
DNAJA1 = REPO.parent / "data_test" / "data_binding" / "DNAJA1_HSPA8"
SPECTRUM = "2D_DNAJA1_0o0.ft"
LINE = f"Processing: {SPECTRUM}"
# Logged by the batch package's configuration manager, not by the batch processor
CONFIG_LINE = "Loaded configuration from"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def batch_run(tmp_path_factory):
    """`python -m lunaNMR batch` with a threshold no peak reaches: a real run, which fails
    fast. Returns (the process, --output)."""
    root = tmp_path_factory.mktemp("batch_log")
    spectra = root / "spectra"
    spectra.mkdir()
    (spectra / SPECTRUM).symlink_to(DNAJA1 / SPECTRUM)
    config = root / "config.json"
    config.write_text(json.dumps({"processing": {"sn_thresholds": {"15N1H": 1000.0}},
                                  "optimization": {"enable_auto_optimization": False}}))
    out = root / "out"
    env = dict(os.environ, PYTHONPATH=str(REPO), QT_QPA_PLATFORM="offscreen")
    process = subprocess.run(
        [sys.executable, "-m", "lunaNMR", "batch", "--config", str(config), str(spectra),
         "--output", str(out)],
        cwd=root, env=env, capture_output=True, text=True, timeout=600)
    return process, out


def test_the_run_reached_the_spectrum(batch_run):
    process, _ = batch_run
    assert process.returncode == 1
    assert "No peaks detected" in process.stdout + process.stderr


def test_a_line_is_printed_once_on_the_console(batch_run):
    process, _ = batch_run
    console = process.stdout + process.stderr
    assert console.count(LINE) == 1


def test_a_line_of_another_batch_module_is_printed_once_on_the_console(batch_run):
    process, _ = batch_run
    console = process.stdout + process.stderr
    assert console.count(CONFIG_LINE) == 1


def test_a_line_is_written_once_in_the_log_file(batch_run):
    _, out = batch_run
    logs = list(out.glob("batch_processing_*.log"))
    assert len(logs) == 1
    assert logs[0].read_text().count(LINE) == 1
