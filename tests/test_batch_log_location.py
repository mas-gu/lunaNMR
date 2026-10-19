# ABOUTME: `lunaNMR batch` writes its log into --output, never into the folder it was run from,
# ABOUTME: and a second processor in one process does not stack handlers on the first's logger.

import logging
import os
from pathlib import Path

import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


@pytest.fixture
def elsewhere(tmp_path, monkeypatch):
    """Run from a folder of its own, so a log written to the working directory shows."""
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    monkeypatch.chdir(cwd)
    return cwd


def _one_spectrum(tmp_path):
    spectra = tmp_path / "spectra"
    spectra.mkdir()
    (spectra / "2D_DNAJA1_0o0.ft").symlink_to(DNAJA1 / "2D_DNAJA1_0o0.ft")
    return spectra


def test_the_log_goes_to_output(tmp_path, elsewhere):
    from lunaNMR.cli import main
    out = tmp_path / "out"
    main(["batch", "--preset", "15N1H", str(_one_spectrum(tmp_path)), "--dry-run",
          "--output", str(out)])
    assert list(elsewhere.glob("*.log")) == []
    assert len(list(out.glob("batch_processing_*.log"))) == 1


def test_two_processors_leave_one_file_handler(tmp_path, elsewhere):
    from lunaNMR.batch_processing.batch_processor import BatchProcessor
    BatchProcessor(output_dir=tmp_path / "a")
    BatchProcessor(output_dir=tmp_path / "b")
    logger = logging.getLogger("lunaNMR.batch_processing.batch_processor")
    files = [h for h in logger.handlers if isinstance(h, logging.FileHandler)]
    assert len(files) == 1
    assert Path(files[0].baseFilename).parent == tmp_path / "b"
    assert len([h for h in logger.handlers if type(h) is logging.StreamHandler]) == 1
