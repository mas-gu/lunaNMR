# ABOUTME: t1t2 and methyl-t2 must label a matrix with no sidecar in ms, the unit the delay
# ABOUTME: parser returns. Labelled seconds, a 711 ms T1 read as 711 s in the JSON summary.

import json
import os
import shutil
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

RUN = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
       / "run_20260825_1405")

pytestmark = pytest.mark.skipif(not RUN.exists(), reason="data_test not present")


def _matrix_without_sidecar(tmp_path, series):
    """A real series matrix as an older run or a hand-built table leaves it: no
    series_metadata.json beside it, so the units come from the CLI's default."""
    target = tmp_path / series
    target.mkdir()
    shutil.copy(RUN / "series" / series / "peak_intensity_matrix.csv", target)
    return target / "peak_intensity_matrix.csv"


def _summary(capsys, argv):
    from lunaNMR.cli import main
    assert main(argv + ["--format", "json"]) == 0
    return json.loads(capsys.readouterr().out)


def test_t1t2_reports_ms_without_a_flag(tmp_path, capsys):
    csv = _matrix_without_sidecar(tmp_path, "600_T1")
    out = tmp_path / "t1"
    summary = _summary(capsys, ["dynamixs", "t1t2", "--exp", "T1", "--input", str(csv),
                                "--out", str(out)])
    assert summary["series_metadata"] is None
    assert summary["time_units"] == "ms"
    # A 600 MHz backbone T1 is a few hundred ms to about a second.
    assert 200 < summary["mean_t2"] < 2000
    fit = json.loads(Path(summary["json_file"]).read_text())
    assert fit["metadata"]["time_units"] == "ms"


def test_methyl_t2_reports_ms_without_a_flag(tmp_path, capsys):
    csv = _matrix_without_sidecar(tmp_path, "600_T2")
    out = tmp_path / "methyl"
    summary = _summary(capsys, ["dynamixs", "methyl-t2", "--input", str(csv),
                                "--out", str(out)])
    assert summary["series_metadata"] is None
    assert summary["time_units"] == "ms"
    fit = json.loads(Path(summary["json_file"]).read_text())
    assert fit["metadata"]["time_units"] == "ms"
