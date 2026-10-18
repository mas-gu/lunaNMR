# ABOUTME: The settings file `kd` and `kd --survey` write for later runs records only what both
# ABOUTME: front ends restore. Real DNAJA1 tidy CSV through the real CLI.

"""
The <prefix>_kd_params.json a run leaves for the next one also recorded protein_conc,
alpha, intensity_value and observables, which `kd` never restores (its flags decide them)
and the page restored over whatever the user had set. The survey recorded observables as
both whatever --observable said. The file now holds the concentrations with their units,
the intensity scales and the survey thresholds, plus what names the input it describes.
"""

import contextlib
import io
import json
import shutil
from pathlib import Path

import pytest

TIDY = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "tidy"
EQUIVALENTS = "0,0.125,0.25,0.5,1,2,3"
STORED = {"input", "points", "concentrations", "conc_units", "intensity_scales",
          "csp_sigma_multiple", "kd_outlier_z", "noise_quantile", "dd_runaway_ratio",
          "ref_max_ratio"}

pytestmark = pytest.mark.skipif(not TIDY.exists(), reason="data_test not present")


@pytest.fixture
def csv(tmp_path):
    path = tmp_path / "DNAJA1_HSPA8_tidy.csv"
    shutil.copy(TIDY / path.name, path)
    return path


def _kd(*args):
    from lunaNMR.cli import main
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
        code = main(["kd", "--format", "json", "--p0", "30", "--alpha", "0.2", *args])
    assert code == 0, out.getvalue()
    summary = json.loads(out.getvalue())
    return summary.get("result", summary)


def test_the_survey_records_only_what_is_restored(csv, tmp_path):
    _kd("--survey", "--observable", "intensity", "--input", str(csv),
        "--out", str(tmp_path / "s"), "--prefix", "DNAJA1",
        "--conc", EQUIVALENTS, "--conc-units", "equivalents")
    written = json.loads((csv.parent / "DNAJA1_kd_params.json").read_text())
    assert set(written) == STORED
    assert written["input"] == csv.name
    assert written["conc_units"] == "equivalents"
    assert "observables" not in written


def test_the_fit_records_only_what_is_restored(csv, tmp_path):
    result = _kd("--observable", "intensity", "--input", str(csv),
                 "--out", str(tmp_path / "f"), "--conc", EQUIVALENTS,
                 "--conc-units", "equivalents")
    written = json.loads(Path(result["params_file"]).read_text())
    assert set(written) == STORED
    assert written["concentrations"] == [0, 0.125, 0.25, 0.5, 1, 2, 3]


def test_the_fit_json_keeps_its_full_record(csv, tmp_path):
    result = _kd("--observable", "intensity", "--input", str(csv),
                 "--out", str(tmp_path / "f"), "--conc", EQUIVALENTS,
                 "--conc-units", "equivalents")
    meta = json.loads(Path(result["json_file"]).read_text())["metadata"]
    assert meta["protein_conc"] == 30.0
    assert meta["alpha"] == 0.2
    assert meta["observables"] == ["intensity"]
    assert meta["intensity_value"] == "height"


def test_units_alone_are_recorded_and_nothing_else_is_invented(csv, tmp_path):
    _kd("--survey", "--input", str(csv), "--out", str(tmp_path / "s"),
        "--prefix", "DNAJA1", "--conc-units", "equivalents")
    written = json.loads((csv.parent / "DNAJA1_kd_params.json").read_text())
    assert written["conc_units"] == "equivalents"
    assert written["concentrations"] is None
