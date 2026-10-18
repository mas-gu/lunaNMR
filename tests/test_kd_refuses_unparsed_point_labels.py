# ABOUTME: `lunaNMR kd` refuses to take concentrations from point labels `series` could not parse.
# ABOUTME: Real DNAJA1 spectra linked under index names, run through `series`, then `kd`.

"""
`series` names a spectrum it cannot read a titration point from after its file stem, and
records `value: null` for it in series_metadata.json. The Kd loader read the trailing
digits of that stem as the point, so `sample_2` became point 2 - and with
`--conc-units equivalents` alone, 2 eq = 100 uM for a spectrum acquired at 3 eq. That is
the mistake the no-concentrations refusal exists to prevent, passed through a name.
"""

import contextlib
import io
import json
import shutil
from pathlib import Path

import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
TIDY = DNAJA1.parent / "tidy"
# Spectrum acquired at each point -> the index name `series` cannot read a point from.
INDEX_NAMES = {"0o0": "sample_0.ft", "0o5": "sample_1.ft", "3o0": "sample_2.ft"}

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


def _main(*argv):
    from lunaNMR.cli import main
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        code = main([*argv, "--format", "json"])
    text = out.getvalue()
    return code, (json.loads(text) if text.strip() else {})


@pytest.fixture(scope="module")
def index_named_series(tmp_path_factory):
    root = tmp_path_factory.mktemp("index_named")
    (root / "spectra").mkdir()
    for point, name in INDEX_NAMES.items():
        (root / "spectra" / name).symlink_to(next(DNAJA1.glob(f"2D_*_{point}.ft")))
    code, _ = _main("series", "--spectra", str(root / "spectra"),
                    "--peaks", str(DNAJA1 / "DNAJA1.txt"), "--mode", "titration",
                    "--out", str(root / "out"))
    assert code == 0
    meta = json.loads((root / "out" / "series_metadata.json").read_text())
    assert meta["n_value_unparsed"] == 3, "the premise: series could not parse these names"
    return root / "out"


def _kd(csv, out, *extra):
    return _main("kd", "--input", str(csv), "--out", str(out), "--p0", "50",
                 "--observable", "intensity", *extra)


@pytest.mark.parametrize("csv", ["series_analysis_tidy.csv",
                                 "comprehensive_peak_tracking.csv",
                                 "peak_intensity_matrix.csv"])
def test_labels_series_marked_unparsed_are_not_concentrations(index_named_series, tmp_path, csv):
    code, summary = _kd(index_named_series / csv, tmp_path / "kd",
                        "--conc-units", "equivalents")
    assert code == 1
    message = summary["error"]["message"]
    for label in ("sample_0", "sample_1", "sample_2"):
        assert label in message
    assert "--conc" in message
    assert not (tmp_path / "kd").exists()


def test_a_survey_is_refused_too(index_named_series, tmp_path):
    code, summary = _kd(index_named_series / "series_analysis_tidy.csv", tmp_path / "s",
                        "--conc-units", "equivalents", "--survey")
    assert code == 1
    assert "sample_2" in summary["error"]["message"]


def test_given_concentrations_fit_those_points(index_named_series, tmp_path):
    code, summary = _kd(index_named_series / "series_analysis_tidy.csv", tmp_path / "kd",
                        "--conc", "0,0.5,3", "--conc-units", "equivalents")
    assert code == 0
    meta = json.loads(Path(summary["json_file"]).read_text())["metadata"]
    assert meta["concentrations"] == [0.0, 25.0, 150.0]


def test_without_the_metadata_a_label_must_be_a_number(index_named_series, tmp_path):
    """The tidy CSV moved away from its series_metadata.json: the label alone decides."""
    csv = tmp_path / "series_analysis_tidy.csv"
    shutil.copy(index_named_series / "series_analysis_tidy.csv", csv)
    code, summary = _kd(csv, tmp_path / "kd", "--conc-units", "equivalents")
    assert code == 1
    assert "sample_1" in summary["error"]["message"]


@pytest.mark.skipif(not TIDY.exists(), reason="data_test not present")
def test_numeric_labels_still_carry_the_concentrations(tmp_path):
    csv = tmp_path / "DNAJA1_HSPA8_tidy.csv"
    shutil.copy(TIDY / csv.name, csv)
    code, summary = _kd(csv, tmp_path / "kd", "--conc-units", "equivalents")
    assert code == 0
    meta = json.loads(Path(summary["json_file"]).read_text())["metadata"]
    assert meta["concentrations"] == [0, 6.25, 12.5, 25, 50, 100, 150]


def test_parsed_series_labels_still_carry_the_concentrations(dnaja1_titration_series, tmp_path):
    code, summary = _kd(dnaja1_titration_series / "series_analysis_tidy.csv", tmp_path / "kd",
                        "--conc-units", "equivalents")
    assert code == 0
    meta = json.loads(Path(summary["json_file"]).read_text())["metadata"]
    assert meta["concentrations"] == [0.0, 25.0, 150.0]
