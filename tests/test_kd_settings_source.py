# ABOUTME: `lunaNMR kd` takes concentrations and their units together, and stored settings only
# ABOUTME: from a file that describes THIS input. Real DNAJA1/ERDJ6 tidy CSVs side by side.

"""
Four ways a Kd run took numbers the caller never gave it, all exiting 0:
- --conc given without --conc-units took the units from a stored file, so absolute
  concentrations were multiplied by P0 (CSP Kd 7219 instead of 70).
- the newest *_kd_params.json in the folder was used whichever dataset wrote it, so ERDJ6
  was fitted with DNAJA1's concentrations and intensity scales;
- a malformed params file was dropped silently, and the valid one behind it with it;
- with no concentrations at all the point labels were fitted as micromolar.
"""

import contextlib
import io
import json
import shutil
from pathlib import Path

import pytest

TIDY = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "tidy"
ABSOLUTE = "0,6.25,12.5,25,50,100,150"
EQUIVALENTS = "0,0.125,0.25,0.5,1,2,3"

pytestmark = pytest.mark.skipif(not TIDY.exists(), reason="data_test not present")


@pytest.fixture
def folder(tmp_path):
    for name in ("DNAJA1_HSPA8_tidy.csv", "ERDJ6_HSPA8_tidy.csv"):
        shutil.copy(TIDY / name, tmp_path / name)
    return tmp_path


def _kd(*args):
    from lunaNMR.cli import main
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        code = main(["kd", "--format", "json", "--p0", "50", "--observable", "intensity", *args])
    text = out.getvalue()
    return code, (json.loads(text) if text.strip() else {})


def _payload(summary):
    return summary.get("result", summary)


def _metadata(summary):
    return json.loads(Path(_payload(summary)["json_file"]).read_text())["metadata"]


def _survey(folder, name, prefix, *extra):
    code, _ = _kd("--survey", "--input", str(folder / name), "--out", str(folder / f"{prefix}_s"),
                  "--prefix", prefix, *extra)
    assert code == 0


def test_given_concentrations_are_never_re_unitised_by_a_stored_file(folder):
    _survey(folder, "DNAJA1_HSPA8_tidy.csv", "DNAJA1", "--conc", EQUIVALENTS,
            "--conc-units", "equivalents")
    code, summary = _kd("--input", str(folder / "DNAJA1_HSPA8_tidy.csv"),
                        "--out", str(folder / "fit"), "--conc", ABSOLUTE)
    assert code == 0
    assert _metadata(summary)["concentrations"] == [0, 6.25, 12.5, 25, 50, 100, 150]


def test_another_datasets_settings_are_not_applied(folder):
    _survey(folder, "DNAJA1_HSPA8_tidy.csv", "DNAJA1", "--conc", ABSOLUTE,
            "--intensity-scale", "1,1,1,1,1,1,4")
    code, summary = _kd("--input", str(folder / "ERDJ6_HSPA8_tidy.csv"),
                        "--out", str(folder / "fit"), "--conc", ABSOLUTE)
    assert code == 0
    assert _payload(summary)["params_source"] is None
    assert _metadata(summary)["intensity_scales"] in (None, [1.0] * 7)
    skipped = _payload(summary)["params_skipped"]
    assert any("DNAJA1_kd_params.json" in s["path"] for s in skipped)


def test_this_inputs_stored_settings_supply_the_pair(folder):
    _survey(folder, "DNAJA1_HSPA8_tidy.csv", "DNAJA1", "--conc", EQUIVALENTS,
            "--conc-units", "equivalents")
    code, summary = _kd("--input", str(folder / "DNAJA1_HSPA8_tidy.csv"),
                        "--out", str(folder / "fit"))
    assert code == 0
    assert _metadata(summary)["concentrations"] == [0, 6.25, 12.5, 25, 50, 100, 150]


def test_a_malformed_settings_file_is_reported_and_the_valid_one_used(folder):
    _survey(folder, "DNAJA1_HSPA8_tidy.csv", "DNAJA1", "--conc", EQUIVALENTS,
            "--conc-units", "equivalents")
    (folder / "kd_params.json").write_text("{ not json")
    code, summary = _kd("--input", str(folder / "DNAJA1_HSPA8_tidy.csv"),
                        "--out", str(folder / "fit"))
    assert code == 0
    assert _payload(summary)["params_source"].endswith("DNAJA1_kd_params.json")
    assert any(s["path"].endswith("kd_params.json") and "kd_params.json" == Path(s["path"]).name
               for s in _payload(summary)["params_skipped"])


def test_no_concentrations_from_anywhere_is_refused(folder):
    code, summary = _kd("--input", str(folder / "ERDJ6_HSPA8_tidy.csv"),
                        "--out", str(folder / "fit"))
    assert code == 1
    assert "--conc" in json.dumps(summary)


def test_point_labels_may_be_declared_as_the_concentrations(folder):
    """--conc-units alone says the spectrum names carry the concentrations."""
    code, summary = _kd("--input", str(folder / "ERDJ6_HSPA8_tidy.csv"),
                        "--out", str(folder / "fit"), "--conc-units", "equivalents")
    assert code == 0
    assert _metadata(summary)["concentrations"] == [0, 6.25, 12.5, 25, 50, 100, 150]


def test_a_survey_reports_the_command_it_was_run_as(folder):
    """`command` is the subcommand path, and a survey is `kd --survey`: the dry-run of
    the same command line already said 'kd'."""
    argv = ("--survey", "--input", str(folder / "DNAJA1_HSPA8_tidy.csv"),
            "--out", str(folder / "s"), "--conc", ABSOLUTE)
    assert _kd(*argv, "--dry-run")[1]["command"] == "kd"
    code, summary = _kd(*argv)
    assert code == 0
    assert summary["command"] == "kd"
    assert summary["survey_file"]


def test_a_refusal_names_the_stored_file_it_passed_over(folder):
    """The refusal said no stored settings were found when one was found and skipped,
    and named it only on stderr: the envelope carries it as the success summary does."""
    _survey(folder, "DNAJA1_HSPA8_tidy.csv", "DNAJA1", "--conc", ABSOLUTE)
    code, summary = _kd("--input", str(folder / "ERDJ6_HSPA8_tidy.csv"),
                        "--out", str(folder / "fit"))
    assert code == 1
    error = summary["error"]
    (skipped,) = error["details"]["params_skipped"]
    assert Path(skipped["path"]).name == "DNAJA1_kd_params.json"
    assert "DNAJA1_HSPA8_tidy.csv" in skipped["reason"]
    assert "were found" not in error["message"]
    assert "DNAJA1_kd_params.json" in error["message"]
    assert skipped["reason"] in error["message"]


def test_a_refusal_with_nothing_stored_says_so(folder):
    code, summary = _kd("--input", str(folder / "ERDJ6_HSPA8_tidy.csv"),
                        "--out", str(folder / "fit"))
    assert code == 1
    assert summary["error"]["details"]["params_skipped"] == []
    assert "No stored settings for this input were found" in summary["error"]["message"]
