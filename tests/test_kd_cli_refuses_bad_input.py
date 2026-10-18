# ABOUTME: `lunaNMR kd` and `export kd` refuse input they cannot fit or render, in the
# ABOUTME: caller's format, instead of exiting 0 on a result made from nothing.

import contextlib
import io
import json
import shutil
from pathlib import Path

import pytest

TIDY = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "tidy"
ABSOLUTE = "0,6.25,12.5,25,50,100,150"

pytestmark = pytest.mark.skipif(not TIDY.exists(), reason="data_test not present")


@pytest.fixture
def tidy(tmp_path):
    dest = tmp_path / "DNAJA1_HSPA8_tidy.csv"
    shutil.copy(TIDY / "DNAJA1_HSPA8_tidy.csv", dest)
    return dest


def _run(*args):
    from lunaNMR.cli import main
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        code = main([*args, "--format", "json"])
    text = out.getvalue()
    return code, (json.loads(text) if text.strip() else {})


@pytest.mark.parametrize("p0", ["0", "-50", "nan"])
def test_a_protein_concentration_that_is_not_positive_is_refused(tidy, tmp_path, p0):
    # P0 = 0 divides fraction_bound by zero, so every CSP fit fails while the intensity
    # fits still count the residue as fitted; a negative P0 returns a 'reliable' Kd.
    code, summary = _run("kd", "--input", str(tidy), "--out", str(tmp_path / "out"),
                         f"--p0={p0}", "--conc", ABSOLUTE)
    assert code == 1
    assert summary["ok"] is False
    assert "P0" in summary["error"]["message"]


def _series_summary():
    """processing_summary.csv of a real series run: a CSV with no titration points."""
    found = sorted((TIDY.parents[1]).glob("**/series_results_*/processing_summary.csv"))
    return found[0] if found else None


@pytest.mark.skipif(_series_summary() is None, reason="data_test not present")
def test_a_csv_that_is_no_titration_is_refused(tmp_path):
    code, summary = _run("kd", "--input", str(_series_summary()), "--out", str(tmp_path / "out"),
                         "--p0", "50", "--observable", "intensity", "--conc-units", "absolute")
    assert code == 1
    assert summary["ok"] is False
    assert "no titration point" in summary["error"]["message"]
    assert not (tmp_path / "out").exists()


def test_an_unknown_export_kind_is_reported_on_stdout_under_json(tidy, tmp_path):
    code, fit = _run("kd", "--input", str(tidy), "--out", str(tmp_path / "fit"), "--p0", "50",
                     "--conc", ABSOLUTE, "--observable", "intensity")
    assert code == 0
    code, summary = _run("export", "kd", "--json", fit["json_file"],
                         "--out", str(tmp_path / "export"), "--kind", "curvs")
    assert code == 1
    assert summary["ok"] is False
    assert "curvs" in summary["error"]["message"]
