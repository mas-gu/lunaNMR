# ABOUTME: The series JSON summary reports the peak tracking the engine actually ran.
# ABOUTME: Titration mode forces cascade with no drift limit; the summary must say so.

"""
`series --mode titration` runs cascade tracking with the absolute drift limit off,
whatever --peak-source says. The dry-run reported the requested source ('reference'), the
real run reported none, and its params echo still showed the drift limit on.
"""

import contextlib
import io
import json
from pathlib import Path

import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"


def _dry_run(mode, *extra):
    from lunaNMR.cli import main
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        code = main(["series", "--spectra", str(DNAJA1), "--peaks", str(DNAJA1 / "DNAJA1.txt"),
                     "--out", "unused", "--mode", mode, "--dry-run", "--format", "json", *extra])
    assert code == 0
    summary = json.loads(out.getvalue())
    return summary.get("result", summary)


@pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")
class TestTheDryRunReportsTheTrackingThatWillRun:

    def test_titration_reports_cascade_without_the_drift_limit(self):
        summary = _dry_run("titration")
        assert summary["peak_source"] == "cascade"
        assert summary["requested_peak_source"] == "reference"
        assert summary["drift_limit"] is False

    def test_time_mode_reports_the_requested_tracking(self):
        summary = _dry_run("time", "--peak-source", "cascade")
        assert summary["peak_source"] == "cascade"
        assert summary["requested_peak_source"] == "cascade"
        assert summary["drift_limit"] is True

    def test_the_drift_limit_applies_only_to_a_tracked_peak_source(self):
        assert _dry_run("time", "--peak-source", "reference")["drift_limit"] is False


def test_the_run_reports_the_tracking_it_ran(dnaja1_titration_series):
    summary = json.loads((dnaja1_titration_series.parent / "summary.json").read_text())
    summary = summary.get("result", summary)
    assert summary["peak_source"] == "cascade"
    assert summary["requested_peak_source"] == "reference"
    assert summary["drift_limit"] is False
