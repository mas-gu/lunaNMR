# ABOUTME: The --format json contract: one parseable object on stdout, success or failure.
# ABOUTME: json.loads(proc.stdout) is the documented usage, and it threw on every error.

import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


def _run(*argv):
    return subprocess.run([sys.executable, "-m", "lunaNMR", *argv],
                          cwd=str(REPO_ROOT), capture_output=True, text=True)


class TestAFailureIsStillParseable:
    """The documented usage is `json.loads(proc.stdout)`. Every failure path wrote the
    message to stderr and nothing at all to stdout, so the parse threw on exactly the
    runs a caller most needs to inspect -- and the exception it raised said nothing
    about what went wrong.
    """

    def test_a_missing_input_still_emits_json(self, tmp_path):
        p = _run("kd", "--input", str(tmp_path / "nope.csv"), "--out", str(tmp_path / "o"),
                 "--p0", "50", "--format", "json")
        assert p.returncode == 1
        payload = json.loads(p.stdout)          # must not raise
        assert payload["ok"] is False
        assert payload["error"]["type"]
        assert "nope.csv" in payload["error"]["message"]

    def test_the_error_names_the_command(self, tmp_path):
        p = _run("dynamixs", "t1t2", "--exp", "T2", "--input", str(tmp_path / "nope.csv"),
                 "--out", str(tmp_path / "o"), "--format", "json")
        assert json.loads(p.stdout)["command"] == "dynamixs t1t2"

    def test_text_mode_is_unchanged(self, tmp_path):
        """Humans keep the stderr line; only --format json gains the object."""
        p = _run("kd", "--input", str(tmp_path / "nope.csv"), "--out", str(tmp_path / "o"),
                 "--p0", "50")
        assert p.returncode == 1
        assert p.stdout == ""
        assert p.stderr.startswith("error:")


class TestASuccessIsLabelled:
    def test_ok_and_schema_version_are_present(self, tmp_path):
        csv = tmp_path / "series_analysis_tidy.csv"
        csv.write_text("spectrum_name,assignment,ppm_x,ppm_y,height,volume\n"
                       "0,K14,7.5,110.0,200.0,2000.0\n1,K14,7.51,110.05,150.0,1500.0\n"
                       "2,K14,7.53,110.2,90.0,900.0\n0,A17,8.0,120.0,100.0,1000.0\n"
                       "1,A17,8.02,120.1,80.0,800.0\n2,A17,8.05,120.3,50.0,500.0\n")
        p = _run("kd", "--input", str(csv), "--out", str(tmp_path / "o"),
                 "--p0", "50", "--format", "json")
        assert p.returncode == 0, p.stderr[-500:]
        payload = json.loads(p.stdout)
        assert payload["ok"] is True
        assert isinstance(payload["schema_version"], int)
        assert payload["command"] == "kd"

    def test_a_dry_run_is_labelled_too(self, tmp_path):
        (tmp_path / "in.csv").write_text("x\n")
        p = _run("kd", "--input", str(tmp_path / "in.csv"), "--out", str(tmp_path / "o"),
                 "--p0", "50", "--dry-run", "--format", "json")
        payload = json.loads(p.stdout)
        assert payload["ok"] is True and payload["dry_run"] is True


class TestNonFiniteNumbersStayValidJson:
    """json.dumps writes bare NaN, which is not JSON: jq, Go and Rust all reject it.
    mean_t2 is NaN when nothing fits -- so the output became unparseable on precisely
    the degenerate runs worth looking at.
    """

    def test_a_run_where_nothing_fits_is_still_parseable(self, tmp_path):
        import numpy as np
        import pandas as pd
        # Flat series: no measurable decay, so every fit is excluded and mean_t2 is NaN.
        times = [8.0, 51.0, 102.0, 204.0]
        rows = [["Peak_Number", "Assignment", "Reference_X", "Reference_Y"]
                + [str(t) for t in times]]
        for i, name in enumerate(("K3", "E4"), start=1):
            rows.append([i, name, 8.0, 120.0] + [900.0] * len(times))
        csv = tmp_path / "peak_intensity_matrix.csv"
        pd.DataFrame(rows).to_csv(csv, header=False, index=False)

        p = _run("dynamixs", "t1t2", "--exp", "T2", "--input", str(csv),
                 "--out", str(tmp_path / "o"), "--format", "json")
        assert "NaN" not in p.stdout, "bare NaN is not valid JSON"
        payload = json.loads(p.stdout)
        assert payload["mean_t2"] is None


class TestStdoutStaysCleanUnderParallelism:
    """`--format json` promises one object on stdout and nothing else. The guarantee is
    made at the file-descriptor level (os.dup2), because worker processes inherit fd 1
    and write past any Python-level redirect.

    That is the cleverest thing in the CLI and nothing exercised it: the existing test
    used capsys, which cannot observe a raw fd-1 write, on a single-process command. A
    real subprocess with real workers is the only way to see it.
    """

    def test_a_parallel_series_run_emits_only_json(self, tmp_path):
        data = REPO_ROOT / "data_example" / "2DHSQC"
        spectra = sorted(data.glob("600_T1_0o*.ft"))
        peaks = data / "600_assi.txt"
        if len(spectra) < 2 or not peaks.exists():
            pytest.skip("data_example/2DHSQC spectra not available")
        p = _run("series", "--spectra", str(data), "--peaks", str(peaks),
                 "--out", str(tmp_path / "out"), "--parallel", "--format", "json")
        assert p.returncode == 0, p.stderr[-800:]
        payload = json.loads(p.stdout)      # the whole of stdout, not a fragment of it
        assert payload["ok"] is True
        assert payload["command"] == "series"
        assert payload["spectra_fitted"] >= 1
        assert p.stderr.strip(), "engine chatter should be on stderr, not discarded"


class TestTheTimeUnitsDefaultsArePinned:
    """The three --time-units defaults differ deliberately: labels-only on t1t2 defaults
    to s, the two that rescale default to ms. Nothing pinned them, so a tidy-up that
    unified them would silently rescale published numbers -- and each flag now defaults
    to None so the sidecar can supply the value, which makes the intended fallback
    invisible in the parser.
    """

    @pytest.mark.parametrize("subcommand,expected", [
        ("_TIME_UNITS_DEFAULT", "s"),
        ("_T1RHO_TIME_UNITS_DEFAULT", "ms"),
        ("_MODELFREE_UNITS_DEFAULT", "ms"),
    ])
    def test_the_documented_fallback(self, subcommand, expected):
        import lunaNMR.cli as cli
        assert getattr(cli, subcommand) == expected

    @pytest.mark.parametrize("argv,dest", [
        (["dynamixs", "t1t2"], "time_units"),
        (["dynamixs", "methyl-t2"], "time_units"),
        (["dynamixs", "t1rho"], "time_units"),
        (["dynamixs", "modelfree"], "f1_t1_units"),
    ])
    def test_the_flag_itself_defaults_to_none(self, argv, dest):
        """None, not the documented value: the handler must be able to tell 'the user
        said nothing' from 'the user asked for this', or the sidecar can never win."""
        import argparse
        from lunaNMR.cli import build_parser
        parser = build_parser()
        for action in parser._actions:
            if isinstance(action, argparse._SubParsersAction):
                sub = action.choices[argv[0]]
                for a2 in sub._actions:
                    if isinstance(a2, argparse._SubParsersAction):
                        target = a2.choices[argv[1]]
                        assert target.get_default(dest) is None
                        return
        pytest.fail("subcommand not found")
