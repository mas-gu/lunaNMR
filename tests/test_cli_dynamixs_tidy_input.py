# ABOUTME: Every dynamixs subcommand reading an intensity matrix refuses a series' tidy table by name.
# ABOUTME: Real JA1 series files; the refusal is a usage error naming the matrix beside it.

import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SERIES = (REPO_ROOT.parent / "data_test" / "data_JA1_relaxation" / "run_20260825_1405"
          / "series")
TIDY = SERIES / "600_T1" / "series_analysis_tidy.csv"
MATRIX = SERIES / "600_T1" / "peak_intensity_matrix.csv"
NOE = REPO_ROOT.parent / "data_test" / "data_JA1_relaxation" / "run_20260825_1405" / "results" / "hetnoe"

pytestmark = pytest.mark.skipif(not TIDY.exists(), reason="data_test not present")


def _commands(out):
    t2 = str(SERIES / "600_T2" / "peak_intensity_matrix.csv")
    noe = ["--f1-noe-sat", str(NOE / "f600_noe_sat.csv"),
           "--f1-noe-unsat", str(NOE / "f600_noe_unsat.csv")]
    rho = ["--peaks", str(MATRIX), "--omega1", "1500", "--carrier", "118", "--theta", "35"]
    return {
        "t1t2": ["dynamixs", "t1t2", "--exp", "T1", "--input", str(TIDY), "--out", out],
        "methyl-t2": ["dynamixs", "methyl-t2", "--input", str(TIDY), "--out", out],
        "t1rho --input": ["dynamixs", "t1rho", "--input", str(TIDY), "--t1", str(MATRIX),
                          *rho, "--out", out],
        "t1rho --t1": ["dynamixs", "t1rho", "--input", str(MATRIX), "--t1", str(TIDY),
                       *rho, "--out", out],
        "modelfree --f1-t1": ["dynamixs", "modelfree", "--f1-t1", str(TIDY), "--f1-t2", t2,
                              *noe, "--out", out],
        "modelfree --f1-t2": ["dynamixs", "modelfree", "--f1-t1", str(MATRIX), "--f1-t2",
                              str(TIDY), *noe, "--out", out],
    }


@pytest.mark.parametrize("case", list(_commands("o")))
@pytest.mark.parametrize("dry_run", [False, True])
def test_a_tidy_table_is_refused_naming_the_matrix(case, dry_run, tmp_path, capsys):
    from lunaNMR.cli import main
    argv = _commands(str(tmp_path / "out"))[case] + (["--dry-run"] if dry_run else [])
    capsys.readouterr()
    assert main(argv + ["--format", "json"]) == 2
    captured = capsys.readouterr()
    error = json.loads(captured.out)["error"]
    assert error["type"] == "UsageError"
    flag = case.split()[-1] if case.split()[-1].startswith("--") else "--input"
    assert error["message"].startswith(f"{flag}: ")
    assert str(MATRIX) in error["message"]
    assert not (tmp_path / "out").exists()


def test_the_text_refusal_reaches_stderr():
    out = subprocess.run([sys.executable, "-m", "lunaNMR", "dynamixs", "t1t2", "--exp", "T1",
                          "--input", str(TIDY), "--out", "unused"],
                         cwd=str(REPO_ROOT), capture_output=True, text=True)
    assert out.returncode == 2
    assert "is a series' tidy table" in out.stderr
    assert "could not convert" not in out.stderr + out.stdout
