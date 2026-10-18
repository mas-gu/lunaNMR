# ABOUTME: A units flag that contradicts the units the delay headers state is refused.
# ABOUTME: Those headers are converted to ms on read, so a seconds flag converted them twice.

"""
`parse_delay_column` reads '..._2o4s' as 2400 ms. A caller who sees seconds in the headers
and truthfully passes `--f1-t1-units s` then had the fitted T1 of ~800 (ms) converted as
800 s: R1 came out 1000x too small, with nothing checking. The sidecar check could not
catch it -- a hand-built matrix has no sidecar -- so the headers are checked the same way.
"""

import json
import os
import re
from pathlib import Path

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

RUN = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
       / "run_20260825_1405")
NOE = RUN / "results" / "hetnoe"

pytestmark = pytest.mark.skipif(not RUN.exists(), reason="data_test not present")

RESERVED = ("Peak_Number", "Assignment", "Reference_X", "Reference_Y")


def _header(ms, repeat, unit):
    """A delay header that states its unit, in the filesystem-safe spelling."""
    marker = "b" if repeat else ""
    if unit == "s":
        return f"T1_{ms / 1000:g}".replace(".", "o") + f"s{marker}"
    return f"T1_{ms:g}ms{marker}"


def _matrix_with_stated_units(tmp_path, series, unit):
    """A real series matrix (bare ms headers, no sidecar) relabelled the way a hand-built
    table would be: every delay header writes out its unit."""
    df = pd.read_csv(RUN / "series" / series / "peak_intensity_matrix.csv")
    renamed = {}
    for col in df.columns:
        if col in RESERVED:
            continue
        m = re.fullmatch(r"(\d+(?:\.\d+)?)(_\d+)?", col)
        renamed[col] = _header(float(m.group(1)), bool(m.group(2)), unit)
    target = tmp_path / f"{series}_{unit}"
    target.mkdir()
    path = target / "peak_intensity_matrix.csv"
    df.rename(columns=renamed).to_csv(path, index=False)
    return path


def _run(argv):
    from lunaNMR.cli import main
    return main(argv)


@pytest.mark.parametrize("header_unit,flag", [("s", "s"), ("ms", "s"), ("s", "us")])
@pytest.mark.parametrize("subcommand", ["t1t2", "methyl-t2"])
def test_a_flag_against_stated_units_is_refused(tmp_path, capsys, subcommand,
                                                header_unit, flag):
    csv = _matrix_with_stated_units(tmp_path, "600_T1", header_unit)
    extra = ["--exp", "T1"] if subcommand == "t1t2" else []
    assert _run(["dynamixs", subcommand, "--input", str(csv), *extra,
                 "--out", str(tmp_path / "o"), "--time-units", flag]) == 1
    err = capsys.readouterr().err
    assert "--time-units" in err and "T1_" in err and "ms" in err


@pytest.mark.parametrize("header_unit", ["s", "ms"])
def test_t1t2_with_the_matching_flag_reports_ms(tmp_path, capsys, header_unit):
    csv = _matrix_with_stated_units(tmp_path, "600_T1", header_unit)
    assert _run(["dynamixs", "t1t2", "--exp", "T1", "--input", str(csv),
                 "--out", str(tmp_path / "o"), "--time-units", "ms",
                 "--format", "json"]) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["time_units"] == "ms"
    assert 200 < summary["mean_t2"] < 2000


def _modelfree(tmp_path, t1_csv, units):
    return ["dynamixs", "modelfree", "--f1-t1", str(t1_csv),
            "--f1-t2", str(RUN / "series" / "600_T2" / "peak_intensity_matrix.csv"),
            "--f1-noe-sat", str(NOE / "f600_noe_sat.csv"),
            "--f1-noe-unsat", str(NOE / "f600_noe_unsat.csv"),
            "--field1-freq", "600", "--out", str(tmp_path / "mf"), "--prefix", "mf",
            "--n-monte-carlo", "2", "--f1-t1-units", units]


def test_modelfree_refuses_a_seconds_flag_on_seconds_headers(tmp_path, capsys):
    """The finding's scenario: headers in seconds, the caller says seconds, and R1 would
    have come out 1000x too small."""
    csv = _matrix_with_stated_units(tmp_path, "600_T1", "s")
    assert _run(_modelfree(tmp_path, csv, "s")) == 1
    assert "--f1-t1-units" in capsys.readouterr().err


def test_modelfree_with_the_matching_flag_gives_physical_rates(tmp_path, capsys):
    csv = _matrix_with_stated_units(tmp_path, "600_T1", "s")
    assert _run(_modelfree(tmp_path, csv, "ms") + ["--format", "json"]) == 0
    summary = json.loads(capsys.readouterr().out)
    basic = pd.read_csv(summary["output_files"]["basic_csv"])
    # A 600 MHz backbone R1 is about 1-2 s^-1; a double conversion gives ~0.001.
    assert 0.5 < basic["R1"].median() < 5


def test_t1rho_refuses_a_seconds_flag_on_seconds_headers(tmp_path, capsys):
    """t1rho rescales by --time-units too, through the same check."""
    csv = _matrix_with_stated_units(tmp_path, "600_T1", "s")
    peaks = next((RUN.parent / "600_WT" / "T1").glob("*assignment*.txt"))
    assert _run(["dynamixs", "t1rho", "--input", str(csv), "--t1", str(csv),
                 "--peaks", str(peaks), "--omega1", "2000", "--carrier", "118",
                 "--theta", "45", "--out", str(tmp_path / "o"),
                 "--time-units", "s"]) == 1
    assert "--time-units" in capsys.readouterr().err
