# ABOUTME: An empty intensity cell in a hetNOE input is an unmeasured residue, not a NaN hetNOE.
# ABOUTME: It is left out of the output and counted apart from a non-positive reference.

"""
The series leaves a bound-pinned point empty. `I_unsat <= 0` is False for NaN, so an
empty cell went through as hetNOE = nan: written, counted in n_residues, and reported
as nothing dropped.
"""

import json
import os
from pathlib import Path

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

NOE = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
       / "run_20260825_1405" / "results" / "hetnoe")

pytestmark = pytest.mark.skipif(not NOE.exists(), reason="data_test not present")


def _blanked(src, dst, row):
    frame = pd.read_csv(src, header=None)
    frame.iloc[row, 1] = float('nan')
    frame.to_csv(dst, header=False, index=False)
    return frame


def test_an_empty_cell_drops_the_residue_and_says_so(tmp_path, capsys):
    from lunaNMR.cli import main
    sat = _blanked(NOE / "f600_noe_sat.csv", tmp_path / "sat.csv", 1)
    unsat = _blanked(NOE / "f600_noe_unsat.csv", tmp_path / "unsat.csv", 0)
    assert main(["dynamixs", "hetnoe", "--sat", str(tmp_path / "sat.csv"),
                 "--unsat", str(tmp_path / "unsat.csv"), "--out", str(tmp_path / "o"),
                 "--prefix", "p", "--format", "json"]) == 0
    summary = json.loads(capsys.readouterr().out)

    blank = {str(sat.iloc[1, 0]), str(unsat.iloc[0, 0])}
    written = pd.read_csv(summary['output'])
    assert not blank & set(written.iloc[:, 0].astype(str))
    assert written.iloc[:, 1].notna().all()
    assert summary['n_dropped_unmeasured'] == 2
    assert summary['n_dropped_nonpositive_ref'] == 0
    assert summary['n_residues'] == summary['n_common'] - 2
