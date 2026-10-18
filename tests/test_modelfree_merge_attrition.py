# ABOUTME: A model-free summary must count the residues lost where R1, R2 and hetNOE merge.
# ABOUTME: Its residue counts started after that intersection, so N/N hid the lost ones.

"""
Model-free needs all three observables, so the merge keeps only residues present in
every input. Who it dropped went to the progress log alone, and the summary's
residue_counts began at the merged set: 55 of 60 fitted residues read as 55/55.
"""

import json
import os
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

RUN = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
       / "run_20260825_1405")
NOE = RUN / "results" / "hetnoe"

pytestmark = pytest.mark.skipif(not RUN.exists(), reason="data_test not present")


def _matrix(name):
    return str(RUN / "series" / name / "peak_intensity_matrix.csv")


def test_the_summary_names_the_residues_the_merge_dropped(tmp_path, capsys):
    from lunaNMR.cli import main
    assert main(["dynamixs", "modelfree", "--f1-t1", _matrix("600_T1"),
                 "--f1-t2", _matrix("600_T2"),
                 "--f1-noe-sat", str(NOE / "f600_noe_sat.csv"),
                 "--f1-noe-unsat", str(NOE / "f600_noe_unsat.csv"),
                 "--field1-freq", "600", "--out", str(tmp_path / "mf"), "--prefix", "mf",
                 "--n-monte-carlo", "2", "--format", "json"]) == 0
    captured = capsys.readouterr()
    summary = json.loads(captured.out)
    counts = summary['residue_counts']
    excluded = summary['excluded_residues']['field1']

    assert excluded, "the real 600 MHz inputs lose residues at the merge"
    assert f"{len(excluded)} residues excluded (not in all datasets)" in captured.err
    assert counts['field1_excluded_at_merge'] == len(excluded)
