# ABOUTME: A model-free summary must count the fits that stopped at the optimiser's iteration
# ABOUTME: cap and name their residues. fit_success is True for those too, so N/N hid them.

"""
The model-free fit stops either on convergence or on its evaluation budget, and reports
`fit_success=True` either way. A fit that ran out of budget is the weaker result, and the
summary gave a caller nothing to tell it apart by. Each fit says which in its own
`at_iteration_cap` column; the summary must agree with it.
"""

import json
import os
from pathlib import Path

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

RUN = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
       / "run_20260825_1405")
NOE = RUN / "results" / "hetnoe"

pytestmark = pytest.mark.skipif(not RUN.exists(), reason="data_test not present")


def _matrix(name):
    return str(RUN / "series" / name / "peak_intensity_matrix.csv")


@pytest.fixture(scope="module")
def modelfree_run(tmp_path_factory):
    from lunaNMR.cli import main
    out = tmp_path_factory.mktemp("mf")
    argv = ["dynamixs", "modelfree", "--f1-t1", _matrix("600_T1"),
            "--f1-t2", _matrix("600_T2"),
            "--f1-noe-sat", str(NOE / "f600_noe_sat.csv"),
            "--f1-noe-unsat", str(NOE / "f600_noe_unsat.csv"),
            "--field1-freq", "600", "--out", str(out), "--prefix", "mf",
            "--n-monte-carlo", "2"]
    return main, argv, out


def test_the_summary_counts_and_names_the_capped_fits(modelfree_run, capsys):
    main, argv, out = modelfree_run
    assert main(argv + ["--format", "json"]) == 0
    summary = json.loads(capsys.readouterr().out)
    basic = pd.read_csv(summary["output_files"]["basic_csv"])
    capped = basic.loc[basic["at_iteration_cap"].astype(bool), "Residue"].astype(str)

    assert summary["n_at_iteration_cap"] == len(capped)
    assert sorted(summary["at_iteration_cap_residues"]) == sorted(capped)
    # Reported, not gated on: which fits count as failed is unchanged.
    assert summary["n_successful"] == summary["n_residues"]


def test_the_text_summary_states_the_count(modelfree_run, capsys):
    main, argv, out = modelfree_run
    assert main(argv) == 0
    text = capsys.readouterr().out
    basic = pd.read_csv(out / "mf_spectral_density_basic.csv")
    n_capped = int(basic["at_iteration_cap"].astype(bool).sum())
    assert f"{n_capped} fits stopped at the iteration cap" in text
