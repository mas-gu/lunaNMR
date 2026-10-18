# ABOUTME: A model-free summary must count the fits that stopped at the optimiser's iteration
# ABOUTME: cap and name their residues. fit_success is True for those too, so N/N hid them.

"""
The Monte-Carlo model-free fit stops either on convergence or after `niter` steps, and
reports `fit_success=True` either way. A fit that ran out of iterations is the weaker
result, and the summary gave a caller nothing to tell it apart by.
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
# The main fit's iteration budget (fit_model_free_individual's niter); `iterations` is the
# last loop index, so a fit that used all of it reports one less.
LAST_ITERATION = 10_000 - 1

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
    capped = basic.loc[basic["iterations"] == LAST_ITERATION, "Residue"].astype(str)

    assert len(capped), "the real 600 MHz fit has residues that use every iteration"
    assert summary["n_at_iteration_cap"] == len(capped)
    assert sorted(summary["at_iteration_cap_residues"]) == sorted(capped)
    # Reported, not gated on: which fits count as failed is unchanged.
    assert summary["n_successful"] == summary["n_residues"]


def test_the_text_summary_states_the_count(modelfree_run, capsys):
    main, argv, out = modelfree_run
    assert main(argv) == 0
    text = capsys.readouterr().out
    basic = pd.read_csv(out / "mf_spectral_density_basic.csv")
    n_capped = int((basic["iterations"] == LAST_ITERATION).sum())
    assert f"{n_capped} fits stopped at the iteration cap" in text
