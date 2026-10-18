# ABOUTME: The methyl-t2 summary must count the fits whose bi-exponential split is undetermined.
# ABOUTME: n_fitted counts every fit that did not raise, collapsed to mono-exp or not.

"""
The fitter flags a residue `bi_exp_unidentifiable` when T2a and T2b cannot be told
apart. The flag reached the per-row files, but the summary said only "77/77 fitted",
which gives a caller no reason to open them.
"""

import json
import os
from pathlib import Path

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

METHYL = (Path(__file__).resolve().parents[2] / "data_test" / "methyl" / "sampleNBD"
          / "series_results_sampleNBD_20260818_124619" / "peak_intensity_matrix.csv")

pytestmark = pytest.mark.skipif(not METHYL.exists(), reason="data_test not present")


def test_the_summary_counts_unidentifiable_fits(tmp_path, capsys):
    from lunaNMR.cli import main
    assert main(["dynamixs", "methyl-t2", "--input", str(METHYL), "--out", str(tmp_path),
                 "--prefix", "m", "--format", "json"]) == 0
    summary = json.loads(capsys.readouterr().out)
    rows = pd.read_csv(summary['results_file'], sep='\t')
    fitted = rows[rows['Success'] == 'Yes']
    expected = int((fitted['bi_exp_unidentifiable'] == 'Yes').sum())
    assert expected > 0
    assert summary['n_bi_exp_unidentifiable'] == expected
