# ABOUTME: Two single-field `dynamixs density` runs sharing an --out keep a plot each.
# ABOUTME: A fixed plot filename let the second run overwrite the first's and report none.

"""
The single-field plotter fell back to one fixed PDF name in the working directory, and
the CLI found it by chdir-ing into --out and diffing the listing. The second run into a
folder overwrote the first run's figure, and listed no plot because the file was not new.
"""

import json
import os
from pathlib import Path

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

BASIC = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
         / "run_20260825_1405" / "modelfree" / "f600" / "f600_spectral_density_basic.csv")

pytestmark = pytest.mark.skipif(not BASIC.exists(), reason="data_test not present")


def _density_input(path):
    """The measured R1/R2/hetNOE of a real 600 MHz model-free run, as density reads them."""
    frame = pd.read_csv(BASIC)
    frame = frame[['Residue', 'R1', 'R1_err', 'R2', 'R2_err', 'hetNOE', 'hetNOE_err']]
    frame.columns = ['Residue', 'R1', 'R1err', 'R2', 'R2err', 'hetNOE', 'hetNOEerr']
    frame.to_csv(path, index=False)


def test_each_run_reports_its_own_plot(tmp_path, capsys):
    from lunaNMR.cli import main
    table = tmp_path / "f600_input.csv"
    _density_input(table)
    out = tmp_path / "dens"
    cwd = os.getcwd()
    plots = {}
    for prefix in ("A", "B"):
        assert main(["dynamixs", "density", "--input", str(table), "--out", str(out),
                     "--field1-freq", "600", "--no-parallel", "--prefix", prefix,
                     "--format", "json"]) == 0
        outputs = json.loads(capsys.readouterr().out)['outputs']
        plots[prefix] = [o for o in outputs if o.endswith('.pdf')]
        assert os.getcwd() == cwd
    assert len(plots['A']) == 1 and len(plots['B']) == 1
    assert plots['A'] != plots['B']
    assert all(Path(p[0]).is_file() for p in plots.values())
