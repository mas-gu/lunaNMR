# ABOUTME: The `<prefix>_dynamixs.json` record holds what a run was given, not the CLI's functions.
# ABOUTME: A function's repr carries its memory address, so two identical runs wrote different records.

import json

import numpy as np
import pandas as pd
import pytest


def _matrix(path):
    delays = [0.0, 0.1, 0.2, 0.4, 0.8, 1.0, 1.5, 2.0]
    header = ["Peak_Number", "Assignment", "Reference_X", "Reference_Y"] + [str(d) for d in delays]
    rows = [[i + 1, f"{i + 3}.0", 8.3, 120.0 + i]
            + [1.0e6 * np.exp(-t / 0.5) for t in delays] for i in range(2)]
    pd.DataFrame(rows, columns=header).to_csv(path, index=False)


@pytest.mark.filterwarnings("ignore:FigureCanvasAgg is non-interactive")
def test_two_identical_t1t2_runs_write_the_same_record(tmp_path):
    from lunaNMR.cli import main
    csv = tmp_path / "peak_intensity_matrix.csv"
    _matrix(csv)
    out = tmp_path / "t1"
    texts = []
    for _ in range(2):
        assert main(["dynamixs", "t1t2", "--input", str(csv), "--out", str(out),
                     "--exp", "T1", "--prefix", "600_T1"]) == 0
        texts.append((out / "600_T1_dynamixs.json").read_text())
    args = json.loads(texts[0])["args"]
    assert "<function" not in texts[0]
    assert "_validate" not in args and "func" not in args
    assert args["exp"] == "T1" and args["prefix"] == "600_T1"
    assert json.loads(texts[0])["args"] == json.loads(texts[1])["args"]
