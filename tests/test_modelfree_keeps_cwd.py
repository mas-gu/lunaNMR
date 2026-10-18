# ABOUTME: `dynamixs modelfree` must not change the process working directory while it runs.
# ABOUTME: In a server running calls on threads, a chdir moves every other call's relative paths.

"""
The working directory is process-wide. The pipeline writes every file under its output
folder by absolute path, so moving the process into that folder for the run only
changes what a concurrent call in the same process resolves its relative paths against.
"""

import os
import threading
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

RUN = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
       / "run_20260825_1405")
NOE = RUN / "results" / "hetnoe"

pytestmark = pytest.mark.skipif(not RUN.exists(), reason="data_test not present")


def test_the_working_directory_holds_still_and_stays_clean(tmp_path, monkeypatch):
    from lunaNMR.cli import main
    here = tmp_path / "caller"
    here.mkdir()
    monkeypatch.chdir(here)
    seen, done = set(), threading.Event()

    def watch():
        while not done.is_set():
            seen.add(os.getcwd())
            done.wait(0.005)

    watcher = threading.Thread(target=watch)
    watcher.start()
    try:
        series = RUN / "series"
        assert main(["dynamixs", "modelfree",
                     "--f1-t1", str(series / "600_T1" / "peak_intensity_matrix.csv"),
                     "--f1-t2", str(series / "600_T2" / "peak_intensity_matrix.csv"),
                     "--f1-noe-sat", str(NOE / "f600_noe_sat.csv"),
                     "--f1-noe-unsat", str(NOE / "f600_noe_unsat.csv"),
                     "--field1-freq", "600", "--out", str(tmp_path / "mf"),
                     "--prefix", "mf", "--n-monte-carlo", "2"]) == 0
    finally:
        done.set()
        watcher.join()
    assert seen == {str(here)}
    assert list(here.iterdir()) == []
