# ABOUTME: `dynamixs t1rho` labels each of its two fits as what it is - T1, and T1rho.
# ABOUTME: Readers go by the label: the GUI conversion finds the T1 fit by its T1 column.

"""
The label only names the output - the time-constant column of the results file, the
JSON's experiment_type, the plot text - but readers go by it: the GUI's T1rho -> T2
conversion looks for a T1/R1 column and found only `T2`.
"""

from pathlib import Path

import pytest

RUN = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
       / "run_20260825_1405")
MATRIX = RUN / "series" / "600_T1" / "peak_intensity_matrix.csv"
PEAKS = RUN.parent / "600_WT" / "T1"


@pytest.mark.skipif(not MATRIX.exists(), reason="data_test not present")
def test_each_fit_is_labelled_by_its_experiment(tmp_path):
    from lunaNMR.cli import main
    peaks = next(PEAKS.glob("*assignment*.txt"))
    # Any decay fits; only the labelling is under test, so one matrix serves as both.
    assert main(["dynamixs", "t1rho", "--input", str(MATRIX), "--t1", str(MATRIX),
                 "--peaks", str(peaks), "--omega1", "2000", "--carrier", "118",
                 "--theta", "45", "--field-freq", "600", "--out", str(tmp_path),
                 "--prefix", "f"]) == 0
    for label in ("T1", "T1rho"):
        header = (tmp_path / f"f_{label}_fit_results.txt").read_text().split("\n")[0]
        assert header.split("\t")[2] == label, header


@pytest.mark.skipif(not MATRIX.exists(), reason="data_test not present")
def test_two_runs_in_one_folder_keep_their_own_fit_json(tmp_path):
    """The fit JSON is named after the run, as its results files are."""
    from lunaNMR.cli import main
    peaks = next(PEAKS.glob("*assignment*.txt"))
    for prefix in ("f600", "f700"):
        assert main(["dynamixs", "t1rho", "--input", str(MATRIX), "--t1", str(MATRIX),
                     "--peaks", str(peaks), "--omega1", "2000", "--carrier", "118",
                     "--theta", "45", "--field-freq", "600", "--out", str(tmp_path),
                     "--prefix", prefix]) == 0
    names = {p.name for p in tmp_path.glob("*_fit_data.json")}
    assert len(names) == 4, names
    assert all(n.startswith(("f600", "f700")) for n in names), names
