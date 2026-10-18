# ABOUTME: A spectrum or peak with no fit leaves its matrix cells EMPTY, never a 0.0 measurement.
# ABOUTME: Real JA1 T1 series with one truncated spectrum; the run names the spectrum it lost.

"""
A spectrum that failed to load or fit, or a peak whose fit failed, was written as 0.0 in
peak_intensity_matrix.csv and peak_volume_matrix.csv. The T1/T2 fitters keep every finite
value, so they fitted it as a decay to zero. The run said ok, and named nothing it lost.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

T1 = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
      / "600_WT" / "T1")
GOOD = ("600_T1_A1_WT_0ms.ft", "600_T1_A1_WT_300ms.ft", "600_T1_A1_WT_900ms.ft")
BROKEN = "600_T1_A1_WT_150ms.ft"

pytestmark = pytest.mark.skipif(not T1.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    """(output folder, exit code, JSON summary) of a series with one truncated spectrum."""
    import contextlib
    import io
    from lunaNMR.cli import main
    root = tmp_path_factory.mktemp("t1bad")
    spectra = root / "spectra"
    spectra.mkdir()
    for name in GOOD:
        (spectra / name).symlink_to(T1 / name)
    (spectra / BROKEN).write_bytes((T1 / BROKEN).read_bytes()[:4096])
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        code = main(["series", "--spectra", str(spectra), "--peaks",
                     str(T1 / "600_T1_assignment.txt"), "--mode", "time",
                     "--peak-source", "reference", "--out", str(root / "out"),
                     "--format", "json"])
    return root / "out", code, json.loads(out.getvalue())


def test_the_lost_spectrum_is_an_empty_column(run):
    out, code, _ = run
    assert code == 0
    for name in ("peak_intensity_matrix.csv", "peak_volume_matrix.csv"):
        matrix = pd.read_csv(out / name)
        assert matrix["150"].isna().all(), name
        assert matrix["300"].notna().any(), name


def test_the_run_names_the_spectrum_it_lost(run):
    _, _, summary = run
    payload = summary.get("result", summary)
    assert payload["failed_spectra"] == [BROKEN]


def test_a_peak_whose_fit_failed_is_left_empty():
    """Failed is the label for a fit that produced nothing; a decayed peak is Poor."""
    from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor
    tracking = pd.DataFrame({
        'Peak_Number': [1, 2], 'Assignment': ['3LysH', '4GluH'],
        'Reference_X': [8.1, 8.2], 'Reference_Y': [120.0, 121.0],
        '0_Height': [5.0, 0.0], '0_Quality': ['Good', 'Failed'],
        '0_Intensity_Pinned': [False, False]})
    matrix = MultiSpectrumProcessor._fitted_matrix(tracking, '_Height')
    assert matrix['0'].tolist()[0] == 5.0
    assert np.isnan(matrix['0'].tolist()[1])
