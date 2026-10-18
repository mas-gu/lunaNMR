# ABOUTME: The tidy series CSV keeps each peak's R-squared and detection sigma as numbers.
# ABOUTME: Real 3-point DNAJA1 titration through `lunaNMR series`, against the tracking table.

"""
Every numeric tidy column but the positions was written to one decimal: R-squared took
eleven values, 0.0 to 1.0, so no gate finer than 0.1 could be applied to it, and a peak at
2.97 sigma read 3.0 and passed the Kd loader's 3-sigma floor. The integer and flag columns
came out as floats, or as True/False, depending on the run.
"""

import numpy as np
import pandas as pd
import pytest


@pytest.fixture(scope="module")
def out(dnaja1_titration_series):
    return dnaja1_titration_series


@pytest.fixture(scope="module")
def tidy(out):
    return pd.read_csv(out / "series_analysis_tidy.csv", dtype=str, keep_default_na=False)


def test_r_squared_matches_the_tracking_table(out, tidy):
    tracking = pd.read_csv(out / "comprehensive_peak_tracking.csv")
    for point, rows in tidy.groupby("spectrum_name"):
        reported = dict(zip(tracking.Assignment, tracking[f"{point}_R_Squared"]))
        for _, row in rows.iterrows():
            if row.assignment in reported:
                assert abs(float(row.r_squared) - reported[row.assignment]) < 6e-4, row.assignment


def test_r_squared_is_not_rounded_to_tenths(tidy):
    r2 = tidy.r_squared.astype(float)
    assert not np.allclose(r2, r2.round(1))


def test_detection_sigma_keeps_its_hundredths(tidy):
    sigma = tidy.detection_sigma.replace("", np.nan).astype(float).dropna()
    assert not np.allclose(sigma, sigma.round(1))


def test_peak_number_is_an_integer(tidy):
    assert tidy.peak_number.str.fullmatch(r"\d+").all()


@pytest.mark.parametrize("flag", ["position_pinned", "intensity_pinned"])
def test_pinned_flags_are_written_as_booleans(tidy, flag):
    assert set(tidy[flag]) <= {"True", "False", ""}
    assert "False" in set(tidy[flag])
