# ABOUTME: A clustered peak is given the same quality label as an isolated peak with its R-squared.
# ABOUTME: Real 3-point DNAJA1 titration through `lunaNMR series`, read from the tidy CSV.

"""
Every member of a successful 2D cluster fit was labelled 'Excellent' if R-squared > 0.9,
else 'Good' - down to R-squared 0.0 - while an isolated peak with the same fit was Poor or
Failed. The Kd loader drops Failed points, so whether a noise-level point was kept depended
on whether it happened to sit in a cluster.
"""

import pandas as pd
import pytest


@pytest.fixture(scope="module")
def tidy(dnaja1_titration_series):
    return pd.read_csv(dnaja1_titration_series / "series_analysis_tidy.csv")


def test_no_fit_below_r_squared_half_is_labelled_fair_or_better(tidy):
    weak = tidy[tidy.r_squared < 0.5]
    assert len(weak)
    mislabelled = weak[weak.quality.isin(["Excellent", "Good", "Fair"])]
    assert mislabelled.empty, mislabelled[["spectrum_name", "assignment", "r_squared", "quality"]]
