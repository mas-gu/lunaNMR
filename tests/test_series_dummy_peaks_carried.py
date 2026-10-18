# ABOUTME: Dummy peaks auto-added on the first spectrum of a series are fitted at every point.
# ABOUTME: Otherwise the first point is deconvolved with a neighbour that no later point has.

"""
Detection adds a dummy for a strong unlisted maximum next to a listed peak, and the dummy
joins that peak's cluster. Dummies are detected on the first spectrum only. With the
reference peak source every later point started again from the bare peak list, so the
dummy existed at point 1 alone: the reference point was fitted with the neighbour modelled
and every other point without it. DNAJA1 gets one dummy, next to 26TyrH and 45PheH.
"""

import contextlib
import io
import warnings
from pathlib import Path

import pandas as pd
import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


def _dummy_points(tidy_path):
    tidy = pd.read_csv(tidy_path)
    dummies = tidy[tidy["assignment"].astype(str).str.startswith("dummy_")]
    return tidy["spectrum_name"].nunique(), dummies.groupby("assignment")["spectrum_name"].nunique()


@pytest.fixture(scope="module")
def reference_series(tmp_path_factory):
    from numba.core.errors import NumbaWarning
    from lunaNMR.cli import main
    out = tmp_path_factory.mktemp("dummies") / "out"
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore", NumbaWarning)
        code = main(["series", "--spectra", str(DNAJA1), "--peaks", str(DNAJA1 / "DNAJA1.txt"),
                     "--out", str(out), "--peak-source", "reference"])
    assert code == 0
    return out


def test_the_reference_peak_source_fits_the_dummy_at_every_point(reference_series):
    n_points, per_dummy = _dummy_points(reference_series / "series_analysis_tidy.csv")
    assert n_points == 7
    assert len(per_dummy) == 1
    assert (per_dummy == n_points).all()


def test_the_dummy_is_fitted_successfully_at_every_point(reference_series):
    tidy = pd.read_csv(reference_series / "series_analysis_tidy.csv")
    dummy = tidy[tidy["assignment"].astype(str).str.startswith("dummy_")]
    assert (dummy["quality"] != "Failed").all()


def test_a_titration_fits_the_dummy_at_every_point(dnaja1_titration_series):
    n_points, per_dummy = _dummy_points(dnaja1_titration_series / "series_analysis_tidy.csv")
    assert len(per_dummy) == 1
    assert (per_dummy == n_points).all()
