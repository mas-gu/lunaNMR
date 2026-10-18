# ABOUTME: peak_detection_matrix.csv reads 1 where a peak stands 3 sigma above the noise, else 0.
# ABOUTME: The floor is the one the Kd loader applies to detection_sigma.

"""
The matrix used to read 1 for every peak a spectrum returned a result for, a fit at 0
sigma included, so it was 0 only where the whole spectrum failed and excluded nothing.
"""

import pandas as pd

from lunaNMR.utils.quality_categories import DETECTION_SIGMA_FLOOR


def _long(series_dir):
    tidy = pd.read_csv(series_dir / "series_analysis_tidy.csv", dtype={"spectrum_name": str})
    matrix = pd.read_csv(series_dir / "peak_detection_matrix.csv")
    points = [c for c in matrix.columns
              if c not in ("Peak_Number", "Assignment", "Reference_X", "Reference_Y")]
    detected = matrix.melt(id_vars="Assignment", value_vars=points,
                           var_name="spectrum_name", value_name="detected")
    return tidy.merge(detected, left_on=["assignment", "spectrum_name"],
                      right_on=["Assignment", "spectrum_name"], validate="one_to_one")


def test_the_floor_is_the_kd_loaders():
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "modules" / "dynamiXs_v2o0"
                           / "dynamiXs_Kd"))
    import kd_input
    assert DETECTION_SIGMA_FLOOR == kd_input.DETECTION_SIGMA_FLOOR == 3.0


def test_a_peak_fitted_at_the_noise_reads_zero(dnaja1_titration_series):
    rows = _long(dnaja1_titration_series)
    weak = rows[rows["detection_sigma"] < 1.0]
    assert len(weak) > 0, "the series has no peak fitted at the noise to check"
    assert (weak["detected"] == 0).all()


def test_a_strong_peak_reads_one(dnaja1_titration_series):
    rows = _long(dnaja1_titration_series)
    strong = rows[rows["detection_sigma"] >= 10 * DETECTION_SIGMA_FLOOR]
    assert len(strong) > 0
    assert (strong["detected"] == 1).all()


def test_every_cell_follows_the_floor(dnaja1_titration_series):
    rows = _long(dnaja1_titration_series)
    expected = (rows["detection_sigma"] >= DETECTION_SIGMA_FLOOR).astype(int)
    assert (rows["detected"] == expected).all()
