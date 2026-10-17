# ABOUTME: A relaxation fit given no output prefix writes to a folder of its own, per call.
# ABOUTME: The default was one fixed path inside the module, shared by every fit of a type.

import sys
from pathlib import Path

import pytest

_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DIR))
sys.path.insert(0, str(_DIR / "dynamiXs_T1_T2"))

MATRIX = (Path(__file__).resolve().parents[4] / "data_test" / "data_JA1_relaxation"
          / "run_20260825_1405" / "series" / "600_T1" / "peak_intensity_matrix.csv")


@pytest.mark.skipif(not MATRIX.exists(), reason="data_test not present")
def test_two_unprefixed_fits_do_not_share_files():
    from dynamiXs_integrated.fitting_wrapper import fit_relaxation_data
    _, first = fit_relaxation_data(str(MATRIX), "T1", field_freq=600.0)
    _, second = fit_relaxation_data(str(MATRIX), "T1", field_freq=600.0)
    assert Path(first).parent != Path(second).parent
    assert not Path(first).resolve().is_relative_to(_DIR)
