# ABOUTME: hetNOE error propagation must use whichever error columns were supplied.
# ABOUTME: Requiring both discarded both, falling back to a fabricated ~2% error.

import sys
from pathlib import Path

import numpy as np
import pytest

_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DIR))

from dynamiXs_integrated.data_converters import calculate_hetnoe_from_intensities

SAT = {"8.TYR": 780.0, "9.LEU": 700.0}
UNSAT = {"8.TYR": 1000.0, "9.LEU": 1000.0}
SAT_ERR = {"8.TYR": 40.0, "9.LEU": 40.0}
UNSAT_ERR = {"8.TYR": 30.0, "9.LEU": 30.0}

# The fabricated fallback: noe * 0.02, ~0.016 at hetNOE 0.78. A realistic per-plane
# floor of ~0.044 relative gives ~0.05, so the fallback is roughly 3x too tight and
# over-weights hetNOE in the model-free fit.
_FABRICATED = 0.02


def _err(sat_err, unsat_err):
    noe = calculate_hetnoe_from_intensities(
        saturated_data=SAT, unsaturated_data=UNSAT,
        saturated_errors=sat_err, unsaturated_errors=unsat_err)
    return noe["8.TYR"]["error"], noe["8.TYR"]["value"]


class TestOneErrorColumnIsBetterThanNone:
    """`if saturated_errors is not None and unsaturated_errors is not None` — an `and`.
    A 3-column --sat paired with a 2-column --unsat threw away the errors the user did
    supply and substituted a flat 2% of the ratio.
    """

    def test_both_supplied_propagates_both(self):
        err, noe = _err(SAT_ERR, UNSAT_ERR)
        expected = noe * np.sqrt((40.0 / 780.0) ** 2 + (30.0 / 1000.0) ** 2)
        assert err == pytest.approx(expected)

    def test_only_saturated_errors_are_still_used(self):
        err, noe = _err(SAT_ERR, None)
        assert err == pytest.approx(noe * (40.0 / 780.0))
        assert err != pytest.approx(noe * _FABRICATED)

    def test_only_unsaturated_errors_are_still_used(self):
        err, noe = _err(None, UNSAT_ERR)
        assert err == pytest.approx(noe * (30.0 / 1000.0))
        assert err != pytest.approx(noe * _FABRICATED)

    def test_neither_falls_back_to_the_estimate(self):
        """With nothing supplied there is nothing to propagate, so the documented
        estimate stands — it is only wrong when real errors were available."""
        err, noe = _err(None, None)
        assert err == pytest.approx(noe * _FABRICATED)

    def test_a_residue_missing_from_a_supplied_error_map_does_not_crash(self):
        err, noe = _err({"9.LEU": 40.0}, None)
        assert np.isfinite(err)
