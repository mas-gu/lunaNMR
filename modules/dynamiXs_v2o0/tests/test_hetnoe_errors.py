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


class TestANegativeNoeHasAPositiveError:
    """A flexible residue's saturated peak is NEGATIVE: on a 600 MHz pair, residue 33 is -47.7
    sigma in the saturated plane and +31.8 sigma in the unsaturated one. Multiplying
    the relative error by the signed hetNOE gave a negative error, which reached the CSV
    and then made the QC plot raise "'yerr' must not contain negative values".
    """

    SAT = {"33": -3614098.0}       # signed sample at residue 33, saturated plane
    UNSAT = {"33": 2603630.3}      # fitted height, unsaturated plane

    def test_the_estimated_error_is_positive(self):
        noe = calculate_hetnoe_from_intensities(saturated_data=self.SAT,
                                                unsaturated_data=self.UNSAT)
        assert noe["33"]["value"] < 0
        assert noe["33"]["error"] == pytest.approx(abs(noe["33"]["value"]) * _FABRICATED)

    def test_the_propagated_error_is_positive(self):
        noe = calculate_hetnoe_from_intensities(
            saturated_data=self.SAT, unsaturated_data=self.UNSAT,
            saturated_errors={"33": 7.58e4}, unsaturated_errors={"33": 7.62e4})
        value = noe["33"]["value"]
        expected = abs(value) * np.sqrt((7.58e4 / 3614098.0) ** 2 + (7.62e4 / 2603630.3) ** 2)
        assert noe["33"]["error"] == pytest.approx(expected)

    def test_the_qc_plot_accepts_it(self, tmp_path):
        import matplotlib
        matplotlib.use('Agg')
        from dynamiXs_integrated.data_converters import plot_hetnoe_vs_residue
        noe = calculate_hetnoe_from_intensities(saturated_data=self.SAT,
                                                unsaturated_data=self.UNSAT)
        out = tmp_path / "hetnoe.pdf"
        plot_hetnoe_vs_residue(noe, str(out))
        assert out.exists() and out.stat().st_size > 0

    def test_the_gui_worker_writes_a_positive_error(self, tmp_path):
        from types import SimpleNamespace
        from workers import SpectraProcessingWorker
        sat, unsat = tmp_path / "sat.csv", tmp_path / "unsat.csv"
        sat.write_text("33,-3614098.0\n")
        unsat.write_text("33,2603630.3\n")
        host = SimpleNamespace(params=SimpleNamespace(output_dir=str(tmp_path)))
        out = SpectraProcessingWorker._calculate_hetnoe_ratio(host, str(sat), str(unsat))
        import pandas as pd
        row = pd.read_csv(out).iloc[0]
        assert row["hetNOE"] < 0 and row["hetNOE_err"] > 0
