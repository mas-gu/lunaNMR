# ABOUTME: Tests for the 1:1 binding isotherm models used by the Kd titration fitter.
# ABOUTME: Covers fraction_bound, CSP/intensity observables, and the CSP magnitude.

import sys
from pathlib import Path

import numpy as np
import pytest

_KD_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_Kd"
sys.path.insert(0, str(_KD_DIR))


class TestFractionBound:
    def test_zero_ligand_is_zero(self):
        from kd_models import fraction_bound
        assert fraction_bound(0.0, 50.0, 10.0) == 0.0

    def test_monotonic_and_bounded(self):
        from kd_models import fraction_bound
        L = np.array([0.0, 5.0, 25.0, 100.0, 1000.0])
        fb = fraction_bound(L, P=50.0, Kd=10.0)
        assert np.all((fb >= 0.0) & (fb <= 1.0))
        assert np.all(np.diff(fb) > 0)          # increases with ligand
        assert fb[-1] > 0.95                     # saturates at high ligand

    def test_quadratic_accounts_for_protein(self):
        # With P comparable to Kd, free ligand < total -> less bound than the
        # ligand-excess hyperbola L/(L+Kd) would predict.
        from kd_models import fraction_bound
        L, P, Kd = 10.0, 50.0, 10.0
        fb = fraction_bound(L, P, Kd)
        hyperbola = L / (L + Kd)
        assert fb < hyperbola


class TestBindingModels:
    def test_csp_model_scales_fraction(self):
        from kd_models import csp_model, fraction_bound
        L, P, Kd, dd_max = 20.0, 50.0, 10.0, 0.3
        assert csp_model(L, dd_max, Kd, P) == pytest.approx(dd_max * fraction_bound(L, P, Kd))

    def test_intensity_decay_with_plateau(self):
        # I = I_inf + (I0 - I_inf)*exp(-L/Kd): anchored at I0, levels off at I_inf.
        from kd_models import intensity_decay
        assert intensity_decay(0.0, I0=2.0, I_inf=0.3, Kd=10.0) == pytest.approx(2.0)
        assert intensity_decay(1e6, I0=2.0, I_inf=0.3, Kd=10.0) == pytest.approx(0.3)  # plateau
        assert intensity_decay(10.0, I0=2.0, I_inf=0.3, Kd=10.0) == pytest.approx(
            0.3 + 1.7 * np.exp(-1.0))


class TestComputeCsp:
    def test_combines_h_and_scaled_n(self):
        from kd_models import compute_csp
        # CSP = sqrt(dH^2 + (alpha*dN)^2)
        assert compute_csp(0.03, 0.0) == pytest.approx(0.03)
        assert compute_csp(0.0, 1.0, alpha=0.14) == pytest.approx(0.14)
        assert compute_csp(0.03, 0.5, alpha=0.14) == pytest.approx(
            np.sqrt(0.03**2 + (0.14 * 0.5) ** 2))


class TestRoundTrip:
    def test_recovers_injected_kd_via_curve_fit(self):
        # Synthetic CSP isotherm -> bounded curve_fit recovers Kd exactly.
        # Bounds (Kd>=0, dd>=0) are essential: unbounded fits land on a
        # degenerate negative-Kd minimum.
        from scipy.optimize import curve_fit
        from kd_models import csp_model
        P, Kd_true, dd_max_true = 50.0, 12.0, 0.25
        L = np.array([0, 5, 10, 20, 40, 80, 160, 320], dtype=float)
        y = csp_model(L, dd_max_true, Kd_true, P)
        (dd_fit, kd_fit), _ = curve_fit(
            lambda L, dd, kd: csp_model(L, dd, kd, P), L, y,
            p0=[0.1, 20.0], bounds=([0, 0], [np.inf, np.inf]))
        assert kd_fit == pytest.approx(Kd_true, rel=1e-4)
        assert dd_fit == pytest.approx(dd_max_true, rel=1e-4)


class TestObservables:
    """Shared Qt-free observable helpers used by both the GUI viewer and the CLI."""

    def _series(self):
        return {"ppm_x": [8.0, 8.1, 8.3], "ppm_y": [120.0, 120.5, 121.0],
                "height": [1000.0, 500.0, 250.0], "volume": [2000.0, 1000.0, 500.0]}

    def test_peak_present_flags_missing_endpoint(self):
        from kd_models import peak_present
        assert peak_present(self._series(), 0, "height")
        # undetected position sentinel (0.0) -> absent
        assert not peak_present({"ppm_x": [0.0], "ppm_y": [0.0], "height": [5.0]},
                                0, "height")

    def test_pair_observable_csp_and_ratio(self):
        import math
        from kd_models import pair_observable
        s = self._series()
        assert pair_observable(s, 0, 2, "csp", alpha=0.14) == pytest.approx(
            math.hypot(0.3, 0.14 * 1.0))
        assert pair_observable(s, 0, 1, "intensity", value="height") == pytest.approx(0.5)

    def test_ref_point_values_marks_absent_as_nan(self):
        import math
        from kd_models import ref_point_values
        fits = [{"residue": "K14", "series": self._series()},
                {"residue": "A17"}]                       # no series -> NaN
        names, vals = ref_point_values(fits, 0, 1, "intensity", value="height")
        assert names == ["K14", "A17"]
        assert vals[0] == pytest.approx(0.5)
        assert math.isnan(vals[1])
