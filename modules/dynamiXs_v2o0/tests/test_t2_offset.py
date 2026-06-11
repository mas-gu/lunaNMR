# ABOUTME: Tests that the T1/T2 fitter recovers a non-zero baseline offset C
# ABOUTME: from synthetic mono-exponential decay data: I(t) = A*exp(-t/tau) + C

import sys
from pathlib import Path

import numpy as np
import pytest

# Make dynamiXs_T1_T2/ importable (no __init__.py in that directory).
_FITTER_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_T1_T2"
sys.path.insert(0, str(_FITTER_DIR))


class TestExpDecaySignature:
    """The exp_decay model accepts an additive offset C."""

    def test_exp_decay_accepts_offset(self):
        from fit_Tx_NMRRE import exp_decay

        x = np.array([0.0, 10.0, 20.0])
        A, t2, C = 100.0, 50.0, 5.0
        y = exp_decay(x, A, t2, C)

        expected = A * np.exp(-x / t2) + C
        np.testing.assert_allclose(y, expected, rtol=1e-12)

    def test_exp_decay_offset_default_zero(self):
        """C must default to 0 so callers passing only (x, A, t2) keep working."""
        from fit_Tx_NMRRE import exp_decay

        x = np.array([0.0, 10.0, 20.0])
        y_no_C = exp_decay(x, 100.0, 50.0)
        y_zero_C = exp_decay(x, 100.0, 50.0, 0.0)

        np.testing.assert_allclose(y_no_C, y_zero_C, rtol=1e-12)


class TestFitRecoversOffset:
    """fit_single_residue recovers a non-zero baseline within fit error."""

    def _synthetic(self, A_true=1000.0, t2_true=50.0, C_true=30.0, noise=2.0, seed=42):
        rng = np.random.default_rng(seed)
        x = np.linspace(0.0, 250.0, 24)
        y_clean = A_true * np.exp(-x / t2_true) + C_true
        y = y_clean + rng.normal(0.0, noise, size=x.size)
        return x, y

    def test_C_present_in_result(self):
        from fit_Tx_NMRRE import fit_single_residue

        x, y = self._synthetic()
        res = fit_single_residue(x, y, "TEST", initial_A=900.0, initial_t2=60.0)

        assert "C" in res, "result dict must include fitted offset 'C'"
        assert "C_err" in res, "result dict must include 'C_err'"

    def test_C_recovered_within_3sigma(self):
        from fit_Tx_NMRRE import fit_single_residue

        C_true = 30.0
        x, y = self._synthetic(C_true=C_true)
        res = fit_single_residue(x, y, "TEST", initial_A=900.0, initial_t2=60.0)

        C_fit = res["C"]
        C_err = res["C_err"]
        assert np.isfinite(C_fit) and np.isfinite(C_err)
        assert abs(C_fit - C_true) < 3.0 * C_err, (
            f"C_fit={C_fit:.3f} not within 3*C_err={3*C_err:.3f} of C_true={C_true}"
        )

    def test_t2_still_recovered_with_offset(self):
        """Adding C must not break T2 recovery on noisy data."""
        from fit_Tx_NMRRE import fit_single_residue

        t2_true = 50.0
        x, y = self._synthetic(t2_true=t2_true)
        res = fit_single_residue(x, y, "TEST", initial_A=900.0, initial_t2=60.0)

        assert abs(res["t2"] - t2_true) < 1.0, (
            f"t2_fit={res['t2']:.3f} not within 1 of t2_true={t2_true}"
        )


class TestMulticoreFitter:
    """Same offset behavior for the parallel fitter."""

    def test_exp_decay_offset_in_multi(self):
        from fitmulti__Tx_NMRRE import exp_decay

        x = np.array([0.0, 10.0, 20.0])
        A, t2, C = 100.0, 50.0, 5.0
        y = exp_decay(x, A, t2, C)

        expected = A * np.exp(-x / t2) + C
        np.testing.assert_allclose(y, expected, rtol=1e-12)


class TestSaveResultsTSV:
    """save_results writes C and C_err columns alongside A/T2."""

    def test_tsv_header_and_row_include_C(self, tmp_path):
        from fit_Tx_NMRRE import save_results

        results = [{
            'residue': 'R1',
            'A': 1000.0, 'A_err': 5.0,
            't2': 50.0, 't2_err': 0.5,
            'C': 30.0, 'C_err': 1.2,
        }]
        out = tmp_path / "out.tsv"
        save_results(results, str(out), experiment_type="T2")

        text = out.read_text()
        header, row = text.strip().splitlines()
        assert header.split("\t") == ["Residue", "A", "T2", "C", "A_err", "T2_err", "C_err"]
        cols = row.split("\t")
        assert cols[0] == "R1"
        assert float(cols[3]) == pytest.approx(30.0)
        assert float(cols[6]) == pytest.approx(1.2)
