# ABOUTME: Tests for the shared-amplitude bi-exponential methyl T2 fitter (fit_methyl_T2.py).
# ABOUTME: Model: I(t) = 0.5*A*[exp(-t/T2a) + exp(-t/T2b)]. 3 free params; T2a >= T2b post-fit.

import json
import sys
from pathlib import Path

import numpy as np
import pytest

_FITTER_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_T1_T2"
sys.path.insert(0, str(_FITTER_DIR))


def _synthetic_methyl(seed=0, A=2000.0, t2_a=120.0, t2_b=15.0,
                      noise=20.0, n=24, t_max=350.0):
    """Build noisy synthetic data following the paper equation."""
    rng = np.random.default_rng(seed)
    x = np.linspace(0.0, t_max, n)
    y = 0.5 * A * (np.exp(-x / t2_a) + np.exp(-x / t2_b))
    y = y + rng.normal(0.0, noise, size=x.size)
    return x, y, dict(A=A, t2_a=t2_a, t2_b=t2_b)


# ---------- equation ----------

class TestBiexpDecay:

    def test_returns_half_a_times_sum_of_components(self):
        from fit_methyl_T2 import biexp_decay

        x = np.array([0.0, 10.0, 50.0])
        y = biexp_decay(x, A=200.0, t2_a=80.0, t2_b=10.0)
        expected = 0.5 * 200.0 * (np.exp(-x / 80.0) + np.exp(-x / 10.0))
        np.testing.assert_allclose(y, expected, rtol=1e-12)

    def test_at_t_zero_equals_A(self):
        """At t=0, both exponentials = 1, so I(0) = 0.5*A*(1+1) = A."""
        from fit_methyl_T2 import biexp_decay

        for A in (1.0, 100.0, 1e6):
            assert biexp_decay(0.0, A, 100.0, 20.0) == pytest.approx(A)


# ---------- fit recovery ----------

class TestFitRecovery:

    def test_recovers_known_params(self):
        from fit_methyl_T2 import fit_single_residue_methyl

        x, y, truth = _synthetic_methyl(seed=1)
        res = fit_single_residue_methyl(x, y, "M1")

        assert abs(res["t2_a"] - truth["t2_a"]) / truth["t2_a"] < 0.10
        assert abs(res["t2_b"] - truth["t2_b"]) / truth["t2_b"] < 0.20
        assert abs(res["A"] - truth["A"]) / truth["A"] < 0.10

    def test_post_fit_ordering_t2a_geq_t2b(self):
        """Even with reversed initial guesses, post-fit t2_a >= t2_b must hold."""
        from fit_methyl_T2 import fit_single_residue_methyl

        x, y, _ = _synthetic_methyl(seed=2)
        res = fit_single_residue_methyl(
            x, y, "M2", initial_t2_a=20.0, initial_t2_b=100.0,
        )
        assert res["t2_a"] >= res["t2_b"]

    def test_eta_equals_ratio(self):
        from fit_methyl_T2 import fit_single_residue_methyl

        x, y, _ = _synthetic_methyl(seed=3)
        res = fit_single_residue_methyl(x, y, "M3")
        assert res["eta"] == pytest.approx(res["t2_a"] / res["t2_b"], rel=1e-9)
        assert res["eta"] >= 1.0

    def test_eta_error_is_finite_and_positive(self):
        """eta_err should be finite and positive for a well-determined fit.

        The exact value depends on the (T2_avg, dT2) covariance and we no longer
        use the naive sqrt-sum-of-squares formula (which ignores Cov(t2_a, t2_b)
        and produces nonsense near degeneracy). Just check sanity here — math
        equivalence is verified by TestReparametrization.test_central_values_match.
        """
        from fit_methyl_T2 import fit_single_residue_methyl

        x, y, _ = _synthetic_methyl(seed=4)
        res = fit_single_residue_methyl(x, y, "M4")
        if np.isfinite(res["eta_err"]):
            assert res["eta_err"] > 0
            assert res["eta_err"] < 100 * res["eta"]  # not absurdly inflated

    def test_data_driven_default_for_amplitude(self):
        """When initial_A is None, default to span(y) without crashing."""
        from fit_methyl_T2 import fit_single_residue_methyl

        x, y, _ = _synthetic_methyl(seed=5)
        res = fit_single_residue_methyl(x, y, "M5")
        assert "A" in res

    def test_result_dict_has_required_keys(self):
        from fit_methyl_T2 import fit_single_residue_methyl

        x, y, _ = _synthetic_methyl(seed=6)
        res = fit_single_residue_methyl(x, y, "M6")
        required = {
            "residue", "A", "t2_a", "t2_b", "eta",
            "A_err", "t2_a_err", "t2_b_err", "eta_err",
            "T2_avg", "dT2", "T2_avg_err", "dT2_err",
            "bi_exp_unidentifiable",
            "x", "y", "result",
        }
        missing = required - set(res.keys())
        assert not missing, f"missing keys: {missing}"

    def test_no_legacy_keys(self):
        """A_a/A_b/C must not appear (paper form has shared A and no offset)."""
        from fit_methyl_T2 import fit_single_residue_methyl

        x, y, _ = _synthetic_methyl(seed=7)
        res = fit_single_residue_methyl(x, y, "M7")
        for legacy in ("A_a", "A_b", "C", "A_a_err", "A_b_err", "C_err"):
            assert legacy not in res, f"legacy key {legacy} should be gone"


# ---------- bootstrap ordering ----------

class TestBootstrapOrdering:

    def test_bootstrap_errs_far_smaller_than_t2_separation(self):
        """Bootstrap stderr must reflect ordering, not be a bimodal-distribution std."""
        from fit_methyl_T2 import fit_single_residue_methyl

        x, y, _ = _synthetic_methyl(seed=99, A=20000, t2_a=120, t2_b=15, noise=50)
        res = fit_single_residue_methyl(
            x, y, "M_BOOT", n_bootstrap=200, error_method="bootstrap"
        )

        separation = abs(res["t2_a"] - res["t2_b"])
        assert res["t2_a_err"] < 0.2 * separation, (
            f"t2_a_err={res['t2_a_err']:.2f} suggests bimodal distribution "
            f"(separation {separation:.2f})"
        )
        assert res["t2_b_err"] < 0.2 * separation, (
            f"t2_b_err={res['t2_b_err']:.2f} suggests bimodal distribution "
            f"(separation {separation:.2f})"
        )


# ---------- bounds ----------

class TestParameterBounds:

    def test_t2_does_not_run_away_on_mono_exp_data(self):
        from fit_methyl_T2 import fit_single_residue_methyl

        rng = np.random.default_rng(0)
        x = np.linspace(0.0, 200.0, 16)
        y = 1000.0 * np.exp(-x / 50.0) + rng.normal(0, 2.0, size=x.size)
        res = fit_single_residue_methyl(x, y, "MONO")

        assert res["t2_a"] < 5.0 * x.max()
        assert res["t2_b"] < 5.0 * x.max()
        assert res["t2_a"] > 0.0
        assert res["t2_b"] > 0.0


class TestReparametrization:
    """Tests that the fit reparameterizes to (A, T2_avg, dT2) with t2_a, t2_b
    derived. Equivalent fit to the original parameterization on well-determined
    data; unidentifiability becomes explicit on degenerate data.
    """

    def test_central_values_match_truth(self):
        """Reparam fit must recover the same (t2_a, t2_b) as before on clean data."""
        from fit_methyl_T2 import fit_single_residue_methyl

        x, y, truth = _synthetic_methyl(seed=42)
        res = fit_single_residue_methyl(x, y, "REPARAM")

        assert abs(res["t2_a"] - truth["t2_a"]) / truth["t2_a"] < 0.10
        assert abs(res["t2_b"] - truth["t2_b"]) / truth["t2_b"] < 0.20

    def test_T2_avg_and_dT2_consistent_with_t2a_t2b(self):
        from fit_methyl_T2 import fit_single_residue_methyl

        x, y, _ = _synthetic_methyl(seed=43)
        res = fit_single_residue_methyl(x, y, "REPARAM_2")

        assert res["T2_avg"] == pytest.approx((res["t2_a"] + res["t2_b"]) / 2.0, rel=1e-9)
        assert res["dT2"] == pytest.approx(res["t2_a"] - res["t2_b"], rel=1e-9)
        assert res["dT2"] >= 0  # slow-first ordering enforced by parameterization

    def test_dT2_lower_bound_at_zero(self):
        """Mono-exp synthetic data must drive dT2 to its 0 lower bound."""
        from fit_methyl_T2 import fit_single_residue_methyl

        rng = np.random.default_rng(7)
        x = np.linspace(0.0, 350.0, 10)
        T2 = 33.0
        y = 350000.0 * np.exp(-x / T2) + rng.normal(0.0, 1500.0, size=x.size)

        res = fit_single_residue_methyl(x, y, "MONO_DT")
        assert res["dT2"] < 1.0  # essentially zero relative to T2_avg ≈ 33


class TestMonoExpCollapse:
    """Behavior on mono-exponential data fed to the bi-exp fitter
    (the 172_cg1_D failure mode: optimizer collapses to t2_a = t2_b).
    """

    @staticmethod
    def _mono_exp_data(seed=0, T2=33.0, A=350000.0, noise=1500.0,
                      n=10, t_max=350.0):
        rng = np.random.default_rng(seed)
        x = np.linspace(0.0, t_max, n)
        y = A * np.exp(-x / T2) + rng.normal(0.0, noise, size=x.size)
        return x, y, T2, A

    def test_mono_exp_data_flagged_unidentifiable(self):
        from fit_methyl_T2 import fit_single_residue_methyl

        x, y, _, _ = self._mono_exp_data(seed=11)
        res = fit_single_residue_methyl(x, y, "MONO_FLAG")

        assert res["bi_exp_unidentifiable"] is True

    def test_mono_exp_T2_avg_well_determined(self):
        """T2_avg should match the true mono-exp T2 to within 10%, with a
        small T2_avg_err — the well-determined direction in the reparam basis.
        """
        from fit_methyl_T2 import fit_single_residue_methyl

        x, y, T2_true, _ = self._mono_exp_data(seed=12)
        res = fit_single_residue_methyl(x, y, "MONO_AVG")

        assert abs(res["T2_avg"] - T2_true) / T2_true < 0.10
        if np.isfinite(res["T2_avg_err"]):
            assert res["T2_avg_err"] / res["T2_avg"] < 0.20

    def test_well_separated_data_not_flagged(self):
        from fit_methyl_T2 import fit_single_residue_methyl

        x, y, _ = _synthetic_methyl(seed=14, A=20000, t2_a=120, t2_b=15, noise=50)
        res = fit_single_residue_methyl(x, y, "WELL_SEPARATED")
        assert res["bi_exp_unidentifiable"] is False

    def test_mono_exp_collapse_suppresses_bound_stderrs(self):
        """When dT2 sits at its 0 lower bound (true mono-exp data), lmfit
        returns a finite-but-meaningless huge stderr from the singular Hessian.
        Suppress to NaN so the viewer doesn't display "0 ± 1,000,000 ms".
        Note: this only applies when dT2 is *at* the bound, not whenever the
        unidentifiable flag fires (small-but-nonzero dT2 has meaningful stderr).
        """
        from fit_methyl_T2 import fit_single_residue_methyl
        import math

        # Noise-free mono-exp data: optimizer converges to dT2 = 0 exactly.
        x = np.linspace(0.0, 350.0, 10)
        T2 = 33.0
        y = 350000.0 * np.exp(-x / T2)
        res = fit_single_residue_methyl(x, y, "MONO_NAN")

        assert res["bi_exp_unidentifiable"] is True
        assert math.isnan(res["dT2_err"]), \
            f"dT2_err should be NaN at bound, got {res['dT2_err']}"
        assert math.isnan(res["t2_a_err"])
        assert math.isnan(res["t2_b_err"])
        assert math.isnan(res["eta_err"])


class TestT2bLowerBoundInvariant:
    """t2_b must never be returned negative — the model internally clips
    when ΔT2 > 2·T2_avg, but the result dict and JSON must stay physical
    so the viewer doesn't plot a growing exponential.
    """

    def test_bound_active_refit_as_mono_exp_matches_visible_amplitude(self):
        """When t2_b hits the lower bound, the bi-exp's A is twice the visible
        amplitude (the 'fast component' carries half but decays before the first
        measurement — a fitting artifact). The fitter should refit as a mono-exp
        in that regime and report A as the actually-visible amplitude.

        Reproduces residue 594_ce_D from real data, where the original fit
        reported A = 2.32e6 with the highest data point at ~1.4e6.
        """
        from fit_methyl_T2 import fit_single_residue_methyl

        x = np.array([2.0, 10.0, 25.0, 37.0, 50.0, 75.0, 100.0, 200.0])
        y = np.array([1.40e6, 0.97e6, 0.98e6, 0.81e6, 0.86e6,
                      0.83e6, 0.75e6, 0.45e6])

        res = fit_single_residue_methyl(x, y, "594_REPRO")
        assert res["bi_exp_unidentifiable"] is True
        # Reported A must be near the visible amplitude (~1.4e6), not 2× it.
        assert res["A"] < 1.7e6, (
            f"A reported as {res['A']:.2e} — should be near visible "
            f"amplitude (~1.4e6), not the bi-exp's 2× artifact"
        )
        # In mono-exp regime t2_a == t2_b == T2_avg.
        assert res["t2_a"] == pytest.approx(res["t2_b"], rel=1e-9)
        assert res["t2_a"] == pytest.approx(res["T2_avg"], rel=1e-9)
        assert res["dT2"] == pytest.approx(0.0, abs=1e-9)
        # T2 should be in the slow-component range (~150 ms here).
        assert res["T2_avg"] > 50.0

    def test_t2_b_never_negative_on_pure_noise(self):
        """Pure-noise data (zero-signal residue) used to come back with
        ΔT2 > 2·T2_avg → derived t2_b < 0, which the viewer plotted as a
        growing exponential. Ensure t2_b stays at the physical lower bound.
        """
        from fit_methyl_T2 import fit_single_residue_methyl, T2_LOWER_BOUND_MS

        rng = np.random.default_rng(0)
        x = np.linspace(0.0, 350.0, 10)
        y = rng.normal(0.0, 200.0, size=x.size)

        res = fit_single_residue_methyl(x, y, "NOISE_232")
        assert res["t2_b"] >= T2_LOWER_BOUND_MS
        assert res["t2_a"] >= res["t2_b"]
        assert res["dT2"] >= 0
        assert res["dT2"] == pytest.approx(res["t2_a"] - res["t2_b"], rel=1e-9)
        # Bound-active case is by definition unidentifiable.
        if res["t2_b"] == T2_LOWER_BOUND_MS:
            assert res["bi_exp_unidentifiable"] is True


# ---------- TSV output ----------

class TestParseDelayColumn:
    """Direct unit-table for parse_delay_column. Each token convention should
    have explicit coverage so a regex regression can't sneak through.
    """

    @pytest.mark.parametrize("col, expected", [
        ("300", 300.0),
        ("50.5", 50.5),
        ("300_2", 300.0),                 # duplicate-measurement suffix
        ("50.5_2", 50.5),
        ("003_T2_ADDA_3ms", 3.0),         # descriptive name + unit
        ("003_T2_ADDA_300ms", 300.0),
        ("600_T1_0o3", 0.3),              # filesystem-safe 'o' for '.'
        ("600_T1_0o0", 0.0),
        ("5s", 5000.0),                   # second → ms conversion
        ("500us", 0.5),                   # microsecond → ms conversion
        ("3MS", 3.0),                     # case-insensitive unit
        ("", None),
        ("Notes", None),
        ("not_a_number", None),
    ])
    def test_parse_table(self, col, expected):
        from fit_methyl_T2 import parse_delay_column

        result = parse_delay_column(col)
        if expected is None:
            assert result is None, f"{col!r} → {result}, expected None"
        else:
            assert result == pytest.approx(expected, rel=1e-9), (
                f"{col!r} → {result}, expected {expected}"
            )


class TestFitCurveInvariants:
    """Whole-pipeline invariants on the dense fit curve written to JSON.
    A growing exponential or non-monotonic curve is always a bug regardless
    of which corner of parameter space the optimizer landed in.
    """

    def test_fit_curve_monotonically_non_increasing_on_real_decay(self, tmp_path):
        """For any clean decay-shaped data, fit_curve['intensity'] must be
        non-increasing — methyl T2 decays, never grows."""
        from fit_methyl_T2 import fit_single_residue_methyl, save_fit_data_json_methyl
        import json

        rng = np.random.default_rng(0)
        x = np.linspace(0.0, 350.0, 12)
        # Exercise several seeds and noise levels to broaden coverage.
        for seed, noise in [(0, 50), (1, 200), (2, 1000)]:
            rng = np.random.default_rng(seed)
            y = 0.5 * 20000 * (np.exp(-x / 120) + np.exp(-x / 15))
            y = y + rng.normal(0, noise, size=x.size)
            res = fit_single_residue_methyl(x, y, f"R{seed}")
            res["success"] = True

            out = tmp_path / f"f_{seed}_methylT2_fit_data.json"
            save_fit_data_json_methyl(
                [res], str(out), time_units="ms", signal_units="Intensity",
                n_bootstrap=0, field_freq=600.0,
            )
            payload = json.loads(out.read_text())
            curve = payload["fits"][0]["fit_curve"]["intensity"]
            for i in range(len(curve) - 1):
                assert curve[i] >= curve[i + 1] - 1e-6, (
                    f"fit_curve grows at index {i} (seed={seed}): "
                    f"{curve[i]:.3e} < {curve[i+1]:.3e}"
                )

    def test_fit_curve_monotonic_on_pure_noise(self):
        """Even on pathological inputs the curve must not grow."""
        from fit_methyl_T2 import fit_single_residue_methyl, save_fit_data_json_methyl

        rng = np.random.default_rng(99)
        x = np.linspace(0.0, 350.0, 10)
        y = rng.normal(0.0, 200.0, size=x.size)
        res = fit_single_residue_methyl(x, y, "NOISE_INV")
        # Build the dense fit_curve directly without writing to disk.
        import numpy as _np
        max_t = float(_np.max(x))
        t_dense = _np.linspace(0.0, max_t * 1.2, 100)
        intensity = 0.5 * res["A"] * (
            _np.exp(-t_dense / res["t2_a"]) + _np.exp(-t_dense / res["t2_b"])
        )
        diffs = _np.diff(intensity)
        assert _np.all(diffs <= 1e-6), (
            f"fit curve grows somewhere; max forward diff = {diffs.max():.3e}"
        )


class TestRunAnalysisWithDuplicatedDelayHeaders:
    """run_methyl_t2_analysis_with_params must accept delay columns like
    '5_1', '5_2', '300_2' (duplicate-measurement suffixes) — same set that
    parse_delay_column already supports.
    """

    def _write_csv(self, tmp_path, header_row, data_rows):
        path = tmp_path / "peak_volume_matrix.csv"
        lines = [",".join(header_row)]
        for row in data_rows:
            lines.append(",".join(str(c) for c in row))
        path.write_text("\n".join(lines) + "\n")
        return path

    def test_runs_with_descriptive_spectrum_headers(self, tmp_path):
        """The lunaNMR series-QC volume matrix uses descriptive headers like
        '003_T2_ADDA_3ms' / '003_T2_ADDA_300ms'. None of these parse as a bare
        float; the loader used to fall through and crash on the Assignment
        column with "could not convert string to float: '27_cg1_D'".
        """
        from fit_methyl_T2 import run_methyl_t2_analysis_with_params

        delays = [3.0, 9.0, 25.0, 50.0, 100.0, 200.0, 300.0]
        header = ["Peak_Number", "Assignment", "Reference_X", "Reference_Y"] + [
            f"003_T2_ADDA_{int(d)}ms" for d in delays
        ]
        rng = np.random.default_rng(0)

        def signal(t):
            return 0.5 * 100000.0 * (np.exp(-t / 80.0) + np.exp(-t / 15.0))
        rows = []
        for k, name in enumerate(["27_cg1_D", "44_cg1_D"]):
            y = [signal(t) + rng.normal(0, 50) for t in delays]
            rows.append([k + 1, name, 0.0, 0.0] + [f"{v:.2f}" for v in y])

        csv_path = self._write_csv(tmp_path, header, rows)
        out_dir = tmp_path / "out"
        out_dir.mkdir()

        result = run_methyl_t2_analysis_with_params({
            "input_csv_file": str(csv_path),
            "output_prefix": str(out_dir / "field1_methylT2"),
            "results_txt_file": str(out_dir / "field1_methylT2_fit_results.txt"),
            "json_folder": str(out_dir / "json"),
            "field_name": "field1",
            "field_freq": 600.0,
            "n_bootstrap": 50,
        })
        assert result["n_fitted"] == 2
        assert result["n_total"] == 2

    def test_runs_with_o_decimal_headers(self, tmp_path):
        """The 2DHSQC example uses 'o' instead of '.' (filesystem-safe) in
        delay headers like '600_T1_0o3' → 0.3 ms.
        """
        from fit_methyl_T2 import run_methyl_t2_analysis_with_params

        delays_str = ["0o0", "0o1", "0o3", "0o5", "1o0", "2o0", "5o0", "10"]
        delays_ms = [0.0, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0]
        header = ["Peak_Number", "Assignment", "Reference_X", "Reference_Y"] + [
            f"600_T1_{s}" for s in delays_str
        ]
        rng = np.random.default_rng(1)

        def signal(t):
            return 0.5 * 100000.0 * (np.exp(-t / 5.0) + np.exp(-t / 1.0))
        rows = []
        for k, name in enumerate(["27_cg1_D"]):
            y = [signal(t) + rng.normal(0, 50) for t in delays_ms]
            rows.append([k + 1, name, 0.0, 0.0] + [f"{v:.2f}" for v in y])

        csv_path = self._write_csv(tmp_path, header, rows)
        out_dir = tmp_path / "out"
        out_dir.mkdir()

        result = run_methyl_t2_analysis_with_params({
            "input_csv_file": str(csv_path),
            "output_prefix": str(out_dir / "f1_methylT2"),
            "results_txt_file": str(out_dir / "f1_methylT2_fit_results.txt"),
            "json_folder": str(out_dir / "json"),
            "field_name": "f1",
            "field_freq": 600.0,
            "n_bootstrap": 50,
        })
        assert result["n_fitted"] == 1

    def test_clear_error_when_no_delay_columns(self, tmp_path):
        """If no header looks like a delay (numeric or numeric_N), raise a
        clear error instead of crashing later on a confusing astype()."""
        from fit_methyl_T2 import run_methyl_t2_analysis_with_params

        header = ["Peak_Number", "Assignment", "Reference_X", "Reference_Y",
                  "Notes", "Comment"]
        rows = [[1, "27_cg1_D", 0.0, 0.0, "x", "y"]]
        csv_path = self._write_csv(tmp_path, header, rows)
        out_dir = tmp_path / "out"
        out_dir.mkdir()

        with pytest.raises(ValueError, match="delay column"):
            run_methyl_t2_analysis_with_params({
                "input_csv_file": str(csv_path),
                "output_prefix": str(out_dir / "field1_methylT2"),
                "results_txt_file": str(out_dir / "field1_methylT2_fit_results.txt"),
                "json_folder": str(out_dir / "json"),
                "field_name": "field1",
                "field_freq": 600.0,
            })


class TestSaveResultsTSV:

    def _result(self, residue="M1"):
        return {
            "residue": residue,
            "A": 2000.0,
            "T2_avg": 67.5, "dT2": 105.0,
            "t2_a": 120.0, "t2_b": 15.0, "eta": 8.0,
            "A_err": 5.0,
            "T2_avg_err": 1.0, "dT2_err": 2.5,
            "t2_a_err": 2.0, "t2_b_err": 0.5, "eta_err": 0.3,
            "bi_exp_unidentifiable": False,
            "success": True,
        }

    def test_header_columns(self, tmp_path):
        from fit_methyl_T2 import save_results_methyl

        out = tmp_path / "methylT2.txt"
        save_results_methyl([self._result()], str(out))

        header = out.read_text().splitlines()[0].split("\t")
        assert header == [
            "Residue", "A", "T2a", "T2b", "T2_avg", "dT2", "eta",
            "A_err", "T2a_err", "T2b_err", "T2_avg_err", "dT2_err", "eta_err",
            "bi_exp_unidentifiable", "Success",
        ]

    def test_row_values(self, tmp_path):
        from fit_methyl_T2 import save_results_methyl

        out = tmp_path / "methylT2.txt"
        save_results_methyl([self._result()], str(out))

        row = out.read_text().splitlines()[1].split("\t")
        assert row[0] == "M1"
        assert float(row[1]) == pytest.approx(2000.0)   # A
        assert float(row[2]) == pytest.approx(120.0)    # T2a
        assert float(row[4]) == pytest.approx(67.5)     # T2_avg
        assert float(row[5]) == pytest.approx(105.0)    # dT2
        assert float(row[6]) == pytest.approx(8.0)      # eta
        assert row[-2] == "No"   # bi_exp_unidentifiable
        assert row[-1] == "Yes"  # Success


# ---------- JSON output ----------

class TestSaveJSONMethyl:

    def _result(self, residue="M1"):
        return {
            "residue": residue,
            "A": 2000.0, "t2_a": 120.0, "t2_b": 15.0, "eta": 8.0,
            "A_err": 5.0, "t2_a_err": 2.0, "t2_b_err": 0.5, "eta_err": 0.3,
            "x": np.array([0.0, 50.0, 100.0]),
            "y": np.array([2000.0, 800.0, 320.0]),
        }

    def test_metadata_experiment_type_methylT2(self, tmp_path):
        from fit_methyl_T2 import save_fit_data_json_methyl

        out = tmp_path / "field1_methylT2_fit_data.json"
        save_fit_data_json_methyl(
            [self._result()], str(out),
            time_units="ms", signal_units="Intensity",
            n_bootstrap=100, field_freq=600.0,
        )
        data = json.loads(out.read_text())
        assert data["metadata"]["experiment_type"] == "methylT2"

    def test_per_residue_includes_all_params(self, tmp_path):
        from fit_methyl_T2 import save_fit_data_json_methyl

        out = tmp_path / "field1_methylT2_fit_data.json"
        save_fit_data_json_methyl(
            [self._result()], str(out),
            time_units="ms", signal_units="Intensity",
            n_bootstrap=100, field_freq=600.0,
        )
        fit = json.loads(out.read_text())["fits"][0]
        for key in ("A", "t2_a", "t2_b", "eta",
                    "A_err", "t2_a_err", "t2_b_err", "eta_err"):
            assert key in fit, f"missing {key}"
        # Legacy keys should not appear
        for legacy in ("A_a", "A_b", "C"):
            assert legacy not in fit, f"legacy key {legacy} should be gone"

    def test_fit_curve_uses_full_dense_range(self, tmp_path):
        from fit_methyl_T2 import save_fit_data_json_methyl

        out = tmp_path / "field1_methylT2_fit_data.json"
        save_fit_data_json_methyl(
            [self._result()], str(out),
            time_units="ms", signal_units="Intensity",
            n_bootstrap=100, field_freq=600.0,
        )
        fit = json.loads(out.read_text())["fits"][0]
        assert fit["fit_curve"]["time"][0] == pytest.approx(0.0)
        assert fit["fit_curve"]["time"][-1] == pytest.approx(100.0 * 1.2, rel=1e-6)
