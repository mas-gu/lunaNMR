# ABOUTME: Tests for the shared-amplitude methyl bi-exp refit helpers in refit.py.
# ABOUTME: Sidecar IO is shared with the mono-exp helpers and tested in test_refit_residue.py.

import json
import sys
from pathlib import Path

import numpy as np
import pytest

_FITTER_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_T1_T2"
sys.path.insert(0, str(_FITTER_DIR))


def _synthetic_methyl_entry(seed=0, with_outlier_at=None):
    """Build a fit_entry/metadata pair matching the new shared-amp schema."""
    rng = np.random.default_rng(seed)
    A, t2_a, t2_b, noise = 2000.0, 120.0, 15.0, 30.0
    x = np.linspace(0.0, 350.0, 24)
    y = 0.5 * A * (np.exp(-x / t2_a) + np.exp(-x / t2_b))
    y = y + rng.normal(0.0, noise, size=x.size)
    if with_outlier_at is not None:
        y[with_outlier_at] += 400.0

    metadata = {
        "experiment_type": "methylT2",
        "time_units": "ms",
        "time_points": x.tolist(),
        "n_bootstrap": 100,
    }
    fit_entry = {
        "residue": "M42",
        "A": 1900.0, "t2_a": 130.0, "t2_b": 14.0, "eta": 9.3,
        "A_err": 50.0, "t2_a_err": 5.0, "t2_b_err": 0.6, "eta_err": 0.5,
        "intensities": y.tolist(),
        "fit_curve": {"time": [0.0, 350.0], "intensity": [2000.0, 0.0]},
    }
    return fit_entry, metadata, dict(A=A, t2_a=t2_a, t2_b=t2_b)


class TestRefitMethyl:

    def test_refit_recovers_clean_params_with_outlier_excluded(self):
        from refit import refit_residue_methyl

        fit_entry, metadata, truth = _synthetic_methyl_entry(seed=11, with_outlier_at=10)
        new_entry = refit_residue_methyl(fit_entry, metadata, excluded_indices=[10])

        assert abs(new_entry["t2_a"] - truth["t2_a"]) / truth["t2_a"] < 0.15
        assert abs(new_entry["t2_b"] - truth["t2_b"]) / truth["t2_b"] < 0.30
        assert new_entry["t2_a"] >= new_entry["t2_b"]

    def test_intensities_preserved(self):
        from refit import refit_residue_methyl

        fit_entry, metadata, _ = _synthetic_methyl_entry(seed=12, with_outlier_at=10)
        original = list(fit_entry["intensities"])

        new_entry = refit_residue_methyl(fit_entry, metadata, excluded_indices=[10])
        assert new_entry["intensities"] == original

    def test_fit_curve_spans_full_time_range(self):
        from refit import refit_residue_methyl

        fit_entry, metadata, _ = _synthetic_methyl_entry(seed=13)
        last_idx = len(metadata["time_points"]) - 1
        new_entry = refit_residue_methyl(fit_entry, metadata, excluded_indices=[last_idx])

        max_t = max(metadata["time_points"])
        assert new_entry["fit_curve"]["time"][0] == pytest.approx(0.0)
        assert new_entry["fit_curve"]["time"][-1] == pytest.approx(max_t * 1.2, rel=1e-6)

    def test_returns_full_methyl_keys(self):
        from refit import refit_residue_methyl

        fit_entry, metadata, _ = _synthetic_methyl_entry(seed=14)
        new_entry = refit_residue_methyl(fit_entry, metadata, excluded_indices=[5])

        for key in ("residue", "A", "t2_a", "t2_b", "eta",
                    "A_err", "t2_a_err", "t2_b_err", "eta_err",
                    "intensities", "fit_curve"):
            assert key in new_entry, f"missing {key}"

    def test_no_legacy_keys(self):
        """A_a, A_b, C must not appear in the new shared-amp output."""
        from refit import refit_residue_methyl

        fit_entry, metadata, _ = _synthetic_methyl_entry(seed=15)
        new_entry = refit_residue_methyl(fit_entry, metadata, excluded_indices=[5])
        for legacy in ("A_a", "A_b", "C", "A_a_err", "A_b_err", "C_err"):
            assert legacy not in new_entry, f"legacy key {legacy} should be gone"

    def test_refit_escapes_bad_warm_start(self):
        """A prior fit anchored at a poor local minimum (e.g. dT2 = 0
        flagged-unidentifiable) used to lock the L-M optimizer in place when
        warm-started. Refit must find a strictly better fit than the prior
        when the data clearly supports one.
        """
        from refit import refit_residue_methyl

        # Data that clearly fits A ≈ 9.4e5, T2 ≈ 75 ms — much better than the
        # prior anchor at A = 7.08e5, T2 = 79.66 (chi² ratio ≈ 18×).
        time_points = [5.0, 25.0, 50.0, 75.0, 200.0]
        intensities = [850000.0, 720000.0, 510000.0, 310000.0, 50000.0]
        fit_entry = {
            "residue": "TRAP",
            "A": 707800.0,
            "T2_avg": 79.66, "dT2": 0.0,
            "t2_a": 79.66, "t2_b": 79.66, "eta": 1.0,
            "A_err": 1.0, "T2_avg_err": 1.0, "dT2_err": float("nan"),
            "t2_a_err": float("nan"), "t2_b_err": float("nan"),
            "eta_err": float("nan"),
            "bi_exp_unidentifiable": True,
            "intensities": intensities,
            "fit_curve": {"time": [0.0, 200.0], "intensity": [707800.0, 0.0]},
        }
        metadata = {
            "experiment_type": "methylT2",
            "time_units": "ms",
            "time_points": time_points,
            "n_bootstrap": 0,
        }
        new = refit_residue_methyl(fit_entry, metadata, excluded_indices=[])

        # The fresh-defaults fit lands at A ≈ 9.4e5 (chi² 18× lower than the
        # warm-start basin). Anything within 15% of that is fine; staying at
        # the prior 7.08e5 means the basin trap survived.
        assert abs(new["A"] - 942000.0) / 942000.0 < 0.15, (
            f"refit stuck in warm-start basin: A={new['A']:.0f}"
        )


class TestUpdateTsvRowMethyl:

    def _make_tsv(self, tmp_path):
        path = tmp_path / "field1_methylT2_fit_results.txt"
        header = ("Residue\tA\tT2a\tT2b\teta\t"
                  "A_err\tT2a_err\tT2b_err\teta_err\tSuccess\n")
        rows = [
            ("M42\t2.000000e+03\t1.200000e+02\t1.500000e+01\t8.000000e+00\t"
             "5.000000e+01\t2.000000e+00\t5.000000e-01\t3.000000e-01\tYes\n"),
            ("M43\t1.800000e+03\t1.000000e+02\t1.800000e+01\t5.555556e+00\t"
             "4.500000e+01\t1.500000e+00\t6.000000e-01\t3.000000e-01\tYes\n"),
        ]
        path.write_text(header + "".join(rows))
        return path

    def _new_fit(self):
        return {
            "A": 2500.0, "t2_a": 150.0, "t2_b": 10.0, "eta": 15.0,
            "A_err": 30.0, "t2_a_err": 1.0, "t2_b_err": 0.3, "eta_err": 0.5,
        }

    def test_replaces_target_row_only(self, tmp_path):
        from refit import update_tsv_row_methyl

        path = self._make_tsv(tmp_path)
        update_tsv_row_methyl(path, "M42", self._new_fit())

        lines = path.read_text().strip().splitlines()
        assert len(lines) == 3
        assert lines[1].split("\t")[0] == "M42"
        assert float(lines[1].split("\t")[2]) == pytest.approx(150.0)  # T2a updated
        assert lines[2].split("\t")[0] == "M43"
        assert float(lines[2].split("\t")[2]) == pytest.approx(100.0)

    def test_preserves_header(self, tmp_path):
        from refit import update_tsv_row_methyl

        path = self._make_tsv(tmp_path)
        original_header = path.read_text().splitlines()[0]
        update_tsv_row_methyl(path, "M42", self._new_fit())
        assert path.read_text().splitlines()[0] == original_header

    def test_unknown_residue_raises(self, tmp_path):
        from refit import update_tsv_row_methyl

        path = self._make_tsv(tmp_path)
        with pytest.raises(KeyError):
            update_tsv_row_methyl(path, "M999", self._new_fit())

    def _make_tsv_with_new_schema(self, tmp_path):
        """TSV with the extended schema (T2_avg, dT2, bi_exp_unidentifiable)."""
        path = tmp_path / "field1_methylT2_fit_results.txt"
        header = (
            "Residue\tA\tT2a\tT2b\tT2_avg\tdT2\teta\t"
            "A_err\tT2a_err\tT2b_err\tT2_avg_err\tdT2_err\teta_err\t"
            "bi_exp_unidentifiable\tSuccess\n"
        )
        rows = [
            ("M42\t2.000000e+03\t1.200000e+02\t1.500000e+01\t6.750000e+01\t"
             "1.050000e+02\t8.000000e+00\t5.000000e+01\t2.000000e+00\t5.000000e-01\t"
             "1.000000e+00\t2.500000e+00\t3.000000e-01\tNo\tYes\n"),
        ]
        path.write_text(header + "".join(rows))
        return path

    def _new_fit_with_extras(self):
        return {
            "A": 2500.0,
            "T2_avg": 80.0, "dT2": 140.0,
            "t2_a": 150.0, "t2_b": 10.0, "eta": 15.0,
            "A_err": 30.0,
            "T2_avg_err": 0.5, "dT2_err": 1.5,
            "t2_a_err": 1.0, "t2_b_err": 0.3, "eta_err": 0.5,
            "bi_exp_unidentifiable": False,
        }

    def test_writes_extended_columns_when_header_has_them(self, tmp_path):
        from refit import update_tsv_row_methyl

        path = self._make_tsv_with_new_schema(tmp_path)
        update_tsv_row_methyl(path, "M42", self._new_fit_with_extras())

        cells = path.read_text().splitlines()[1].split("\t")
        assert cells[0] == "M42"
        assert float(cells[2]) == pytest.approx(150.0)   # T2a
        assert float(cells[4]) == pytest.approx(80.0)    # T2_avg
        assert float(cells[5]) == pytest.approx(140.0)   # dT2
        assert cells[-2] == "No"   # bi_exp_unidentifiable
        assert cells[-1] == "Yes"  # Success

    def test_unidentifiable_flag_writes_yes(self, tmp_path):
        from refit import update_tsv_row_methyl

        path = self._make_tsv_with_new_schema(tmp_path)
        new_fit = self._new_fit_with_extras()
        new_fit["bi_exp_unidentifiable"] = True
        update_tsv_row_methyl(path, "M42", new_fit)

        cells = path.read_text().splitlines()[1].split("\t")
        assert cells[-2] == "Yes"  # bi_exp_unidentifiable
        assert cells[-1] == "Yes"  # Success
