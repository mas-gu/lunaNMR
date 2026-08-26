# ABOUTME: Tests for headless refit-with-outlier-rejection logic in dynamiXs.
# ABOUTME: Covers sidecar (.outliers.json) IO, single-residue refit, and surgical JSON/TSV updates.

import json
import sys
from pathlib import Path

import numpy as np
import pytest

# Make dynamiXs_T1_T2/ importable (no __init__.py in that directory).
_FITTER_DIR = Path(__file__).resolve().parent.parent / "dynamiXs_T1_T2"
sys.path.insert(0, str(_FITTER_DIR))


# ---------- helpers ----------

def _synthetic_fit_entry(seed=0, with_outlier_at=None):
    """Build a fit_entry/metadata pair for one residue, optionally with one outlier."""
    rng = np.random.default_rng(seed)
    # C_true = 0: the fitter fixes the baseline at zero, so a synthetic offset would
    # only measure that convention's bias rather than what refit is for.
    A_true, t2_true, C_true, noise = 1000.0, 50.0, 0.0, 2.0
    x = np.linspace(0.0, 250.0, 24)
    y = A_true * np.exp(-x / t2_true) + C_true + rng.normal(0.0, noise, size=x.size)
    if with_outlier_at is not None:
        y[with_outlier_at] += 500.0  # large positive outlier

    metadata = {
        "experiment_type": "T2",
        "time_units": "ms",
        "time_points": x.tolist(),
        "n_bootstrap": 100,
    }
    fit_entry = {
        "residue": "142",
        "A": 950.0,
        "t2": 55.0,
        "C": 0.0,
        "A_err": 5.0,
        "t2_err": 1.0,
        "C_err": 1.0,
        "intensities": y.tolist(),
        "fit_curve": {"time": [0.0, 250.0], "intensity": [1000.0, 0.0]},
    }
    return fit_entry, metadata, (A_true, t2_true, C_true)


# ---------- sidecar IO ----------

class TestSidecarIO:

    def test_missing_sidecar_returns_empty(self, tmp_path):
        from refit import load_outliers

        json_path = tmp_path / "field1_T2_fit_data.json"
        json_path.write_text("{}")  # main JSON exists but no sidecar

        assert load_outliers(json_path) == {}

    def test_save_then_load_roundtrip(self, tmp_path):
        from refit import load_outliers, save_outliers

        json_path = tmp_path / "field1_T2_fit_data.json"
        json_path.write_text("{}")
        exclusions = {"142": [2, 5], "143": [], "200": [0]}

        save_outliers(json_path, exclusions)
        loaded = load_outliers(json_path)

        assert loaded == exclusions

    def test_sidecar_path_is_alongside_json(self, tmp_path):
        from refit import sidecar_path_for

        json_path = tmp_path / "field1_T2_fit_data.json"
        sc = sidecar_path_for(json_path)

        assert sc.parent == json_path.parent
        assert sc.name == "field1_T2_fit_data.outliers.json"

    def test_save_overwrites_previous(self, tmp_path):
        from refit import load_outliers, save_outliers

        json_path = tmp_path / "field1_T2_fit_data.json"
        json_path.write_text("{}")

        save_outliers(json_path, {"142": [0, 1]})
        save_outliers(json_path, {"142": [2]})

        assert load_outliers(json_path) == {"142": [2]}


# ---------- refit ----------

class TestRefit:

    def test_refit_with_one_outlier_recovers_clean_params(self):
        """With a single bad point excluded, T2 must come back close to truth."""
        from refit import refit_residue

        fit_entry, metadata, (A_true, t2_true, C_true) = _synthetic_fit_entry(
            seed=1, with_outlier_at=10
        )

        # Confirm the outlier *does* hurt the prior fit (sanity that the test is real).
        # (We pretend the prior fit was poor; refit should still converge cleanly.)
        new_entry = refit_residue(fit_entry, metadata, excluded_indices=[10])

        assert abs(new_entry["t2"] - t2_true) < 1.5, (
            f"T2 not recovered: got {new_entry['t2']}, expected ~{t2_true}"
        )
        assert new_entry["C"] == 0.0        # baseline is fixed, never warm-started
        assert abs(new_entry["A"] - A_true) < 50.0

    def test_refit_preserves_full_intensities(self):
        """The raw intensities array must NOT be mutated by exclusion."""
        from refit import refit_residue

        fit_entry, metadata, _ = _synthetic_fit_entry(seed=2, with_outlier_at=10)
        original_intensities = list(fit_entry["intensities"])

        new_entry = refit_residue(fit_entry, metadata, excluded_indices=[10])

        assert new_entry["intensities"] == original_intensities, (
            "intensities field must remain the raw data; exclusions live in the sidecar"
        )

    def test_refit_curve_spans_full_time_range(self):
        """Even if the last data point is excluded, the dense fit_curve must still
        span 0..max(original time_points)*1.2 so plots look the same."""
        from refit import refit_residue

        fit_entry, metadata, _ = _synthetic_fit_entry(seed=3)
        last_idx = len(metadata["time_points"]) - 1

        new_entry = refit_residue(fit_entry, metadata, excluded_indices=[last_idx])

        max_t = max(metadata["time_points"])
        fit_t = new_entry["fit_curve"]["time"]
        assert min(fit_t) == pytest.approx(0.0)
        assert max(fit_t) == pytest.approx(max_t * 1.2, rel=1e-6)

    def test_refit_empty_exclusion_close_to_prior(self):
        """An empty exclusion list re-fits the original data; result should be near prior."""
        from refit import refit_residue

        fit_entry, metadata, (A_true, t2_true, C_true) = _synthetic_fit_entry(seed=4)

        new_entry = refit_residue(fit_entry, metadata, excluded_indices=[])

        # No outlier, so result should match the truth within noise
        assert abs(new_entry["t2"] - t2_true) < 1.0
        assert abs(new_entry["C"] - C_true) < 4.0

    def test_refit_returns_fit_entry_with_required_keys(self):
        from refit import refit_residue

        fit_entry, metadata, _ = _synthetic_fit_entry(seed=5)
        new_entry = refit_residue(fit_entry, metadata, excluded_indices=[5])

        for key in ("residue", "A", "t2", "C", "A_err", "t2_err", "C_err",
                    "intensities", "fit_curve"):
            assert key in new_entry, f"missing key: {key}"
        assert new_entry["residue"] == fit_entry["residue"]

    def test_refit_escapes_bad_warm_start_anchor(self):
        """The mono-exp refit warm-starts from the prior fit's params. If the
        prior was a bad fit (e.g. amplitude or t2 wildly off), the optimizer
        must still find the correct minimum on clean data — analogous to the
        bug fixed for refit_residue_methyl.
        """
        from refit import refit_residue

        # Generate clean data with a clear-best mono-exp solution.
        rng = np.random.default_rng(101)
        A_true, t2_true, C_true = 1000.0, 50.0, 0.0   # baseline is fixed at zero
        x = np.linspace(0.0, 250.0, 24)
        y = A_true * np.exp(-x / t2_true) + C_true + rng.normal(0.0, 2.0, size=x.size)

        # Prior fit anchored at wildly wrong values (bad warm-start).
        bad_fit_entry = {
            "residue": "BAD_ANCHOR",
            "A": 100000.0,   # 100× too large
            "t2": 5.0,       # 10× too small
            "C": -200.0,     # wrong sign (must not be warm-started back in)
            "A_err": 1.0, "t2_err": 1.0, "C_err": 1.0,
            "intensities": y.tolist(),
            "fit_curve": {"time": [0.0, 250.0], "intensity": [100000.0, 0.0]},
        }
        metadata = {
            "experiment_type": "T2",
            "time_units": "ms",
            "time_points": x.tolist(),
            "n_bootstrap": 0,
        }
        new = refit_residue(bad_fit_entry, metadata, excluded_indices=[])

        assert abs(new["A"] - A_true) / A_true < 0.10, (
            f"refit stuck near bad warm-start: A={new['A']:.0f}"
        )
        assert abs(new["t2"] - t2_true) / t2_true < 0.10, (
            f"refit stuck near bad warm-start: t2={new['t2']:.2f}"
        )


# ---------- JSON surgical update ----------

class TestUpdateJsonFitEntry:

    def _make_json(self, tmp_path):
        data = {
            "metadata": {"experiment_type": "T2", "time_units": "ms",
                         "time_points": [0.0, 50.0, 100.0], "n_bootstrap": 1000,
                         "n_residues": 2, "field_freq": 600.0,
                         "signal_units": "Intensity"},
            "fits": [
                {"residue": "142", "A": 1.0, "t2": 50.0, "C": 0.0,
                 "A_err": 0.1, "t2_err": 1.0, "C_err": 0.1,
                 "intensities": [1.0, 0.4, 0.1],
                 "fit_curve": {"time": [0.0, 100.0], "intensity": [1.0, 0.1]}},
                {"residue": "143", "A": 2.0, "t2": 60.0, "C": 0.5,
                 "A_err": 0.2, "t2_err": 2.0, "C_err": 0.2,
                 "intensities": [2.0, 0.8, 0.6],
                 "fit_curve": {"time": [0.0, 100.0], "intensity": [2.0, 0.6]}},
            ],
        }
        path = tmp_path / "field1_T2_fit_data.json"
        path.write_text(json.dumps(data))
        return path

    def test_replaces_target_residue_only(self, tmp_path):
        from refit import update_json_fit_entry

        path = self._make_json(tmp_path)
        new_entry = {"residue": "142", "A": 9.9, "t2": 99.0, "C": 9.0,
                     "A_err": 0.5, "t2_err": 0.5, "C_err": 0.5,
                     "intensities": [1.0, 0.4, 0.1],
                     "fit_curve": {"time": [0.0, 100.0], "intensity": [9.9, 9.0]}}

        update_json_fit_entry(path, "142", new_entry)

        out = json.loads(path.read_text())
        residues = {f["residue"]: f for f in out["fits"]}
        assert residues["142"]["t2"] == 99.0  # replaced
        assert residues["143"]["t2"] == 60.0  # untouched

    def test_preserves_metadata(self, tmp_path):
        from refit import update_json_fit_entry

        path = self._make_json(tmp_path)
        original_meta = json.loads(path.read_text())["metadata"]

        new_entry = {"residue": "142", "A": 9.9, "t2": 99.0, "C": 9.0,
                     "A_err": 0.5, "t2_err": 0.5, "C_err": 0.5,
                     "intensities": [1.0, 0.4, 0.1],
                     "fit_curve": {"time": [0.0, 100.0], "intensity": [9.9, 9.0]}}
        update_json_fit_entry(path, "142", new_entry)

        assert json.loads(path.read_text())["metadata"] == original_meta

    def test_unknown_residue_raises(self, tmp_path):
        from refit import update_json_fit_entry

        path = self._make_json(tmp_path)
        with pytest.raises(KeyError):
            update_json_fit_entry(path, "999", {"residue": "999"})


# ---------- TSV surgical update ----------

class TestUpdateTsvRow:

    def _make_tsv(self, tmp_path, experiment_type="T2"):
        path = tmp_path / "fit_results.txt"
        header = f"Residue\tA\t{experiment_type}\tC\tA_err\t{experiment_type}_err\tC_err\n"
        rows = [
            "142\t1.000000e+00\t5.000000e+01\t0.000000e+00\t1.000000e-01\t1.000000e+00\t1.000000e-01\n",
            "143\t2.000000e+00\t6.000000e+01\t5.000000e-01\t2.000000e-01\t2.000000e+00\t2.000000e-01\n",
        ]
        path.write_text(header + "".join(rows))
        return path

    def test_replaces_target_row_only(self, tmp_path):
        from refit import update_tsv_row

        path = self._make_tsv(tmp_path)
        new_fit = {"A": 9.9, "t2": 99.0, "C": 9.0,
                   "A_err": 0.5, "t2_err": 0.5, "C_err": 0.5}

        update_tsv_row(path, "142", new_fit, experiment_type="T2")

        lines = path.read_text().strip().splitlines()
        assert len(lines) == 3  # header + 2 rows
        assert lines[0].startswith("Residue\t")
        # row 142 replaced
        assert lines[1].split("\t")[0] == "142"
        assert float(lines[1].split("\t")[2]) == pytest.approx(99.0)
        # row 143 untouched
        assert lines[2].split("\t")[0] == "143"
        assert float(lines[2].split("\t")[2]) == pytest.approx(60.0)

    def test_preserves_header(self, tmp_path):
        from refit import update_tsv_row

        path = self._make_tsv(tmp_path)
        original_header = path.read_text().splitlines()[0]

        update_tsv_row(path, "142",
                       {"A": 9.9, "t2": 99.0, "C": 9.0,
                        "A_err": 0.5, "t2_err": 0.5, "C_err": 0.5},
                       experiment_type="T2")

        assert path.read_text().splitlines()[0] == original_header

    def test_unknown_residue_raises(self, tmp_path):
        from refit import update_tsv_row

        path = self._make_tsv(tmp_path)
        with pytest.raises(KeyError):
            update_tsv_row(path, "999",
                           {"A": 1.0, "t2": 1.0, "C": 0.0,
                            "A_err": 0.1, "t2_err": 0.1, "C_err": 0.1},
                           experiment_type="T2")
