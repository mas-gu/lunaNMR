# ABOUTME: The model-free outputs must be readable by one consumer across single- and dual-field runs.
# ABOUTME: Covers the per-field column aliases, the plot filename, and the results path the run reports.

import sys
from pathlib import Path

import pandas as pd
import pytest

_MOD = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_MOD))
sys.path.insert(0, str(_MOD / "dynamiXs_density_functions"))

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")


def _single_field_frame():
    """A single-field results frame, shaped like ZZ_multi_density_087.analyze_csv output."""
    return pd.DataFrame({
        "Residue": ["3LysH", "8TyrH"],
        "R1": [1.39, 1.42], "R1_err": [0.05, 0.06],
        "R2": [7.85, 8.10], "R2_err": [0.30, 0.31],
        "hetNOE": [0.78, 0.80], "hetNOE_err": [0.05, 0.05],
        "J0": [2.5e-9, 2.6e-9], "J0_err": [1e-10, 1e-10],
        "JwN": [2.5e-10, 2.6e-10], "JwN_err": [1e-11, 1e-11],
        "JwH_087": [1.1e-11, 1.2e-11], "JwH_087_err": [1e-12, 1e-12],
        "S2": [0.69, 0.71], "S2_err": [0.04, 0.04],
        "te": [21.5, 15.1], "te_err": [6.0, 5.0],
        "Rex": [0.0, 2.1], "Rex_err": [0.3, 0.5],
        "Rex_95CI": ["[0. 0.7]", "[1.1 3.4]"],
        "fit_success": [True, True],
        "chi2": [3.5e-4, 2.7e-4],
    })


class TestFieldAliases:
    """A reader written against the dual-field names must also work on single-field output."""

    def test_per_field_quantities_gain_f1_aliases(self):
        from dynamiXs_integrated.integrated_analysis import add_field1_aliases
        out = add_field1_aliases(_single_field_frame())
        for name in ("R1_f1", "R2_f1", "hetNOE_f1", "J0_f1", "JwN_f1", "JwH_087_f1", "Rex_f1"):
            assert name in out.columns, f"missing alias {name}"

    def test_err_and_ci_suffixes_keep_the_dual_field_spelling(self):
        from dynamiXs_integrated.integrated_analysis import add_field1_aliases
        out = add_field1_aliases(_single_field_frame())
        # dual-field writes Rex_f1_err / Rex_f1_95CI, not Rex_err_f1
        for name in ("R1_f1_err", "R2_f1_err", "hetNOE_f1_err", "J0_f1_err",
                     "JwN_f1_err", "JwH_087_f1_err", "Rex_f1_err", "Rex_f1_95CI"):
            assert name in out.columns, f"missing alias {name}"

    def test_aliases_carry_the_same_values(self):
        from dynamiXs_integrated.integrated_analysis import add_field1_aliases
        df = _single_field_frame()
        out = add_field1_aliases(df)
        assert list(out["Rex_f1"]) == list(df["Rex"])
        assert list(out["R1_f1_err"]) == list(df["R1_err"])

    def test_field_independent_columns_are_not_aliased(self):
        from dynamiXs_integrated.integrated_analysis import add_field1_aliases
        out = add_field1_aliases(_single_field_frame())
        for name in ("S2", "te", "Residue", "fit_success", "chi2"):
            assert f"{name}_f1" not in out.columns, f"{name} is not a per-field quantity"

    def test_original_columns_survive(self):
        """Additive only — anything already reading `Rex` keeps working."""
        from dynamiXs_integrated.integrated_analysis import add_field1_aliases
        df = _single_field_frame()
        out = add_field1_aliases(df)
        for name in df.columns:
            assert name in out.columns


class TestPlotFilename:
    """Single-field plots must land where the caller asks, not on a hardcoded name."""

    def test_plot_results_honours_plot_filename(self, tmp_path):
        from ZZ_multi_density_087 import ReducedSpectralDensityAnalysis
        target = tmp_path / "single600_spectral_density_plots.pdf"
        analyzer = ReducedSpectralDensityAnalysis(spectrometer_frequency=600)
        analyzer.plot_results(_single_field_frame(), save_plots=True,
                              plot_filename=str(target))
        assert target.is_file(), "plot_filename was ignored"

    def test_default_name_is_unchanged_when_no_filename_given(self, tmp_path, monkeypatch):
        from ZZ_multi_density_087 import ReducedSpectralDensityAnalysis
        monkeypatch.chdir(tmp_path)
        analyzer = ReducedSpectralDensityAnalysis(spectrometer_frequency=600)
        analyzer.plot_results(_single_field_frame(), save_plots=True)
        assert (tmp_path / "ZZ_rsdm_multicore_analysis_results_087.pdf").is_file()


class TestReportedResultsPath:
    """The path the run logs as saved must be a file that exists."""

    def test_save_final_results_returns_an_existing_file(self, tmp_path):
        from dynamiXs_integrated.integrated_analysis import IntegratedAnalysisPipeline
        basic = tmp_path / "single600_spectral_density_basic.csv"
        _single_field_frame().to_csv(basic, index=False)
        results = {"output_files": {"basic_csv": str(basic),
                                    "detailed_csv": str(tmp_path / "d.csv"),
                                    "plots_pdf": str(tmp_path / "p.pdf")}}
        returned = IntegratedAnalysisPipeline.__new__(IntegratedAnalysisPipeline)._save_final_results(results)
        assert Path(returned).is_file(), f"run reports a file that was never written: {returned}"
