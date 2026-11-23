# ABOUTME: DynamiXs and Kd order residues by one rule: by residue number, a name without a
# ABOUTME: number last. Real viewers, tables, merges and fits, each fed names with and without one.

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DIR))
sys.path.insert(0, str(_DIR / "dynamiXs_T1_T2"))
sys.path.insert(0, str(_DIR / "dynamiXs_Kd"))

# Listed out of order; the two side-chain names carry no residue number.
NAMES = ("TrpSide", "10ValH", "3LysH", "ArgSide")
IN_ORDER = ["3LysH", "10ValH", "ArgSide", "TrpSide"]
DELAYS_MS = [10, 20, 40, 60, 80, 120, 160, 200]


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


def _silenced(module, monkeypatch):
    for box in ("show_info", "show_warning", "show_error"):
        monkeypatch.setattr(module, box, lambda *a, **k: None)
    return module


def _rates(values):
    return {name: {"value": value, "error": 0.02 * value} for name, value in values.items()}


class TestOneFunction:
    def test_kd_orders_with_the_function_the_relaxation_modules_use(self):
        import kd_models
        import residue_order
        assert kd_models.residue_sort_key is residue_order.residue_sort_key

    def test_names_without_a_number_come_last_alphabetically(self):
        from residue_order import residue_sort_key
        assert sorted(NAMES, key=residue_sort_key) == IN_ORDER

    @pytest.mark.parametrize("folder, module", [("dynamiXs_Kd", "kd_models"),
                                                ("dynamiXs_integrated", "data_converters"),
                                                ("dynamiXs_integrated", "data_validation")])
    def test_a_module_imported_by_bare_name_finds_the_function(self, folder, module):
        code = ("import sys; "
                "sys.path[:] = [p for p in sys.path if 'dynamiXs' not in p] "
                f"+ [{str(_DIR / folder)!r}]; import {module}")
        run = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                             cwd=str(_DIR.parent))
        assert run.returncode == 0, run.stderr


class TestAxisPositions:
    def test_a_numbered_name_sits_at_its_number(self):
        from residue_order import residue_axis_positions
        assert residue_axis_positions(["10ValH", "3LysH", 45.0, "A17"]) == [10, 3, 45, 17]

    def test_names_without_a_number_follow_the_highest_number_in_order(self):
        from residue_order import residue_axis_positions
        assert residue_axis_positions(NAMES) == [12, 10, 3, 11]

    def test_nothing_numbered_still_gives_every_name_a_position(self):
        from residue_order import residue_axis_positions
        assert residue_axis_positions(["TrpSide", "ArgSide"]) == [2, 1]


@pytest.fixture
def t1rho_fit_json(tmp_path):
    """The fit JSON of a real T1rho fit of NAMES."""
    import fit_Tx_NMRRE
    rng = np.random.default_rng(3)
    rows = []
    for number, name in enumerate(NAMES, start=1):
        row = {"Assignment": name}
        for delay in DELAYS_MS:
            row[str(delay)] = 1.0e6 * number * np.exp(-delay / (60.0 + 10 * number)) \
                * (1 + 0.01 * rng.standard_normal())
        rows.append(row)
    matrix = tmp_path / "t1rho_matrix.csv"
    pd.DataFrame(rows).to_csv(matrix, index=False)
    result = fit_Tx_NMRRE.run_analysis_with_params({
        "input_csv_file": str(matrix),
        "output_prefix": str(tmp_path / "run_field1_T1rho"),
        "results_txt_file": str(tmp_path / "run_field1_T1rho_fit_results.txt"),
        "experiment_type": "T1rho", "time_units": "ms", "signal_units": "Intensity",
        "initial_A": 1.0e6, "initial_t2": 80.0, "error_method": "analytical",
        "n_bootstrap": 100, "n_plots_per_figure": 20,
        "json_folder": str(tmp_path / "json"), "field_name": "run_field1",
        "field_freq": 600.0})
    assert result["n_fitted"] == len(NAMES)
    return result["json_file"]


class TestFitViewer:
    @pytest.fixture
    def viewer(self, app, t1rho_fit_json, monkeypatch):
        import visualization.fit_viewer as module
        viewer = _silenced(module, monkeypatch).FitViewer(json_files=[str(t1rho_fit_json)])
        yield viewer
        viewer.close()

    def test_the_residue_list_ends_with_the_names_without_a_number(self, viewer):
        listed = [viewer.residue_list.item(i).text() for i in range(viewer.residue_list.count())]
        assert listed == IN_ORDER

    def test_the_plots_of_several_residues_are_drawn_in_that_order(self, viewer):
        viewer.t1_checkbox.setChecked(False)
        viewer.t2_checkbox.setChecked(False)
        viewer._generate_plots(list(NAMES))
        assert [ax.get_title() for ax in viewer.axes] == [
            f"Residue {name} - T1rho" for name in IN_ORDER]


def test_the_methyl_viewer_lists_names_without_a_number_last(app, tmp_path, monkeypatch):
    from fit_methyl_T2 import run_methyl_t2_analysis_with_params
    import visualization.methyl_t2_fit_viewer as module
    delays = np.array([2.0, 4.0, 8.0, 12.0, 18.0, 25.0, 35.0, 50.0, 75.0, 100.0,
                       150.0, 200.0, 250.0, 300.0])
    rows = []
    for seed, name in enumerate(NAMES):
        rng = np.random.default_rng(seed)
        y = 0.5 * 2.0e6 * (np.exp(-delays / 80.0) + np.exp(-delays / 15.0))
        rows.append({"Assignment": name,
                     **{f"{d:g}": v for d, v in zip(delays, y + rng.normal(0.0, 4.0e3, y.size))}})
    pd.DataFrame(rows).to_csv(tmp_path / "matrix.csv", index=False)
    result = run_methyl_t2_analysis_with_params({
        "input_csv_file": str(tmp_path / "matrix.csv"), "output_prefix": str(tmp_path / "p"),
        "results_txt_file": str(tmp_path / "p_fit_results.txt"),
        "json_folder": str(tmp_path / "json"), "field_name": "p_field1", "field_freq": 600.0,
        "time_units": "ms"})
    viewer = _silenced(module, monkeypatch).MethylT2FitViewer(json_file=result["json_file"])
    try:
        assert viewer.available_residues == IN_ORDER
    finally:
        viewer.close()


class TestResultsViewers:
    """A name without a number is plotted after the highest residue number."""

    def test_the_model_free_viewer_plots_it_after_the_numbered_residues(
            self, app, tmp_path, monkeypatch):
        import visualization.results_viewer as module
        pd.DataFrame({"Residue": list(NAMES), "R1": [1.1, 1.2, 1.3, 1.4], "R1_err": 0.02,
                      "S2": [0.80, 0.82, 0.84, 0.86], "S2_err": 0.01}).to_csv(
            tmp_path / "results.csv", index=False)
        viewer = _silenced(module, monkeypatch).ResultsViewer(
            results_file=str(tmp_path / "results.csv"))
        try:
            viewer.error_checkbox.setChecked(False)
            viewer.parameter_buttons["R1"].click()
            (ax,) = viewer.figure.axes
            (points,) = ax.lines
            assert dict(zip(points.get_xdata(), points.get_ydata())) == {
                12: 1.1, 10: 1.2, 3: 1.3, 11: 1.4}
        finally:
            viewer.close()

    def test_the_t1t2_viewer_plots_it_after_the_numbered_residues(
            self, app, tmp_path, monkeypatch):
        import visualization.t1t2_results_viewer as module
        lines = ["Residue\tA\tT1\tC\tA_err\tT1_err\tC_err\tSuccess\tWindowRatio"]
        values = dict(zip(NAMES, (700.0, 710.0, 720.0, 730.0)))
        for name, value in values.items():
            lines.append(f"{name}\t1.000000e+06\t{value:.6e}\t0.000000e+00\t1.000000e+03\t"
                         f"{value * 0.02:.6e}\t0.000000e+00\tYes\t3.000")
        (tmp_path / "T1_fit_results.txt").write_text("\n".join(lines) + "\n")
        viewer = _silenced(module, monkeypatch).T1T2ResultsViewer(session_results={
            "field1_t1": {"results_file": str(tmp_path / "T1_fit_results.txt"),
                          "frequency_mhz": 600.0, "experiment_type": "T1", "time_units": "ms"}})
        try:
            (ax,) = viewer.figure.axes
            plotted = {}
            for bars in ax.containers:
                plotted.update(zip(bars.lines[0].get_xdata(), bars.lines[0].get_ydata()))
            assert plotted == pytest.approx({12: 700.0, 10: 710.0, 3: 720.0, 11: 730.0})
        finally:
            viewer.close()


class TestHetNOE:
    def test_a_name_without_a_number_is_drawn_after_the_highest_residue_number(self):
        from matplotlib.figure import Figure
        from dynamiXs_integrated.data_converters import draw_hetnoe_vs_residue
        values = {"TrpSide": 0.5, "70ValH": 0.8, "3LysH": 0.6, "ArgSide": 0.4}
        ax = Figure().add_subplot()
        draw_hetnoe_vs_residue(ax, _rates(values))
        points = ax.containers[0].lines[0]
        assert dict(zip(points.get_xdata(), points.get_ydata())) == {
            3: 0.6, 70: 0.8, 71: 0.4, 72: 0.5}

    def test_the_table_order_puts_names_without_a_number_last(self):
        from dynamiXs_integrated.data_converters import in_sequence_order
        ordered = in_sequence_order(_rates(dict.fromkeys(NAMES, 0.7)))
        assert [name for name, _ in ordered] == IN_ORDER


class TestModelFreeInputs:
    def test_the_density_input_file_ends_with_the_names_without_a_number(self, tmp_path):
        from dynamiXs_integrated.data_converters import create_spectral_density_input_csv
        r1 = _rates(dict.fromkeys(NAMES, 1.4))
        r2 = _rates(dict.fromkeys(NAMES, 11.0))
        noe = _rates(dict.fromkeys(NAMES, 0.8))
        create_spectral_density_input_csv(r1, r2, noe, str(tmp_path / "input.csv"))
        assert list(pd.read_csv(tmp_path / "input.csv")["Residue"]) == IN_ORDER

    def test_the_residues_dropped_at_the_merge_end_with_the_names_without_a_number(self):
        from dynamiXs_integrated.data_validation import merge_datasets_by_residue
        everyone = _rates(dict.fromkeys(NAMES + ("8TyrH",), 1.4))
        merged, excluded = merge_datasets_by_residue(everyone, everyone, _rates({"8TyrH": 0.8}))
        assert list(merged) == ["8TyrH"]
        assert excluded == IN_ORDER


class TestDualFieldDensity:
    @pytest.fixture
    def analyzer(self):
        from dynamiXs_density_functions.ZZ_multi_2fields_density import (
            DualFieldSpectralDensityAnalysis)
        return DualFieldSpectralDensityAnalysis(field1_freq=600.0, field2_freq=700.0,
                                                rNH=1.015e-10, csaN=-172.0e-6,
                                                omega_h_factor=1.0)

    @pytest.fixture
    def results(self, analyzer, tmp_path):
        """A real dual-field fit of NAMES, the two input tables listing them out of order."""
        for field, (r1, r2, noe) in {"f1": (1.60, 12.5, 0.78), "f2": (1.35, 14.6, 0.80)}.items():
            pd.DataFrame({"Residue": list(NAMES),
                          "R1": [r1 + 0.01 * i for i in range(len(NAMES))], "R1err": 0.03,
                          "R2": [r2 + 0.1 * i for i in range(len(NAMES))], "R2err": 0.2,
                          "hetNOE": noe, "hetNOEerr": 0.03}).to_csv(
                tmp_path / f"{field}.csv", index=False)
        return analyzer.analyze_dual_field_csv(
            str(tmp_path / "f1.csv"), str(tmp_path / "f2.csv"), use_monte_carlo_errors=False,
            use_multiprocessing=False)

    def test_the_results_rows_end_with_the_names_without_a_number(self, results):
        assert list(results["Residue"]) == IN_ORDER
        assert list(results["R1_f1"]) == pytest.approx([1.62, 1.61, 1.63, 1.60])

    def test_the_plot_draws_them_after_the_highest_residue_number(
            self, analyzer, results, tmp_path, monkeypatch):
        import dynamiXs_density_functions.ZZ_multi_2fields_density as module
        saved = []
        monkeypatch.setattr(module.Figure, "savefig",
                            lambda figure, *a, **k: saved.append(figure))
        analyzer.plot_dual_field_results(results, plot_filename=str(tmp_path / "plots.pdf"))

        (figure,) = saved
        points = figure.axes[0].containers[0].lines[0]
        assert dict(zip(points.get_xdata(), points.get_ydata())) == pytest.approx(
            {3: 1.62, 10: 1.61, 11: 1.63, 12: 1.60})


def test_the_kd_viewer_records_excluded_residues_in_that_order(app, tmp_path):
    from kd_models import csp_model
    from kd_fit import run_kd_analysis_with_params
    from visualization.kd_titration_fit_viewer import open_kd_titration_viewer
    points = [0.0, 10.0, 25.0, 60.0, 150.0, 300.0]
    rows = []
    for number, name in enumerate(NAMES):
        shift = csp_model(np.array(points), 0.2 + 0.05 * number, 15.0, 50.0)
        rows += [(str(p), name, 8.0 + d, 120.0, 1000.0 - 2 * p, 2000.0)
                 for p, d in zip(points, shift)]
    pd.DataFrame(rows, columns=["spectrum_name", "assignment", "ppm_x", "ppm_y", "height",
                                "volume"]).to_csv(tmp_path / "series_analysis_tidy.csv",
                                                  index=False)
    fit_json = run_kd_analysis_with_params({
        "input_csv_file": str(tmp_path / "series_analysis_tidy.csv"),
        "output_dir": str(tmp_path), "output_prefix": "kd", "concentrations": points,
        "protein_conc": 50.0, "alpha": 0.14, "observables": ["csp"],
        "n_bootstrap": 0})["json_file"]
    viewer = open_kd_titration_viewer(parent=None, json_file=fit_json)
    try:
        viewer.edit_mode = True
        for name in NAMES:
            artist = object()
            viewer._pick_registry[id(artist)] = (name, "RESIDUE", None)
            viewer._on_pick(type("Pick", (), {"artist": artist})())
        assert json.loads(Path(fit_json).read_text())["excluded_residues"]["csp"] == IN_ORDER
    finally:
        viewer.deleteLater()
