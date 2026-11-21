# ABOUTME: The Fit Viewer shows, refits and hands to Inspect Peak the T1rho fits it loads.
# ABOUTME: A real T1rho fit of a synthetic decay matrix, opened in the real viewer.

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DIR))
sys.path.insert(0, str(_DIR / "dynamiXs_T1_T2"))

DELAYS_MS = [10, 20, 40, 60, 80, 120, 160, 200]
RESIDUES = ("3LysH", "8TyrH", "10ValH")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def t1rho_run(tmp_path):
    """The fit JSON and the results file of a T1rho fit of three residues."""
    import fit_Tx_NMRRE
    rng = np.random.default_rng(3)
    rows = []
    for number, name in enumerate(RESIDUES, start=1):
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
    assert result["n_fitted"] == len(RESIDUES)
    return Path(result["json_file"]), Path(result["results_file"])


@pytest.fixture
def viewer(app, t1rho_run, monkeypatch):
    import visualization.fit_viewer as fit_viewer
    for box in ("show_info", "show_warning", "show_error"):
        monkeypatch.setattr(fit_viewer, box, lambda *a, **k: None)
    viewer = fit_viewer.FitViewer(json_files=[str(t1rho_run[0])],
                                  source_series_map={"field1_t1rho": "spinlock_series"})
    yield viewer
    viewer.close()


def _row(results_file, residue):
    for line in Path(results_file).read_text().splitlines()[1:]:
        if line.split("\t", 1)[0] == residue:
            return line
    raise KeyError(residue)


def test_the_t1rho_fit_of_the_selected_residue_is_drawn(viewer):
    assert viewer.t1rho_checkbox.isChecked()
    viewer.t1_checkbox.setChecked(False)
    viewer.t2_checkbox.setChecked(False)
    viewer.residue_list.setCurrentRow(1)

    titles = [ax.get_title() for ax in viewer.axes]
    assert titles == ["Residue 8TyrH - T1rho"]
    (ax,) = viewer.axes
    assert ax.lines, "the fitted curve is not drawn"
    assert "T1rho (600.0 MHz)" in " ".join(t.get_text() for t in ax.texts)


def test_a_viewer_without_t1rho_fits_draws_no_t1rho_panel(app, t1rho_run, monkeypatch, tmp_path):
    import visualization.fit_viewer as fit_viewer
    for box in ("show_info", "show_warning", "show_error"):
        monkeypatch.setattr(fit_viewer, box, lambda *a, **k: None)
    t1_json = tmp_path / "json" / "run_field1_T1_fit_data.json"
    t1_json.write_text(t1rho_run[0].read_text())
    viewer = fit_viewer.FitViewer(json_files=[str(t1_json)])
    try:
        viewer.residue_list.setCurrentRow(0)
        assert not viewer.t1rho_checkbox.isChecked()
        assert [ax.get_title() for ax in viewer.axes][0] == "Residue 3LysH - T1"
        assert not any("T1rho" in ax.get_title() for ax in viewer.axes)
    finally:
        viewer.close()


def test_refit_rewrites_the_t1rho_results_row(viewer, t1rho_run):
    results_file = t1rho_run[1]
    viewer.residue_list.setCurrentRow(0)
    before = _row(results_file, "3LysH")
    viewer.exclusions.setdefault("field1_T1rho", {})["3LysH"] = [0, 1]

    viewer._on_refit_clicked()

    after = _row(results_file, "3LysH")
    assert after != before
    fit = next(f for f in viewer.data["field1_T1rho"]["fits"] if f["residue"] == "3LysH")
    assert float(after.split("\t")[2]) == pytest.approx(fit["t2"], rel=1e-6)


def test_inspect_peak_hands_over_the_t1rho_fits(viewer):
    requests = []
    viewer.inspect_peak_requested.connect(lambda *args: requests.append(args))
    viewer.residue_list.setCurrentRow(2)

    viewer._on_inspect_peak_clicked()

    ((residue, series_name, fits),) = requests
    assert (residue, series_name) == ("10ValH", "spinlock_series")
    assert set(fits) == set(RESIDUES)
    assert fits["10ValH"]["experiment_type"] == "T1rho"
    assert fits["10ValH"]["t_value"] == pytest.approx(90.0, rel=0.1)
