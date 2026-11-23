# ABOUTME: The model-free results viewer returns to its blank message when nothing is
# ABOUTME: selected, and offers every quantity a results file holds, J(wH) and tc included.

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DIR))


def _results_csv(path, jwh_column):
    """A single-field model-free results table whose J at the proton frequency is in
    `jwh_column` ('JwH' or 'JwH_087')."""
    n = 12
    rng = np.random.default_rng(0)
    frame = pd.DataFrame({
        "Residue": [f"{i}AlaH" for i in range(3, 3 + n)],
        "R1": 1.4 + 0.05 * rng.standard_normal(n), "R1_err": 0.02,
        "R2": 11.0 + 0.5 * rng.standard_normal(n), "R2_err": 0.2,
        "hetNOE": 0.8 + 0.02 * rng.standard_normal(n), "hetNOE_err": 0.02,
        "J0": 3.0 + 0.1 * rng.standard_normal(n), "J0_err": 0.05,
        "JwN": 0.3 + 0.01 * rng.standard_normal(n), "JwN_err": 0.01,
        jwh_column: 0.005 + 0.0005 * rng.standard_normal(n), f"{jwh_column}_err": 0.0005,
        "S2": 0.85 + 0.02 * rng.standard_normal(n), "S2_err": 0.01,
        "tc": 8.2, "te": 40.0 + 5 * rng.standard_normal(n), "te_err": 5.0,
        "Rex": 0.5, "Rex_err": 0.1,
    })
    frame.to_csv(path, index=False)
    return str(path)


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


def _dual_field_results_csv(path):
    """A dual-field model-free results table: R1, R2 and Rex per field, one S2 and te."""
    n = 12
    rng = np.random.default_rng(1)
    frame = pd.DataFrame({
        "Residue": [f"{i}AlaH" for i in range(3, 3 + n)],
        "R1_f1": 1.6 + 0.05 * rng.standard_normal(n), "R1_f1_err": 0.02,
        "R1_f2": 1.35 + 0.05 * rng.standard_normal(n), "R1_f2_err": 0.02,
        "R2_f1": 10.5 + 0.5 * rng.standard_normal(n), "R2_f1_err": 0.2,
        "R2_f2": 11.8 + 0.5 * rng.standard_normal(n), "R2_f2_err": 0.2,
        "S2": 0.85 + 0.02 * rng.standard_normal(n), "S2_err": 0.01,
        "tc": 8.2, "te": 40.0 + 5 * rng.standard_normal(n), "te_err": 5.0,
        "Rex_f1": 0.5, "Rex_f1_err": 0.1, "Rex_f2": 0.7, "Rex_f2_err": 0.1,
    })
    frame.to_csv(path, index=False)
    return str(path)


def _viewer_module(monkeypatch):
    """The results viewer module with its message boxes silenced."""
    import visualization.results_viewer as module
    for box in ("show_info", "show_warning", "show_error"):
        monkeypatch.setattr(module, box, lambda *a, **k: None)
    return module


@pytest.fixture
def dual_field_viewer(app, tmp_path, monkeypatch):
    viewer = _viewer_module(monkeypatch).ResultsViewer(
        results_file=_dual_field_results_csv(tmp_path / "dual.csv"))
    yield viewer
    viewer.close()


@pytest.fixture
def open_viewer(app, tmp_path, monkeypatch):
    module = _viewer_module(monkeypatch)
    made = []

    def make(jwh_column):
        made.append(module.ResultsViewer(
            results_file=_results_csv(tmp_path / f"{jwh_column}.csv", jwh_column)))
        return made[-1]
    yield make
    for viewer in made:
        viewer.close()


def _texts(viewer):
    return [t.get_text() for ax in viewer.figure.axes for t in ax.texts]


def test_deselecting_the_last_parameter_shows_the_blank_message(open_viewer):
    viewer = open_viewer("JwH_087")
    viewer.parameter_buttons["R1"].click()
    assert [ax.get_ylabel() for ax in viewer.figure.axes] == ["R₁ (s⁻¹)"]

    viewer.parameter_buttons["R1"].click()

    assert _texts(viewer) == ["Select a parameter above to display results"]
    (ax,) = viewer.figure.axes
    assert not ax.lines and not ax.collections


@pytest.mark.parametrize("control", ["error_checkbox", "field2_radio", "overlay_radio"])
def test_a_display_option_changed_with_nothing_selected_keeps_the_blank_message(
        dual_field_viewer, control):
    viewer = dual_field_viewer
    viewer.parameter_buttons["R1"].click()
    assert [ax.get_ylabel() for ax in viewer.figure.axes] == ["R₁ (s⁻¹)"]
    viewer.parameter_buttons["R1"].click()

    getattr(viewer, control).click()

    assert not any(button.isChecked() for button in viewer.parameter_buttons.values())
    assert _texts(viewer) == ["Select a parameter above to display results"]
    (ax,) = viewer.figure.axes
    assert not ax.lines and not ax.collections


def test_tc_has_a_button_and_plots(open_viewer):
    viewer = open_viewer("JwH_087")
    assert viewer.parameter_buttons["tc"].isEnabled()

    viewer.parameter_buttons["tc"].click()

    assert [ax.get_ylabel() for ax in viewer.figure.axes] == ["τc (ns)"]


@pytest.mark.parametrize("jwh_column, offered, not_offered", [
    ("JwH", "J(wH)", "J(0.87wH)"),
    ("JwH_087", "J(0.87wH)", "J(wH)"),
])
def test_the_proton_frequency_density_the_file_holds_is_offered(open_viewer, jwh_column,
                                                                offered, not_offered):
    viewer = open_viewer(jwh_column)
    assert viewer.parameter_buttons[offered].isEnabled()
    assert not viewer.parameter_buttons[not_offered].isEnabled()

    viewer.parameter_buttons[offered].click()

    (ax,) = viewer.figure.axes
    assert ax.get_ylabel() == viewer._get_ylabel(offered)
    assert ax.lines or ax.containers


def test_show_all_draws_the_quantities_the_file_holds(open_viewer):
    viewer = open_viewer("JwH")

    viewer.show_all_btn.click()

    titles = [ax.get_title() for ax in viewer.figure.axes]
    assert titles == ["R1", "R2", "hetNOE", "J(0)", "J(wN)", "J(wH)", "S2", "te", "Rex", "tc"]
    assert not any("not found" in text for text in _texts(viewer))
