# ABOUTME: Every DynamiXs page that takes a series offers the project's series to drag, the same way.
# ABOUTME: Real series (DNAJA1 titration, JA1 hetNOE); one series panel shared by all the pages.

"""
The T1/T2, model-free and Kd pages each carried their own copy of the series list, and
the hetNOE and methyl pages had a drop target with no list to drag from. One panel on
BasePage now serves them all. A dropped series reaches each page as the file that page
reads: the tidy CSV for T1/T2 and Kd, the intensity matrix for methyl T2 and hetNOE.
"""

import json
import os
import sys
import types
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

DATA = Path(__file__).resolve().parents[2] / "data_test"
DNAJA1 = DATA / "data_binding" / "DNAJA1_HSPA8"
HETNOE = DATA / "data_JA1_relaxation" / "run_20260825_1405"

pytestmark = pytest.mark.skipif(not (DNAJA1.exists() and HETNOE.exists()),
                                reason="data_test not present")


@pytest.fixture(scope="module")
def series(tmp_path_factory):
    """{name: series output folder} for a titration and a time series, run for real."""
    from lunaNMR.cli import main
    root = tmp_path_factory.mktemp("series")
    runs = {
        "DNAJA1": ["--spectra", str(DNAJA1), "--peaks", str(DNAJA1 / "DNAJA1.txt"),
                   "--mode", "titration", "--peak-source", "cascade"],
        "hetNOE600": ["--spectra", str(HETNOE / "hetnoe_planes" / "600"),
                      "--peaks", str(HETNOE.parent / "600_WT" / "hetNOE" /
                                     "600_hetNOE_assignment.txt"),
                      "--mode", "time", "--peak-source", "reference"],
    }
    for name, args in runs.items():
        assert main(["series", *args, "--out", str(root / name)]) == 0
    return {name: root / name for name in runs}


@pytest.fixture
def window(series):
    from lunaNMR.utils.project_manager import batch_results_from_dict
    saved = {name: batch_results_from_dict(json.loads((out / "batch_results.json").read_text()))
             for name, out in series.items()}
    return types.SimpleNamespace(saved_series=saved, current_project_path=None)


@pytest.fixture
def dialog(window):
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication(sys.argv)
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    d = DynamiXsDialog(parent=None, main_window=window)
    yield d
    d.deleteLater()


def _listed(page):
    page._populate_series_list()
    widget = page.series_list_widget
    return {widget.item(i).text() for i in range(widget.count())}


@pytest.mark.parametrize("page_name", ["t1t2_page", "integrated_page", "methyl_t2_page",
                                       "hetnoe_page"])
def test_a_page_lists_every_series(dialog, page_name):
    assert _listed(getattr(dialog, page_name)) == {"DNAJA1", "hetNOE600"}


def test_the_kd_page_lists_titrations_only(window):
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication(sys.argv)
    from lunaNMR.gui.dialogs.kd_titration_dialog import KdTitrationDialog
    d = KdTitrationDialog(parent=None, main_window=window)
    try:
        assert _listed(d.kd_page) == {"DNAJA1"}
    finally:
        d.deleteLater()


def test_a_series_in_the_project_resolves_to_its_bundled_csv(dialog, window, tmp_path):
    """A series opened from the project reads the bundle's own copy first: its run
    folder may have moved."""
    bundled = tmp_path / "p.lunaNMR" / "series_results" / "DNAJA1" / "series_analysis_tidy.csv"
    bundled.parent.mkdir(parents=True)
    bundled.write_text("spectrum_name,assignment\n")
    window.current_project_path = str(tmp_path / "p.lunaNMR")
    window.saved_series["DNAJA1"].bundled_folder = str(bundled.parent)   # as loaded
    found = {s["name"]: s["csv_path"] for s in dialog.t1t2_page._get_available_series()}
    assert found["DNAJA1"] == str(bundled)


def test_a_series_dropped_on_hetnoe_brings_its_two_planes(dialog, series):
    page = dialog.hetnoe_page
    page._on_series_dropped("hetnoe", "hetNOE600",
                            str(series["hetNOE600"] / "series_analysis_tidy.csv"))
    assert page.matrix_file == str(series["hetNOE600"] / "peak_intensity_matrix.csv")
    assert page.sat_combo.currentText() == "hetnoe_600_saturated"
    assert page.unsat_combo.currentText() == "hetnoe_600_unsaturated"


def test_a_series_dropped_on_methyl_t2_brings_its_intensity_matrix(dialog, series):
    """The methyl fitter reads the matrix format; a tidy CSV is not one it can parse."""
    page = dialog.methyl_t2_page
    page._on_series_dropped("methyl_t2", "hetNOE600",
                            str(series["hetNOE600"] / "series_analysis_tidy.csv"))
    assert page.input_file == str(series["hetNOE600"] / "peak_intensity_matrix.csv")


def test_a_series_in_the_project_drops_its_bundled_matrix(dialog, window, tmp_path):
    """The bundle's copy first, as for the tidy CSV: the run folder may have moved."""
    bundled = tmp_path / "p.lunaNMR" / "series_results" / "hetNOE600" / "peak_intensity_matrix.csv"
    bundled.parent.mkdir(parents=True)
    bundled.write_text("Assignment,hetnoe_600_saturated,hetnoe_600_unsaturated\n3LysH,1.0,2.0\n")
    window.current_project_path = str(tmp_path / "p.lunaNMR")
    window.saved_series["hetNOE600"].bundled_folder = str(bundled.parent)   # as loaded
    page = dialog.methyl_t2_page
    page._on_series_dropped("methyl_t2", "hetNOE600", "")
    assert page.input_file == str(bundled)
