# ABOUTME: The DynamiXs hetNOE page: I_sat / I_unsat per residue, from a series or two files.
# ABOUTME: Real JA1 600 MHz hetNOE data; checked against `dynamixs hetnoe` and a plain ratio.

"""
hetNOE had no page: it existed only as a model-free input. The page computes it the way
the CLI does (calculate_hetnoe_from_intensities), from a hetNOE series' intensity matrix
- both planes as columns - or from the CLI's two intensity files, and shows the table
and the hetNOE-vs-residue plot. It is a named DynamiXs analysis, saved with a project.
"""

import os
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

RUN = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
       / "run_20260825_1405")
MATRIX = RUN / "series" / "600_hetNOE" / "peak_intensity_matrix.csv"
SAT = RUN / "results" / "hetnoe" / "f600_noe_sat.csv"
UNSAT = RUN / "results" / "hetnoe" / "f600_noe_unsat.csv"

pytestmark = pytest.mark.skipif(not MATRIX.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def dialog(app):
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    d = DynamiXsDialog(parent=None, main_window=types.SimpleNamespace())
    yield d
    d.deleteLater()


def _run_matrix(page, out, prefix="f600"):
    page.mode_combo.setCurrentIndex(0)
    page._set_matrix(str(MATRIX))
    page.output_dir = str(out)
    page.prefix_edit.setText(prefix)
    page._run()
    return Path(out) / f"{prefix}_hetnoe.csv"


def test_the_plot_leaves_the_backend_alone(app, tmp_path):
    """The page draws in-process; switching the GUI's backend to Agg would break every
    other plot window."""
    import matplotlib
    from lunaNMR.cli import _add_modules_path
    _add_modules_path("dynamiXs_v2o0")
    matplotlib.use("QtAgg")
    from dynamiXs_integrated.data_converters import plot_hetnoe_vs_residue
    plot_hetnoe_vs_residue({"3LysH": {"value": 0.8, "error": 0.05}}, str(tmp_path / "p.pdf"))
    assert (tmp_path / "p.pdf").stat().st_size > 0
    assert matplotlib.get_backend().lower() == "qtagg"


def test_a_series_matrix_offers_its_saturated_and_unsaturated_planes(dialog):
    page = dialog.hetnoe_page
    page._set_matrix(str(MATRIX))
    assert page.sat_combo.currentText() == "hetnoe_600_saturated"
    assert page.unsat_combo.currentText() == "hetnoe_600_unsaturated"


def test_the_hetnoe_is_the_plane_ratio(dialog, tmp_path):
    out = _run_matrix(dialog.hetnoe_page, tmp_path)
    got = pd.read_csv(out).set_index("Residue")["hetNOE"]
    m = pd.read_csv(MATRIX).set_index("Assignment")
    ratio = (m["hetnoe_600_saturated"] / m["hetnoe_600_unsaturated"]).dropna()
    ratio = ratio[m.loc[ratio.index, "hetnoe_600_unsaturated"] > 0]
    assert set(got.index) == set(ratio.index.astype(str))
    np.testing.assert_allclose(got.loc[ratio.index.astype(str)].to_numpy(), ratio.to_numpy())


def test_two_files_give_what_the_cli_gives(dialog, tmp_path):
    from lunaNMR.cli import main
    assert main(["dynamixs", "hetnoe", "--sat", str(SAT), "--unsat", str(UNSAT),
                 "--out", str(tmp_path / "cli"), "--prefix", "f600"]) == 0
    page = dialog.hetnoe_page
    page.mode_combo.setCurrentIndex(1)
    page._set_sat_file(str(SAT))
    page._set_unsat_file(str(UNSAT))
    page.output_dir = str(tmp_path / "gui")
    page.prefix_edit.setText("f600")
    page._run()
    gui = pd.read_csv(tmp_path / "gui" / "f600_hetnoe.csv").set_index("Residue")
    cli = pd.read_csv(tmp_path / "cli" / "f600_hetnoe.csv").set_index("Residue")
    pd.testing.assert_frame_equal(gui.sort_index(), cli.sort_index())


def test_the_results_show_as_a_table_and_a_plot(dialog, tmp_path):
    page = dialog.hetnoe_page
    out = _run_matrix(page, tmp_path)
    assert page.results_table.rowCount() == len(pd.read_csv(out))
    assert page.figure.axes and page.figure.axes[0].lines or page.figure.axes[0].collections


def test_a_saved_state_reopens_its_results(app, dialog, tmp_path):
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    out = _run_matrix(dialog.hetnoe_page, tmp_path)
    state = dialog.hetnoe_page.get_session_state()
    other = DynamiXsDialog(parent=None, main_window=types.SimpleNamespace())
    try:
        other.hetnoe_page.restore_session_state(state)
        assert other.hetnoe_page.results_table.rowCount() == len(pd.read_csv(out))
        assert other.hetnoe_page.sat_combo.currentText() == "hetnoe_600_saturated"
    finally:
        other.deleteLater()


def test_a_project_saves_and_reopens_the_analysis(app, tmp_path):
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    from lunaNMR.utils.project_manager import ProjectManager
    window = types.SimpleNamespace(integrator=types.SimpleNamespace(peak_list=None))
    first = DynamiXsDialog(parent=None, main_window=window)
    try:
        first.show_hetnoe()
        _run_matrix(first.hetnoe_page, tmp_path / "run")
        name = first.ensure_current_saved()
        assert name and name.endswith("hetNOE")
        assert ProjectManager(window).save_project(tmp_path / "p.lunaNMR")[0]
    finally:
        first.deleteLater()
    loaded = types.SimpleNamespace(integrator=types.SimpleNamespace(peak_list=None))
    assert ProjectManager(loaded).load_project(tmp_path / "p.lunaNMR")
    second = DynamiXsDialog(parent=None, main_window=loaded)
    try:
        second.open_analysis(loaded.dynamixs_analyses[name], name)
        page = second.hetnoe_page
        assert second.stack.currentWidget() is page
        bundled = Path(page.last_results_file).resolve()
        assert bundled.is_relative_to((tmp_path / "p.lunaNMR").resolve()), bundled
        assert page.results_table.rowCount() > 0
    finally:
        second.deleteLater()
