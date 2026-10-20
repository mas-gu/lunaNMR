# ABOUTME: The DynamiXs main menu lists the project's saved runs and opens or deletes them, as the Kd
# ABOUTME: page does its fits. Real JA1 data -> CLI t1t2 + hetnoe -> project build -> loaded dialog.

import os
import shutil
import sys
import types
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

RUN = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
       / "run_20260825_1405")
NOE = RUN / "results" / "hetnoe"
MATRIX = RUN / "series" / "600_T1" / "peak_intensity_matrix.csv"

pytestmark = pytest.mark.skipif(not RUN.exists(), reason="data_test not present")

HETNOE = "f600_hetNOE"


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture(scope="module")
def bundle(tmp_path_factory):
    from lunaNMR.cli import main
    run = tmp_path_factory.mktemp("run")
    assert main(["dynamixs", "t1t2", "--input", str(MATRIX), "--out", str(run / "fits" / "600_T1"),
                 "--exp", "T1", "--field-freq", "600", "--prefix", "600_T1"]) == 0
    assert main(["dynamixs", "hetnoe", "--sat", str(NOE / "f600_noe_sat.csv"),
                 "--unsat", str(NOE / "f600_noe_unsat.csv"),
                 "--out", str(run / "hetnoe" / "f600"), "--prefix", "f600"]) == 0
    out = run.parent / "saved_runs.lunaNMR"
    assert main(["project", "build", str(run), "--out", str(out)]) == 0
    shutil.rmtree(run)
    return out


def _loaded(bundle):
    from lunaNMR.utils.project_manager import ProjectManager
    window = types.SimpleNamespace(integrator=types.SimpleNamespace(peak_list=None))
    assert ProjectManager(window).load_project(bundle)
    return window


@pytest.fixture
def dialog(app, bundle):
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    window = _loaded(bundle)
    d = DynamiXsDialog(parent=None, main_window=window)
    yield d
    d.deleteLater()


def _listed(dialog):
    from PySide6.QtCore import Qt
    lst = dialog.saved_runs_list
    return [lst.item(i).data(Qt.UserRole) for i in range(lst.count())]


def _select(dialog, name):
    from PySide6.QtCore import Qt
    lst = dialog.saved_runs_list
    item = next(lst.item(i) for i in range(lst.count()) if lst.item(i).data(Qt.UserRole) == name)
    lst.setCurrentItem(item)
    return item


def test_the_main_menu_lists_the_projects_saved_runs(dialog):
    names = sorted(dialog.main_window.dynamixs_analyses)
    assert len(names) == 2 and HETNOE in names
    assert _listed(dialog) == names
    assert not dialog.saved_runs_list.isHidden()
    assert dialog.no_saved_runs_label.isHidden()


def test_each_listed_run_names_its_analysis_type(dialog):
    texts = [dialog.saved_runs_list.item(i).text() for i in range(dialog.saved_runs_list.count())]
    assert any(HETNOE in t and "hetNOE" in t for t in texts)
    assert any("T1/T2" in t for t in texts)


def test_opening_a_listed_run_puts_it_on_its_page(dialog):
    _select(dialog, HETNOE)
    dialog._open_selected_run()
    assert dialog.stack.currentWidget() is dialog.hetnoe_page
    assert dialog._page_has_results(dialog.get_state(), 'hetnoe')


def test_a_double_click_opens_the_run(dialog):
    item = _select(dialog, HETNOE)
    dialog.saved_runs_list.itemDoubleClicked.emit(item)
    assert dialog.stack.currentWidget() is dialog.hetnoe_page


def test_deleting_a_listed_run_removes_it_from_the_session(dialog, monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.Yes))
    _select(dialog, HETNOE)
    dialog._delete_selected_run()
    assert HETNOE not in dialog.main_window.dynamixs_analyses
    assert HETNOE not in _listed(dialog)


def test_a_declined_delete_keeps_the_run(dialog, monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.No))
    _select(dialog, HETNOE)
    dialog._delete_selected_run()
    assert HETNOE in dialog.main_window.dynamixs_analyses


def test_a_run_saved_from_a_page_is_listed_back_at_the_menu(dialog):
    _select(dialog, HETNOE)
    dialog._open_selected_run()
    dialog._store_run('hetnoe', 'f600_hetNOE_copy')
    dialog.show_main_menu()
    assert 'f600_hetNOE_copy' in _listed(dialog)


def test_the_run_buttons_wait_for_a_selection(dialog):
    assert dialog.saved_runs_list.currentItem() is None
    assert not dialog.open_run_button.isEnabled()
    assert not dialog.delete_run_button.isEnabled()
    _select(dialog, HETNOE)
    assert dialog.open_run_button.isEnabled()
    assert dialog.delete_run_button.isEnabled()


def test_a_project_without_runs_says_so(app):
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    window = types.SimpleNamespace(integrator=types.SimpleNamespace(peak_list=None))
    d = DynamiXsDialog(parent=None, main_window=window)
    try:
        assert d.saved_runs_list.count() == 0
        assert d.saved_runs_list.isHidden()
        assert not d.no_saved_runs_label.isHidden()
    finally:
        d.deleteLater()
