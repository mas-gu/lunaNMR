# ABOUTME: A DynamiXs page saves its run to the project under a name it suggests from its own
# ABOUTME: dropped series, never one already taken; a clash asks Replace / Save as / Cancel.

import os
import sys
import types

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def window():
    return types.SimpleNamespace(dynamixs_analyses=None)


@pytest.fixture
def dialog(app, window):
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    d = DynamiXsDialog(parent=None, main_window=window)
    yield d
    d.deleteLater()


def _with_t1t2_run(dialog, tmp_path, series=("600_T1", "600_T2")):
    """The T1/T2 page, active, with a stored fit and the series dropped on its slots."""
    page = dialog.t1t2_page
    dialog.stack.setCurrentWidget(page)
    page.output_dir = str(tmp_path)
    page.fitted_experiments = {"field1_t1"}
    page.source_series_map = dict(zip(("field1_t1", "field1_t2"), series))
    return page


def _with_model_free_run(dialog, tmp_path):
    page = dialog.integrated_page
    dialog.stack.setCurrentWidget(page)
    page.output_dir = str(tmp_path)
    page.session_results = {"analysis_complete": True, "results_prefix": "mf"}
    page.source_series_map = {"field1_t1": "600_T1", "field1_t2": "600_T2"}
    return page


class TestTheSuggestedName:

    def test_joins_the_pages_differing_series_before_its_type(self, dialog, tmp_path):
        _with_t1t2_run(dialog, tmp_path)
        assert dialog._suggested_analysis_name("t1t2") == "600_T1+600_T2_T1T2"

    def test_names_one_series_once(self, dialog, tmp_path):
        _with_t1t2_run(dialog, tmp_path, series=("600_T1", "600_T1"))
        assert dialog._suggested_analysis_name("t1t2") == "600_T1_T1T2"

    def test_is_never_a_name_already_taken(self, dialog, window, tmp_path):
        _with_t1t2_run(dialog, tmp_path)
        window.dynamixs_analyses = {"600_T1+600_T2_T1T2": {}, "600_T1+600_T2_T1T2_2": {}}
        assert dialog._suggested_analysis_name("t1t2") == "600_T1+600_T2_T1T2_3"

    @pytest.mark.parametrize("page_key, attr, label", [
        ("methyl_t2", "methyl_t2_page", "methylT2"),
        ("hetnoe", "hetnoe_page", "hetNOE"),
        ("integrated", "integrated_page", "modelFree")])
    def test_comes_from_each_pages_own_series(self, dialog, page_key, attr, label):
        dialog.t1t2_page.source_series_map = {"field1_t1": "other"}
        getattr(dialog, attr).source_series_map = {"slot": "series_a"}
        assert dialog._suggested_analysis_name(page_key) == f"series_a_{label}"


class TestSavingARun:

    def test_stores_it_under_the_chosen_name(self, dialog, window, tmp_path):
        _with_t1t2_run(dialog, tmp_path)
        assert dialog._store_run("t1t2", "mine") == "mine"
        assert window.dynamixs_analyses["mine"]["state"]["t1t2"]["output_dir"] == str(tmp_path)

    def test_t1t2_and_model_free_keep_their_own_names(self, dialog, window, tmp_path):
        _with_t1t2_run(dialog, tmp_path)
        dialog._store_run("t1t2", "relax")
        _with_model_free_run(dialog, tmp_path)
        assert dialog.ensure_current_saved() == "600_T1+600_T2_modelFree"
        assert set(window.dynamixs_analyses) == {"relax", "600_T1+600_T2_modelFree"}

    def test_the_project_save_upserts_under_the_name_the_page_saved_with(self, dialog, window,
                                                                        tmp_path):
        _with_t1t2_run(dialog, tmp_path)
        dialog._store_run("t1t2", "relax")
        assert dialog.ensure_current_saved() == "relax"
        assert list(window.dynamixs_analyses) == ["relax"]

    def test_a_reopened_analysis_is_saved_back_under_its_name(self, dialog, window, tmp_path):
        _with_t1t2_run(dialog, tmp_path)
        dialog._store_run("t1t2", "relax")
        entry = window.dynamixs_analyses["relax"]
        from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
        other = DynamiXsDialog(parent=None, main_window=window)
        try:
            other.open_analysis(entry, "relax")
            assert other.ensure_current_saved() == "relax"
        finally:
            other.deleteLater()

    def test_a_page_with_no_run_saves_nothing(self, dialog, window):
        dialog.stack.setCurrentWidget(dialog.methyl_t2_page)
        assert dialog.save_run_to_project(dialog.methyl_t2_page) is None
        assert not window.dynamixs_analyses

    def test_every_page_with_a_run_has_the_button(self, dialog):
        for page in (dialog.t1t2_page, dialog.integrated_page, dialog.methyl_t2_page,
                     dialog.hetnoe_page, dialog.spectral_page):
            assert page.save_run_btn.text() == "Save run to project"


def _answer_modals(app, *answers):
    """Answer the next modal dialogs, in order, as a user would: each answer is a callable
    given the active modal widget. Nothing closes a modal offscreen."""
    from PySide6.QtCore import QTimer
    pending = list(answers)

    def answer_next():
        widget = app.activeModalWidget()
        if widget is None:
            QTimer.singleShot(20, answer_next)
            return
        pending.pop(0)(widget)
        if pending:
            QTimer.singleShot(20, answer_next)
    QTimer.singleShot(20, answer_next)
    return pending


def test_the_real_prompts_name_the_run_and_ask_on_a_clash(app, dialog, window, tmp_path):
    from PySide6.QtWidgets import QInputDialog, QMessageBox
    page = _with_t1t2_run(dialog, tmp_path)
    window.dynamixs_analyses = {"taken": {"state": {}, "file_refs": {}}}
    offered = []

    def type_taken(box):
        assert isinstance(box, QInputDialog)
        offered.append(box.textValue())
        box.setTextValue("taken")
        box.accept()

    def replace(box):
        assert isinstance(box, QMessageBox)
        offered.append(box.text())
        next(b for b in box.buttons() if b.text() == "Replace").click()

    pending = _answer_modals(app, type_taken, replace)
    assert dialog.save_run_to_project(page) == "taken"
    assert not pending
    assert offered[0] == "600_T1+600_T2_T1T2"
    assert "'taken' already exists" in offered[1]
    assert window.dynamixs_analyses["taken"]["state"]["t1t2"]["output_dir"] == str(tmp_path)


def test_cancelling_the_name_saves_nothing(app, dialog, window, tmp_path):
    page = _with_t1t2_run(dialog, tmp_path)
    pending = _answer_modals(app, lambda box: box.reject())
    assert dialog.save_run_to_project(page) is None
    assert not pending
    assert not window.dynamixs_analyses


def test_the_kd_page_asks_through_the_same_prompt():
    import lunaNMR.gui.dialogs.dynamixs_dialog  # noqa: F401  puts the pages on the path
    import kd_titration_page
    import gui_components
    assert kd_titration_page.prompt_name_collision is gui_components.prompt_name_collision
