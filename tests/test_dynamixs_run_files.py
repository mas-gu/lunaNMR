# ABOUTME: Each DynamiXs GUI run names its files after its own prefix, its viewer opens only
# ABOUTME: those, and a project bundles an analysis with only the files its own runs wrote.

import json
import os
import sys
import types
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

DATA = Path(__file__).resolve().parents[2] / "data_test"
SERIES = DATA / "data_JA1_relaxation" / "run_20260825_1405" / "series"
MATRIX_A = SERIES / "600_T1" / "peak_intensity_matrix.csv"
MATRIX_B = SERIES / "700_T1" / "peak_intensity_matrix.csv"
METHYL_MATRIX = (DATA / "methyl" / "sampleNBD" / "series_results_sampleNBD_20260818_124619"
                 / "peak_intensity_matrix.csv")

pytestmark = pytest.mark.skipif(not MATRIX_A.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def quiet_viewers(monkeypatch):
    """The fit viewer reports each load in a message box; nothing closes one offscreen."""
    import visualization.fit_viewer as fit_viewer
    for box in ("show_info", "show_warning", "show_error"):
        monkeypatch.setattr(fit_viewer, box, lambda *a, **k: None)


def _dialog(window=None):
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    return DynamiXsDialog(parent=None, main_window=window or types.SimpleNamespace())


def _wait(app, page):
    if page.worker is not None:
        page.worker.wait()
    app.processEvents()


def _fit_t1(app, page, matrix, prefix):
    page.field1_t1_file = str(matrix)
    page.analysis_combo.setCurrentIndex(0)          # Field 1 - T1
    page.prefix_edit.setText(prefix)
    page._start_analysis()
    _wait(app, page)


@pytest.fixture(scope="module")
def two_t1t2_runs(app, tmp_path_factory):
    """Two runs of the real T1/T2 page into one output folder, prefixes A then B, each in
    its own dialog; returns (out, state of A's dialog, state of B's dialog, B's page)."""
    out = tmp_path_factory.mktemp("shared_out")
    states, dialogs = [], []
    for prefix, matrix in (("A", MATRIX_A), ("B", MATRIX_B)):
        d = _dialog()
        d.stack.setCurrentWidget(d.t1t2_page)
        d.t1t2_page.output_dir = str(out)
        _fit_t1(app, d.t1t2_page, matrix, prefix)
        states.append(d.get_state())
        dialogs.append(d)
    yield out, states[0], states[1], dialogs[1].t1t2_page
    for d in dialogs:
        d.deleteLater()


class TestTheT1T2PageNamesItsFilesAfterItsPrefix:

    def test_the_fit_json_carries_the_prefix_before_the_field(self, two_t1t2_runs):
        out = two_t1t2_runs[0]
        names = {p.name for p in (out / "json").iterdir()}
        assert {"A_field1_T1_fit_data.json", "B_field1_T1_fit_data.json"} <= names

    def test_a_second_run_leaves_the_first_runs_json_as_it_was(self, two_t1t2_runs):
        out = two_t1t2_runs[0]
        a = json.loads((out / "json" / "A_field1_T1_fit_data.json").read_text())
        b = json.loads((out / "json" / "B_field1_T1_fit_data.json").read_text())
        assert [f["t2"] for f in a["fits"]] != [f["t2"] for f in b["fits"]]

    def test_the_results_file_starts_with_the_prefix(self, two_t1t2_runs):
        out = two_t1t2_runs[0]
        assert (out / "A_field1_T1_fit_results.txt").exists()

    def test_the_run_records_its_file_prefix(self, two_t1t2_runs):
        _, state_a, _, _ = two_t1t2_runs
        assert state_a["t1t2"]["file_prefixes"] == ["A_field1_T1"]


class TestTheFitViewerOpensTheRunsOwnFits:

    def test_the_page_opens_only_its_own_runs_json(self, two_t1t2_runs, quiet_viewers):
        page = two_t1t2_runs[3]
        page._open_fit_viewer()
        try:
            assert [p.name for p in page.fit_viewer.json_paths.values()] == [
                "B_field1_T1_fit_data.json"]
        finally:
            page.fit_viewer.close()

    def test_a_prefixed_json_is_shown_as_its_field_and_experiment(self, app, two_t1t2_runs,
                                                                  quiet_viewers):
        from visualization.fit_viewer import FitViewer
        out = two_t1t2_runs[0]
        viewer = FitViewer(json_files=[str(out / "json" / "A_field1_T1_fit_data.json")])
        try:
            assert list(viewer.data) == ["field1_T1"]
            assert viewer.field1_radio.isEnabled()
        finally:
            viewer.close()


class TestAProjectBundlesEachRunWithItsOwnFiles:

    @pytest.fixture
    def bundle(self, two_t1t2_runs, tmp_path):
        from lunaNMR.utils.project_manager import ProjectManager
        _, state_a, state_b, _ = two_t1t2_runs
        mw = types.SimpleNamespace(dynamixs_dialog=None, dynamixs_analyses={
            "runA": {"state": json.loads(json.dumps(state_a)), "file_refs": {}},
            "runB": {"state": json.loads(json.dumps(state_b)), "file_refs": {}}})
        bundle = tmp_path / "proj.lunaNMR"
        bundle.mkdir()
        assert ProjectManager(mw)._save_dynamixs_state(bundle)
        return bundle / "dynamixs" / "analyses"

    def _files(self, analysis_dir):
        root = analysis_dir / "results" / "t1t2"
        return sorted(str(p.relative_to(root)) for p in root.rglob("*") if p.is_file())

    def test_each_analysis_holds_only_its_runs_files(self, bundle):
        assert self._files(bundle / "runA") == ["A_field1_T1_fit_results.txt",
                                                "json/A_field1_T1_fit_data.json"]
        assert self._files(bundle / "runB") == ["B_field1_T1_fit_results.txt",
                                                "json/B_field1_T1_fit_data.json"]


def test_a_shorter_prefix_does_not_claim_a_longer_ones_files(tmp_path):
    from lunaNMR.utils.project_manager import ProjectManager
    out = tmp_path / "out"
    (out / "json").mkdir(parents=True)
    for name in ("a_field1_T1_fit_results.txt", "a_b_field1_T1_fit_results.txt"):
        (out / name).write_text("x")
    for name in ("a_field1_T1_fit_data.json", "a_b_field1_T1_fit_data.json"):
        (out / "json" / name).write_text("{}")
    page = lambda p: {"t1t2": {"output_dir": str(out), "file_prefixes": [p]}}
    mw = types.SimpleNamespace(dynamixs_dialog=None, dynamixs_analyses={
        "short": {"state": page("a_field1_T1"), "file_refs": {}},
        "long": {"state": page("a_b_field1_T1"), "file_refs": {}}})
    bundle = tmp_path / "proj.lunaNMR"
    bundle.mkdir()
    ProjectManager(mw)._save_dynamixs_state(bundle)
    short = bundle / "dynamixs" / "analyses" / "short" / "results" / "t1t2"
    assert sorted(p.name for p in short.rglob("*") if p.is_file()) == [
        "a_field1_T1_fit_data.json", "a_field1_T1_fit_results.txt"]


@pytest.mark.skipif(not METHYL_MATRIX.exists(), reason="data_test not present")
class TestTheMethylPageNamesItsFilesAfterItsPrefix:

    @pytest.fixture(scope="class")
    def run(self, app, tmp_path_factory):
        d = _dialog()
        page = d.methyl_t2_page
        page.output_dir = str(tmp_path_factory.mktemp("methyl"))
        page._set_input_file(str(METHYL_MATRIX))
        page.prefix_edit.setText("NBD")
        page._start_analysis()
        _wait(app, page)
        yield page
        d.deleteLater()

    def test_the_prefix_is_editable(self, run):
        from PySide6.QtWidgets import QLineEdit
        assert isinstance(run.prefix_edit, QLineEdit)

    def test_the_files_carry_the_prefix(self, run):
        out = Path(run.output_dir)
        assert (out / "json" / "NBD_field1_methylT2_fit_data.json").exists()
        assert (out / "NBD_field1_methylT2_fit_results.txt").exists()

    def test_the_run_records_its_file_prefix_and_json(self, run):
        state = run.get_session_state()
        assert state["file_prefixes"] == ["NBD_field1_methylT2"]
        assert Path(state["last_json_file"]).name == "NBD_field1_methylT2_fit_data.json"

    def test_its_viewer_opens_only_the_runs_own_json(self, run, monkeypatch):
        import shutil
        import visualization.methyl_t2_fit_viewer as methyl_viewer
        for box in ("show_info", "show_warning", "show_error"):
            monkeypatch.setattr(methyl_viewer, box, lambda *a, **k: None)
        json_dir = Path(run.output_dir) / "json"
        shutil.copy(json_dir / "NBD_field1_methylT2_fit_data.json",
                    json_dir / "other_field1_methylT2_fit_data.json")
        run._open_viewer()
        try:
            assert [p.name for p in run._viewer.json_paths.values()] == [
                "NBD_field1_methylT2_fit_data.json"]
        finally:
            run._viewer.close()


class TestAStateRecordingNoRunFilesBundlesAndOpensNothing:

    def test_a_page_recording_no_file_prefixes_bundles_none_of_its_folder(self, tmp_path,
                                                                          caplog):
        import logging
        from lunaNMR.utils.project_manager import ProjectManager
        out = tmp_path / "out"
        (out / "json").mkdir(parents=True)
        (out / "a_field1_T1_fit_results.txt").write_text("x")
        (out / "json" / "a_field1_T1_fit_data.json").write_text("{}")
        mw = types.SimpleNamespace(dynamixs_dialog=None, dynamixs_analyses={
            "run": {"state": {"t1t2": {"output_dir": str(out)}}, "file_refs": {}}})
        bundle = tmp_path / "proj.lunaNMR"
        bundle.mkdir()
        with caplog.at_level(logging.WARNING, logger="lunaNMR.utils.project_manager"):
            ProjectManager(mw)._save_dynamixs_state(bundle)
        assert not (bundle / "dynamixs" / "analyses" / "run" / "results" / "t1t2").exists()
        said = [r.getMessage() for r in caplog.records if "t1t2" in r.getMessage()]
        assert len(said) == 1 and "file prefixes" in said[0]

    def test_a_t1t2_session_recording_no_json_file_opens_no_viewer(self, app, two_t1t2_runs,
                                                                  quiet_viewers):
        _, _, state_b, _ = two_t1t2_runs
        state = json.loads(json.dumps(state_b["t1t2"]))
        for field in state["session_results"].values():
            for entry in (field or {}).values():
                if isinstance(entry, dict):
                    entry.pop("json_file", None)
        d = _dialog()
        try:
            page = d.t1t2_page
            page.restore_session_state(state)
            page._open_fit_viewer()
            assert getattr(page, "fit_viewer", None) is None
            assert "no fit JSON" in page.results_text.toPlainText().splitlines()[-1]
        finally:
            d.deleteLater()

    def test_a_modelfree_session_recording_no_json_files_opens_no_viewer(self, app,
                                                                         two_t1t2_runs,
                                                                         quiet_viewers):
        out = two_t1t2_runs[0]
        d = _dialog()
        try:
            page = d.integrated_page
            page.output_dir = str(out)
            page.session_results = {"analysis_complete": True, "results_prefix": "A",
                                    "output_dir": str(out), "json_folder": str(out / "json")}
            page._open_fit_viewer()
            assert getattr(page, "fit_viewer", None) is None
            assert "no fit JSON" in page.progress_text.toPlainText().splitlines()[-1]
        finally:
            d.deleteLater()
