# ABOUTME: The DynamiXs relaxation pages fit what the CLI fits and say on the page what it reports.
# ABOUTME: Real JA1 600 MHz T1 series and the sampleNBD methyl series, driven through the real pages.

import os
import sys
import types
from pathlib import Path

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

RUN = (Path(__file__).resolve().parents[2] / "data_test" / "data_JA1_relaxation"
       / "run_20260825_1405")
T1_SERIES = RUN / "series" / "600_T1"
T1_MATRIX = T1_SERIES / "peak_intensity_matrix.csv"

pytestmark = pytest.mark.skipif(not T1_MATRIX.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


def _dialog(window=None):
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    return DynamiXsDialog(parent=None, main_window=window or types.SimpleNamespace())


@pytest.fixture
def t1t2_page(app, tmp_path):
    d = _dialog()
    page = d.t1t2_page
    page.output_dir = str(tmp_path / "out")
    os.makedirs(page.output_dir)
    yield page
    d.deleteLater()


def _wait(app, page):
    worker = page.worker
    if worker is not None:
        worker.wait()
    app.processEvents()


def _fit_t1(app, page, matrix):
    page.field1_t1_file = str(matrix)
    page.analysis_combo.setCurrentIndex(0)          # Field 1 - T1
    page._start_analysis()
    _wait(app, page)
    return page.results_text.toPlainText()


def _matrix(tmp_path, columns=None, rename=None):
    df = pd.read_csv(T1_MATRIX)
    if columns is not None:
        df = df[["Peak_Number", "Assignment", "Reference_X", "Reference_Y", *columns]]
    if rename:
        df = df.rename(columns=rename)
    path = tmp_path / "matrix.csv"
    df.to_csv(path, index=False)
    return path


class TestTheT1T2PageReportsTheFit:

    def test_a_column_with_no_delay_is_named_on_the_page(self, app, t1t2_page, tmp_path):
        log = _fit_t1(app, t1t2_page, _matrix(tmp_path, rename={"2400": "last_plane"}))
        assert "last_plane" in log
        assert "no parseable delay" in log

    def test_the_fitted_count_is_given_out_of_the_total(self, app, t1t2_page, tmp_path):
        log = _fit_t1(app, t1t2_page, _matrix(tmp_path))
        n_total = len(pd.read_csv(T1_MATRIX))
        assert f"/{n_total} residues fitted" in log

    def test_a_fit_where_nothing_fitted_is_not_reported_as_a_success(self, app, t1t2_page,
                                                                      tmp_path):
        log = _fit_t1(app, t1t2_page, _matrix(tmp_path, columns=["0", "50"]))
        assert "completed successfully" not in log
        assert "No residue fitted" in log
        assert "field1_t1" not in t1t2_page.fitted_experiments


class TestAFileChosenForASlotReplacesTheSeriesThere:
    """A slot that once held a dropped series is fitted from the file chosen after it,
    not from the series still in memory."""

    @pytest.fixture
    def page(self, app, tmp_path, series_from_tidy):
        br, _ = series_from_tidy(T1_SERIES / "series_analysis_tidy.csv")
        d = _dialog(types.SimpleNamespace(saved_series={"t1_series": br},
                                          current_project_path=None))
        page = d.t1t2_page
        page.output_dir = str(tmp_path / "out")
        os.makedirs(page.output_dir)
        page._on_series_dropped("field1_t1", "t1_series",
                                str(T1_SERIES / "series_analysis_tidy.csv"))
        yield page
        d.deleteLater()

    def _run(self, app, page):
        page.analysis_combo.setCurrentIndex(0)
        page._start_analysis()
        _wait(app, page)
        return page.results_text.toPlainText()

    def test_the_dropped_series_is_fitted_from_memory(self, app, page):
        assert "Using data from series 't1_series' in memory" in self._run(app, page)

    def test_browsing_a_file_replaces_it(self, app, page, tmp_path, monkeypatch):
        import dynamiXs_gui
        matrix = _matrix(tmp_path)
        monkeypatch.setattr(dynamiXs_gui, "open_file_dialog", lambda *a, **k: str(matrix))
        page._browse_data_file("field1_t1")
        assert "in memory" not in self._run(app, page)
        assert page.field1_t1_file == str(matrix)

    def test_auto_loading_a_folder_replaces_it(self, app, page, tmp_path, monkeypatch):
        import dynamiXs_gui
        folder = tmp_path / "folder"
        folder.mkdir()
        pd.read_csv(T1_MATRIX).to_csv(folder / "T1_matrix.csv", index=False)
        monkeypatch.setattr(dynamiXs_gui, "open_directory_dialog", lambda *a, **k: str(folder))
        page._auto_load_folder("field1")
        assert "in memory" not in self._run(app, page)


METHYL_MATRIX = (Path(__file__).resolve().parents[2] / "data_test" / "methyl" / "sampleNBD"
                 / "series_results_sampleNBD_20260818_124619" / "peak_intensity_matrix.csv")


@pytest.mark.skipif(not METHYL_MATRIX.exists(), reason="data_test not present")
class TestTheMethylPageReportsTheFit:

    @pytest.fixture(scope="class")
    def run(self, app, tmp_path_factory):
        """The real page's fit of the methyl series, and the result the worker returned."""
        d = _dialog()
        page = d.methyl_t2_page
        page.output_dir = str(tmp_path_factory.mktemp("methyl"))
        page._set_input_file(str(METHYL_MATRIX))
        results = []
        page._start_analysis()
        page.worker.finished.connect(results.append)
        _wait(app, page)
        yield page.log_text.toPlainText(), results[0]
        d.deleteLater()

    def test_the_log_counts_the_unidentifiable_fits(self, run):
        log, result = run
        n = result["n_bi_exp_unidentifiable"]
        assert n > 0
        assert (f"{result['n_fitted']}/{result['n_total']} residues fitted, "
                f"{n} of them bi-exp unidentifiable") in log


class TestTheModelFreePageReportsTheFit:

    @pytest.fixture(scope="class")
    def run(self, app, tmp_path_factory):
        """The real page's single-field model-free fit of the 600 MHz data."""
        d = _dialog()
        page = d.integrated_page
        noe = RUN / "results" / "hetnoe"
        page.field1_t1_file = str(T1_MATRIX)
        page.field1_t2_file = str(RUN / "series" / "600_T2" / "peak_intensity_matrix.csv")
        page.field1_noe_sat_file = str(noe / "f600_noe_sat.csv")
        page.field1_noe_unsat_file = str(noe / "f600_noe_unsat.csv")
        page.field1_freq_spin.setValue(600.0)
        next(radio for radio, value in page.method_radios if value == "single_jwh").setChecked(True)
        page.mc_spin.setValue(2)
        page.output_dir = str(tmp_path_factory.mktemp("modelfree"))
        results = []
        # The page opens both viewers when the fit finishes, and each reports its load
        # in a message box: answered here, since nothing closes a modal offscreen.
        import visualization.fit_viewer as fit_viewer
        import visualization.results_viewer as results_viewer
        with pytest.MonkeyPatch.context() as mp:
            for viewer in (fit_viewer, results_viewer):
                for box in ("show_info", "show_warning", "show_error"):
                    mp.setattr(viewer, box, lambda *a, **k: None)
            page._start_analysis()
            page.worker.finished.connect(results.append)
            _wait(app, page)
        yield page.progress_text.toPlainText(), results[0], page
        for viewer in ("results_viewer", "fit_viewer"):
            if getattr(page, viewer, None) is not None:
                getattr(page, viewer).close()
        d.deleteLater()

    def test_the_fit_jsons_carry_the_runs_prefix(self, run):
        out = Path(run[2].output_dir)
        assert {p.name for p in (out / "json").iterdir()} == {
            "integrated_analysis_field1_T1_fit_data.json",
            "integrated_analysis_field1_T2_fit_data.json"}

    def test_its_fit_viewer_shows_both_fits_of_the_run(self, run):
        assert set(run[2].fit_viewer.data) == {"field1_T1", "field1_T2"}

    def test_the_log_gives_the_fitted_count_out_of_the_total(self, run):
        log, result, _ = run
        assert (f"Model-free complete: {result['n_successful']}/{result['n_residues']} "
                f"residues") in log

    def test_the_log_counts_the_fits_at_the_iteration_cap(self, run):
        log, result, _ = run
        assert (f"{result['n_at_iteration_cap']} fits stopped at the iteration cap "
                f"without converging") in log

    def test_the_log_names_the_residues_the_merge_dropped(self, run):
        log, result, _ = run
        excluded = result["excluded_residues"]["field1"]
        assert excluded
        assert (f"Field 1: {len(excluded)} residues excluded at the merge (not in all "
                f"datasets): {', '.join(excluded)}") in log


class TestASeriesDroppedOnTheModelFreePage:

    @pytest.fixture
    def page(self, app, series_from_tidy):
        br, _ = series_from_tidy(T1_SERIES / "series_analysis_tidy.csv")
        br.metadata["output_folder"] = str(T1_SERIES)
        d = _dialog(types.SimpleNamespace(saved_series={"t1_series": br},
                                          current_project_path=None))
        yield d.integrated_page
        d.deleteLater()

    def test_a_relaxation_slot_takes_the_series_intensity_matrix(self, page):
        page._on_series_dropped("field1_t1", "t1_series",
                                str(T1_SERIES / "series_analysis_tidy.csv"))
        assert page.field1_t1_file == str(T1_MATRIX)

    def test_a_hetnoe_slot_is_not_given_the_tidy_csv(self, page):
        page._on_series_dropped("field1_noe_sat", "t1_series",
                                str(T1_SERIES / "series_analysis_tidy.csv"))
        assert page.field1_noe_sat_file is None
        assert "hetNOE" in page.progress_text.toPlainText()


def test_the_analysis_name_prompt_suggests_a_name_the_project_does_not_hold(app, monkeypatch):
    import lunaNMR.gui.dialogs.dynamixs_dialog as dynamixs_dialog
    suggested = []

    def answer(parent, title, label, text=""):     # the name box, answered as offered
        suggested.append(text)
        return text, True
    monkeypatch.setattr(dynamixs_dialog.QInputDialog, "getText", staticmethod(answer))
    d = _dialog(types.SimpleNamespace(dynamixs_analyses={"Analysis_1": {}, "Analysis_2": {}}))
    try:
        assert d._prompt_analysis_name("T1/T2 Fitting") == "Analysis_3"
        assert suggested == ["Analysis_3"]
    finally:
        d.deleteLater()
