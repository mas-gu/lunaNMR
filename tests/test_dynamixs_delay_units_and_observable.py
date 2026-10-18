# ABOUTME: The DynamiXs relaxation pages state their delay units (refused against the series'
# ABOUTME: record or unit-stating headers) and fit the intensity the user chose: height or volume.

import json
import os
import shutil
import sys
import types
from pathlib import Path

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

DATA = Path(__file__).resolve().parents[2] / "data_test"
RUN = DATA / "data_JA1_relaxation" / "run_20260825_1405"
T1_SERIES = RUN / "series" / "600_T1"
T1_MATRIX = T1_SERIES / "peak_intensity_matrix.csv"
NOE = RUN / "results" / "hetnoe"

pytestmark = pytest.mark.skipif(not T1_MATRIX.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def errors(monkeypatch):
    """Every error box a page opens, answered and recorded: none closes offscreen."""
    import lunaNMR.gui.dialogs.dynamixs_dialog  # noqa: F401  puts the pages on the path
    import dynamiXs_gui
    import methyl_t2_page
    shown = []
    for module in (dynamiXs_gui, methyl_t2_page):
        monkeypatch.setattr(module, "show_error", lambda parent, title, text: shown.append(text))
    return shown


def _dialog(window=None):
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    return DynamiXsDialog(parent=None, main_window=window or types.SimpleNamespace())


def _wait(app, page):
    if page.worker is not None:
        page.worker.wait()
    app.processEvents()


def _series_folder(tmp_path, units="ms", seconds=False):
    """The 600 MHz T1 series under a neutral name, its sidecar recording `units`; with
    `seconds`, its matrix delay headers written in seconds."""
    folder = tmp_path / "series_a"
    folder.mkdir()
    for name in ("series_analysis_tidy.csv", "peak_volume_matrix.csv"):
        os.symlink(T1_SERIES / name, folder / name)
    matrix = pd.read_csv(T1_MATRIX)
    if seconds:
        def in_seconds(column):          # '300' -> '0.3', a repeat '300_2' -> '0.3_2'
            delay, _, repeat = column.partition("_")
            return f"{float(delay) / 1000:g}" + (f"_{repeat}" if repeat else "")
        matrix = matrix.rename(columns={c: in_seconds(c) for c in matrix.columns
                                        if c.split("_")[0].isdigit()})
    matrix.to_csv(folder / "peak_intensity_matrix.csv", index=False)
    meta = json.loads((T1_SERIES / "series_metadata.json").read_text())
    meta["value_units"] = units
    (folder / "series_metadata.json").write_text(json.dumps(meta))
    return folder


@pytest.fixture
def series_window(series_from_tidy):
    """A main window holding the 600 MHz T1 series, built from its tidy CSV."""
    def build(folder):
        br, _ = series_from_tidy(T1_SERIES / "series_analysis_tidy.csv")
        br.metadata["output_folder"] = str(folder)
        return types.SimpleNamespace(saved_series={"series_a": br}, current_project_path=None)
    return build


# ---------- delay units ----------

@pytest.mark.parametrize("page_name", ["t1t2_page", "methyl_t2_page", "integrated_page"])
def test_each_relaxation_page_has_one_delay_units_selector_defaulting_to_ms(app, page_name):
    d = _dialog()
    try:
        combo = getattr(d, page_name).units_combo
        assert [combo.itemText(i) for i in range(combo.count())] == ["ms", "s", "us"]
        assert combo.currentText() == "ms"
    finally:
        d.deleteLater()


class TestADroppedSeriesGivesItsRecordedUnits:

    def test_on_the_t1t2_page(self, app, tmp_path, series_window):
        folder = _series_folder(tmp_path, units="s")
        d = _dialog(series_window(folder))
        try:
            page = d.t1t2_page
            page._on_series_dropped("field1_t1", "series_a",
                                    str(folder / "series_analysis_tidy.csv"))
            assert page.units_combo.currentText() == "s"
            assert "series_metadata.json" in page.results_text.toPlainText()
        finally:
            d.deleteLater()

    def test_on_the_methyl_page(self, app, tmp_path, series_window):
        folder = _series_folder(tmp_path, units="s")
        d = _dialog(series_window(folder))
        try:
            page = d.methyl_t2_page
            page._on_series_dropped("methyl_t2", "series_a", "")
            assert page.units_combo.currentText() == "s"
            assert "series_metadata.json" in page.log_text.toPlainText()
        finally:
            d.deleteLater()

    def test_on_the_model_free_page(self, app, tmp_path, series_window):
        folder = _series_folder(tmp_path, units="s")
        d = _dialog(series_window(folder))
        try:
            page = d.integrated_page
            page._on_series_dropped("field1_t1", "series_a", "")
            assert page.units_combo.currentText() == "s"
            assert "series_metadata.json" in page.progress_text.toPlainText()
        finally:
            d.deleteLater()


class TestUnitsThatContradictTheRecordAreRefused:

    def test_on_the_t1t2_page(self, app, tmp_path, errors):
        folder = _series_folder(tmp_path, units="ms")
        d = _dialog()
        try:
            page = d.t1t2_page
            page.output_dir = str(tmp_path)
            page.field1_t1_file = str(folder / "peak_intensity_matrix.csv")
            page.analysis_combo.setCurrentIndex(0)
            page.units_combo.setCurrentText("s")
            page._start_analysis()
            assert page.worker is None
            assert "records the delays as 'ms'" in errors[0]
        finally:
            d.deleteLater()

    def test_against_headers_that_state_their_units(self, app, tmp_path, errors):
        matrix = pd.read_csv(T1_MATRIX)
        matrix = matrix.rename(columns={c: f"T1_{c}ms" for c in matrix.columns if c.isdigit()})
        path = tmp_path / "stated.csv"
        matrix.to_csv(path, index=False)
        d = _dialog()
        try:
            page = d.t1t2_page
            page.output_dir = str(tmp_path)
            page.field1_t1_file = str(path)
            page.analysis_combo.setCurrentIndex(0)
            page.units_combo.setCurrentText("s")
            page._start_analysis()
            assert page.worker is None
            assert "state their own units" in errors[0]
        finally:
            d.deleteLater()

    def test_on_the_methyl_page(self, app, tmp_path, errors):
        folder = _series_folder(tmp_path, units="ms")
        d = _dialog()
        try:
            page = d.methyl_t2_page
            page.output_dir = str(tmp_path)
            page._set_input_file(str(folder / "peak_intensity_matrix.csv"))
            page.units_combo.setCurrentText("us")
            page._start_analysis()
            assert getattr(page, "worker", None) is None
            assert "records the delays as 'ms'" in errors[0]
        finally:
            d.deleteLater()

    def test_on_the_model_free_page(self, app, tmp_path, errors):
        d = _dialog()
        try:
            page = d.integrated_page
            page.field1_t1_file = str(T1_MATRIX)
            page.field1_t2_file = str(RUN / "series" / "600_T2" / "peak_intensity_matrix.csv")
            page.field1_noe_sat_file = str(NOE / "f600_noe_sat.csv")
            page.field1_noe_unsat_file = str(NOE / "f600_noe_unsat.csv")
            next(radio for radio, value in page.method_radios
                 if value == "single_jwh").setChecked(True)
            page.output_dir = str(tmp_path)
            page.units_combo.setCurrentText("s")
            page._start_analysis()
            assert getattr(page, "worker", None) is None
            assert "records the delays as 'ms'" in errors[0]
        finally:
            d.deleteLater()


def test_a_t1t2_series_in_seconds_is_fitted_in_seconds(app, tmp_path):
    folder = _series_folder(tmp_path, units="s", seconds=True)
    d = _dialog()
    try:
        page = d.t1t2_page
        page.output_dir = str(tmp_path / "out")
        os.makedirs(page.output_dir)
        page.field1_t1_file = str(folder / "peak_intensity_matrix.csv")
        page.analysis_combo.setCurrentIndex(0)
        page.units_combo.setCurrentText("s")
        page._start_analysis()
        _wait(app, page)
        (json_file,) = (Path(page.output_dir) / "json").glob("*_fit_data.json")
        data = json.loads(json_file.read_text())
        assert data["metadata"]["time_units"] == "s"
        t1 = pd.Series([f["t2"] for f in data["fits"]]).median()
        assert 0.2 < t1 < 2.0             # a JA1 600 MHz T1, in seconds
    finally:
        d.deleteLater()


def test_the_model_free_page_passes_its_units_to_every_series(app, tmp_path):
    import visualization.fit_viewer as fit_viewer
    import visualization.results_viewer as results_viewer
    d = _dialog()
    try:
        page = d.integrated_page
        page.field1_t1_file = str(T1_MATRIX)
        page.field1_t2_file = str(RUN / "series" / "600_T2" / "peak_intensity_matrix.csv")
        page.field1_noe_sat_file = str(NOE / "f600_noe_sat.csv")
        page.field1_noe_unsat_file = str(NOE / "f600_noe_unsat.csv")
        next(radio for radio, value in page.method_radios if value == "single_jwh").setChecked(True)
        page.mc_spin.setValue(2)
        page.output_dir = str(tmp_path)
        with pytest.MonkeyPatch.context() as mp:
            for viewer in (fit_viewer, results_viewer):
                for box in ("show_info", "show_warning", "show_error"):
                    mp.setattr(viewer, box, lambda *a, **k: None)
            page._start_analysis()
            params = page.worker.params
            _wait(app, page)
        assert params.time_units == "ms"
        for viewer in ("results_viewer", "fit_viewer"):
            if getattr(page, viewer, None) is not None:
                getattr(page, viewer).close()
    finally:
        d.deleteLater()


def test_the_pipeline_takes_the_pages_units_for_each_series():
    from workers import IntegratedAnalysisParams, pipeline_units
    assert pipeline_units(IntegratedAnalysisParams(time_units="s")) == {
        "field1_t1_units": "s", "field1_t2_units": "s",
        "field2_t1_units": "s", "field2_t2_units": "s"}


# ---------- height or volume ----------

@pytest.mark.parametrize("page_name", ["t1t2_page", "methyl_t2_page", "integrated_page",
                                       "hetnoe_page"])
def test_each_page_says_what_it_fits_on_and_defaults_to_height(app, page_name):
    d = _dialog()
    try:
        combo = getattr(d, page_name).observable_combo
        assert [combo.itemText(i) for i in range(combo.count())] == ["Height", "Volume"]
        assert combo.currentText() == "Height"
    finally:
        d.deleteLater()


class TestADroppedSeriesGivesTheMatrixOfTheChosenObservable:

    @pytest.mark.parametrize("observable, matrix", [("Height", "peak_intensity_matrix.csv"),
                                                    ("Volume", "peak_volume_matrix.csv")])
    def test_on_the_methyl_page(self, app, tmp_path, series_window, observable, matrix):
        folder = _series_folder(tmp_path)
        d = _dialog(series_window(folder))
        try:
            page = d.methyl_t2_page
            page.observable_combo.setCurrentText(observable)
            page._on_series_dropped("methyl_t2", "series_a", "")
            assert page.input_file == str(folder / matrix)
        finally:
            d.deleteLater()

    def test_a_choice_made_after_the_drop_takes_the_other_matrix(self, app, tmp_path,
                                                                 series_window):
        folder = _series_folder(tmp_path)
        d = _dialog(series_window(folder))
        try:
            page = d.methyl_t2_page
            page._on_series_dropped("methyl_t2", "series_a", "")
            page.observable_combo.setCurrentText("Volume")
            assert page.input_file == str(folder / "peak_volume_matrix.csv")
        finally:
            d.deleteLater()

    def test_on_the_model_free_page(self, app, tmp_path, series_window):
        folder = _series_folder(tmp_path)
        d = _dialog(series_window(folder))
        try:
            page = d.integrated_page
            page.observable_combo.setCurrentText("Volume")
            page._on_series_dropped("field1_t1", "series_a", "")
            assert page.field1_t1_file == str(folder / "peak_volume_matrix.csv")
        finally:
            d.deleteLater()

    def test_on_the_hetnoe_page(self, app, tmp_path, series_window):
        folder = _series_folder(tmp_path)
        d = _dialog(series_window(folder))
        try:
            page = d.hetnoe_page
            page.observable_combo.setCurrentText("Volume")
            page._on_series_dropped("hetnoe_matrix", "series_a", "")
            assert page.matrix_file == str(folder / "peak_volume_matrix.csv")
        finally:
            d.deleteLater()


@pytest.mark.parametrize("observable, column", [("Height", "height"), ("Volume", "volume")])
def test_the_t1t2_page_fits_the_chosen_observable_of_a_series_in_memory(
        app, tmp_path, series_window, observable, column):
    folder = _series_folder(tmp_path)
    d = _dialog(series_window(folder))
    try:
        page = d.t1t2_page
        page.output_dir = str(tmp_path / "out")
        os.makedirs(page.output_dir)
        page._on_series_dropped("field1_t1", "series_a", str(folder / "series_analysis_tidy.csv"))
        page.observable_combo.setCurrentText(observable)
        page.analysis_combo.setCurrentIndex(0)
        page._start_analysis()
        _wait(app, page)
        assert "in memory" in page.results_text.toPlainText()
        (json_file,) = (Path(page.output_dir) / "json").glob("*_fit_data.json")
        fit = json.loads(json_file.read_text())["fits"][0]
        tidy = pd.read_csv(T1_SERIES / "series_analysis_tidy.csv", dtype={"spectrum_name": str})
        expected = tidy[tidy["assignment"] == fit["residue"]][column]
        assert max(fit["intensities"]) == pytest.approx(expected.max())
    finally:
        d.deleteLater()


def test_a_project_carries_a_series_volume_matrix(tmp_path, series_from_tidy):
    from lunaNMR.utils.project_manager import ProjectManager
    br, _ = series_from_tidy(T1_SERIES / "series_analysis_tidy.csv")
    br.metadata["output_folder"] = str(T1_SERIES)
    dest = tmp_path / "bundled"
    dest.mkdir()
    ProjectManager(types.SimpleNamespace())._copy_series_csvs(br, dest)
    assert (dest / "peak_volume_matrix.csv").read_text() == (
        T1_SERIES / "peak_volume_matrix.csv").read_text()
