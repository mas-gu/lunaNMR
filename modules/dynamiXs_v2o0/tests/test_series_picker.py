# ABOUTME: Tests the series picker of the real T1/T2 page: the draggable list of the project's
# ABOUTME: series, its drag payload, the drop targets, and what a dropped series assigns.

import os
import sys
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DIR))

DATA = Path(__file__).resolve().parents[4] / "data_test"
SERIES = DATA / "data_JA1_relaxation" / "run_20260825_1405" / "series"
T1_TIDY = SERIES / "600_T1" / "series_analysis_tidy.csv"
T2_TIDY = SERIES / "600_T2" / "series_analysis_tidy.csv"


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


def _series(folder=None, bundled_folder=None):
    """A series as LunaNMR holds it: run into `folder` this session, or opened from a
    project whose copy of its files sits in `bundled_folder`."""
    from dynamiXs_gui import _ensure_lunaNMR_path
    _ensure_lunaNMR_path()
    from lunaNMR.processors.multi_spectrum_processor import BatchResults
    br = BatchResults()
    if folder:
        br.metadata['output_folder'] = str(folder)
    if bundled_folder:
        br.bundled_folder = str(bundled_folder)
    return br


@pytest.fixture
def window(app):
    """The standalone DynamiXs window, the main window its pages read the series from."""
    from dynamiXs_gui import DynamiXsMainWindow
    w = DynamiXsMainWindow()
    w.saved_series = {}
    yield w
    w.deleteLater()


def _item(name, csv_path=None):
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QListWidgetItem
    item = QListWidgetItem(name)
    item.setData(Qt.UserRole, name)
    if csv_path is not None:
        item.setData(Qt.UserRole + 1, csv_path)
    return item


def _drag_payload(text):
    from PySide6.QtCore import QMimeData
    mime = QMimeData()
    mime.setText(text)
    return mime


def _drag_enter(label, text):
    """Hand `label` a drag entering it with `text`; whether it accepted."""
    from PySide6.QtCore import QPoint, Qt
    from PySide6.QtGui import QDragEnterEvent
    mime = _drag_payload(text)
    event = QDragEnterEvent(QPoint(1, 1), Qt.CopyAction, mime, Qt.LeftButton, Qt.NoModifier)
    event.ignore()
    label.dragEnterEvent(event)
    return event.isAccepted()


def _drop(label, mime):
    """Drop `mime` on `label`; the (field, series, csv path) it emitted, if any."""
    from PySide6.QtCore import QPointF, Qt
    from PySide6.QtGui import QDropEvent
    emitted = []
    label.series_dropped.connect(lambda *args: emitted.append(args))
    event = QDropEvent(QPointF(1, 1), Qt.CopyAction, mime, Qt.LeftButton, Qt.NoModifier)
    label.dropEvent(event)
    return emitted


class TestDraggableSeriesList:

    @pytest.mark.skipif(not T1_TIDY.exists(), reason="data_test not present")
    def test_the_page_lists_each_series_with_its_name_and_tidy_csv(self, window):
        from PySide6.QtCore import Qt
        window.saved_series = {"600_T1": _series(SERIES / "600_T1")}
        page = window.t1t2_page
        page._populate_series_list()
        item = page.series_list_widget.item(0)
        assert item.data(Qt.UserRole) == "600_T1"
        assert item.data(Qt.UserRole + 1) == str(T1_TIDY)

    def test_the_drag_payload_is_series_name_path(self, app):
        from dynamiXs_gui import DraggableSeriesList
        lst = DraggableSeriesList()
        item = _item("T1_series", "/path/to/series_analysis_tidy.csv")
        lst.addItem(item)
        mime = lst.mimeData([item])
        assert mime.hasText()
        assert mime.text() == "series:T1_series:/path/to/series_analysis_tidy.csv"

    def test_a_series_with_no_csv_path_still_yields_a_payload(self, app):
        from dynamiXs_gui import DraggableSeriesList
        lst = DraggableSeriesList()
        item = _item("S")
        lst.addItem(item)
        assert lst.mimeData([item]).text() == "series:S:"


class TestDropTargetLabel:

    def test_it_keeps_its_field_name(self, app):
        from dynamiXs_gui import DropTargetLabel
        assert DropTargetLabel("field1_t1").field_name == "field1_t1"

    def test_it_accepts_a_series_drag(self, app):
        from dynamiXs_gui import DropTargetLabel
        assert _drag_enter(DropTargetLabel("field1_t1"), "series:T1_series:/path/to/csv")

    def test_it_refuses_anything_else(self, app):
        from dynamiXs_gui import DropTargetLabel
        assert not _drag_enter(DropTargetLabel("field1_t1"), "file:///path/to/file.csv")

    def test_a_drop_emits_its_field_the_series_and_its_path(self, app):
        from dynamiXs_gui import DropTargetLabel
        emitted = _drop(DropTargetLabel("field1_t1"),
                        _drag_payload("series:T1_series:/path/to/csv"))
        assert emitted == [("field1_t1", "T1_series", "/path/to/csv")]

    def test_a_path_with_a_drive_colon_survives_the_drop(self, app):
        from dynamiXs_gui import DropTargetLabel
        emitted = _drop(DropTargetLabel("field1_t1"),
                        _drag_payload("series:T1_series:C:/Users/data/series.csv"))
        assert emitted == [("field1_t1", "T1_series", "C:/Users/data/series.csv")]


@pytest.mark.skipif(not T1_TIDY.exists(), reason="data_test not present")
class TestSeriesDropOnThePage:

    def test_dragging_a_listed_series_onto_a_slot_assigns_its_tidy_csv(self, window):
        window.saved_series = {"600_T1": _series(SERIES / "600_T1")}
        page = window.t1t2_page
        page._populate_series_list()
        item = page.series_list_widget.item(0)
        _drop(page.field1_t1_display, page.series_list_widget.mimeData([item]))
        assert page.field1_t1_file == str(T1_TIDY)
        assert page.source_series_map["field1_t1"] == "600_T1"
        assert page.field1_t1_display.text() == "📊 600_T1"

    def test_each_slot_records_the_series_dropped_on_it(self, window):
        page = window.t1t2_page
        page._on_series_dropped("field1_t1", "600_T1", str(T1_TIDY))
        page._on_series_dropped("field1_t2", "600_T2", str(T2_TIDY))
        assert page.source_series_map == {"field1_t1": "600_T1", "field1_t2": "600_T2"}
        assert page.field1_t2_file == str(T2_TIDY)


class TestGetAvailableSeries:

    def test_a_window_holding_no_series_offers_none(self, window):
        del window.saved_series
        assert window.t1t2_page._get_available_series() == []

    @pytest.mark.skipif(not T1_TIDY.exists(), reason="data_test not present")
    def test_every_saved_series_is_offered(self, window):
        window.saved_series = {"600_T1": _series(SERIES / "600_T1"),
                               "600_T2": _series(SERIES / "600_T2")}
        assert window.t1t2_page._get_available_series() == [
            {"name": "600_T1", "csv_path": str(T1_TIDY)},
            {"name": "600_T2", "csv_path": str(T2_TIDY)}]

    @pytest.mark.skipif(not T1_TIDY.exists(), reason="data_test not present")
    def test_a_series_opened_from_a_project_offers_its_bundled_tidy_csv(self, window,
                                                                         tmp_path):
        import shutil
        bundled = tmp_path / "proj.lunaNMR" / "series_results" / "600_T1"
        bundled.mkdir(parents=True)
        shutil.copy(T1_TIDY, bundled / "series_analysis_tidy.csv")
        window.saved_series = {"600_T1": _series(SERIES / "600_T1", bundled_folder=bundled)}
        assert window.t1t2_page._get_available_series() == [
            {"name": "600_T1", "csv_path": str(bundled / "series_analysis_tidy.csv")}]


class TestSeriesListPopulation:

    def test_with_no_series_the_page_says_so_and_hides_the_list(self, window):
        page = window.t1t2_page
        page._populate_series_list()
        assert not page.no_series_label.isHidden()
        assert page.series_list_widget.isHidden()

    def test_with_series_the_page_shows_the_list(self, window):
        window.saved_series = {"600_T1": _series(SERIES / "600_T1")}
        page = window.t1t2_page
        page._populate_series_list()
        assert page.no_series_label.isHidden()
        assert not page.series_list_widget.isHidden()
        assert page.series_list_widget.count() == 1


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
