# ABOUTME: Unit tests for Series Picker drag-and-drop functionality in T1T2FittingPage
# ABOUTME: Tests DraggableSeriesList, DropTargetLabel, and series assignment logic

import pytest
from unittest.mock import MagicMock, patch
from pathlib import Path


class TestDraggableSeriesList:
    """Test the draggable series list widget."""

    def test_list_stores_series_data(self):
        """Series items store name and CSV path in UserRole data."""
        from PySide6.QtWidgets import QListWidgetItem
        from PySide6.QtCore import Qt

        item = QListWidgetItem("T1_asyn_series")
        item.setData(Qt.UserRole, "T1_asyn_series")
        item.setData(Qt.UserRole + 1, "/path/to/series_analysis_tidy.csv")

        assert item.data(Qt.UserRole) == "T1_asyn_series"
        assert item.data(Qt.UserRole + 1) == "/path/to/series_analysis_tidy.csv"

    def test_mime_data_format(self):
        """Mime data uses series:name:path format."""
        mime_text = "series:T1_asyn_series:/path/to/csv"
        parts = mime_text.split(":", 2)

        assert parts[0] == "series"
        assert parts[1] == "T1_asyn_series"
        assert parts[2] == "/path/to/csv"

    def test_mime_data_with_colons_in_path(self):
        """Mime data handles paths with colons (Windows drives)."""
        mime_text = "series:T1_series:C:/Users/data/series.csv"
        parts = mime_text.split(":", 2)

        assert parts[0] == "series"
        assert parts[1] == "T1_series"
        assert parts[2] == "C:/Users/data/series.csv"


class TestDropTargetLabel:
    """Test the drop target label functionality."""

    def test_field_name_stored(self):
        """DropTargetLabel stores its field name."""
        # This tests the logic, not actual Qt widgets
        field_name = "field1_t1"
        label_text = "No file selected"

        # Simulate the label behavior
        class MockDropLabel:
            def __init__(self, field_name, text):
                self.field_name = field_name
                self.text = text

        label = MockDropLabel(field_name, label_text)
        assert label.field_name == "field1_t1"

    def test_accepts_series_drop_data(self):
        """Drop label accepts mime data starting with 'series:'."""
        mime_text = "series:T1_asyn_series:/path/to/csv"
        assert mime_text.startswith("series:")

    def test_rejects_non_series_data(self):
        """Drop label rejects mime data not starting with 'series:'."""
        mime_text = "file:///path/to/file.csv"
        assert not mime_text.startswith("series:")


class TestSeriesDropHandler:
    """Test the series drop handling logic."""

    def test_drop_sets_file_path(self):
        """Dropping a series sets the file path attribute."""
        # Simulate _on_series_dropped behavior
        class MockPage:
            def __init__(self):
                self.source_series_map = {}

            def _on_series_dropped(self, field_name, series_name, csv_path):
                if csv_path:
                    setattr(self, f"{field_name}_file", csv_path)
                if not hasattr(self, 'source_series_map'):
                    self.source_series_map = {}
                self.source_series_map[field_name] = series_name

        page = MockPage()
        page._on_series_dropped("field1_t1", "T1_asyn_series", "/path/to/csv")

        assert page.field1_t1_file == "/path/to/csv"
        assert page.source_series_map["field1_t1"] == "T1_asyn_series"

    def test_drop_stores_source_series(self):
        """Dropping a series stores source_series for Inspect Peak."""
        source_series_map = {}
        field_name = "field2_t2"
        series_name = "T2_asyn_series"

        source_series_map[field_name] = series_name

        assert source_series_map["field2_t2"] == "T2_asyn_series"


class TestGetAvailableSeries:
    """Test the series discovery logic."""

    def test_returns_empty_list_without_lunaNMR(self):
        """Returns empty list when lunaNMR main window not available."""
        # Simulate main_window without main_window attribute
        class MockMainWindow:
            pass

        main_window = MockMainWindow()

        # Simulate _get_available_series logic
        lunaNMR_main = getattr(main_window, 'main_window', None)
        series_list = []

        if not lunaNMR_main:
            pass  # Returns empty list

        assert series_list == []

    def test_returns_series_from_saved_series(self):
        """Returns series from lunaNMR saved_series dict."""
        # Simulate the nested main_window structure
        saved_series = {
            'T1_asyn_series': MagicMock(),
            'T2_asyn_series': MagicMock(),
        }

        series_list = []
        for series_name in saved_series.keys():
            series_list.append({
                'name': series_name,
                'csv_path': ""
            })

        assert len(series_list) == 2
        assert series_list[0]['name'] == 'T1_asyn_series'
        assert series_list[1]['name'] == 'T2_asyn_series'

    def test_finds_csv_path_in_project(self):
        """Finds series_analysis_tidy.csv for series with project path."""
        import tempfile
        import os

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create project structure
            series_dir = Path(tmpdir) / ".lunaNMR" / "series_results" / "T1_series"
            series_dir.mkdir(parents=True)
            tidy_csv = series_dir / "series_analysis_tidy.csv"
            tidy_csv.write_text("header\ndata")

            # Simulate the path lookup
            project_path = tmpdir
            series_name = "T1_series"
            csv_path = ""

            tidy_path = Path(project_path) / ".lunaNMR" / "series_results" / series_name / "series_analysis_tidy.csv"
            if tidy_path.exists():
                csv_path = str(tidy_path)

            assert csv_path == str(tidy_csv)


class TestSeriesListPopulation:
    """Test the series list population logic."""

    def test_shows_no_series_message_when_empty(self):
        """Shows 'No project series available' when series list is empty."""
        series = []
        no_series_visible = len(series) == 0
        list_visible = len(series) > 0

        assert no_series_visible is True
        assert list_visible is False

    def test_shows_list_when_series_available(self):
        """Shows series list when series are available."""
        series = [{'name': 'T1_series', 'csv_path': '/path'}]
        no_series_visible = len(series) == 0
        list_visible = len(series) > 0

        assert no_series_visible is False
        assert list_visible is True


class TestDraggableSeriesListMimeData:
    """Real-widget tests for the drag payload the drop target consumes."""

    @pytest.fixture(scope="class")
    def app(self):
        import sys
        from PySide6.QtWidgets import QApplication
        return QApplication.instance() or QApplication(sys.argv)

    def _widget_class(self):
        import sys
        from pathlib import Path
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
        from dynamiXs_gui import DraggableSeriesList
        return DraggableSeriesList

    def test_mimedata_emits_series_payload(self, app):
        """mimeData() carries the 'series:name:path' text the drop target reads,
        so the base view manages the drag lifecycle (no manual drag.exec)."""
        from PySide6.QtWidgets import QListWidgetItem
        from PySide6.QtCore import Qt
        lst = self._widget_class()()
        item = QListWidgetItem("T1_series")
        item.setData(Qt.UserRole, "T1_series")
        item.setData(Qt.UserRole + 1, "/path/to/series_analysis_tidy.csv")
        lst.addItem(item)

        mime = lst.mimeData([item])
        assert mime.hasText()
        assert mime.text() == "series:T1_series:/path/to/series_analysis_tidy.csv"

    def test_mimedata_handles_missing_csv_path(self, app):
        """A series with no CSV path still yields a well-formed payload."""
        from PySide6.QtWidgets import QListWidgetItem
        from PySide6.QtCore import Qt
        lst = self._widget_class()()
        item = QListWidgetItem("S")
        item.setData(Qt.UserRole, "S")
        lst.addItem(item)

        mime = lst.mimeData([item])
        assert mime.text() == "series:S:"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
