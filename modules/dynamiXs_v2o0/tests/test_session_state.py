# ABOUTME: Tests that the real T1/T2 page and the standalone DynamiXs window save and restore
# ABOUTME: the page's field assignments (T1 Data, T2 Data) and the series dropped on them.

import os
import sys
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_DIR))

SERIES = (Path(__file__).resolve().parents[4] / "data_test" / "data_JA1_relaxation"
          / "run_20260825_1405" / "series")
T1_TIDY = SERIES / "600_T1" / "series_analysis_tidy.csv"
T2_TIDY = SERIES / "600_T2" / "series_analysis_tidy.csv"

pytestmark = pytest.mark.skipif(not T1_TIDY.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def windows(app):
    """A factory of standalone DynamiXs windows, each closed afterwards."""
    from dynamiXs_gui import DynamiXsMainWindow
    made = []

    def make():
        made.append(DynamiXsMainWindow())
        return made[-1]
    yield make
    for w in made:
        w.deleteLater()


def _assigned(window):
    """The window's T1/T2 page with series dropped on Field 1's T1 and T2 slots."""
    page = window.t1t2_page
    page._on_series_dropped("field1_t1", "600_T1", str(T1_TIDY))
    page._on_series_dropped("field1_t2", "600_T2", str(T2_TIDY))
    return page


class TestT1T2FittingPageSessionState:

    def test_get_session_state_includes_source_series_map(self, windows):
        state = _assigned(windows()).get_session_state()
        assert state['source_series_map'] == {'field1_t1': '600_T1', 'field1_t2': '600_T2'}
        assert state['field1_t1_file'] == str(T1_TIDY)
        assert state['field1_t2_file'] == str(T2_TIDY)

    def test_restore_session_state_restores_source_series_map(self, windows):
        page = windows().t1t2_page
        page.restore_session_state({
            'field1_t1_file': str(T1_TIDY), 'field1_t2_file': str(T2_TIDY),
            'source_series_map': {'field1_t1': '600_T1', 'field1_t2': '600_T2'}})
        assert page.source_series_map == {'field1_t1': '600_T1', 'field1_t2': '600_T2'}
        assert page.field1_t1_file == str(T1_TIDY)
        assert page.field1_t1_display.text() == "📊 600_T1"
        assert page.field1_t2_display.text() == "📊 600_T2"

    def test_session_state_round_trip(self, windows):
        original = _assigned(windows())
        restored = windows().t1t2_page
        restored.restore_session_state(original.get_session_state())
        assert restored.source_series_map == original.source_series_map
        assert restored.field1_t1_file == original.field1_t1_file
        assert restored.field1_t2_file == original.field1_t2_file
        assert restored.field1_t1_display.text() == "📊 600_T1"


class TestDynamiXsMainWindowState:

    def test_get_state_carries_the_pages_series_and_files(self, windows):
        window = windows()
        _assigned(window)
        assert window.get_state()['t1t2']['source_series_map'] == {
            'field1_t1': '600_T1', 'field1_t2': '600_T2'}
        assert window.get_file_refs()['t1t2']['field1_t1_file'] == str(T1_TIDY)

    def test_restore_state_hands_them_back_to_the_page(self, windows):
        saved = windows()
        _assigned(saved)
        window = windows()
        window.restore_state(saved.get_state(), saved.get_file_refs())
        page = window.t1t2_page
        assert page.source_series_map == {'field1_t1': '600_T1', 'field1_t2': '600_T2'}
        assert page.field1_t1_file == str(T1_TIDY)
        assert page.field1_t1_display.text() == "📊 600_T1"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
