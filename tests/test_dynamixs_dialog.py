# ABOUTME: Tests the real DynamiXs dialog's analysis naming, metadata, source-series detection
# ABOUTME: and the normalisation of peak data read back from JSON for Inspect Peak.

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
    return types.SimpleNamespace(dynamixs_analyses=None, saved_series={})


@pytest.fixture
def dialog(app, window):
    from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
    d = DynamiXsDialog(parent=None, main_window=window)
    yield d
    d.deleteLater()


def _with_t1t2_run(dialog, tmp_path, series):
    """The T1/T2 page, active, with a stored fit and `series` dropped on its T1 slot."""
    page = dialog.t1t2_page
    dialog.stack.setCurrentWidget(page)
    page.output_dir = str(tmp_path)
    page.fitted_experiments = {"field1_t1"}
    page.source_series_map = {"field1_t1": series}
    return page


class TestAnalysisNaming:

    def test_the_analysis_name_is_stored_in_the_state(self, dialog, tmp_path):
        _with_t1t2_run(dialog, tmp_path, "600_T1")
        dialog.ensure_current_saved()
        assert dialog.get_state()['analysis_name'] == "600_T1_T1T2"

    def test_a_restored_state_names_the_analysis(self, dialog):
        dialog.set_state({'analysis_name': 'T1_relaxation_2024'})
        assert dialog.current_analysis_name == 'T1_relaxation_2024'

    def test_invalid_filename_characters_are_replaced_in_the_name(self, dialog, tmp_path):
        _with_t1t2_run(dialog, tmp_path, 'T1:data/test<2>')
        name = dialog._suggested_analysis_name("t1t2")
        assert name == "T1_data_test_2__T1T2"
        assert not any(c in name for c in '<>:"/\\|?*')

    def test_with_no_series_the_suggestion_is_the_analysis_type(self, dialog, tmp_path):
        _with_t1t2_run(dialog, tmp_path, None)
        assert dialog._suggested_analysis_name("t1t2") == "T1T2"


class TestAnalysisMetadata:

    def test_metadata_includes_the_name_type_and_creation_time(self, dialog, tmp_path):
        _with_t1t2_run(dialog, tmp_path, "600_T1")
        dialog.ensure_current_saved()
        metadata = dialog.get_analysis_metadata()
        assert metadata['analysis_name'] == "600_T1_T1T2"
        assert metadata['analysis_type'] == 't1t2'
        assert metadata['created_at']

    def test_metadata_includes_the_source_series(self, dialog):
        dialog.stack.setCurrentWidget(dialog.integrated_page)
        dialog.set_source_series('T1_asyn_series')
        metadata = dialog.get_analysis_metadata()
        assert metadata['source_series'] == 'T1_asyn_series'
        assert metadata['analysis_type'] == 'model_free'

    def test_metadata_travels_in_the_state(self, dialog):
        dialog.set_analysis_metadata({'analysis_name': 'Test_Analysis',
                                      'source_series': 'Source_Series_1',
                                      'created_at': '2024-01-15T10:30:00'})
        state = dialog.get_state()
        assert state['analysis_name'] == 'Test_Analysis'
        assert state['analysis_metadata']['source_series'] == 'Source_Series_1'
        assert state['analysis_metadata']['created_at'] == '2024-01-15T10:30:00'


class TestSourceSeriesDetection:

    def test_detect_from_series_results_path(self, dialog):
        path = "/project/.lunaNMR/series_results/T1_asyn_series/series_analysis_tidy.csv"
        assert dialog.detect_source_series_from_path(path) == 'T1_asyn_series'

    def test_detect_from_saved_series_match(self, dialog, window):
        window.saved_series = {'T1_relax_2024': object(), 'T2_relax_2024': object()}
        assert dialog.detect_source_series_from_path(
            "/data/T1_relax_2024/results.csv") == 'T1_relax_2024'

    def test_no_detection_for_unrelated_path(self, dialog):
        assert dialog.detect_source_series_from_path("/some/random/data.csv") is None


class TestPeakDataNormalization:

    def test_float_values_remain_floats(self, dialog):
        peak = dialog._normalize_peak_data({'pos_f1': 8.5, 'pos_f2': 120.5,
                                            'intensity': 1000.0})
        assert peak == {'pos_f1': 8.5, 'pos_f2': 120.5, 'intensity': 1000.0}

    def test_single_element_lists_become_floats(self, dialog):
        peak = dialog._normalize_peak_data({'pos_f1': [8.5], 'pos_f2': [120.5],
                                            'intensity': [1000.0]})
        assert peak == {'pos_f1': 8.5, 'pos_f2': 120.5, 'intensity': 1000.0}

    def test_nested_all_peaks_are_normalized(self, dialog):
        peak = dialog._normalize_peak_data({
            'assignment': 'A.8.LEU.H',
            'all_peaks': [{'pos_f1': [120.0], 'pos_f2': [8.5], 'intensity': [500.0]},
                          {'pos_f1': [121.0], 'pos_f2': [8.6], 'intensity': [600.0]}]})
        assert peak['assignment'] == 'A.8.LEU.H'
        assert peak['all_peaks'][0]['pos_f1'] == 120.0
        assert peak['all_peaks'][1]['intensity'] == 600.0

    def test_to_float_handles_none_and_empty(self, dialog):
        assert dialog._to_float(None) == 0.0
        assert dialog._to_float([]) == 0.0
        assert dialog._to_float([5.5]) == 5.5
        assert dialog._to_float(3.14) == 3.14
