# ABOUTME: Where a series' spectra are on disk: one rule, and every view that opens a spectrum follows it.
# ABOUTME: Real DNAJA1 series in a bundle, opened through the main window, the series and DynamiXs dialogs.

"""
A series result records each spectrum by its file name and the series records the
folder they were read from. The overlay viewer is opened from three places and the
Spectrum Browser opens the Spectrum Viewer: all four must find the same file.
"""

import json
import os
import sys

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

SERIES = "DNAJA1"


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture(scope="module")
def stored(dnaja1_series_bundle):
    return json.loads((dnaja1_series_bundle / "series_results" / SERIES
                       / "batch_results.json").read_text())


@pytest.fixture(scope="module")
def spectra_on_disk(stored):
    """{spectrum name: its file in the folder the series was read from}."""
    folder = stored['metadata']['data_folder']
    files = {name: os.path.join(folder, name) for name in stored['results']}
    assert len(files) == 3 and all(os.path.exists(path) for path in files.values())
    return files


class TestTheRule:
    def test_a_recorded_file_name_is_looked_for_in_the_series_data_folder(self, stored, spectra_on_disk):
        from lunaNMR.gui.dialogs.spectrum_files import series_data_folder, series_spectrum_file
        from lunaNMR.utils.project_manager import batch_results_from_dict
        batch = batch_results_from_dict(stored)
        folder = series_data_folder(batch)
        assert folder == stored['metadata']['data_folder']
        for name, result in batch.results.items():
            assert not os.path.isabs(result['spectrum_file'])
            assert series_spectrum_file(folder, name, result) == spectra_on_disk[name]

    def test_a_series_held_as_a_dict_records_its_folder_the_same_way(self, stored):
        from lunaNMR.gui.dialogs.spectrum_files import series_data_folder
        assert series_data_folder(stored) == stored['metadata']['data_folder']

    def test_a_series_that_records_no_folder_falls_back_on_the_one_given(self, stored, tmp_path):
        from lunaNMR.gui.dialogs.spectrum_files import series_data_folder
        from lunaNMR.utils.project_manager import batch_results_from_dict
        batch = batch_results_from_dict(stored)
        del batch.metadata['data_folder']
        assert series_data_folder(batch) is None
        assert series_data_folder(batch, str(tmp_path)) == str(tmp_path)

    def test_an_absolute_recorded_path_is_kept(self, stored, spectra_on_disk, tmp_path):
        from lunaNMR.gui.dialogs.spectrum_files import series_spectrum_file
        name, path = next(iter(spectra_on_disk.items()))
        assert series_spectrum_file(str(tmp_path), name, {'spectrum_file': path}) == path

    def test_with_no_folder_known_the_recorded_name_is_all_there_is(self, stored):
        from lunaNMR.gui.dialogs.spectrum_files import series_spectrum_file
        name, result = next(iter(stored['results'].items()))
        assert series_spectrum_file(None, name, result) == result['spectrum_file']


@pytest.fixture
def window(app, dnaja1_series_bundle):
    """A main window that opened the project bundle."""
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    set_ps2d_config('15N')
    w = LunaNMRMainWindow()
    assert w._get_project_manager().load_project(dnaja1_series_bundle)[0]
    yield w
    w.close()
    w.deleteLater()
    set_ps2d_config('15N')


def _files_of(overlay):
    return {spec['name']: spec['file_path'] for spec in overlay.spectra}


class TestEveryOpenerFollowsIt:
    def test_the_main_window_s_overlay_viewer(self, window, spectra_on_disk):
        from lunaNMR.gui.dialogs.multi_spectrum_viewer_dialog import MultiSpectrumViewerDialog
        window.open_multi_spectrum_viewer()
        [overlay] = window.findChildren(MultiSpectrumViewerDialog)
        assert _files_of(overlay) == spectra_on_disk

    def test_the_series_dialog_s_overlay_viewer(self, window, spectra_on_disk):
        from lunaNMR.gui.dialogs.multi_spectrum_viewer_dialog import MultiSpectrumViewerDialog
        from lunaNMR.gui.dialogs.series_integration_dialog import SeriesIntegrationDialog
        dialog = SeriesIntegrationDialog(main_window=window)
        try:
            dialog.batch_results = window.batch_results
            dialog._open_overlay_viewer()
            [overlay] = dialog.findChildren(MultiSpectrumViewerDialog)
            assert _files_of(overlay) == spectra_on_disk
            log = dialog.log_text.toPlainText()
            assert f"Data folder for spectrum paths: {os.path.dirname(next(iter(spectra_on_disk.values())))}" in log
            for name, path in spectra_on_disk.items():
                assert f"Constructed full path for {name}: {path}" in log
        finally:
            dialog.deleteLater()

    def test_the_series_dialog_says_when_it_knows_no_folder(self, window):
        from lunaNMR.gui.dialogs.multi_spectrum_viewer_dialog import MultiSpectrumViewerDialog
        from lunaNMR.gui.dialogs.series_integration_dialog import SeriesIntegrationDialog
        import copy
        dialog = SeriesIntegrationDialog(main_window=window)
        try:
            unplaced = copy.copy(window.batch_results)
            unplaced.metadata = {k: v for k, v in unplaced.metadata.items() if k != 'data_folder'}
            window.current_nmr_folder = None
            dialog.batch_results = unplaced
            dialog._open_overlay_viewer()
            [overlay] = dialog.findChildren(MultiSpectrumViewerDialog)
            assert _files_of(overlay) == {name: name for name in unplaced.results}
            log = dialog.log_text.toPlainText()
            for name in unplaced.results:
                assert f"WARNING: No data folder - cannot construct path for {name}" in log
        finally:
            dialog.deleteLater()

    def test_the_dynamixs_dialog_s_overlay_viewer(self, window, dnaja1_series_bundle, spectra_on_disk):
        from lunaNMR.gui.dialogs.dynamixs_dialog import DynamiXsDialog
        window.current_project_path = dnaja1_series_bundle
        dialog = DynamiXsDialog(parent=None, main_window=window)
        try:
            listed = dialog._convert_batch_results_to_list(window.batch_results)
            assert {entry['spectrum_name']: entry['spectrum_file'] for entry in listed} == spectra_on_disk
        finally:
            dialog.deleteLater()

    def test_the_spectrum_browser_s_viewer(self, window, spectra_on_disk):
        from PySide6.QtWidgets import QMenu
        from lunaNMR.gui.dialogs.spectrum_browser_dialog import SpectrumBrowserDialog
        from lunaNMR.gui.dialogs.spectrum_viewer_dialog import SpectrumViewerDialog
        next(action for menu in window.menuBar().findChildren(QMenu) for action in menu.actions()
             if action.text() == "Browse Individual Spectra...").trigger()
        [browser] = window.findChildren(SpectrumBrowserDialog)
        name, path = next(iter(spectra_on_disk.items()))
        browser.selected_spectrum = name
        browser.open_spectrum_viewer()
        [viewer] = browser.findChildren(SpectrumViewerDialog)
        assert viewer.spectrum_file_path == path
