# ABOUTME: Tests that the dialogs which list or browse for spectra offer the spectrum files
# ABOUTME: NMRFileManager supports - its extensions and Bruker processed-data names - not a list of their own.

import os
import sys
import types

import pytest

from PySide6.QtWidgets import QApplication


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication(sys.argv)


def _supported_extensions():
    from lunaNMR.utils.file_manager import NMRFileManager
    return list(NMRFileManager().supported_nmr_formats)


def _supported_patterns():
    from lunaNMR.utils.file_manager import NMRFileManager
    return {f"*.{ext}" for ext in _supported_extensions()} | set(NMRFileManager.BRUKER_PDATA_NAMES)


def _patterns(file_filter, entry="NMR Files"):
    """The patterns of one named entry of a QFileDialog filter string."""
    for part in file_filter.split(";;"):
        name, _, patterns = part.partition(" (")
        if name == entry:
            return set(patterns.rstrip(")").split())
    raise AssertionError(f"no '{entry}' entry in: {file_filter}")


@pytest.fixture
def asked_filters(monkeypatch):
    """The filter of every file dialog opened, each dialog cancelled at once."""
    from PySide6.QtWidgets import QFileDialog
    asked = []

    def cancelled(parent=None, caption='', directory='', file_filter='', *args, **kwargs):
        asked.append(file_filter)
        return '', ''

    def cancelled_many(parent=None, caption='', directory='', file_filter='', *args, **kwargs):
        asked.append(file_filter)
        return [], ''
    monkeypatch.setattr(QFileDialog, 'getOpenFileName', cancelled)
    monkeypatch.setattr(QFileDialog, 'getOpenFileNames', cancelled_many)
    return asked


def test_the_series_dialog_lists_every_supported_spectrum_of_its_folder(app, tmp_path):
    from lunaNMR.gui.dialogs.series_integration_dialog import SeriesIntegrationDialog
    for ext in _supported_extensions():
        (tmp_path / f"point_1.{ext}").write_bytes(b"")
    (tmp_path / "notes.txt").write_text("")
    dialog = SeriesIntegrationDialog(
        main_window=types.SimpleNamespace(current_nmr_folder=str(tmp_path)))
    try:
        listed = sorted(os.path.basename(f) for f in dialog.nmr_files)
        assert listed == sorted(f"point_1.{ext}" for ext in _supported_extensions())
    finally:
        dialog.deleteLater()


def test_the_load_dialog_browses_for_every_supported_spectrum(app, asked_filters):
    from lunaNMR.gui.dialogs.data_loading_dialog import DataLoadingDialog
    dialog = DataLoadingDialog()
    try:
        dialog.browse_nmr()
        assert _patterns(asked_filters[-1]) == _supported_patterns()
    finally:
        dialog.deleteLater()


def test_the_missing_files_dialog_browses_for_every_supported_spectrum(app, asked_filters):
    from lunaNMR.gui.dialogs.missing_files_dialog import MissingFilesDialog
    dialog = MissingFilesDialog()
    try:
        for ext in _supported_extensions():
            dialog.browse_for_file(f"/gone/spectrum.{ext}")
            assert _patterns(asked_filters[-1]) == _supported_patterns(), ext
    finally:
        dialog.deleteLater()


def test_the_inspector_browses_for_every_supported_spectrum(app, asked_filters):
    from lunaNMR.gui.dialogs.spectral_inspector import SpectralInspector
    inspector = SpectralInspector(parent=None)
    try:
        inspector._library_panel._on_load_spectrum()
        assert _patterns(asked_filters[-1]) == _supported_patterns()
    finally:
        inspector.deleteLater()
