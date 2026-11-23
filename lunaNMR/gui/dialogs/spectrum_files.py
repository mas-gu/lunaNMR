# ABOUTME: The spectrum files the dialogs list and browse for, read from NMRFileManager (the supported
# ABOUTME: extensions, the QFileDialog filter entry that offers them), and where a series' spectra are on disk.

import os
from typing import List, Optional

from lunaNMR.utils.file_manager import NMRFileManager


def spectrum_extensions() -> List[str]:
    """The spectrum file extensions the application supports, without the dot."""
    return list(NMRFileManager().supported_nmr_formats)


def spectrum_file_filter() -> str:
    """The QFileDialog filter entry offering every spectrum file the application
    supports: its extensions, and the Bruker processed-data names, which carry none."""
    patterns = [f"*.{ext}" for ext in spectrum_extensions()]
    patterns += NMRFileManager.BRUKER_PDATA_NAMES
    return f"NMR Files ({' '.join(patterns)})"


def series_data_folder(batch_results, fallback_folder: Optional[str] = None) -> Optional[str]:
    """The folder a series' spectra are in: the one the series records, else
    `fallback_folder`. The series is a BatchResults, or its contents as a dict."""
    if isinstance(batch_results, dict):
        metadata = batch_results.get('metadata', {})
    else:
        metadata = getattr(batch_results, 'metadata', None) or {}
    return metadata.get('data_folder') or fallback_folder


def series_spectrum_file(data_folder: Optional[str], spectrum_name: str, result: dict) -> str:
    """Where a spectrum of a series is on disk: the path its result records when that is
    absolute, else its name in the series' data folder (series_data_folder). With no
    folder known it is what the result records."""
    recorded = result.get('spectrum_file', '')
    if os.path.isabs(recorded) or not data_folder:
        return recorded
    return os.path.join(data_folder, spectrum_name)
