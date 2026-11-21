# ABOUTME: The spectrum files the dialogs list and browse for, read from NMRFileManager:
# ABOUTME: the supported extensions, and the QFileDialog filter entry that offers them.

from typing import List

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
