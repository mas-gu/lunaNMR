# ABOUTME: Tests that the file manager recognises the NMRPipe 1D extension .ft1.
# ABOUTME: .ft1 was absent from supported_nmr_formats, so 1D spectra were invisible to discovery.

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from lunaNMR.utils.file_manager import NMRFileManager


class TestSupportedFormats:
    def test_nmrpipe_1d_extension_is_supported(self):
        """NMRPipe writes 1D spectra as .ft1; without it a folder of 1D data reads
        as empty while .ft and .ft2 are found."""
        assert 'ft1' in NMRFileManager().supported_nmr_formats

    def test_existing_formats_are_retained(self):
        formats = NMRFileManager().supported_nmr_formats
        for extension in ('ft', 'ser', 'ft2', 'ft3', 'pipe', 'ucsf'):
            assert extension in formats
