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


class TestOneExtensionList:
    """Three lists of spectrum extensions existed and disagreed: cli._discover_spectra
    carried a stale default tuple, spectra_check had its own, and only
    NMRFileManager.supported_nmr_formats had 'ft1'. The `series` path was safe because
    it passes the manager's list explicitly -- but `diagnose` used the third list and
    CLI_AGENT.md mirrored the first, so NMRPipe 1D spectra were invisible to both.
    """

    def test_spectra_check_uses_the_managers_list(self):
        from lunaNMR.validation.spectra_check import SPECTRUM_EXTENSIONS
        from lunaNMR.utils.file_manager import NMRFileManager
        assert set(SPECTRUM_EXTENSIONS) == set(NMRFileManager().supported_nmr_formats)

    def test_discover_spectra_takes_the_default_from_the_loader(self, tmp_path):
        """The literal default that used to sit here had fallen an extension behind,
        so a caller taking it silently missed NMRPipe 1D spectra."""
        from lunaNMR.cli import _discover_spectra
        from lunaNMR.utils.file_manager import NMRFileManager
        for ext in NMRFileManager().supported_nmr_formats:
            (tmp_path / f"s_300ms.{ext}").write_bytes(b"")
        found = {Path(f).suffix.lstrip('.') for f in _discover_spectra(str(tmp_path))}
        assert found == set(NMRFileManager().supported_nmr_formats)

    def test_ft1_is_discoverable(self, tmp_path):
        from lunaNMR.cli import _discover_spectra
        from lunaNMR.utils.file_manager import NMRFileManager
        (tmp_path / "a_300ms.ft1").write_bytes(b"")
        (tmp_path / "b_300ms.ft").write_bytes(b"")
        found = _discover_spectra(str(tmp_path), NMRFileManager().supported_nmr_formats)
        assert sorted(Path(f).suffix for f in found) == ['.ft', '.ft1']
