# ABOUTME: A series refuses a peak list in which two peaks share an assignment, naming them.
# ABOUTME: Every series join is keyed by assignment, so one peak's values were copied onto the other.

"""
Clusters, the cascade seed, the matrices and the tidy CSV all look a peak up by its
assignment. Two peaks under one label came out as one peak's heights written on both
rows, with exit 0 and no warning.
"""

from pathlib import Path

import pytest

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


@pytest.fixture
def duplicated_list(tmp_path):
    """The real DNAJA1 peak list with its 4GluH row written twice."""
    lines = (DNAJA1 / "DNAJA1.txt").read_text().splitlines()
    row = next(line for line in lines if line.startswith("4GluH,"))
    peaks = tmp_path / "DNAJA1_dup.txt"
    peaks.write_text("\n".join(lines + [row]) + "\n")
    return peaks


def test_the_cli_refuses_it_and_names_the_duplicate(duplicated_list, tmp_path, capsys):
    from lunaNMR.cli import main
    code = main(["series", "--spectra", str(DNAJA1), "--peaks", str(duplicated_list),
                 "--out", str(tmp_path / "o"), "--mode", "titration"])
    assert code == 1
    err = capsys.readouterr().err
    assert "4GluH" in err and "duplicate" in err.lower()
    assert not (tmp_path / "o").exists()


def test_the_series_engine_refuses_it_too(duplicated_list, tmp_path):
    """The GUI series dialog calls the engine directly, so the engine must refuse."""
    from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor
    from lunaNMR.utils.file_manager import NMRFileManager
    from lunaNMR.cli import _default_series_params
    peaks = NMRFileManager().load_peak_list(str(duplicated_list))
    spectra = sorted(str(p) for p in DNAJA1.glob("*.ft"))
    with pytest.raises(ValueError, match="4GluH"):
        MultiSpectrumProcessor(_default_series_params()).process_nmr_series(
            spectra, peaks, str(tmp_path / "o"), series_mode="titration")
