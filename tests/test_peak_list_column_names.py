# ABOUTME: A peak list whose columns carry the alternative names the file manager maps
# ABOUTME: (1H / 15N / label, x / y, f2 / f1 ...) is accepted and read under the standard names.

from pathlib import Path

import pandas as pd
import pytest

from lunaNMR.utils.file_manager import NMRFileManager

PEAKS = Path(__file__).resolve().parents[1] / "data_example" / "2DHSQC" / "600_peaks.txt"

pytestmark = pytest.mark.skipif(not PEAKS.exists(), reason="data_example/2DHSQC not present")


def _renamed(tmp_path, columns, name="peaks.csv", sep=","):
    listed = pd.read_csv(PEAKS, skipinitialspace=True)
    listed.columns = columns
    path = tmp_path / name
    listed.to_csv(path, index=False, sep=sep)
    return path


@pytest.mark.parametrize("columns", [("label", "1H", "15N"), ("assign", "x", "y"),
                                     ("Assignment", "f2", "f1"), ("assignment", "ppm_x", "ppm_y")])
def test_alternative_column_names_validate_and_load(tmp_path, columns):
    path = _renamed(tmp_path, list(columns))
    manager = NMRFileManager()
    valid, message = manager.validate_peak_file(str(path))
    assert valid, message
    loaded = manager.load_peak_list(str(path))
    standard = manager.load_peak_list(str(PEAKS))
    for column in ("Assignment", "Position_X", "Position_Y"):
        assert list(loaded[column]) == list(standard[column]), column


def test_a_tab_separated_list_with_alternative_names_loads(tmp_path):
    path = _renamed(tmp_path, ["label", "1H", "15N"], name="peaks.txt", sep="\t")
    loaded = NMRFileManager().load_peak_list(str(path))
    assert {"Assignment", "Position_X", "Position_Y"} <= set(loaded.columns)


def test_a_list_without_positions_is_still_refused(tmp_path):
    path = _renamed(tmp_path, ["Assignment", "first", "second"])
    valid, message = NMRFileManager().validate_peak_file(str(path))
    assert not valid
    assert "Position_X" in message and "Position_Y" in message
    with pytest.raises(RuntimeError):
        NMRFileManager().load_peak_list(str(path))


def test_validation_and_loading_read_the_file_the_same_way(tmp_path):
    """One reader behind both: a list validation accepts is a list loading returns."""
    path = _renamed(tmp_path, ["label", "1H", "15N"], name="peaks.txt", sep="\t")
    manager = NMRFileManager()
    valid, message = manager.validate_peak_file(str(path))
    assert valid, message
    assert f"{len(manager.load_peak_list(str(path)))} peaks" in message
