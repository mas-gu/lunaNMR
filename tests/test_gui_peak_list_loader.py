# ABOUTME: The integrator loads a peak list (the GUI's path) through the file manager, as the
# ABOUTME: CLI does: a tab-separated list reads the same in both.

import contextlib
import io
from pathlib import Path

import pandas as pd
import pytest

PEAKS = Path(__file__).resolve().parents[1] / "data_example" / "2DHSQC" / "600_peaks.txt"

pytestmark = pytest.mark.skipif(not PEAKS.exists(), reason="data_example/2DHSQC not present")


def _through_the_gui_path(path):
    from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
    integrator = EnhancedVoigtIntegrator()
    with contextlib.redirect_stdout(io.StringIO()):
        assert integrator.load_peak_list_file(str(path)) is True
    return integrator.peak_list


def _through_the_cli_path(path):
    from lunaNMR.utils.file_manager import NMRFileManager
    return NMRFileManager().load_peak_list(str(path))


def test_a_tab_separated_list_has_its_columns(tmp_path):
    tabbed = tmp_path / "peaks.txt"
    pd.read_csv(PEAKS, skipinitialspace=True).to_csv(tabbed, sep="\t", index=False)
    loaded = _through_the_gui_path(tabbed)
    assert {"Assignment", "Position_X", "Position_Y"} <= set(loaded.columns)
    assert len(loaded) == len(_through_the_cli_path(tabbed))


def test_both_paths_read_the_same_list():
    gui, cli = _through_the_gui_path(PEAKS), _through_the_cli_path(PEAKS)
    for column in ("Assignment", "Position_X", "Position_Y"):
        assert list(gui[column]) == list(cli[column]), column
