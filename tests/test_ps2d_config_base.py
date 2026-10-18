# ABOUTME: Every fit starts from the base PS2D config - nucleus defaults plus the user's own settings.
# ABOUTME: The radii a run's adaptive search chose stay for refits afterwards but never reach the next run.

"""
The adaptive search writes the radii and overlap thresholds it chose into the global PS2D
config, and nothing put them back. The next run in the same process started from them:
two identical DNAJA1 series run back to back gave 207 rows, then 204 - the smaller radii
no longer reached the peak a dummy is auto-added beside.

The BASE is what the nucleus selector and the expert panel set. A fresh fit (no clusters
or statistics inherited from a reference spectrum) restores it before it starts.
"""

import os
import sys
import types
from pathlib import Path

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
POINTS = ("2D_DNAJA1_0o0.ft", "2D_DNAJA1_HSPA8_0o125.ft", "2D_DNAJA1_HSPA8_0o25.ft")


@pytest.fixture(autouse=True)
def fresh_config():
    from lunaNMR.core.ps2d_config import set_ps2d_config
    set_ps2d_config('15N')
    yield
    set_ps2d_config('15N')


def test_restoring_the_base_undoes_what_a_run_tuned():
    from lunaNMR.core.ps2d_config import get_ps2d_config, restore_base_ps2d_config
    config = get_ps2d_config()
    config.radF1, config.overlap_threshold_x = 0.158, 0.0258
    restore_base_ps2d_config()
    assert (config.radF1, config.overlap_threshold_x) == (0.4, 0.04)


def test_the_base_keeps_the_users_own_settings():
    from lunaNMR.core.ps2d_config import (get_ps2d_config, mark_ps2d_config_as_base,
                                          restore_base_ps2d_config)
    config = get_ps2d_config()
    config.radF1 = 0.3
    mark_ps2d_config_as_base()
    config.radF1 = 0.158
    restore_base_ps2d_config()
    assert config.radF1 == 0.3


def test_restoring_keeps_the_object_every_caller_holds():
    from lunaNMR.core.ps2d_config import get_ps2d_config, restore_base_ps2d_config
    held = get_ps2d_config()
    restore_base_ps2d_config()
    assert get_ps2d_config() is held


def test_the_expert_panel_sets_the_base():
    from PySide6.QtWidgets import QApplication, QDoubleSpinBox, QSpinBox
    QApplication.instance() or QApplication(sys.argv)
    from lunaNMR.core.ps2d_config import get_ps2d_config, restore_base_ps2d_config
    from lunaNMR.gui.main_window import LunaNMRMainWindow

    def spin(cls, value, hi):
        box = cls()
        box.setRange(0, hi)
        if cls is QDoubleSpinBox:
            box.setDecimals(4)
        box.setValue(value)
        return box
    host = types.SimpleNamespace(
        expert_radF1_spin=spin(QDoubleSpinBox, 0.3, 2), expert_radF2_spin=spin(QDoubleSpinBox, 0.03, 1),
        expert_max_iter_spin=spin(QSpinBox, 500, 10000), expert_max_cluster_spin=spin(QSpinBox, 8, 100),
        expert_overlap_x_spin=spin(QDoubleSpinBox, 0.03, 1), expert_overlap_y_spin=spin(QDoubleSpinBox, 0.3, 2),
        update_status=lambda msg: None, update_spectrum_plot=lambda: None)
    LunaNMRMainWindow._on_ps2d_params_apply(host)
    get_ps2d_config().radF1 = 0.158
    restore_base_ps2d_config()
    assert get_ps2d_config().radF1 == 0.3


@pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")
def test_two_identical_runs_give_the_same_series(tmp_path):
    from lunaNMR.cli import main
    (tmp_path / "spectra").mkdir()
    for name in POINTS:
        (tmp_path / "spectra" / name).symlink_to(DNAJA1 / name)
    runs = []
    for out in ("first", "second"):
        assert main(["series", "--spectra", str(tmp_path / "spectra"),
                     "--peaks", str(DNAJA1 / "DNAJA1.txt"), "--mode", "titration",
                     "--peak-source", "cascade", "--out", str(tmp_path / out)]) == 0
        runs.append(pd.read_csv(tmp_path / out / "series_analysis_tidy.csv"))
    pd.testing.assert_frame_equal(runs[0], runs[1])
