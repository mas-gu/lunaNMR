# ABOUTME: The Expert "Show Ellipses" box draws each peak's PS2D fit and selector ellipses on the main spectrum.
# ABOUTME: Real LunaNMRMainWindow offscreen on the DNAJA1 spectrum, a 10-peak list cut from its peak list.

import os
import sys
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
SPECTRUM = "2D_DNAJA1_0o0.ft"

pytestmark = pytest.mark.skipif(not DNAJA1.exists(), reason="data_test not present")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def window(app, tmp_path, monkeypatch):
    """A window showing the spectrum and a 10-peak list, its Expert Mode dialog built."""
    from PySide6.QtWidgets import QDialog
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    set_ps2d_config('15N')
    w = LunaNMRMainWindow()
    lines = (DNAJA1 / "DNAJA1.txt").read_text().splitlines()
    peaks = tmp_path / "ten.txt"
    peaks.write_text("\n".join(lines[:11]) + "\n")
    w.load_data(str(DNAJA1 / SPECTRUM), str(peaks))
    monkeypatch.setattr(QDialog, "exec", lambda self: 0)
    w.on_expert_mode()
    yield w
    w.close()
    w.deleteLater()
    set_ps2d_config('15N')


def _ellipses(window):
    """(width, height) of every ellipse drawn on the main spectrum, in ppm."""
    from matplotlib.patches import Ellipse
    return sorted((round(patch.width, 6), round(patch.height, 6))
                  for patch in window.spectrum_plotter.ax.patches if isinstance(patch, Ellipse))


def test_no_ellipse_is_drawn_until_the_box_is_ticked(window):
    assert _ellipses(window) == []


def test_ticking_the_box_draws_the_fit_and_selector_ellipse_of_every_peak(window):
    from lunaNMR.core.ps2d_config import get_ps2d_config
    config = get_ps2d_config()

    window.expert_show_ellipses_check.setChecked(True)

    fit = (round(2 * config.radF2, 6), round(2 * config.radF1, 6))
    selector = (round(2 * config.radF2_selector, 6), round(2 * config.radF1_selector, 6))
    assert _ellipses(window) == sorted([fit] * 10 + [selector] * 10)


def test_unticking_the_box_removes_them(window):
    window.expert_show_ellipses_check.setChecked(True)

    window.expert_show_ellipses_check.setChecked(False)

    assert _ellipses(window) == []
