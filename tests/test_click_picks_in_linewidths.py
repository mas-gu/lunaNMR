# ABOUTME: A click in a spectrum view selects the peak nearest in linewidths of the nucleus, not in raw ppm.
# ABOUTME: Real main window on DNAJA1 0o0 with its peak list, and the overlay viewer on a real 13C methyl series.

"""
In raw ppm a 15N ppm counts as a 1H ppm, though a 15N linewidth is ten times a 1H one.
A click a fifth of a 15N linewidth beside 25AlaH, on DNAJA1, then selected 50GlnH, which
sits 0.064 ppm away in 1H - one and a half 1H linewidths. The overlay viewer measured
in fixed 15N radii, so on a 13C spectrum a click half a 1H linewidth beside 278_cg2_D
selected 27_cg2_D, 0.1 ppm away in 13C - two thirds of a 13C linewidth.

The reach of a click is unchanged: any click that reached a peak in ppm still does.
"""

import os
import sys
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

DNAJA1 = Path(__file__).resolve().parents[2] / "data_test" / "data_binding" / "DNAJA1_HSPA8"
SPECTRUM = "2D_DNAJA1_0o0.ft"


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture
def window(app):
    if not DNAJA1.exists():
        pytest.skip("data_test not present")
    from lunaNMR.core.ps2d_config import set_ps2d_config
    from lunaNMR.gui.main_window import LunaNMRMainWindow
    set_ps2d_config('15N')
    w = LunaNMRMainWindow()
    w.load_nmr_file(str(DNAJA1 / SPECTRUM))
    w.load_peak_file(str(DNAJA1 / "DNAJA1.txt"))
    yield w
    app.processEvents()
    w.close()
    w.deleteLater()
    set_ps2d_config('15N')


def _listed(window, assignment):
    row = window.integrator.peak_list.set_index('Assignment').loc[assignment]
    return float(row['Position_X']), float(row['Position_Y'])


def _middle_click(window, x, y):
    """Select the peak a middle-click at (x, y) lands on; returns its assignment."""
    window._on_peak_select_for_move(x, y)
    selected = window.selected_peak_info
    window.selected_peak_info = None
    return selected and str(selected['assignment'])


class TestTheMainWindow:
    def test_a_click_a_fifth_of_a_15n_linewidth_beside_a_peak_selects_it(self, window):
        x, y = _listed(window, '25AlaH')
        assert _middle_click(window, x, y - 0.075) == '25AlaH'

    @pytest.mark.parametrize("dx, dy", [(0.14, 0.0), (-0.14, 0.0), (0.0, 0.14), (0.0, -0.14)])
    def test_a_click_that_reached_a_peak_in_ppm_still_does(self, window, dx, dy):
        x, y = _listed(window, '36AspH')
        assert _middle_click(window, x + dx, y + dy) == '36AspH'

    @pytest.mark.parametrize("dx", [0.16, -0.16])
    def test_a_click_beyond_the_reach_selects_nothing(self, window, dx):
        x, y = _listed(window, '36AspH')
        assert _middle_click(window, x + dx, y) is None


@pytest.fixture
def overlay(app, methyl_series, open_overlay_viewer):
    """The overlay viewer on the methyl series while the global config is 15N, as after
    the main window loaded a 15N spectrum."""
    from lunaNMR.core.ps2d_config import set_ps2d_config
    set_ps2d_config('15N')
    viewer = open_overlay_viewer(methyl_series)
    yield viewer
    viewer.deleteLater()


def _overlay_middle_click(viewer, x, y):
    viewer._on_overlay_middle_click(x, y)
    selected = viewer._edit_selected_peak
    viewer._edit_selected_peak = None
    return selected and str(selected['assignment'])


def _fitted(viewer, assignment):
    return viewer._peak_center(viewer._find_peak_obj(0, assignment))


class TestTheOverlayViewer:
    def test_a_click_half_a_1h_linewidth_beside_a_13c_peak_selects_it(self, overlay):
        assert overlay.spectra[0]['name'] == 'T2_3ms.ft' and overlay.spectra[0]['visible']
        x, y = _fitted(overlay, '278_cg2_D')
        assert _overlay_middle_click(overlay, x + 0.02, y) == '278_cg2_D'

    def test_a_click_that_reached_a_peak_in_ppm_still_does(self, overlay):
        x, y = _fitted(overlay, '278_cg2_D')
        assert _overlay_middle_click(overlay, x - 0.07, y) == '278_cg2_D'
