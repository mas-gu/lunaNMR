# ABOUTME: Tests the overlay viewer drawing one peak list per visible spectrum.
# ABOUTME: Each spectrum's peaks render in its own color, gated by visible + show_peaks.

import sys
from pathlib import Path

import pytest

_SRC = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_SRC))

from PySide6.QtWidgets import QApplication


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication(sys.argv)


def _drawn(axes):
    """(colour name, point count) of every marker set drawn on the real axes. Colours
    are compared without alpha, which the markers set for themselves."""
    from matplotlib.colors import to_rgb
    names = {to_rgb(c): c for c in ('red', 'green')}
    return [(names[tuple(coll.get_facecolor()[0][:3])], len(coll.get_offsets()))
            for coll in axes.collections]


def _viewer():
    from lunaNMR.gui.dialogs.multi_spectrum_viewer_dialog import MultiSpectrumViewerDialog
    o = MultiSpectrumViewerDialog.__new__(MultiSpectrumViewerDialog)
    o.spectra = [
        {'name': 'A', 'color': 'red', 'visible': True, 'show_peaks': True,
         'fitted_peaks': [{'ppm_x': 8.0, 'ppm_y': 120.0, 'assignment': 'R1'}]},
        {'name': 'B', 'color': 'green', 'visible': True, 'show_peaks': True,
         'fitted_peaks': [{'ppm_x': 7.0, 'ppm_y': 118.0, 'assignment': 'R2'}]},
    ]
    o.overlay_selected_spectrum_idx = 0
    o.edit_mode = False
    # A real matplotlib Axes in place of the dialog's plot widget, which this test
    # does not build: the markers are drawn on it, and read back from it.
    import types
    from matplotlib.figure import Figure
    o.plot_widget = types.SimpleNamespace(axes=Figure().add_subplot())
    return o


def test_each_visible_spectrum_draws_its_list_in_its_color(app):
    o = _viewer()
    o._plot_overlay_peak_markers()
    assert _drawn(o.plot_widget.axes) == [('red', 1), ('green', 1)]


def test_per_spectrum_toggle_hides_one_list(app):
    o = _viewer()
    o.spectra[1]['show_peaks'] = False
    o._plot_overlay_peak_markers()
    assert _drawn(o.plot_widget.axes) == [('red', 1)]   # only spectrum A


def test_hidden_spectrum_draws_no_peaks(app):
    o = _viewer()
    o.spectra[0]['visible'] = False
    o._plot_overlay_peak_markers()
    assert _drawn(o.plot_widget.axes) == [('green', 1)]  # only visible spectrum B
