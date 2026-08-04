# ABOUTME: Tests the 1D spectrum plot - drawing peaks, including ones absent from a spectrum.
# ABOUTME: Runs headless via QT_QPA_PLATFORM=offscreen; drives the real widget, not mocks.

import os
import sys
from pathlib import Path

import numpy as np
import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

_MODULE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_MODULE_DIR))

pytest.importorskip('PySide6')


@pytest.fixture(scope='module')
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


@pytest.fixture
def plotter(app):
    from oned_loader import from_arrays
    from oned_plotter import OneDPlotter

    ppm = np.linspace(9.0, 7.0, 400)
    data = 1000.0 * np.exp(-((ppm - 8.2) ** 2) / (2 * 0.005 ** 2))

    widget = OneDPlotter()
    widget.set_spectrum(from_arrays(data, ppm, noise=1.0))
    return widget


def test_a_peak_missing_from_this_spectrum_is_still_drawn(plotter):
    """A picked peak that finds no maximum here has no height. It is drawn at the
    baseline rather than taking the whole plot down with it."""
    plotter.set_peaks([{'ppm': 8.2, 'height': None, 'assignment': 'A'}])

    plotter.figure.canvas.draw()

    label = plotter.axes.texts[0]
    assert label.xy == (8.2, 0)


def test_a_measured_peak_is_labelled_at_its_height(plotter):
    plotter.set_peaks([{'ppm': 8.2, 'height': 1000.0, 'assignment': 'A'}])

    plotter.figure.canvas.draw()

    assert plotter.axes.texts[0].xy == (8.2, 1000.0)
