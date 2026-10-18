# ABOUTME: Both peak navigators colour a row as the spectrum marks that peak: by its quality label.
# ABOUTME: Real DNAJA1 series results; real PeakNavigator, PeakNavigatorTable and SpectrumPlotter offscreen.
"""
The navigators coloured rows from R-squared at one 0.7 threshold while the spectrum
markers coloured by the fit's quality label, so an Excellent peak was lime on the
spectrum and the navigator's green beside it, and a Poor one red against orange.
"""

import os
import sys

import pandas as pd
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture(scope="module")
def fitted_peaks(dnaja1_titration_series):
    """The last titration point's fitted peaks, as the GUI holds them."""
    tidy = pd.read_csv(dnaja1_titration_series / "series_analysis_tidy.csv",
                       dtype={'spectrum_name': str})
    last = tidy[tidy['spectrum_name'] == tidy['spectrum_name'].iloc[-1]]
    peaks = [{k: (None if pd.isna(v) else v) for k, v in row.items()}
             for row in last.to_dict('records')]
    assert len({p['quality'] for p in peaks}) >= 3, "too few quality labels to compare"
    return peaks


def _marker_rgb(app, peaks):
    from lunaNMR.gui.components.spectrum_plotter import SpectrumPlotter
    plotter = SpectrumPlotter()
    plotter._plot_detected_peaks(peaks)
    faces = plotter.ax.collections[-1].get_facecolors()
    assert len(faces) == len(peaks)
    return [tuple(round(c * 255) for c in face[:3]) for face in faces]


def _row_rgb(table, row):
    colour = table.item(row, 0).background().color()
    return colour.red(), colour.green(), colour.blue()


def _lightened(rgb):
    from PySide6.QtGui import QColor
    colour = QColor(*rgb).lighter(180)
    return colour.red(), colour.green(), colour.blue()


def test_the_overlay_navigator_colours_rows_as_the_markers(app, fitted_peaks):
    from lunaNMR.gui.components.peak_navigator_table import PeakNavigatorTable
    table = PeakNavigatorTable()
    table.load_peaks(fitted_peaks)

    markers = _marker_rgb(app, fitted_peaks)
    mismatched = [(p['assignment'], p['quality']) for row, (p, rgb) in
                  enumerate(zip(fitted_peaks, markers))
                  if _row_rgb(table.table, row) != _lightened(rgb)]
    assert not mismatched, f"rows coloured unlike their markers: {mismatched[:10]}"


def test_the_main_navigator_colours_rows_as_the_markers(app, fitted_peaks):
    from lunaNMR.gui.components.peak_navigator import PeakNavigator
    navigator = PeakNavigator()
    navigator.load_detected_peaks(fitted_peaks)
    navigator.peak_type_combo.setCurrentText("Detected Peaks")
    navigator.refresh_peak_list()

    markers = _marker_rgb(app, fitted_peaks)
    assert navigator.table.rowCount() == len(fitted_peaks)
    mismatched = [(p['assignment'], p['quality']) for row, (p, rgb) in
                  enumerate(zip(fitted_peaks, markers))
                  if _row_rgb(navigator.table, row) != _lightened(rgb)]
    assert not mismatched, f"rows coloured unlike their markers: {mismatched[:10]}"


def test_the_main_navigator_keeps_the_label_through_a_fit_update(app, fitted_peaks):
    """After a fit the main window updates rows by assignment; the label goes with it."""
    from lunaNMR.gui.components.peak_navigator import PeakNavigator
    navigator = PeakNavigator()
    navigator.load_detected_peaks([{k: v for k, v in p.items() if k != 'quality'}
                                   for p in fitted_peaks])
    navigator.update_heights_from_results(fitted_peaks)
    navigator.peak_type_combo.setCurrentText("Detected Peaks")
    navigator.refresh_peak_list()

    markers = _marker_rgb(app, fitted_peaks)
    mismatched = [(p['assignment'], p['quality']) for row, (p, rgb) in
                  enumerate(zip(fitted_peaks, markers))
                  if _row_rgb(navigator.table, row) != _lightened(rgb)]
    assert not mismatched, f"rows coloured unlike their markers: {mismatched[:10]}"
