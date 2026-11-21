# ABOUTME: Tests the Results Browser on a real series: positions, linewidths, height, volume and
# ABOUTME: detection shown per peak and per spectrum must be those the series wrote to its tracking CSV.

import json
import math
import sys

import pandas as pd
import pytest

from PySide6.QtWidgets import QApplication


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication(sys.argv)


@pytest.fixture(scope="module")
def series(dnaja1_titration_series):
    """The stored series, its tracking table by assignment, and {spectrum: point label}."""
    from lunaNMR.utils.project_manager import batch_results_from_dict
    folder = dnaja1_titration_series
    batch = batch_results_from_dict(json.loads((folder / "batch_results.json").read_text()))
    tracking = pd.read_csv(folder / "comprehensive_peak_tracking.csv",
                           dtype={'Assignment': str}).set_index('Assignment')
    columns = json.loads((folder / "series_metadata.json").read_text())['columns']
    return batch, tracking, {c['spectrum']: str(c['column']) for c in columns}


@pytest.fixture
def browser(app, series):
    from lunaNMR.gui.dialogs.results_browser_dialog import ResultsBrowserDialog
    dialog = ResultsBrowserDialog(batch_results=series[0])
    yield dialog
    dialog.deleteLater()


def _cells(table, row):
    return [table.item(row, col).text() for col in range(table.columnCount())]


def _select_spectrum(browser, name):
    browser.spectrum_combo.setCurrentText(name)
    browser.on_spectrum_selected(browser.spectrum_combo.currentIndex())


def _select_peak(browser, assignment):
    browser.peak_combo.setCurrentText(assignment)
    browser.on_peak_selected(browser.peak_combo.currentIndex())


SPECTRUM = '2D_DNAJA1_HSPA8_3o0.ft'


class TestDetailsOfASpectrum:
    def test_each_peak_shows_the_position_and_linewidth_the_series_wrote(self, browser, series):
        _batch, tracking, label_of = series
        label = label_of[SPECTRUM]
        _select_spectrum(browser, SPECTRUM)
        table = browser.details_table
        assert table.rowCount() == len(tracking)
        for row in range(table.rowCount()):
            assignment, pos_f2, pos_f1, _height, _volume, lw_f2, lw_f1, _r2 = _cells(table, row)
            written = tracking.loc[assignment]
            assert math.isclose(float(pos_f2), written[f'{label}_Position_X'], abs_tol=1e-4)
            assert math.isclose(float(pos_f1), written[f'{label}_Position_Y'], abs_tol=1e-3)
            assert math.isclose(float(lw_f2), written[f'{label}_LW_X'], abs_tol=1e-4)
            assert math.isclose(float(lw_f1), written[f'{label}_LW_Y'], abs_tol=1e-3)

    def test_intensity_is_the_height_and_volume_the_volume(self, browser, series):
        _batch, tracking, label_of = series
        label = label_of[SPECTRUM]
        _select_spectrum(browser, SPECTRUM)
        table = browser.details_table
        for row in range(table.rowCount()):
            cells = _cells(table, row)
            written = tracking.loc[cells[0]]
            assert math.isclose(float(cells[3]), written[f'{label}_Height'], rel_tol=1e-2, abs_tol=1.0)
            assert math.isclose(float(cells[4]), written[f'{label}_Volume'], rel_tol=1e-2, abs_tol=1.0)

    def test_the_detection_rate_is_the_share_of_peaks_the_series_marked_detected(self, browser, series):
        _batch, tracking, label_of = series
        detected = tracking[f'{label_of[SPECTRUM]}_Detected']
        assert 0 < detected.sum() < len(detected)
        _select_spectrum(browser, SPECTRUM)
        assert browser.spectrum_peaks_label.text() == f"Peaks: {int(detected.sum())}/{len(detected)}"
        assert browser.spectrum_detection_label.text() == f"Detection Rate: {detected.mean():.1%}"


class TestOverview:
    def test_the_average_detection_rate_is_the_mean_over_the_spectra(self, browser, series):
        _batch, tracking, label_of = series
        rates = [tracking[f'{label}_Detected'].mean() for label in label_of.values()]
        assert browser.avg_detection_label.text() == f"Avg Detection Rate: {sum(rates) / len(rates):.1%}"


class TestTrackingOfAPeak:
    def test_each_point_shows_the_position_height_and_volume_the_series_wrote(self, browser, series):
        _batch, tracking, label_of = series
        _select_peak(browser, '42GlyH')
        table = browser.tracking_table
        assert table.rowCount() == len(label_of)
        written = tracking.loc['42GlyH']
        for row in range(table.rowCount()):
            name, pos_f2, pos_f1, height, volume, _r2, _status = _cells(table, row)
            label = label_of[name]
            assert math.isclose(float(pos_f2), written[f'{label}_Position_X'], abs_tol=1e-4)
            assert math.isclose(float(pos_f1), written[f'{label}_Position_Y'], abs_tol=1e-3)
            assert math.isclose(float(height), written[f'{label}_Height'], rel_tol=1e-2)
            assert math.isclose(float(volume), written[f'{label}_Volume'], rel_tol=1e-2)


class TestEvolutionOfAPeak:
    @pytest.mark.parametrize("shown, column", [("Intensity", "Height"), ("Volume", "Volume"),
                                               ("Position F2", "Position_X"),
                                               ("Position F1", "Position_Y")])
    def test_the_plot_follows_what_the_series_wrote(self, browser, series, shown, column):
        _batch, tracking, label_of = series
        browser.evolution_peak_combo.setCurrentText('42GlyH')
        browser.property_combo.setCurrentText(shown)
        browser.on_evolution_peak_selected()
        plotted = browser.evolution_chart.figure.axes[0].lines[0].get_ydata()
        written = [tracking.loc['42GlyH', f'{label_of[name]}_{column}'] for name in sorted(label_of)]
        assert len(plotted) == len(written)
        for value, expected in zip(plotted, written):
            assert math.isclose(value, expected, rel_tol=1e-4)
