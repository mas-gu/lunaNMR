# ABOUTME: Tests the series results table - spectra against peaks, value toggle, copy and CSV export.
# ABOUTME: Runs headless via QT_QPA_PLATFORM=offscreen and drives the real widgets.

import os
import sys
from pathlib import Path

import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

_MODULE_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_MODULE_DIR))

pytest.importorskip('PySide6')


@pytest.fixture(scope='module')
def app():
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def _rows():
    """Two peaks over three spectra, with one point where a peak was not found."""
    def row(point, spectrum, assignment, height, area, ppm, matched=True):
        return {'point': point, 'spectrum': spectrum, 'assignment': assignment,
                'height': height, 'area': area, 'ppm': ppm, 'matched': matched,
                'fit_area': None if area is None else area * 1.1,
                'r_squared': None if area is None else 0.995,
                'fwhm': None if area is None else 0.0026}

    return [
        row(0, 's1', 'A', 10.0, 1.0, 8.1847),
        row(0, 's1', 'B', 20.0, 2.0, 8.1752),
        row(1, 's2', 'A', 30.0, 3.0, 8.1848),
        row(1, 's2', 'B', 40.0, 4.0, 8.1753),
        row(2, 's3', 'A', 50.0, 5.0, 8.1849),
        row(2, 's3', 'B', None, None, None, matched=False),
    ]


@pytest.fixture
def dialog(app):
    from oned_series_table import SeriesTableDialog
    return SeriesTableDialog(_rows())


class TestLayout:
    def test_rows_are_spectra_and_columns_are_peaks(self, dialog):
        assert dialog.table.rowCount() == 3
        assert dialog.table.columnCount() == 2
        assert [dialog.table.horizontalHeaderItem(c).text()
                for c in range(2)] == ['A', 'B']
        assert [dialog.table.verticalHeaderItem(r).text()
                for r in range(3)] == ['s1', 's2', 's3']

    def test_shows_intensity_first(self, dialog):
        assert dialog.table.item(0, 0).text() == '1.000e+01'

    def test_missing_measurement_is_blank_not_zero(self, dialog):
        assert dialog.table.item(2, 1).text() == ''


class TestValueToggle:
    def test_switching_to_ppm(self, dialog):
        dialog.set_value('ppm')
        assert dialog.table.item(0, 0).text() == '8.1847'
        assert dialog.table.item(1, 1).text() == '8.1753'

    def test_switching_to_area(self, dialog):
        dialog.set_value('area')
        assert dialog.table.item(2, 0).text() == '5.000e+00'

    def test_switching_back_to_intensity(self, dialog):
        dialog.set_value('ppm')
        dialog.set_value('height')
        assert dialog.table.item(0, 1).text() == '2.000e+01'

    def test_switching_to_the_fitted_area(self, dialog):
        dialog.set_value('fit_area')
        assert dialog.table.item(0, 0).text() == '1.100e+00'

    def test_fit_quality_reads_as_a_plain_fraction(self, dialog):
        """R2 and linewidth in scientific notation are unreadable."""
        dialog.set_value('r_squared')
        assert dialog.table.item(0, 0).text() == '0.9950'

    def test_linewidth_reads_as_a_plain_number(self, dialog):
        dialog.set_value('fwhm')
        assert dialog.table.item(0, 0).text() == '0.0026'

    def test_every_observable_is_offered(self, dialog):
        from oned_series_table import VALUES
        keys = [key for key, _label in VALUES]
        assert keys == ['height', 'area', 'fit_area', 'r_squared', 'fwhm', 'ppm']
        assert dialog.value_box.count() == len(VALUES)

    def test_toggle_widget_drives_the_table(self, dialog):
        index = dialog.value_box.findText('Position (ppm)')
        dialog.value_box.setCurrentIndex(index)
        assert dialog.table.item(0, 0).text() == '8.1847'


class TestCopy:
    def test_copies_the_whole_table_as_tab_separated_text(self, dialog):
        # only the final newline is stripped - a bare strip() would eat the trailing
        # tab that marks the last row's empty cell
        text = dialog.table_as_text()
        lines = text.rstrip('\n').split('\n')
        assert lines[0].split('\t') == ['spectrum', 'A', 'B']
        assert lines[1].split('\t') == ['s1', '1.000e+01', '2.000e+01']
        assert lines[3].split('\t') == ['s3', '5.000e+01', '']

    def test_copy_puts_it_on_the_clipboard(self, dialog, app):
        dialog.copy_all()
        assert app.clipboard().text().startswith('spectrum\tA\tB')

    def test_copy_follows_the_selected_value(self, dialog):
        dialog.set_value('ppm')
        assert '8.1847' in dialog.table_as_text()

    def test_copy_selection_returns_only_the_chosen_cells(self, dialog):
        dialog.table.setRangeSelected(
            _range(dialog, 0, 0, 1, 0), True)
        text = dialog.selection_as_text()
        assert text.strip().split('\n') == ['1.000e+01', '3.000e+01']


def _range(dialog, top, left, bottom, right):
    from PySide6.QtWidgets import QTableWidgetSelectionRange
    return QTableWidgetSelectionRange(top, left, bottom, right)


class TestCsvExport:
    def test_writes_the_displayed_value(self, dialog, tmp_path, monkeypatch):
        from PySide6.QtWidgets import QMessageBox
        monkeypatch.setattr(QMessageBox, 'information', lambda *a, **k: None)

        out = tmp_path / "series.csv"
        dialog.export_csv(path=str(out))
        lines = out.read_text().strip().splitlines()
        assert lines[0] == 'spectrum,A,B'
        assert lines[1].startswith('s1,10.0')

    def test_exporting_ppm_writes_positions(self, dialog, tmp_path, monkeypatch):
        from PySide6.QtWidgets import QMessageBox
        monkeypatch.setattr(QMessageBox, 'information', lambda *a, **k: None)

        dialog.set_value('ppm')
        out = tmp_path / "ppm.csv"
        dialog.export_csv(path=str(out))
        assert '8.1847' in out.read_text()

    def test_missing_value_stays_an_empty_cell(self, dialog, tmp_path, monkeypatch):
        from PySide6.QtWidgets import QMessageBox
        monkeypatch.setattr(QMessageBox, 'information', lambda *a, **k: None)

        out = tmp_path / "gap.csv"
        dialog.export_csv(path=str(out))
        assert out.read_text().strip().splitlines()[-1] == 's3,50.0,'
