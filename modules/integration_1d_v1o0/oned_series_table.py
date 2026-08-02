# ABOUTME: Results table for a 1D series run - spectra down, peaks across, value toggle.
# ABOUTME: Copies as tab-separated text for pasting into a spreadsheet, and exports CSV.

from pathlib import Path

from PySide6.QtGui import QGuiApplication, QKeySequence, QShortcut
from PySide6.QtWidgets import (QAbstractItemView, QComboBox, QDialog, QFileDialog,
                               QHBoxLayout, QLabel, QMessageBox, QPushButton,
                               QTableWidget, QTableWidgetItem, QVBoxLayout)

from oned_series import write_series_csv

# Selectable observables: the key in a series row, and how it reads in the toggle.
VALUES = (('height', 'Intensity (height)'),
          ('area', 'Area'),
          ('ppm', 'Position (ppm)'))

# Positions want fixed decimals; intensities and areas span orders of magnitude.
PPM_FORMAT = '{:.4f}'
MAGNITUDE_FORMAT = '{:.3e}'


class SeriesTableDialog(QDialog):
    """Series results as a spectrum-by-peak grid.

    A missing measurement is left blank rather than shown as zero: a peak that could
    not be found is not a measurement of nothing.
    """

    def __init__(self, rows, parent=None, default_value='height'):
        super().__init__(parent)
        self.setWindowTitle("Series results")
        self.resize(760, 560)

        self.rows = list(rows)
        self.value = default_value

        self.assignments = []
        self.points = []
        for row in self.rows:
            if row['assignment'] not in self.assignments:
                self.assignments.append(row['assignment'])
            if row['point'] not in self.points:
                self.points.append(row['point'])
        self.points.sort()

        self._build_ui()
        self._fill()

    # ---------------------------------------------------------------- ui

    def _build_ui(self):
        self.value_box = QComboBox()
        self.value_box.addItems([label for _key, label in VALUES])
        self.value_box.setCurrentIndex(
            [key for key, _label in VALUES].index(self.value))
        self.value_box.currentIndexChanged.connect(self._value_box_changed)

        self.table = QTableWidget()
        self.table.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)

        copy_button = QPushButton("Copy all")
        copy_button.clicked.connect(self.copy_all)

        export_button = QPushButton("Integrate series → CSV")
        export_button.clicked.connect(self.export_csv)

        close_button = QPushButton("Close")
        close_button.clicked.connect(self.accept)

        QShortcut(QKeySequence.Copy, self.table, self.copy_selection)

        controls = QHBoxLayout()
        controls.addWidget(QLabel("Show:"))
        controls.addWidget(self.value_box)
        controls.addStretch()
        controls.addWidget(copy_button)
        controls.addWidget(export_button)
        controls.addWidget(close_button)

        self.hint = QLabel("Select cells and press Ctrl+C, or use Copy all, "
                           "to paste into a spreadsheet.")
        self.hint.setStyleSheet("color: #666; font-size: 10px;")

        layout = QVBoxLayout(self)
        layout.addLayout(controls)
        layout.addWidget(self.table)
        layout.addWidget(self.hint)

    # ------------------------------------------------------------- values

    def _value_box_changed(self, index):
        self.set_value(VALUES[index][0])

    def set_value(self, key):
        """Show `key` - one of the row keys in VALUES - in every cell."""
        self.value = key
        index = [k for k, _label in VALUES].index(key)
        if self.value_box.currentIndex() != index:
            self.value_box.blockSignals(True)
            self.value_box.setCurrentIndex(index)
            self.value_box.blockSignals(False)
        self._fill()

    def _format(self, value):
        if value is None:
            return ''
        template = PPM_FORMAT if self.value == 'ppm' else MAGNITUDE_FORMAT
        return template.format(value)

    def _lookup(self):
        """{(point, assignment): value} for the currently selected observable."""
        return {(row['point'], row['assignment']): row.get(self.value)
                for row in self.rows}

    def _spectrum_names(self):
        names = {row['point']: row['spectrum'] for row in self.rows}
        return [names.get(point, str(point)) for point in self.points]

    def _fill(self):
        values = self._lookup()

        self.table.setRowCount(len(self.points))
        self.table.setColumnCount(len(self.assignments))
        self.table.setHorizontalHeaderLabels([str(a) for a in self.assignments])
        self.table.setVerticalHeaderLabels(self._spectrum_names())

        for r, point in enumerate(self.points):
            for c, assignment in enumerate(self.assignments):
                self.table.setItem(r, c, QTableWidgetItem(
                    self._format(values.get((point, assignment)))))

        self.table.resizeColumnsToContents()

    # -------------------------------------------------------------- copy

    def table_as_text(self):
        """The whole grid as tab-separated text, headers included."""
        lines = ['\t'.join(['spectrum'] + [str(a) for a in self.assignments])]

        for r, name in enumerate(self._spectrum_names()):
            cells = [self.table.item(r, c).text() if self.table.item(r, c) else ''
                     for c in range(self.table.columnCount())]
            lines.append('\t'.join([name] + cells))

        return '\n'.join(lines) + '\n'

    def selection_as_text(self):
        """Only the selected cells, laid out as they appear."""
        ranges = self.table.selectedRanges()
        if not ranges:
            return self.table_as_text()

        lines = []
        for selected in ranges:
            for r in range(selected.topRow(), selected.bottomRow() + 1):
                cells = []
                for c in range(selected.leftColumn(), selected.rightColumn() + 1):
                    item = self.table.item(r, c)
                    cells.append(item.text() if item else '')
                lines.append('\t'.join(cells))

        return '\n'.join(lines) + '\n'

    def copy_all(self):
        QGuiApplication.clipboard().setText(self.table_as_text())

    def copy_selection(self):
        QGuiApplication.clipboard().setText(self.selection_as_text())

    # ------------------------------------------------------------ export

    def export_csv(self, _=None, path=None):
        path = path or QFileDialog.getSaveFileName(
            self, "Save series table", str(Path.cwd() / "series_results.csv"),
            "CSV (*.csv)")[0]
        if not path:
            return None

        write_series_csv(self.rows, path, value=self.value)

        QMessageBox.information(self, "Series results", f"Saved to:\n{path}")

        return path
