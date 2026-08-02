# ABOUTME: 1D integration page - load a series, pick peaks interactively, integrate all points.
# ABOUTME: Wires the plotter's selections to the series extraction and writes the results CSV.

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QAbstractItemView, QComboBox, QFileDialog, QHBoxLayout,
                               QHeaderView, QLabel, QListWidget, QMessageBox,
                               QPushButton, QSplitter, QTableWidget, QTableWidgetItem,
                               QVBoxLayout, QWidget)

from oned_loader import load_spectrum
from oned_plotter import OneDPlotter
from oned_series import (integrate_series, load_series, peak_at_position,
                         write_series_csv)

SPECTRUM_PATTERNS = ('*.ft1', '*.ft', '*.dat')


class OneDIntegrationPage(QWidget):
    """Pick peaks on one spectrum, integrate them across the whole series."""

    def __init__(self, parent=None):
        super().__init__(parent)

        self.paths = []
        self.spectrum = None
        self.peaks = []
        self.table_rows = []

        self._build_ui()

    # ---------------------------------------------------------------- ui

    def _build_ui(self):
        self.plotter = OneDPlotter()
        self.plotter.peaks_selected.connect(self.add_peaks)
        self.plotter.position_clicked.connect(self.add_peak_at)

        self.spectrum_list = QListWidget()
        self.spectrum_list.currentRowChanged.connect(self.show_spectrum)

        self.peak_table = QTableWidget(0, 3)
        self.peak_table.setHorizontalHeaderLabels(['Peak', 'ppm', 'Height'])
        self.peak_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.peak_table.setSelectionBehavior(QAbstractItemView.SelectRows)

        load_button = QPushButton("Load spectra…")
        load_button.clicked.connect(self.load_folder)

        remove_button = QPushButton("Remove peak")
        remove_button.clicked.connect(self.remove_selected_peak)

        clear_button = QPushButton("Clear peaks")
        clear_button.clicked.connect(self.clear_peaks)

        self.observable = QComboBox()
        self.observable.addItems(['Intensity (height)', 'Area (region sum)'])

        self.integrate_button = QPushButton("Integrate series → CSV")
        self.integrate_button.clicked.connect(self.integrate)
        self.integrate_button.setEnabled(False)

        self.status = QLabel("Load a folder of 1D spectra to begin.")
        self.status.setWordWrap(True)

        left = QVBoxLayout()
        left.addWidget(load_button)
        left.addWidget(QLabel("Spectra"))
        left.addWidget(self.spectrum_list)
        left.addWidget(QLabel("Selected peaks"))
        left.addWidget(self.peak_table)
        peak_buttons = QHBoxLayout()
        peak_buttons.addWidget(remove_button)
        peak_buttons.addWidget(clear_button)
        left.addLayout(peak_buttons)
        left.addWidget(QLabel("Observable"))
        left.addWidget(self.observable)
        left.addWidget(self.integrate_button)

        left_panel = QWidget()
        left_panel.setLayout(left)
        left_panel.setMaximumWidth(340)

        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(left_panel)
        splitter.addWidget(self.plotter)
        splitter.setStretchFactor(1, 1)

        layout = QVBoxLayout(self)
        layout.addWidget(splitter)
        layout.addWidget(self.status)

    # ------------------------------------------------------------- loading

    def load_folder(self, _=None, folder=None):
        folder = folder or QFileDialog.getExistingDirectory(self, "Folder of 1D spectra")
        if not folder:
            return

        paths = []
        for pattern in SPECTRUM_PATTERNS:
            paths.extend(sorted(Path(folder).glob(pattern)))

        if not paths:
            self.status.setText(f"No 1D spectra found in {folder}")
            return

        self.paths = paths
        self.spectrum_list.clear()
        self.spectrum_list.addItems([p.name for p in paths])
        self.spectrum_list.setCurrentRow(0)
        self.status.setText(f"Loaded {len(paths)} spectra. "
                            "Drag a box over a region to detect peaks, or switch to "
                            "click mode to place one.")

    def show_spectrum(self, row):
        if not (0 <= row < len(self.paths)):
            return
        self.spectrum = load_spectrum(self.paths[row])
        self.plotter.set_spectrum(self.spectrum, keep_view=True)
        self.plotter.set_peaks(self.peaks)

    # ---------------------------------------------------------- peak picking

    def add_peaks(self, found):
        for peak in found:
            self._append_peak(peak['ppm'], peak.get('height'))
        self._refresh_peaks()

    def add_peak_at(self, ppm):
        if self.spectrum is None:
            return
        peak = peak_at_position(self.spectrum, ppm)
        self._append_peak(peak['ppm'], peak.get('height'))
        self._refresh_peaks()

    def _append_peak(self, ppm, height):
        # a second pick on the same resonance is a mis-click, not a new peak
        if any(abs(p['position'] - ppm) < 1e-6 for p in self.peaks):
            return
        self.peaks.append({'assignment': f"peak_{len(self.peaks) + 1}",
                           'position': float(ppm), 'height': height})

    def remove_selected_peak(self):
        row = self.peak_table.currentRow()
        if 0 <= row < len(self.peaks):
            self.peaks.pop(row)
            self._refresh_peaks()

    def clear_peaks(self):
        self.peaks = []
        self._refresh_peaks()

    def _refresh_peaks(self):
        self.peak_table.setRowCount(len(self.peaks))
        for row, peak in enumerate(self.peaks):
            height = peak.get('height')
            for col, text in enumerate((peak['assignment'],
                                        f"{peak['position']:.4f}",
                                        '' if height is None else f"{height:.3e}")):
                self.peak_table.setItem(row, col, QTableWidgetItem(text))

        self.plotter.set_peaks([{'ppm': p['position'], 'height': p.get('height'),
                                 'assignment': p['assignment']} for p in self.peaks])
        self.integrate_button.setEnabled(bool(self.peaks) and bool(self.paths))
        self.status.setText(f"{len(self.peaks)} peak(s) selected across {len(self.paths)} spectra.")

    # ----------------------------------------------------------- integration

    def integrate(self, _=None, out_path=None):
        if not (self.peaks and self.paths):
            return None

        value = 'height' if self.observable.currentIndex() == 0 else 'area'

        out_path = out_path or QFileDialog.getSaveFileName(
            self, "Save series table", str(Path(self.paths[0]).parent / "series_intensity.csv"),
            "CSV (*.csv)")[0]
        if not out_path:
            return None

        self.status.setText(f"Integrating {len(self.peaks)} peak(s) over {len(self.paths)} spectra…")

        spectra = load_series(self.paths)
        self.table_rows = integrate_series(spectra, [dict(p) for p in self.peaks],
                                           names=[p.stem for p in self.paths])
        write_series_csv(self.table_rows, out_path, value=value)

        self.status.setText(f"Wrote {out_path}")
        QMessageBox.information(self, "Series integration",
                                f"Integrated {len(self.peaks)} peak(s) over "
                                f"{len(self.paths)} spectra.\n\nSaved to:\n{out_path}")
        return out_path
