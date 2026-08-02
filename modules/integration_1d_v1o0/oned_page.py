# ABOUTME: 1D integration page - load a series, pick peaks interactively, integrate all points.
# ABOUTME: Wires the plotter's selections to the series extraction and writes the results CSV.

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QKeySequence, QShortcut
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QComboBox, QFileDialog,
                               QHBoxLayout, QHeaderView, QLabel, QListWidget,
                               QListWidgetItem, QMessageBox, QPushButton, QSplitter, QTableWidget, QTableWidgetItem, QVBoxLayout,
                               QWidget)

from oned_loader import load_spectrum
from oned_plotter import OneDPlotter
from oned_series import (integrate_series, load_series, locate_peaks,
                         peak_at_position)
from oned_series_table import SeriesTableDialog

SPECTRUM_PATTERNS = ('*.ft1', '*.ft', '*.ft2', '*.ft3', '*.dat')

# NMRPipe writes 1D spectra as .ft1, so it leads the filter. .ft2/.ft3 are offered
# because they are picked by habit; a 2D file is rejected on load with a message
# rather than being hidden from the dialog.
SPECTRUM_FILTER = ("1D spectra (*.ft1 *.ft *.ft2 *.ft3 *.dat);;"
                   "NMRPipe 1D (*.ft1);;All files (*)")


class OneDIntegrationPage(QWidget):
    """Pick peaks on one spectrum, integrate them across the whole series."""

    # Appended to a spectrum's list entry as the run reaches it, so the propagation
    # is visible while it happens rather than only in the finished CSV.
    MEASURED_MARK = '✓'
    MISSING_MARK = '✗'

    def __init__(self, parent=None):
        super().__init__(parent)

        self.paths = []
        self.spectrum = None
        self.peaks = []
        self.table_rows = []

        # Peak picking is the undoable state. Each entry is a full snapshot of the
        # peak list: it is a handful of small dicts, so storing the whole thing is
        # simpler and less error-prone than reversing individual edits, and a box
        # drag that adds several peaks undoes as the single action it looked like.
        self._undo_stack = []
        self._redo_stack = []

        self._build_ui()

    # ---------------------------------------------------------------- ui

    def _build_ui(self):
        self.plotter = OneDPlotter()
        self.plotter.peaks_selected.connect(self.add_peaks)
        self.plotter.position_clicked.connect(self.add_peak_at)

        self.spectrum_list = QListWidget()
        self.spectrum_list.currentRowChanged.connect(self.show_spectrum)
        self.spectrum_list.itemChanged.connect(lambda _item: self._update_run_state())

        check_all_button = QPushButton("All")
        check_all_button.clicked.connect(self.check_all)
        check_none_button = QPushButton("None")
        check_none_button.clicked.connect(self.check_none)

        self.peak_table = QTableWidget(0, 3)
        self.peak_table.setHorizontalHeaderLabels(['Peak', 'ppm', 'Height'])
        self.peak_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.peak_table.setSelectionBehavior(QAbstractItemView.SelectRows)

        load_button = QPushButton("Load spectra…")
        load_button.clicked.connect(self.load_spectra)

        load_folder_button = QPushButton("Load folder…")
        load_folder_button.clicked.connect(self.load_folder)

        remove_button = QPushButton("Remove peak")
        remove_button.clicked.connect(self.remove_selected_peak)

        clear_button = QPushButton("Clear peaks")
        clear_button.clicked.connect(self.clear_peaks)

        self.observable = QComboBox()
        self.observable.addItems(['Intensity (height)', 'Area (region sum)'])

        self.integrate_button = QPushButton("Integrate series")
        self.integrate_button.clicked.connect(self.integrate)
        self.integrate_button.setEnabled(False)

        self.status = QLabel("Load a folder of 1D spectra to begin.")
        self.status.setWordWrap(True)

        # Every standard binding Qt lists for undo/redo is bound, not just the first:
        # the platform resolves Cmd+Z and Cmd+Shift+Z on macOS against Ctrl+Z and
        # Ctrl+Y elsewhere, and which one comes first depends on the running platform
        # theme. WindowShortcut so they fire wherever focus sits in the module.
        for standard, handler in ((QKeySequence.Undo, self.undo),
                                  (QKeySequence.Redo, self.redo)):
            for sequence in QKeySequence.keyBindings(standard):
                shortcut = QShortcut(sequence, self)
                shortcut.setContext(Qt.WindowShortcut)
                shortcut.activated.connect(handler)

        load_buttons = QHBoxLayout()
        load_buttons.addWidget(load_button)
        load_buttons.addWidget(load_folder_button)

        left = QVBoxLayout()
        left.addLayout(load_buttons)
        spectra_header = QHBoxLayout()
        spectra_header.addWidget(QLabel("Spectra (checked are integrated)"))
        spectra_header.addStretch()
        spectra_header.addWidget(check_all_button)
        spectra_header.addWidget(check_none_button)

        left.addLayout(spectra_header)
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

    def load_spectra(self, _=None, paths=None):
        """Load an explicit list of spectrum files, chosen in the file dialog."""
        if paths is None:
            paths, _filter = QFileDialog.getOpenFileNames(
                self, "Select 1D spectra", '', SPECTRUM_FILTER)
        if not paths:
            return

        self._set_paths([Path(p) for p in paths])

    def load_folder(self, _=None, folder=None):
        """Load every spectrum in a folder, in name order."""
        folder = folder or QFileDialog.getExistingDirectory(self, "Folder of 1D spectra")
        if not folder:
            return

        paths = []
        for pattern in SPECTRUM_PATTERNS:
            paths.extend(sorted(Path(folder).glob(pattern)))

        if not paths:
            self.status.setText(f"No 1D spectra found in {folder}")
            return

        self._set_paths(sorted(paths))

    def _set_paths(self, paths):
        """Keep the files that load as 1D, reporting the ones that do not.

        A folder often holds 2D spectra beside the 1D ones, and picking one should
        say so rather than aborting the whole selection.
        """
        readable, rejected = [], []
        for path in paths:
            try:
                load_spectrum(path)
                readable.append(path)
            except Exception as exc:
                rejected.append(f"{Path(path).name}: {exc}")

        self.paths = readable
        self.spectrum_list.blockSignals(True)
        self.spectrum_list.clear()
        for path in readable:
            item = QListWidgetItem(path.name)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Checked)
            self.spectrum_list.addItem(item)
        self.spectrum_list.blockSignals(False)

        if not readable:
            self._refresh_peaks()
            self.status.setText("Could not read any 1D spectra from the selection. "
                                + (rejected[0] if rejected else ''))
            return

        self.spectrum_list.setCurrentRow(0)

        # A newly opened set gets the full sweep. Browsing within a set keeps the view,
        # so a peak can be followed at zoom, but inheriting the previous set's zoom (or
        # matplotlib's empty-axes default on the very first load) shows the wrong region.
        self.plotter.reset_view()

        message = (f"Loaded {len(readable)} spectra. Drag a box over a region to "
                   "detect peaks, or switch to click mode to place one.")
        if rejected:
            message += f"  Skipped {len(rejected)} unreadable file(s)."
        self.status.setText(message)

    def show_spectrum(self, row):
        if not (0 <= row < len(self.paths)):
            return
        self.spectrum = load_spectrum(self.paths[row])
        self.plotter.set_spectrum(self.spectrum, keep_view=True)
        self._refresh_peaks()

    # ------------------------------------------------------------- undo/redo

    def _record(self):
        """Snapshot the peak list before changing it."""
        self._undo_stack.append([dict(peak) for peak in self.peaks])
        self._redo_stack.clear()

    def undo(self):
        if not self._undo_stack:
            return
        self._redo_stack.append([dict(peak) for peak in self.peaks])
        self.peaks = self._undo_stack.pop()
        self._refresh_peaks()

    def redo(self):
        if not self._redo_stack:
            return
        self._undo_stack.append([dict(peak) for peak in self.peaks])
        self.peaks = self._redo_stack.pop()
        self._refresh_peaks()

    # ---------------------------------------------------------- peak picking

    def add_peaks(self, found):
        self._record()
        for peak in found:
            self._append_peak(peak['ppm'], peak.get('height'))
        self._refresh_peaks()

    def add_peak_at(self, ppm):
        if self.spectrum is None:
            return
        self._record()
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
            self._record()
            self.peaks.pop(row)
            self._refresh_peaks()

    def clear_peaks(self):
        if not self.peaks:
            return
        self._record()
        self.peaks = []
        self._refresh_peaks()

    # ------------------------------------------------------- spectrum selection

    def checked_paths(self):
        """Spectra the user has left ticked, in list order."""
        return [self.paths[i] for i in range(self.spectrum_list.count())
                if self.spectrum_list.item(i).checkState() == Qt.Checked]

    def _set_all_checked(self, checked):
        self.spectrum_list.blockSignals(True)
        for i in range(self.spectrum_list.count()):
            self.spectrum_list.item(i).setCheckState(Qt.Checked if checked else Qt.Unchecked)
        self.spectrum_list.blockSignals(False)
        self._update_run_state()

    def check_all(self):
        self._set_all_checked(True)

    def check_none(self):
        self._set_all_checked(False)

    def _base_name(self, text):
        """The spectrum name without whichever status mark was appended to it."""
        return text.rstrip().rstrip(self.MEASURED_MARK + self.MISSING_MARK).rstrip()

    def _clear_marks(self):
        """Drop the status marks so a new run does not show the previous one's."""
        self.spectrum_list.blockSignals(True)
        for i in range(self.spectrum_list.count()):
            item = self.spectrum_list.item(i)
            item.setText(self._base_name(item.text()))
            item.setToolTip('')
        self.spectrum_list.blockSignals(False)

    def _mark_measured(self, row, rows):
        """Annotate a spectrum in the list once it has been measured.

        This is the propagation readout: the list fills in top to bottom as the run
        proceeds, and a peak that could not be matched is marked so a gap in the
        series is visible rather than buried in the CSV.
        """
        item = self.spectrum_list.item(row)
        if item is None:
            return

        unmatched = [r for r in rows if not r['matched']]
        mark = self.MEASURED_MARK if not unmatched else self.MISSING_MARK

        self.spectrum_list.blockSignals(True)
        item.setText(f"{self._base_name(item.text())}  {mark}")
        item.setToolTip('  '.join(
            f"{r['assignment']}: {'—' if r['height'] is None else format(r['height'], '.3e')}"
            for r in rows))
        self.spectrum_list.blockSignals(False)

    def _update_run_state(self):
        self.integrate_button.setEnabled(bool(self.peaks) and bool(self.checked_paths()))

    def _located_peaks(self):
        """Each picked peak as it appears in the spectrum currently on screen.

        The reference position is where the user picked; a resonance drifts through a
        series, so a marker frozen at the reference slides off the peak as you browse.
        Located here with the same matching the series run uses, so the display agrees
        with what will be measured. Falls back to the reference where nothing matched.
        """
        if self.spectrum is None or not self.peaks:
            return [dict(peak) for peak in self.peaks]

        located = []
        for peak, match in zip(self.peaks, locate_peaks(self.spectrum, self.peaks)):
            located.append({
                'assignment': peak['assignment'],
                'position': peak['position'],
                'ppm': match['ppm'] if match['matched'] else peak['position'],
                'height': match['height'],
                'matched': match['matched'],
            })
        return located

    def _plotter_peaks(self):
        """The page keys a peak by its reference `position`; the plotter draws at `ppm`.
        Every hand-off goes through here so the two shapes cannot drift apart."""
        return [{'ppm': peak.get('ppm', peak['position']), 'height': peak.get('height'),
                 'assignment': peak['assignment']} for peak in self._located_peaks()]

    def _refresh_peaks(self):
        located = self._located_peaks()

        self.peak_table.setRowCount(len(located))
        for row, peak in enumerate(located):
            height = peak.get('height')
            ppm = peak.get('ppm', peak['position'])
            marker = '' if peak.get('matched', True) else ' ?'
            for col, text in enumerate((peak['assignment'],
                                        f"{ppm:.4f}{marker}",
                                        '' if height is None else f"{height:.3e}")):
                self.peak_table.setItem(row, col, QTableWidgetItem(text))

        self.plotter.set_peaks([{'ppm': p.get('ppm', p['position']),
                                 'height': p.get('height'),
                                 'assignment': p['assignment']} for p in located])
        self._update_run_state()
        self.status.setText(f"{len(self.peaks)} peak(s) selected across "
                            f"{len(self.checked_paths())} checked spectra.")

    # ----------------------------------------------------------- integration

    def integrate(self, _=None, show_table=True):
        """Measure every checked spectrum, then show the results for review and export."""
        selected = self.checked_paths()
        if not (self.peaks and selected):
            return None

        self._clear_marks()
        rows_by_path = [self.paths.index(p) for p in selected]

        def on_point(point, rows):
            self._mark_measured(rows_by_path[point], rows)
            self.status.setText(f"Integrating… {point + 1}/{len(selected)} "
                                f"({selected[point].name})")
            QApplication.processEvents()      # let the list fill in as the run proceeds

        spectra = load_series(selected)
        self.table_rows = integrate_series(spectra, [dict(p) for p in self.peaks],
                                           names=[p.stem for p in selected],
                                           progress=on_point)

        missing = sum(1 for r in self.table_rows if not r['matched'])
        self.status.setText(f"Integrated {len(self.peaks)} peak(s) over "
                            f"{len(selected)} spectra."
                            + (f"  {missing} unmatched point(s)." if missing else ''))

        if show_table:
            self.show_results_table()

        return self.table_rows

    def show_results_table(self):
        """Open the results grid: spectra down, peaks across, with export and copy."""
        if not self.table_rows:
            return None

        default = 'height' if self.observable.currentIndex() == 0 else 'area'
        self.results_dialog = SeriesTableDialog(self.table_rows, parent=self,
                                                default_value=default)
        self.results_dialog.show()

        return self.results_dialog
