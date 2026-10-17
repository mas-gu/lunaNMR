# ABOUTME: dynamiXs page for the heteronuclear NOE: I_sat / I_unsat per residue, table and plot.
# ABOUTME: Input is a hetNOE series' intensity matrix (both planes as columns) or two intensity files.

import os

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFrame, QHBoxLayout, QLabel, QFileDialog, QDoubleSpinBox, QComboBox, QLineEdit,
    QPlainTextEdit, QStackedWidget, QTableWidget, QTableWidgetItem, QHeaderView,
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from constants import (
    PRIMARY_TEXT, SECONDARY_TEXT, FRAME_BG_COLOR,
    SPACING_XS, SPACING_SM, FONT_SIZE_SMALL,
)
from gui_components import (
    create_primary_button, create_secondary_button, create_label,
    create_header_label, create_v_layout, create_h_layout, get_font, show_error,
    open_directory_dialog,
)
from dynamiXs_gui import BasePage, DropTargetLabel
from dynamiXs_integrated.data_converters import (
    calculate_hetnoe_from_intensities, draw_hetnoe_vs_residue, in_sequence_order,
    intensities_from_matrix, parse_intensity_csv, plot_hetnoe_vs_residue, read_hetnoe_csv,
    write_hetnoe_csv,
)

MATRIX_MODE, FILES_MODE = 0, 1
MATRIX_ID_COLUMNS = ('Peak_Number', 'Assignment', 'Reference_X', 'Reference_Y')


class HetNOEPage(BasePage):
    """Heteronuclear NOE per residue, computed as `dynamixs hetnoe` computes it.

    One field at a time. From a series matrix the saturated and unsaturated planes are
    two of its point columns, picked by name; from two files, each is the CLI's
    `residue,intensity[,error]`. A matrix carries no intensity errors, so its hetNOE
    error is the flat 2% estimate - the files are the route for measured ones.
    """

    def __init__(self, main_window):
        super().__init__(main_window, "Heteronuclear NOE Analysis")

        # State
        self.matrix_file = None
        self.sat_file = None
        self.unsat_file = None
        self.output_dir = None
        self.source_series_map = {}  # for project save/restore parity
        self.last_results_file = None

        self._setup_content()

    # ---------- UI construction ----------

    def _setup_content(self):
        layout = self.content_layout

        desc = create_label(
            "hetNOE = I_sat / I_unsat per residue. A residue with no positive unsaturated "
            "intensity has no hetNOE. A flexible residue's is negative - real data, kept."
        )
        desc.setStyleSheet(f"color: {SECONDARY_TEXT};")
        desc.setWordWrap(True)
        layout.addWidget(desc)

        # ---- Input ----
        input_section = self._make_section("hetNOE input data")
        body = input_section.layout()

        self.mode_combo = QComboBox()
        self.mode_combo.addItems(["Series intensity matrix", "Two intensity files"])
        self.mode_combo.currentIndexChanged.connect(self._on_mode_changed)
        body.addWidget(self._make_field_row("Input", self.mode_combo))

        self.input_stack = QStackedWidget()
        self.input_stack.addWidget(self._matrix_inputs())
        self.input_stack.addWidget(self._file_inputs())
        body.addWidget(self.input_stack)
        layout.addWidget(input_section)

        # ---- Output ----
        out_section = self._make_section("Output")
        outb = out_section.layout()
        out_row = QHBoxLayout()
        out_row.setSpacing(SPACING_SM)
        self.outdir_label = QLabel("No output directory selected")
        self.outdir_label.setStyleSheet(
            f"color: {SECONDARY_TEXT}; background-color: {FRAME_BG_COLOR}; "
            f"padding: 4px 8px; border-radius: 4px;"
        )
        out_row.addWidget(self.outdir_label, stretch=1)
        out_row.addWidget(create_secondary_button("Choose…", clicked=self._choose_output_dir,
                                                  width=110))
        outb.addLayout(out_row)
        self.prefix_edit = QLineEdit("field1")
        outb.addWidget(self._make_field_row("Results prefix", self.prefix_edit))
        self.field_freq_spin = QDoubleSpinBox()
        self.field_freq_spin.setRange(100.0, 1500.0)
        self.field_freq_spin.setDecimals(2)
        self.field_freq_spin.setValue(600.0)
        outb.addWidget(self._make_field_row("Field frequency (MHz)", self.field_freq_spin))
        layout.addWidget(out_section)

        run_row = QHBoxLayout()
        run_row.setSpacing(SPACING_SM)
        self.run_btn = create_primary_button("Compute hetNOE", clicked=self._run, width=200)
        run_row.addWidget(self.run_btn)
        run_row.addStretch()
        layout.addLayout(run_row)

        # ---- Results ----
        results_section = self._make_section("Results")
        rb = results_section.layout()
        results_row = QHBoxLayout()
        results_row.setSpacing(SPACING_SM)
        self.results_table = QTableWidget(0, 3)
        self.results_table.setHorizontalHeaderLabels(["Residue", "hetNOE", "± error"])
        self.results_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.results_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.results_table.setMinimumWidth(260)
        self.results_table.setMinimumHeight(340)
        results_row.addWidget(self.results_table)
        self.figure = Figure(figsize=(6, 3))
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setMinimumHeight(340)
        results_row.addWidget(self.canvas, stretch=1)
        rb.addLayout(results_row)
        layout.addWidget(results_section, stretch=1)

        self.log_text = QPlainTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setMaximumBlockCount(500)
        self.log_text.setMaximumHeight(110)
        self.log_text.setStyleSheet(
            f"QPlainTextEdit {{ color: {PRIMARY_TEXT}; background-color: {FRAME_BG_COLOR}; "
            f"border-radius: 4px; padding: 4px; }}"
        )
        layout.addWidget(self.log_text)

    def _matrix_inputs(self) -> QFrame:
        frame = QFrame()
        box = create_v_layout(SPACING_XS)
        frame.setLayout(box)
        row = QHBoxLayout()
        row.setSpacing(SPACING_SM)
        self.matrix_drop = DropTargetLabel("hetnoe", "Drop a hetNOE series here, or browse")
        self.matrix_drop.setMinimumHeight(36)
        self.matrix_drop.series_dropped.connect(self._on_series_dropped)
        row.addWidget(self.matrix_drop, stretch=1)
        row.addWidget(create_secondary_button("Browse CSV…", clicked=self._browse_matrix,
                                              width=130))
        box.addLayout(row)
        self.sat_combo = QComboBox()
        self.unsat_combo = QComboBox()
        box.addWidget(self._make_field_row("Saturated plane", self.sat_combo))
        box.addWidget(self._make_field_row("Unsaturated plane", self.unsat_combo))
        hint = create_label("A matrix carries no intensity errors: each hetNOE error is a flat "
                            "2% estimate. For measured errors, use two intensity files.")
        hint.setStyleSheet(f"color: {SECONDARY_TEXT};")
        hint.setFont(get_font(FONT_SIZE_SMALL))
        hint.setWordWrap(True)
        box.addWidget(hint)
        return frame

    def _file_inputs(self) -> QFrame:
        frame = QFrame()
        box = create_v_layout(SPACING_XS)
        frame.setLayout(box)
        self.sat_label = create_label("No saturated file", font_size=FONT_SIZE_SMALL)
        self.unsat_label = create_label("No unsaturated file", font_size=FONT_SIZE_SMALL)
        for label, browse in ((self.sat_label, self._browse_sat),
                              (self.unsat_label, self._browse_unsat)):
            row = QHBoxLayout()
            row.setSpacing(SPACING_SM)
            label.setStyleSheet(f"color: {SECONDARY_TEXT};")
            row.addWidget(label, stretch=1)
            row.addWidget(create_secondary_button("Browse…", clicked=browse, width=110))
            box.addLayout(row)
        hint = create_label("residue,intensity[,error] - as `dynamixs hetnoe --sat/--unsat`.")
        hint.setStyleSheet(f"color: {SECONDARY_TEXT};")
        hint.setFont(get_font(FONT_SIZE_SMALL))
        box.addWidget(hint)
        return frame

    def _make_section(self, title: str) -> QFrame:
        frame = QFrame()
        frame.setProperty("class", "card")
        frame.setStyleSheet(
            f"QFrame[class='card'] {{ background-color: {FRAME_BG_COLOR}; "
            f"border-radius: 8px; }}"
        )
        body = create_v_layout(SPACING_XS, (SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM))
        frame.setLayout(body)
        body.addWidget(create_header_label(title, level=3))
        return frame

    def _make_field_row(self, label_text: str, widget):
        row = QFrame()
        row_layout = create_h_layout(SPACING_SM)
        row.setLayout(row_layout)
        row_layout.addWidget(create_label(label_text), stretch=1)
        widget.setMaximumWidth(260)
        row_layout.addWidget(widget)
        return row

    # ---------- input selection ----------

    def _on_mode_changed(self, index: int):
        self.input_stack.setCurrentIndex(index)

    def _on_series_dropped(self, field_name: str, series_name: str, csv_path: str):
        self.source_series_map[field_name] = series_name
        matrix = os.path.join(os.path.dirname(csv_path), "peak_intensity_matrix.csv") \
            if csv_path else None
        if matrix and os.path.exists(matrix):
            self._set_matrix(matrix)
            self.matrix_drop.setText(f"Series '{series_name}'  →  {os.path.basename(matrix)}")
        else:
            self.matrix_drop.setText(f"Series '{series_name}'  (no intensity matrix)")

    def _browse_matrix(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select a hetNOE series intensity matrix", self._start_dir(),
            "CSV files (*.csv);;All files (*)")
        if path:
            self._set_matrix(path)

    def _set_matrix(self, path, sat_column=None, unsat_column=None):
        """Load a matrix's point columns into both plane pickers; preselect by name."""
        import pandas as pd
        self.matrix_file = path
        self.matrix_drop.setText(os.path.basename(path))
        columns = [c for c in pd.read_csv(path, nrows=0).columns if c not in MATRIX_ID_COLUMNS]
        for combo in (self.sat_combo, self.unsat_combo):
            combo.clear()
            combo.addItems(columns)
        unsat = unsat_column or next((c for c in columns if 'unsat' in c.lower()), None)
        sat = sat_column or next((c for c in columns if 'sat' in c.lower()
                                  and 'unsat' not in c.lower()), None)
        if sat:
            self.sat_combo.setCurrentText(sat)
        if unsat:
            self.unsat_combo.setCurrentText(unsat)
        self._log(f"Matrix: {path}")

    def _browse_sat(self):
        path = self._pick_file("Select the saturated intensities")
        if path:
            self._set_sat_file(path)

    def _browse_unsat(self):
        path = self._pick_file("Select the unsaturated intensities")
        if path:
            self._set_unsat_file(path)

    def _set_sat_file(self, path):
        self.sat_file = path
        self.sat_label.setText(os.path.basename(path))

    def _set_unsat_file(self, path):
        self.unsat_file = path
        self.unsat_label.setText(os.path.basename(path))

    def _pick_file(self, title):
        path, _ = QFileDialog.getOpenFileName(self, title, self._start_dir(),
                                              "CSV files (*.csv);;All files (*)")
        return path

    def _choose_output_dir(self):
        folder = open_directory_dialog(self, "Select output directory", self._start_dir())
        if folder:
            self.output_dir = folder
            self.outdir_label.setText(folder)

    def _start_dir(self) -> str:
        for path in (self.output_dir, self.matrix_file, self.sat_file):
            if path and os.path.exists(path):
                return path if os.path.isdir(path) else os.path.dirname(path)
        return getattr(self.main_window, "current_dir", "") or ""

    # ---------- run ----------

    def _input_names(self) -> str:
        """The current inputs, by file name, for a message."""
        if self.mode_combo.currentIndex() == MATRIX_MODE:
            return os.path.basename(self.matrix_file or '')
        return " and ".join(os.path.basename(p or '') for p in (self.sat_file, self.unsat_file))

    def _inputs(self):
        """(sat intensities, unsat intensities, sat errors, unsat errors), or a reason
        the inputs are not ready."""
        if self.mode_combo.currentIndex() == MATRIX_MODE:
            if not self.matrix_file:
                return "Choose a hetNOE series intensity matrix."
            sat, unsat = self.sat_combo.currentText(), self.unsat_combo.currentText()
            if not sat or not unsat or sat == unsat:
                return "Pick two different planes: the saturated and the unsaturated."
            return (intensities_from_matrix(self.matrix_file, sat),
                    intensities_from_matrix(self.matrix_file, unsat), None, None)
        if not (self.sat_file and self.unsat_file):
            return "Choose both the saturated and the unsaturated intensity files."
        sat_i, sat_e = parse_intensity_csv(self.sat_file)
        unsat_i, unsat_e = parse_intensity_csv(self.unsat_file)
        return sat_i, unsat_i, sat_e, unsat_e

    def _run(self):
        if not self.output_dir:
            show_error(self, "hetNOE", "Choose an output directory.")
            return
        prefix = self.prefix_edit.text().strip() or "field1"
        try:
            inputs = self._inputs()
            if isinstance(inputs, str):          # inputs not chosen yet: say which
                show_error(self, "hetNOE", inputs)
                return
            sat_i, unsat_i, sat_e, unsat_e = inputs
            noe = calculate_hetnoe_from_intensities(
                saturated_data=sat_i, unsaturated_data=unsat_i,
                saturated_errors=sat_e, unsaturated_errors=unsat_e)
            os.makedirs(self.output_dir, exist_ok=True)
            results = os.path.join(self.output_dir, f"{prefix}_hetnoe.csv")
            write_hetnoe_csv(noe, results)
            plot_hetnoe_vs_residue(noe, os.path.join(self.output_dir, f"{prefix}_hetnoe.pdf"),
                                   title=prefix)
        except Exception as error:
            # A file that is not an intensity table, a plane missing from the matrix, an
            # unwritable folder: said, and the results shown stay the previous ones.
            message = f"Could not compute hetNOE from {self._input_names()}: {error}"
            self._log(message)
            show_error(self, "hetNOE", message)
            return
        self.last_results_file = results
        self._show(noe)
        self._log(f"hetNOE: {len(noe)} residues ({len(sat_i)} saturated, "
                  f"{len(unsat_i)} unsaturated) -> {results}")

    def _show(self, noe):
        """Fill the table and the plot with `noe` ({residue: {'value', 'error'}})."""
        self.results_table.setRowCount(len(noe))
        for row, (residue, rec) in enumerate(in_sequence_order(noe)):
            for column, text in enumerate((str(residue), f"{rec['value']:.3f}",
                                           f"{rec['error']:.3f}")):
                item = QTableWidgetItem(text)
                item.setTextAlignment(Qt.AlignCenter)
                self.results_table.setItem(row, column, item)
        self.figure.clear()
        if noe:
            draw_hetnoe_vs_residue(self.figure.add_subplot(), noe,
                                   self.prefix_edit.text().strip() or None)
            self.figure.tight_layout()
        self.canvas.draw_idle()

    # ---------- project state ----------

    def get_session_state(self) -> dict:
        """Serializable state for project save."""
        return {
            'mode': self.mode_combo.currentIndex(),
            'matrix_file': self.matrix_file,
            'sat_column': self.sat_combo.currentText() or None,
            'unsat_column': self.unsat_combo.currentText() or None,
            'sat_file': self.sat_file,
            'unsat_file': self.unsat_file,
            'output_dir': self.output_dir,
            'prefix': self.prefix_edit.text(),
            'field_freq': self.field_freq_spin.value(),
            'source_series_map': dict(self.source_series_map),
            'last_results_file': self.last_results_file,
        }

    def restore_session_state(self, state: dict):
        """Restore state produced by get_session_state(), and the results it points at."""
        if not state:
            return
        self.mode_combo.setCurrentIndex(state.get('mode', MATRIX_MODE))
        self.source_series_map = dict(state.get('source_series_map', {}))
        matrix = state.get('matrix_file')
        if matrix and os.path.exists(matrix):
            self._set_matrix(matrix, state.get('sat_column'), state.get('unsat_column'))
        else:
            self.matrix_file = matrix
        for path, setter in ((state.get('sat_file'), self._set_sat_file),
                             (state.get('unsat_file'), self._set_unsat_file)):
            if path:
                setter(path)
        self.output_dir = state.get('output_dir')
        if self.output_dir:
            self.outdir_label.setText(self.output_dir)
        self.prefix_edit.setText(state.get('prefix', 'field1'))
        self.field_freq_spin.setValue(state.get('field_freq', 600.0))
        self.last_results_file = state.get('last_results_file')
        if self.last_results_file and os.path.exists(self.last_results_file):
            self._show(read_hetnoe_csv(self.last_results_file))

    def set_file_refs(self, refs: dict):
        """Apply saved input file references (file_refs() of a saved analysis)."""
        if refs.get('matrix_file') and os.path.exists(refs['matrix_file']):
            self._set_matrix(refs['matrix_file'])
        if refs.get('sat_file'):
            self._set_sat_file(refs['sat_file'])
        if refs.get('unsat_file'):
            self._set_unsat_file(refs['unsat_file'])

    def file_refs(self) -> dict:
        """The input files this analysis was computed from."""
        return {key: getattr(self, key) for key in ('matrix_file', 'sat_file', 'unsat_file')
                if getattr(self, key)}

    # ---------- logging ----------

    def _log(self, msg: str):
        self.log_text.appendPlainText(msg)
