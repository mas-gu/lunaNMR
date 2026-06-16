# ABOUTME: dynamiXs page for methyl bi-exp T2 fitting analysis (I(t) = A_a*exp(-t/T2a) + A_b*exp(-t/T2b) + C).
# ABOUTME: Mirrors the T1/T2 page visual pattern but exposes only initial guesses relevant to bi-exp fitting.

import os
from pathlib import Path

import pandas as pd
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFrame, QHBoxLayout, QLabel, QFileDialog, QDoubleSpinBox, QSpinBox,
    QComboBox, QProgressBar, QPlainTextEdit,
)

from constants import (
    BG_COLOR, PANEL_BG_COLOR, FRAME_BG_COLOR,
    PRIMARY_TEXT, SECONDARY_TEXT, SECONDARY_BUTTON_BORDER,
    SPACING_XS, SPACING_SM, SPACING_MD,
    BUTTON_CORNER_RADIUS, BUTTON_HEIGHT_DIALOG,
    FONT_SIZE_BODY, FONT_SIZE_SMALL,
)
from gui_components import (
    create_primary_button, create_secondary_button, create_label,
    create_header_label, create_v_layout, create_h_layout,
    get_font, show_info, show_error, show_warning, open_directory_dialog,
)


# Imported from dynamiXs_gui to avoid duplication of widget classes.
from dynamiXs_gui import BasePage, DropTargetLabel
from workers import MethylT2FittingParams, MethylT2FittingWorker


class MethylT2FittingPage(BasePage):
    """Page for bi-exponential methyl T2 fitting.

    Streamlined version of the T1/T2 page: single experiment, single series,
    bi-exp-specific initial guesses (A_a fraction, T2_a, T2_b), and one Run button.
    """

    def __init__(self, main_window):
        super().__init__(main_window, "T2 Methyl Fitting Analysis (bi-exponential)")

        # State
        self.input_file = None
        self.output_dir = None
        self.source_series_map = {}  # for project save/restore parity
        self.last_json_folder = None
        self.last_results_file = None

        self._setup_content()

    # ---------- UI construction ----------

    def _setup_content(self):
        layout = self.content_layout

        # Description
        desc = create_label(
            "Methyl bi-exponential T2 fit (shared amplitude):\n"
            "    I(t) = ½ · A · [ exp(-t / T2a) + exp(-t / T2b) ]\n"
            "Post-fit ordering enforces T2a ≥ T2b; ratio η = T2a / T2b is reported."
        )
        desc.setStyleSheet(f"color: {SECONDARY_TEXT};")
        desc.setWordWrap(True)
        layout.addWidget(desc)

        # ---- Input data row ----
        input_section = self._make_section("Methyl T2 input data")
        section_body = input_section.layout()

        file_row = QHBoxLayout()
        file_row.setSpacing(SPACING_SM)

        self.file_drop = DropTargetLabel("methyl_t2", "Drop a series here, or browse for a CSV")
        self.file_drop.setMinimumHeight(36)
        self.file_drop.series_dropped.connect(self._on_series_dropped)
        file_row.addWidget(self.file_drop, stretch=1)

        browse_btn = create_secondary_button("Browse CSV…", clicked=self._browse_input_file, width=130)
        file_row.addWidget(browse_btn)

        section_body.addLayout(file_row)
        layout.addWidget(input_section)

        # ---- Initial guesses ----
        guess_section = self._make_section("Initial guesses (bi-exp)")
        gb = guess_section.layout()

        gb.addWidget(self._make_field_row("Initial T2a (slow, ms)",
                                          self._make_float_spin("t2a_spin", 0.1, 1e5, 100.0, decimals=2, step=10.0)))
        gb.addWidget(self._make_field_row("Initial T2b (fast, ms)",
                                          self._make_float_spin("t2b_spin", 0.1, 1e5, 20.0, decimals=2, step=5.0)))

        hint = create_label(
            "A defaults to span(y) = max(y) − min(y). "
            "The 0.5 weighting between the two components is fixed (paper convention)."
        )
        hint.setStyleSheet(f"color: {SECONDARY_TEXT};")
        hint.setFont(get_font(FONT_SIZE_SMALL))
        hint.setWordWrap(True)
        gb.addWidget(hint)

        layout.addWidget(guess_section)

        # ---- Fit options ----
        opt_section = self._make_section("Fit options")
        ob = opt_section.layout()

        ob.addWidget(self._make_field_row("Field frequency (MHz)",
                                          self._make_float_spin("field_freq_spin", 100.0, 1500.0, 600.0, decimals=2, step=10.0)))

        n_boot = QSpinBox()
        n_boot.setRange(50, 100000)
        n_boot.setValue(1000)
        n_boot.setSingleStep(100)
        self.boot_spin = n_boot
        ob.addWidget(self._make_field_row("Bootstrap iterations", n_boot))

        err = QComboBox()
        err.addItems(["Analytical (covariance)", "Bootstrap"])
        self.error_combo = err
        ob.addWidget(self._make_field_row("Error method", err))

        layout.addWidget(opt_section)

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
        out_btn = create_secondary_button("Choose…", clicked=self._choose_output_dir, width=110)
        out_row.addWidget(out_btn)
        outb.addLayout(out_row)

        prefix_row = QHBoxLayout()
        prefix_row.addWidget(QLabel("Results prefix:"))
        self.prefix_edit = QLabel("field1_methylT2")
        self.prefix_edit.setStyleSheet(
            f"color: {PRIMARY_TEXT}; background-color: {FRAME_BG_COLOR}; "
            f"padding: 4px 8px; border-radius: 4px;"
        )
        prefix_row.addWidget(self.prefix_edit, stretch=1)
        outb.addLayout(prefix_row)

        layout.addWidget(out_section)

        # ---- Run + viewer buttons ----
        run_row = QHBoxLayout()
        run_row.setSpacing(SPACING_SM)
        self.run_btn = create_primary_button("Run methyl T2 fit", clicked=self._start_analysis, width=200)
        run_row.addWidget(self.run_btn)
        self.viewer_btn = create_secondary_button("Open Methyl T2 Viewer",
                                                  clicked=self._open_viewer, width=200)
        self.viewer_btn.setEnabled(False)
        run_row.addWidget(self.viewer_btn)
        self.summary_btn = create_secondary_button("Plot T2 vs Residue",
                                                   clicked=self._open_viewer_with_summary, width=200)
        self.summary_btn.setEnabled(False)
        self.summary_btn.setToolTip(
            "Open the Methyl T2 viewer with the T2a, T2b, and η across-residues "
            "panels pre-enabled."
        )
        run_row.addWidget(self.summary_btn)
        self.open_prev_btn = create_secondary_button("Open previous dataset…",
                                                     clicked=self._open_previous_dataset,
                                                     width=200)
        self.open_prev_btn.setToolTip(
            "Load a folder of *_methylT2_fit_data.json from a prior run\n"
            "into the viewer without re-running the fit."
        )
        run_row.addWidget(self.open_prev_btn)
        run_row.addStretch()
        layout.addLayout(run_row)

        # ---- Progress + log ----
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.hide()
        layout.addWidget(self.progress_bar)

        self.log_text = QPlainTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setMaximumBlockCount(500)
        self.log_text.setMinimumHeight(140)
        self.log_text.setStyleSheet(
            f"QPlainTextEdit {{ color: {PRIMARY_TEXT}; background-color: {FRAME_BG_COLOR}; "
            f"border-radius: 4px; padding: 4px; }}"
        )
        layout.addWidget(self.log_text, stretch=1)

    # ---------- helpers for compact field rows ----------

    def _make_section(self, title: str) -> QFrame:
        frame = QFrame()
        frame.setProperty("class", "card")
        frame.setStyleSheet(
            f"QFrame[class='card'] {{ background-color: {FRAME_BG_COLOR}; "
            f"border-radius: 8px; }}"
        )
        body = create_v_layout(SPACING_XS, (SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM))
        frame.setLayout(body)

        title_label = create_header_label(title, level=3)
        body.addWidget(title_label)
        return frame

    def _make_field_row(self, label_text: str, widget):
        row = QFrame()
        row_layout = create_h_layout(SPACING_SM)
        row.setLayout(row_layout)
        row_layout.addWidget(create_label(label_text), stretch=1)
        widget.setMaximumWidth(180)
        row_layout.addWidget(widget)
        return row

    def _make_float_spin(self, attr_name, lo, hi, default, decimals=2, step=1.0):
        spin = QDoubleSpinBox()
        spin.setRange(lo, hi)
        spin.setDecimals(decimals)
        spin.setSingleStep(step)
        spin.setValue(default)
        setattr(self, attr_name, spin)
        return spin

    # ---------- input selection ----------

    def _on_series_dropped(self, field_name: str, series_name: str, csv_path: str):
        self.source_series_map[field_name] = series_name
        if csv_path and os.path.exists(csv_path):
            self._set_input_file(csv_path)
            self.file_drop.setText(f"Series '{series_name}'  →  {os.path.basename(csv_path)}")
        else:
            self.file_drop.setText(f"Series '{series_name}'  (no CSV path)")

    def _browse_input_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select methyl T2 CSV",
            self.main_window.current_dir if hasattr(self.main_window, "current_dir") else "",
            "CSV files (*.csv);;All files (*)"
        )
        if path:
            self._set_input_file(path)
            self.file_drop.setText(os.path.basename(path))

    def _set_input_file(self, path):
        self.input_file = path
        self._log(f"Input set: {path}")

    def _choose_output_dir(self):
        initial = self._suggested_output_start_dir()
        folder = open_directory_dialog(self, "Select output directory", initial)
        if folder:
            self.output_dir = folder
            self.outdir_label.setText(folder)

    def _suggested_output_start_dir(self) -> str:
        """Open the output picker next to where the user is already working:
        previously-chosen output dir, else the input CSV's folder, else the
        main window's current_dir.
        """
        if self.output_dir and os.path.isdir(self.output_dir):
            return self.output_dir
        if getattr(self, "input_file", None) and os.path.exists(self.input_file):
            return os.path.dirname(self.input_file)
        return self.main_window.current_dir if hasattr(self.main_window, "current_dir") else ""

    # ---------- run ----------

    def _start_analysis(self):
        if not self.input_file or not os.path.exists(self.input_file):
            show_error(self, "Input required", "Select a CSV file or drop a series first.")
            return
        if not self.output_dir:
            self._choose_output_dir()
            if not self.output_dir:
                return

        params = MethylT2FittingParams(
            input_file=self.input_file,
            output_dir=self.output_dir,
            results_prefix="field1_methylT2",
            json_folder=os.path.join(self.output_dir, "json"),
            field_name="field1",
            field_freq=float(self.field_freq_spin.value()),
            initial_A=None,  # data-driven: span(y)
            initial_t2_a=float(self.t2a_spin.value()),
            initial_t2_b=float(self.t2b_spin.value()),
            n_bootstrap=int(self.boot_spin.value()),
            error_method="bootstrap" if self.error_combo.currentIndex() == 1 else "analytical",
        )

        self.progress_bar.setValue(0)
        self.progress_bar.show()
        self.run_btn.setEnabled(False)
        self.viewer_btn.setEnabled(False)
        self._log("Starting methyl T2 fit...")

        self.worker = MethylT2FittingWorker(params)
        self.worker.progress.connect(lambda msg: self._log(msg))
        self.worker.progress_value.connect(self.progress_bar.setValue)
        self.worker.finished.connect(self._on_finished)
        self.worker.error.connect(self._on_error)
        self.worker.start()

    def _on_finished(self, result: dict):
        self.run_btn.setEnabled(True)
        self.progress_bar.setValue(100)
        n_fitted = result.get("n_fitted", 0)
        self._log(f"Fit complete: {n_fitted} residues, results -> {result.get('results_file')}")
        self.last_json_folder = result.get("json_folder") or os.path.dirname(result.get("json_file") or "")
        self.last_results_file = result.get("results_file")
        if n_fitted:
            self.viewer_btn.setEnabled(True)
            self.summary_btn.setEnabled(True)
            show_info(self, "Fit complete", f"Fitted {n_fitted} residues. Open the viewer to inspect.")
        else:
            show_warning(self, "No fits succeeded", "All residues failed to fit. Check the log.")

    def _on_error(self, error_msg: str):
        self.run_btn.setEnabled(True)
        self._log(f"ERROR: {error_msg}")
        show_error(self, "Fit failed", str(error_msg))

    # ---------- viewer ----------

    def _open_viewer(self):
        if not self.last_json_folder:
            show_warning(self, "No results", "Run a fit first.")
            return
        try:
            from visualization.methyl_t2_fit_viewer import open_methyl_t2_fit_viewer
            self._viewer = open_methyl_t2_fit_viewer(parent=self, json_folder=self.last_json_folder)
        except Exception as e:
            show_error(self, "Viewer error", f"Could not open viewer:\n{e}")

    def _open_viewer_with_summary(self):
        """Open the across-residues summary plot (no per-residue fit panel)."""
        if not self.last_json_folder:
            show_warning(self, "No results", "Run a fit first.")
            return
        try:
            from visualization.methyl_t2_fit_viewer import open_methyl_t2_summary_viewer
            self._viewer = open_methyl_t2_summary_viewer(
                parent=self, json_folder=self.last_json_folder
            )
        except Exception as e:
            show_error(self, "Viewer error", f"Could not open viewer:\n{e}")

    def _open_previous_dataset(self):
        """Pick a folder of *_methylT2_fit_data.json from a prior run and
        open the viewer pointed at it — no need to re-run the fit.
        """
        initial = self.last_json_folder
        if not initial:
            initial = self.output_dir or (
                self.main_window.current_dir
                if hasattr(self.main_window, "current_dir") else ""
            )
        folder = open_directory_dialog(
            self, "Select folder containing *_methylT2_fit_data.json", initial
        )
        if not folder:
            return
        json_matches = list(Path(folder).glob("*_methylT2_fit_data.json"))
        if not json_matches:
            show_warning(
                self, "No methyl T2 fits found",
                f"No '*_methylT2_fit_data.json' files in:\n{folder}"
            )
            return
        try:
            from visualization.methyl_t2_fit_viewer import open_methyl_t2_fit_viewer
            viewer = open_methyl_t2_fit_viewer(parent=self, json_folder=folder)
        except Exception as e:
            show_error(self, "Viewer error", f"Could not open viewer:\n{e}")
            return
        # Only mutate page state once the viewer is successfully open; otherwise
        # a failed open would leave last_json_folder pointing at a bad path and
        # the subsequent buttons enabled.
        self._viewer = viewer
        self.last_json_folder = folder
        self.viewer_btn.setEnabled(True)
        self.summary_btn.setEnabled(True)
        self._log(f"Opened previous dataset: {folder} "
                  f"({len(json_matches)} field(s) found)")

    # ---------- session state (project save/load) ----------

    def get_session_state(self) -> dict:
        """Serializable state for project save."""
        return {
            'input_file': self.input_file,
            'output_dir': self.output_dir,
            'source_series_map': dict(self.source_series_map),
            'last_json_folder': self.last_json_folder,
            'last_results_file': self.last_results_file,
            'initial_t2_a': self.t2a_spin.value(),
            'initial_t2_b': self.t2b_spin.value(),
            'field_freq': self.field_freq_spin.value(),
            'bootstrap_iter': self.boot_spin.value(),
            'error_method': self.error_combo.currentIndex(),
        }

    def restore_session_state(self, state: dict):
        """Restore state produced by get_session_state()."""
        if not state:
            return
        self.input_file = state.get('input_file')
        self.output_dir = state.get('output_dir')
        self.source_series_map = dict(state.get('source_series_map', {}))
        self.last_json_folder = state.get('last_json_folder')
        self.last_results_file = state.get('last_results_file')

        self.t2a_spin.setValue(state.get('initial_t2_a', 100.0))
        self.t2b_spin.setValue(state.get('initial_t2_b', 20.0))
        self.field_freq_spin.setValue(state.get('field_freq', 600.0))
        self.boot_spin.setValue(state.get('bootstrap_iter', 1000))
        self.error_combo.setCurrentIndex(state.get('error_method', 0))

        if self.input_file:
            self.file_drop.setText(os.path.basename(self.input_file))
        if self.output_dir:
            self.outdir_label.setText(self.output_dir)
        if self.last_json_folder:
            self.viewer_btn.setEnabled(True)
            self.summary_btn.setEnabled(True)

    # ---------- logging ----------

    def _log(self, msg: str):
        self.log_text.appendPlainText(msg)
