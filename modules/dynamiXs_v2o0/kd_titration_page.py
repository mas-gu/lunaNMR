# ABOUTME: dynamiXs page for Kd / titration binding fits (CSP and intensity ratio).
# ABOUTME: 1:1 quadratic isotherm, per-residue + global Kd; mirrors the methyl-T2 page.

import os
import sys
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFrame, QHBoxLayout, QVBoxLayout, QLabel, QFileDialog, QDoubleSpinBox, QSpinBox,
    QComboBox, QFormLayout, QProgressBar, QPlainTextEdit, QScrollArea, QWidget,
    QListWidget, QListWidgetItem, QMessageBox,
)

from constants import (
    FRAME_BG_COLOR, PRIMARY_TEXT, SECONDARY_TEXT,
    SPACING_XS, SPACING_SM, FONT_SIZE_SMALL,
)
from gui_components import (
    create_primary_button, create_secondary_button, create_label,
    create_header_label, create_v_layout, create_h_layout,
    get_font, show_info, show_error, show_warning, open_directory_dialog,
)

from dynamiXs_gui import BasePage, DropTargetLabel, DraggableSeriesList
from workers import KdTitrationFittingParams, KdTitrationFittingWorker

# Make the Kd analysis package importable (bare-name imports inside it).
_KD_DIR = os.path.join(os.path.dirname(__file__), "dynamiXs_Kd")
if _KD_DIR not in sys.path:
    sys.path.insert(0, _KD_DIR)


class KdTitrationPage(BasePage):
    """Page for 1:1 Kd titration fitting from a LunaNMR series result."""

    def __init__(self, main_window):
        super().__init__(main_window, "Kd / Titration Analysis (CSP + intensity)")
        self.input_file = None
        self.output_dir = None
        self.detected_points = []
        self.last_json_folder = None
        self.last_json_file = None
        self.worker = None
        # Series identity — every series writes the SAME generic filename
        # (series_analysis_tidy.csv), so the saved-fit name must come from the series
        # name (drag) or the series folder, not the file stem.
        self.series_name = None            # set when a series is dropped
        self._current_analysis_name = None  # set when a saved fit is opened (keeps its name)
        self._setup_content()

    # ---------- UI ----------

    def _setup_content(self):
        # Wrap everything in a scroll area so the page (7+ concentration rows + all
        # parameters + log) is reachable without squishing on small windows.
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)
        inner = QWidget()
        layout = QVBoxLayout(inner)
        layout.setSpacing(SPACING_SM)
        scroll.setWidget(inner)
        self.content_layout.addWidget(scroll)

        # Input
        sec = self._make_section("Titration input data (series_analysis_tidy.csv or comprehensive_peak_tracking.csv)")
        body = sec.layout()

        # Titration-mode series runs from this project — drag one onto the box below.
        series_hdr = create_label("Titration-mode series runs — drag one onto the box below:")
        series_hdr.setStyleSheet(f"color: {SECONDARY_TEXT};")
        series_hdr.setFont(get_font(FONT_SIZE_SMALL))
        body.addWidget(series_hdr)
        self.series_list_widget = DraggableSeriesList()
        self.series_list_widget.setMaximumHeight(110)
        body.addWidget(self.series_list_widget)
        self.no_series_label = create_label(
            "No titration-mode series in this project — run a series in Titration mode, "
            "or use Browse CSV…")
        self.no_series_label.setStyleSheet(f"color: {SECONDARY_TEXT};")
        self.no_series_label.setFont(get_font(FONT_SIZE_SMALL))
        self.no_series_label.setWordWrap(True)
        body.addWidget(self.no_series_label)

        row = QHBoxLayout()
        row.setSpacing(SPACING_SM)
        self.file_drop = DropTargetLabel("kd", "Drop a series here, or browse for a CSV")
        self.file_drop.setMinimumHeight(36)
        self.file_drop.series_dropped.connect(self._on_series_dropped)
        row.addWidget(self.file_drop, stretch=1)
        row.addWidget(create_secondary_button("Browse CSV…", clicked=self._browse_input_file, width=130))
        body.addLayout(row)
        self.points_label = create_label("No file loaded.")
        self.points_label.setStyleSheet(f"color: {SECONDARY_TEXT};")
        self.points_label.setFont(get_font(FONT_SIZE_SMALL))
        body.addWidget(self.points_label)
        layout.addWidget(sec)

        # Saved Kd fits in this project — reopen / delete (management tool).
        msec = self._make_section("Saved Kd fits in this project")
        mbody = msec.layout()
        self.saved_fits_list = QListWidget()
        self.saved_fits_list.setMaximumHeight(110)
        self.saved_fits_list.itemDoubleClicked.connect(
            lambda _i: self._open_selected_saved_fit())
        mbody.addWidget(self.saved_fits_list)
        self.no_saved_fits_label = create_label(
            "No saved Kd fits yet — run a fit; it is saved with the project automatically "
            "(named after its series).")
        self.no_saved_fits_label.setStyleSheet(f"color: {SECONDARY_TEXT};")
        self.no_saved_fits_label.setFont(get_font(FONT_SIZE_SMALL))
        self.no_saved_fits_label.setWordWrap(True)
        mbody.addWidget(self.no_saved_fits_label)
        mrow = QHBoxLayout()
        mrow.setSpacing(SPACING_SM)
        mrow.addWidget(create_secondary_button(
            "Open selected fit", clicked=self._open_selected_saved_fit, width=160))
        mrow.addWidget(create_secondary_button(
            "Delete selected", clicked=self._delete_selected_saved_fit, width=140))
        mrow.addStretch()
        mbody.addLayout(mrow)
        layout.addWidget(msec)

        # Binding inputs
        sec2 = self._make_section("Binding parameters")
        b2 = sec2.layout()

        # Per-point ligand concentration + intensity scale (populated on CSV load).
        conc_hdr = create_label(
            "Per titration point: ligand concentration [L], and an optional intensity "
            "scale (× height/volume only, e.g. ×2 or ×0.5 to correct a spectrum acquired "
            "with a different number of scans). Positions/CSP are never scaled.")
        conc_hdr.setStyleSheet(f"color: {SECONDARY_TEXT};")
        conc_hdr.setWordWrap(True)
        b2.addWidget(conc_hdr)
        self.conc_form_frame = QFrame()
        self.conc_form = QFormLayout(self.conc_form_frame)
        self.conc_form.setContentsMargins(SPACING_SM, 0, SPACING_SM, 0)
        b2.addWidget(self.conc_form_frame)
        self.conc_spins = []
        self.scale_spins = []

        b2.addWidget(self._make_field_row("Protein concentration [P]₀",
                                          self._make_float_spin("p0_spin", 0.001, 1e6, 50.0, decimals=3, step=10.0)))
        b2.addWidget(self._make_field_row("CSP ¹⁵N scaling α",
                                          self._make_float_spin("alpha_spin", 0.01, 1.0, 0.14, decimals=3, step=0.01)))
        obs = QComboBox()
        obs.addItems(["CSP + Intensity", "CSP only", "Intensity only"])
        self.obs_combo = obs
        b2.addWidget(self._make_field_row("Observable", obs))
        iv = QComboBox()
        iv.addItems(["height", "volume"])
        self.intvalue_combo = iv
        b2.addWidget(self._make_field_row("Intensity from", iv))
        nb = QSpinBox()
        nb.setRange(0, 100000)
        nb.setValue(0)
        nb.setSingleStep(100)
        self.boot_spin = nb
        b2.addWidget(self._make_field_row("Bootstrap iterations (0 = covariance errors)", nb))
        layout.addWidget(sec2)

        # Output
        sec3 = self._make_section("Output")
        b3 = sec3.layout()
        orow = QHBoxLayout()
        orow.setSpacing(SPACING_SM)
        self.outdir_label = QLabel("No output directory selected")
        self.outdir_label.setStyleSheet(
            f"color: {SECONDARY_TEXT}; background-color: {FRAME_BG_COLOR}; padding: 4px 8px; border-radius: 4px;")
        orow.addWidget(self.outdir_label, stretch=1)
        orow.addWidget(create_secondary_button("Choose…", clicked=self._choose_output_dir, width=110))
        b3.addLayout(orow)
        layout.addWidget(sec3)

        # Buttons
        rrow = QHBoxLayout()
        rrow.setSpacing(SPACING_SM)
        self.run_btn = create_primary_button("Run Kd fit", clicked=self._start_analysis, width=180)
        rrow.addWidget(self.run_btn)
        self.viewer_btn = create_secondary_button("Open Kd Viewer", clicked=self._open_viewer, width=180)
        self.viewer_btn.setEnabled(False)
        rrow.addWidget(self.viewer_btn)
        self.open_prev_btn = create_secondary_button("Open previous dataset…",
                                                     clicked=self._open_previous_dataset, width=200)
        rrow.addWidget(self.open_prev_btn)
        self.save_to_project_btn = create_secondary_button(
            "Save analysis to project", clicked=self._save_analysis_to_project, width=200)
        self.save_to_project_btn.setEnabled(False)
        rrow.addWidget(self.save_to_project_btn)
        rrow.addStretch()
        layout.addLayout(rrow)

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.hide()
        layout.addWidget(self.progress_bar)

        self.log_text = QPlainTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setMaximumBlockCount(500)
        self.log_text.setMinimumHeight(120)
        self.log_text.setStyleSheet(
            f"QPlainTextEdit {{ color: {PRIMARY_TEXT}; background-color: {FRAME_BG_COLOR}; "
            f"border-radius: 4px; padding: 4px; }}")
        layout.addWidget(self.log_text, stretch=1)

        self._populate_series_list()
        self._populate_saved_fits()

    def showEvent(self, event):
        # Refresh the lists each time the page is shown (a series may have been run,
        # or a project with saved fits loaded, after the dialog was created).
        super().showEvent(event)
        self._populate_series_list()
        self._populate_saved_fits()

    # ---------- series picker (drag source) ----------

    def _get_available_series(self):
        """Titration-mode series runs in this project, as {'name', 'csv_path'} dicts.
        Filters saved_series to series_mode == 'titration' (set when the series ran)."""
        lunaNMR_main = self._lunanmr_main_window()
        saved_series = getattr(lunaNMR_main, 'saved_series', {}) or {}
        project_path = getattr(lunaNMR_main, 'current_project_path', None)
        out = []
        for name, br in saved_series.items():
            if getattr(br, 'series_mode', 'time') != 'titration':
                continue
            csv_path = ""
            meta = getattr(br, 'metadata', None) or {}
            if meta.get('csv_path'):
                csv_path = meta['csv_path']
            elif meta.get('output_folder'):
                csv_path = os.path.join(meta['output_folder'], "series_analysis_tidy.csv")
            elif project_path:
                csv_path = str(Path(project_path) / ".lunaNMR" / "series_results"
                               / name / "series_analysis_tidy.csv")
            out.append({'name': name, 'csv_path': csv_path})
        return out

    def _populate_series_list(self):
        """Fill the draggable series list; show a hint when there are none."""
        if not hasattr(self, 'series_list_widget'):
            return
        self.series_list_widget.clear()
        series = self._get_available_series()
        self.no_series_label.setVisible(not series)
        self.series_list_widget.setVisible(bool(series))
        for s in series:
            item = QListWidgetItem(s['name'])
            item.setData(Qt.UserRole, s['name'])
            item.setData(Qt.UserRole + 1, s['csv_path'])
            if s['csv_path']:
                item.setToolTip(f"CSV: {s['csv_path']}")
            self.series_list_widget.addItem(item)

    # ---------- UI helpers ----------

    def _make_section(self, title):
        frame = QFrame()
        frame.setStyleSheet(f"QFrame {{ background-color: {FRAME_BG_COLOR}; border-radius: 8px; }}")
        body = create_v_layout(SPACING_XS, (SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM))
        frame.setLayout(body)
        body.addWidget(create_header_label(title, level=3))
        return frame

    def _make_field_row(self, label_text, widget, wide=False):
        row = QFrame()
        rl = create_h_layout(SPACING_SM)
        row.setLayout(rl)
        rl.addWidget(create_label(label_text), stretch=1)
        widget.setMaximumWidth(360 if wide else 180)
        rl.addWidget(widget)
        return row

    def _make_float_spin(self, attr, lo, hi, default, decimals=2, step=1.0):
        spin = QDoubleSpinBox()
        spin.setRange(lo, hi)
        spin.setDecimals(decimals)
        spin.setSingleStep(step)
        spin.setValue(default)
        setattr(self, attr, spin)
        return spin

    # ---------- input ----------

    def _on_series_dropped(self, field_name, series_name, csv_path):
        if csv_path and os.path.exists(csv_path):
            self._set_input_file(csv_path)
            self.series_name = series_name        # identity for the saved-fit name
            self.file_drop.setText(f"Series '{series_name}'  →  {os.path.basename(csv_path)}")

    def _browse_input_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select titration CSV",
            self.main_window.current_dir if hasattr(self.main_window, "current_dir") else "",
            "CSV files (*.csv);;All files (*)")
        if path:
            self._set_input_file(path)
            self.file_drop.setText(os.path.basename(path))

    def _set_input_file(self, path):
        self._log(f"Input set: {path}")
        try:
            from kd_input import load_titration
            points, residues = load_titration(path)
            self.input_file = path        # only commit after a successful load
            # New input = new analysis identity (a drop re-sets series_name after this).
            self.series_name = None
            self._current_analysis_name = None
            # Default the output to where the series data lives, so results land next to
            # the fitted data without the user having to navigate (still overridable).
            self._set_output_dir(os.path.dirname(path))
            self.detected_points = points
            # CSP needs per-point positions; intensity matrices have none.
            import math
            has_positions = any(
                any(not (isinstance(v, float) and math.isnan(v)) for v in r.get('ppm_x', []))
                for r in residues.values())
            note = "" if has_positions else "  —  intensity only (no positions → CSP unavailable)"
            self.points_label.setText(
                f"Detected {len(points)} points  ({len(residues)} residues){note}")
            if not has_positions:
                self.obs_combo.setCurrentIndex(2)  # Intensity only
            if not points:
                self.points_label.setText(
                    f"Detected 0 points — is this a titration series CSV? ({len(residues)} residues)")
            self._build_conc_rows(points)
        except Exception as e:
            # Don't leave the PREVIOUS series committed — otherwise a later "Run Kd fit"
            # silently fits the old series ("shows previous series").
            self.input_file = None
            self.series_name = None
            self._current_analysis_name = None
            self.detected_points = []
            self.points_label.setText(f"Could not read titration points: {e}")

    def _build_conc_rows(self, points):
        """Per titration point: a concentration field (prefilled with the point
        label) and an intensity-scale field (default 1.0)."""
        while self.conc_form.rowCount():
            self.conc_form.removeRow(0)
        self.conc_form.setVerticalSpacing(SPACING_XS)
        self.conc_spins = []
        self.scale_spins = []
        for p in points:
            conc = QDoubleSpinBox()
            conc.setRange(0.0, 1e9)
            conc.setDecimals(4)
            conc.setSingleStep(1.0)
            conc.setValue(float(p))
            conc.setMinimumHeight(26)
            conc.setMaximumWidth(140)
            scale = QDoubleSpinBox()
            scale.setRange(1e-6, 1e9)
            scale.setDecimals(3)
            scale.setSingleStep(0.5)
            scale.setValue(1.0)
            scale.setMinimumHeight(26)
            scale.setMaximumWidth(100)
            cell = QWidget()
            h = QHBoxLayout(cell)
            h.setContentsMargins(0, 0, 0, 0)
            h.setSpacing(SPACING_SM)
            h.addWidget(QLabel("[L]"))
            h.addWidget(conc)
            h.addWidget(QLabel("× int"))
            h.addWidget(scale)
            h.addStretch()
            self.conc_form.addRow(create_label(f"Point  {p}"), cell)
            self.conc_spins.append(conc)
            self.scale_spins.append(scale)
        # Let the form take its natural (tall) height inside the page scroll area.
        self.conc_form_frame.setMinimumHeight(self.conc_form.sizeHint().height())

    def _set_output_dir(self, folder):
        self.output_dir = folder
        self.outdir_label.setText(folder or "No output directory selected")

    def _choose_output_dir(self):
        initial = self.output_dir or (os.path.dirname(self.input_file) if self.input_file else "")
        folder = open_directory_dialog(self, "Select output directory", initial)
        if folder:
            self._set_output_dir(folder)

    # ---------- run ----------

    def _parse_concentrations(self):
        if not self.conc_spins:
            return None
        return [s.value() for s in self.conc_spins]

    def _parse_intensity_scales(self):
        if not self.scale_spins:
            return None
        scales = [s.value() for s in self.scale_spins]
        # only pass them through when the user actually changed something
        return scales if any(abs(s - 1.0) > 1e-9 for s in scales) else None

    def _start_analysis(self):
        if not self.input_file or not os.path.exists(self.input_file):
            show_error(self, "Input required", "Select a CSV file first.")
            return
        concs = self._parse_concentrations()
        if not concs:
            show_error(self, "No concentrations", "Load a CSV so the titration points appear, then set [L].")
            return
        if not self.output_dir:
            self._choose_output_dir()
            if not self.output_dir:
                return

        idx = self.obs_combo.currentIndex()
        observables = {0: ["csp", "intensity"], 1: ["csp"], 2: ["intensity"]}[idx]

        params = KdTitrationFittingParams(
            input_file=self.input_file,
            output_dir=self.output_dir,
            output_prefix="kd",
            concentrations=concs,
            intensity_scales=self._parse_intensity_scales(),
            protein_conc=float(self.p0_spin.value()),
            alpha=float(self.alpha_spin.value()),
            observables=observables,
            intensity_value=self.intvalue_combo.currentText(),
            n_bootstrap=int(self.boot_spin.value()),
        )

        self.progress_bar.setValue(0)
        self.progress_bar.show()
        self.run_btn.setEnabled(False)
        self.viewer_btn.setEnabled(False)
        self._log("Starting Kd fit...")

        self.worker = KdTitrationFittingWorker(params)
        self.worker.progress.connect(self._log)
        self.worker.progress_value.connect(self.progress_bar.setValue)
        self.worker.finished.connect(self._on_finished)
        self.worker.error.connect(self._on_error)
        self.worker.start()

    def _cleanup_worker(self):
        if self.worker is not None:
            self.worker.wait()
            self.worker.deleteLater()
            self.worker = None

    def _on_finished(self, result):
        self.run_btn.setEnabled(True)
        self.progress_bar.setValue(100)
        self.progress_bar.hide()
        self._cleanup_worker()
        n = result.get("n_fitted", 0)
        self.last_json_file = result.get("json_file")
        self.last_json_folder = os.path.dirname(result.get("json_file") or "") or result.get("output_dir")
        self._log(f"Fit complete: {n} residues -> {result.get('json_file')}")
        if n:
            self.viewer_btn.setEnabled(True)
            self.save_to_project_btn.setEnabled(True)
            show_info(self, "Fit complete", f"Fitted {n} residues. Open the viewer to inspect.")
        else:
            show_warning(self, "No fits succeeded", "All residues failed. Check the log.")

    def _on_error(self, error_msg):
        self.run_btn.setEnabled(True)
        self.progress_bar.hide()
        self._cleanup_worker()
        self._log(f"ERROR: {error_msg}")
        show_error(self, "Fit failed", str(error_msg))

    # ---------- viewer ----------

    def _open_viewer(self):
        if not self.last_json_file:
            show_warning(self, "No results", "Run a fit first.")
            return
        try:
            from visualization.kd_titration_fit_viewer import open_kd_titration_viewer
            existing = getattr(self, "_viewer", None)
            if existing is not None and existing.isVisible():
                existing.close()
            self._viewer = open_kd_titration_viewer(parent=self, json_file=self.last_json_file)
        except Exception as e:
            show_error(self, "Viewer error", f"Could not open viewer:\n{e}")

    def _open_previous_dataset(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select *_kd_fit_data.json",
            self.output_dir or "", "JSON files (*_kd_fit_data.json);;All files (*)")
        if not path:
            return
        try:
            from visualization.kd_titration_fit_viewer import open_kd_titration_viewer
            self._viewer = open_kd_titration_viewer(parent=self, json_file=path)
            self.last_json_file = path
            self.viewer_btn.setEnabled(True)
            self.save_to_project_btn.setEnabled(True)
        except Exception as e:
            show_error(self, "Viewer error", f"Could not open viewer:\n{e}")

    # Generic per-series output filenames — never use these as the saved-fit name.
    _GENERIC_STEMS = ("series_analysis_tidy", "comprehensive_peak_tracking")

    def _analysis_base_name(self):
        """A meaningful identity for the saved fit: the opened fit's name (preserved on
        re-save), else the dropped series name, else derived from the input path —
        falling back to the SERIES FOLDER when the file is a generic series output
        (e.g. .../series_results_ERDj6/series_analysis_tidy.csv -> 'ERDj6')."""
        if self._current_analysis_name:
            return self._current_analysis_name
        if self.series_name:
            return self.series_name
        if not self.input_file:
            return "analysis"
        p = Path(self.input_file)
        if p.stem in self._GENERIC_STEMS:
            folder = p.parent.name
            for pre in ("series_results_", "series_results"):
                if folder.startswith(pre):
                    folder = folder[len(pre):]
                    break
            return folder or p.stem
        return p.stem

    def _unique_analysis_name(self, base, existing):
        """A name not already in `existing`; suffix _2, _3, … on collision."""
        base = base or "analysis"
        if base not in existing:
            return base
        i = 2
        while f"{base}_{i}" in existing:
            i += 1
        return f"{base}_{i}"

    def _lunanmr_main_window(self):
        """The real LunaNMR main window. In the app the page's main_window is the
        hosting dialog, which holds the main window as .main_window (BasePage
        convention); in isolation the main window is passed directly."""
        mw = self.main_window
        return getattr(mw, 'main_window', None) or mw

    def _store_current_analysis(self, overwrite=False):
        """Snapshot the current fit + params into main_window.kd_analyses, named after
        the source series. overwrite=True upserts that name (auto-save on project save);
        otherwise a new snapshot is suffixed on collision. Returns the name, or None if
        there is no current fit."""
        if not self.last_json_file or not os.path.exists(self.last_json_file):
            return None
        import json
        fit_data = json.loads(Path(self.last_json_file).read_text())
        mw = self._lunanmr_main_window()
        if getattr(mw, 'kd_analyses', None) is None:
            mw.kd_analyses = {}
        base = self._analysis_base_name()
        name = base if overwrite else self._unique_analysis_name(base, mw.kd_analyses)
        meta = self.get_session_state()
        meta['input_basename'] = os.path.basename(self.input_file) if self.input_file else None
        viewer = getattr(self, "_viewer", None)
        if viewer is not None and viewer.isVisible():
            try:
                meta['comparison'] = {
                    'ref': viewer.ref_point_combo.currentIndex(),
                    'cmp': viewer.cmp_point_combo.currentIndex(),
                    'observable': viewer.obs_combo.currentIndex()}
            except Exception:
                pass
        mw.kd_analyses[name] = {'fit_data': fit_data, 'meta': meta}
        return name

    def _save_analysis_to_project(self):
        name = self._store_current_analysis()       # explicit snapshot (suffixed)
        if name is None:
            show_warning(self, "No results", "Run or open a fit first.")
            return
        self._populate_saved_fits()
        self._log(f"Saved analysis '{name}' to project (save the project to persist).")
        show_info(self, "Saved to project",
                  f"Analysis '{name}' added. Save the project (File ▸ Save) to write it to disk.")

    def ensure_current_saved(self):
        """Auto-capture the fit currently on screen under its source-series name
        (upsert), so saving the project always keeps the visible fit. No-op if none."""
        return self._store_current_analysis(overwrite=True)

    # ---------- saved-fits management (reopen / delete) ----------

    def _populate_saved_fits(self):
        """List the Kd fits saved in this project (main_window.kd_analyses)."""
        if not hasattr(self, 'saved_fits_list'):
            return
        self.saved_fits_list.clear()
        mw = self._lunanmr_main_window()
        names = sorted((getattr(mw, 'kd_analyses', {}) or {}).keys())
        self.no_saved_fits_label.setVisible(not names)
        self.saved_fits_list.setVisible(bool(names))
        for n in names:
            self.saved_fits_list.addItem(n)

    def _open_selected_saved_fit(self):
        item = self.saved_fits_list.currentItem()
        if item is None:
            show_warning(self, "No fit selected", "Select a saved Kd fit to open.")
            return
        mw = self._lunanmr_main_window()
        entry = (getattr(mw, 'kd_analyses', {}) or {}).get(item.text())
        if entry is not None:
            self.open_saved_analysis(entry, name=item.text())

    def _delete_selected_saved_fit(self):
        item = self.saved_fits_list.currentItem()
        if item is None:
            return
        name = item.text()
        if QMessageBox.question(
                self, "Delete saved fit",
                f"Remove the saved Kd fit '{name}' from this project?\n"
                "(applies when you next save the project)") != QMessageBox.Yes:
            return
        mw = self._lunanmr_main_window()
        (getattr(mw, 'kd_analyses', {}) or {}).pop(name, None)
        self._populate_saved_fits()

    def open_saved_analysis(self, entry, name=None):
        """Load a saved analysis (from the project) back into the page + viewer:
        restore its parameters, point the viewer at its fit JSON, and reapply the
        comparison-view selections. `name` is the analysis's key, remembered so a
        re-save upserts the same entry. Returns the viewer, or None on failure."""
        if not isinstance(entry, dict):
            return None
        self._current_analysis_name = name    # keep its name on re-save
        path = entry.get('fit_data_path')
        if not path or not os.path.exists(path):
            # in-memory only (saved this session, project not written yet)
            import json, tempfile
            tmp = tempfile.NamedTemporaryFile(
                'w', suffix='_kd_fit_data.json', delete=False)
            json.dump(entry.get('fit_data', {}), tmp)
            tmp.close()
            path = tmp.name
        meta = entry.get('meta') or {}
        if meta:
            self.restore_session_state(meta)
        self.last_json_file = path
        self.viewer_btn.setEnabled(True)
        self.save_to_project_btn.setEnabled(True)
        from visualization.kd_titration_fit_viewer import open_kd_titration_viewer
        existing = getattr(self, "_viewer", None)
        if existing is not None and existing.isVisible():
            existing.close()
        self._viewer = open_kd_titration_viewer(parent=self, json_file=path)
        cmp = meta.get('comparison')
        if cmp and self._viewer is not None:
            self._viewer.obs_combo.setCurrentIndex(int(cmp.get('observable', 0)))
            self._viewer.view_combo.setCurrentIndex(3)        # "Ref vs point (bars)"
            self._viewer.ref_point_combo.setCurrentIndex(int(cmp.get('ref', 0)))
            self._viewer.cmp_point_combo.setCurrentIndex(int(cmp.get('cmp', 1)))
        return self._viewer

    # ---------- session state (project save/load) ----------

    def get_session_state(self) -> dict:
        """Serializable state for project save."""
        return {
            'input_file': self.input_file,
            'output_dir': self.output_dir,
            'detected_points': list(self.detected_points),
            'concentrations': [s.value() for s in self.conc_spins],
            'intensity_scales': [s.value() for s in self.scale_spins],
            'p0': self.p0_spin.value(),
            'alpha': self.alpha_spin.value(),
            'observable': self.obs_combo.currentIndex(),
            'intensity_value': self.intvalue_combo.currentText(),
            'bootstrap_iter': self.boot_spin.value(),
            'last_json_file': self.last_json_file,
            'last_json_folder': self.last_json_folder,
        }

    def restore_session_state(self, state: dict):
        """Restore state produced by get_session_state()."""
        if not state:
            return
        self.input_file = state.get('input_file')
        self.output_dir = state.get('output_dir')
        self.detected_points = list(state.get('detected_points', []))
        self.last_json_file = state.get('last_json_file')
        self.last_json_folder = state.get('last_json_folder')

        # Rebuild the per-point rows before restoring their values.
        self._build_conc_rows(self.detected_points)
        for spin, val in zip(self.conc_spins, state.get('concentrations', [])):
            spin.setValue(val)
        for spin, val in zip(self.scale_spins, state.get('intensity_scales', [])):
            spin.setValue(val)

        self.p0_spin.setValue(state.get('p0', 50.0))
        self.alpha_spin.setValue(state.get('alpha', 0.14))
        self.obs_combo.setCurrentIndex(state.get('observable', 0))
        self.intvalue_combo.setCurrentText(state.get('intensity_value', 'height'))
        self.boot_spin.setValue(state.get('bootstrap_iter', 0))

        if self.input_file:
            self.file_drop.setText(os.path.basename(self.input_file))
            self.points_label.setText(
                f"Detected {len(self.detected_points)} points")
        if self.output_dir:
            self.outdir_label.setText(self.output_dir)
        if self.last_json_file:
            self.viewer_btn.setEnabled(True)

    def _log(self, msg):
        self.log_text.appendPlainText(str(msg))
