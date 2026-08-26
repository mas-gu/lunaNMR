# ABOUTME: dynamiXs page for Kd / titration binding fits (CSP and intensity ratio).
# ABOUTME: 1:1 quadratic isotherm, per-residue + global Kd; mirrors the methyl-T2 page.

import os
import sys
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFrame, QHBoxLayout, QVBoxLayout, QLabel, QFileDialog, QDoubleSpinBox, QSpinBox,
    QComboBox, QFormLayout, QProgressBar, QPlainTextEdit, QScrollArea, QWidget,
    QListWidget, QListWidgetItem, QMessageBox, QInputDialog,
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


def resolve_series_csv(name, meta, project_path):
    """Path to a series' series_analysis_tidy.csv for the Kd input.

    Prefers the copy bundled inside the project (so a moved/renamed source run
    folder still resolves), then the recorded csv_path, then output_folder.
    Returns the first candidate that exists on disk; if none exist, returns the
    recorded csv_path (or best-effort candidate) for display.
    """
    candidates = []
    if project_path:
        candidates.append(os.path.join(str(project_path), "series_results",
                                       name, "series_analysis_tidy.csv"))
    if meta.get('csv_path'):
        candidates.append(meta['csv_path'])
    if meta.get('output_folder'):
        candidates.append(os.path.join(meta['output_folder'], "series_analysis_tidy.csv"))
    for c in candidates:
        if c and os.path.exists(c):
            return c
    return meta.get('csv_path') or (candidates[0] if candidates else "")


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
            "Load for refit", clicked=self._refit_selected_saved_fit, width=140))
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
        cu = QComboBox()
        cu.addItems(["absolute", "equivalents of [P]0"])
        self.conc_units_combo = cu
        b2.addWidget(self._make_field_row("Concentration units", cu))
        nb = QSpinBox()
        nb.setRange(0, 100000)
        nb.setValue(0)
        nb.setSingleStep(100)
        self.boot_spin = nb
        b2.addWidget(self._make_field_row("Bootstrap iterations (0 = covariance errors)", nb))
        prow = QHBoxLayout()
        prow.setSpacing(SPACING_SM)
        prow.addWidget(create_secondary_button("Load params (JSON)…", clicked=self._load_params_dialog))
        prow.addWidget(create_secondary_button("Save params (JSON)…", clicked=self._save_params_dialog))
        prow.addStretch()
        b2.addLayout(prow)
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
            meta = getattr(br, 'metadata', None) or {}
            csv_path = resolve_series_csv(name, meta, project_path)
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
            # Default the output to a 'kd_analysis' subfolder of where the series data
            # lives (the peak-integration folder), so Kd results are grouped together
            # next to the fitted data without cluttering it (still overridable).
            self._set_output_dir(os.path.join(os.path.dirname(path), "kd_analysis"))
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
            self._maybe_apply_params_json(path, points)
            if not has_positions:
                self.obs_combo.setCurrentIndex(2)  # no positions -> intensity only wins
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

    # ---------- importable binding parameters ----------

    def _apply_params(self, params, n_points):
        """Populate the Binding-parameters panel from a normalized params dict. Per-point
        lists ([L], intensity scales) are applied only when their length matches the
        detected points, so a mismatched file can't silently corrupt the run."""
        from kd_params import observables_to_combo_index
        concs = params.get('concentrations')
        if concs and len(concs) == n_points:
            for spin, val in zip(self.conc_spins, concs):
                spin.setValue(float(val))
        scales = params.get('intensity_scales')
        if scales and len(scales) == n_points:
            for spin, val in zip(self.scale_spins, scales):
                spin.setValue(float(val))
        self.p0_spin.setValue(float(params.get('protein_conc', 50.0)))
        self.alpha_spin.setValue(float(params.get('alpha', 0.14)))
        self.obs_combo.setCurrentIndex(observables_to_combo_index(params.get('observables')))
        self.intvalue_combo.setCurrentText(str(params.get('intensity_value', 'height')))
        self.conc_units_combo.setCurrentIndex(
            1 if str(params.get('conc_units', 'absolute')).lower().startswith('eq') else 0)
        self.boot_spin.setValue(int(params.get('n_bootstrap', 0)))

    def _maybe_apply_params_json(self, csv_path, points):
        """Auto-detect a sibling binding-parameters JSON (or a prior fit JSON) next to the
        input CSV and populate the panel. Values stay editable; absent → manual entry."""
        try:
            from kd_params import find_params_source, load_params
            src = find_params_source(csv_path)
            if not src:
                return
            self._apply_params(load_params(src), n_points=len(points))
            self._log(f"Loaded binding parameters from {os.path.basename(src)}")
        except Exception as e:
            self._log(f"Could not read binding-parameters JSON: {e}")

    def _gather_params_dict(self):
        """The current panel values as a params dict (for saving a params JSON)."""
        from kd_params import combo_index_to_observables
        return {
            'points': list(self.detected_points),
            'concentrations': self._parse_concentrations(),
            'intensity_scales': [s.value() for s in self.scale_spins] or None,
            'protein_conc': float(self.p0_spin.value()),
            'alpha': float(self.alpha_spin.value()),
            'observables': combo_index_to_observables(self.obs_combo.currentIndex()),
            'intensity_value': self.intvalue_combo.currentText(),
            'conc_units': ('equivalents' if self.conc_units_combo.currentIndex() == 1
                           else 'absolute'),
            'n_bootstrap': int(self.boot_spin.value()),
        }

    def _load_params_dialog(self):
        if not self.conc_spins:
            show_error(self, "Load a CSV first",
                       "Load a titration CSV so the points appear, then load parameters.")
            return
        start = self.output_dir or (os.path.dirname(self.input_file) if self.input_file else "")
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Kd binding parameters", start,
            "Params JSON (*_kd_params.json *_kd_fit_data.json *.json);;All files (*)")
        if not path:
            return
        try:
            from kd_params import load_params
            self._apply_params(load_params(path), n_points=len(self.conc_spins))
            self._log(f"Loaded binding parameters from {os.path.basename(path)}")
        except Exception as e:
            show_error(self, "Could not load parameters", str(e))

    def _save_params_dialog(self):
        if not self.conc_spins:
            show_error(self, "Nothing to save",
                       "Load a titration CSV and set the parameters first.")
            return
        from kd_params import dump_params_json, PARAMS_SUFFIX
        folder = self.output_dir or (os.path.dirname(self.input_file) if self.input_file else "")
        default = os.path.join(folder, f"{self._analysis_base_name()}{PARAMS_SUFFIX}")
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Kd binding parameters", default, "Params JSON (*.json)")
        if not path:
            return
        try:
            dump_params_json(path, self._gather_params_dict())
            self._log(f"Saved binding parameters to {os.path.basename(path)}")
        except Exception as e:
            show_error(self, "Could not save parameters", str(e))

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

        from kd_params import combo_index_to_observables
        observables = combo_index_to_observables(self.obs_combo.currentIndex())

        params = KdTitrationFittingParams(
            input_file=self.input_file,
            output_dir=self.output_dir,
            output_prefix=self._analysis_base_name(),
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
        warnings_ = result.get("quality_warnings") or []
        for w in warnings_:
            self._log(w)
        if n:
            self.viewer_btn.setEnabled(True)
            self.save_to_project_btn.setEnabled(True)
            if warnings_:
                # A dataset the titration cannot answer must not be reported as a plain
                # success — the number would look ordinary and mean nothing.
                show_warning(self, "Fit complete, but the result is not usable",
                             "\n\n".join(warnings_))
            else:
                show_info(self, "Fit complete",
                          f"Fitted {n} residues. Open the viewer to inspect.")
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

    def _json_browse_dir(self):
        """Where fit JSONs live: the data/ subfolder of the output, when it exists.

        last_json_folder is already the JSON's own directory, so it is used as-is;
        only an output_dir needs the data/ suffix appended.
        """
        if self.last_json_folder and os.path.isdir(self.last_json_folder):
            return self.last_json_folder
        base = self.output_dir or ""
        nested = os.path.join(base, "data") if base else ""
        return nested if nested and os.path.isdir(nested) else base

    def _open_previous_dataset(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select *_kd_fit_data.json",
            self._json_browse_dir(), "JSON files (*_kd_fit_data.json);;All files (*)")
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
        (e.g. .../series_results_<name>/series_analysis_tidy.csv -> '<name>')."""
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

    def _store_current_analysis(self, overwrite=False, name=None):
        """Snapshot the current fit + params into main_window.kd_analyses. `name` stores
        under that exact name (upsert); otherwise the source-series base name is used,
        upserted when overwrite=True or suffixed on collision. Returns the stored name,
        or None if there is no current fit."""
        if not self.last_json_file or not os.path.exists(self.last_json_file):
            return None
        import json
        fit_data = json.loads(Path(self.last_json_file).read_text())
        mw = self._lunanmr_main_window()
        if getattr(mw, 'kd_analyses', None) is None:
            mw.kd_analyses = {}
        base = self._analysis_base_name()
        if name is None:
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
        if not self.last_json_file or not os.path.exists(self.last_json_file):
            show_warning(self, "No results", "Run or open a fit first.")
            return
        mw = self._lunanmr_main_window()
        existing = getattr(mw, 'kd_analyses', None) or {}
        name = self._analysis_base_name()
        if name in existing:
            name = self._prompt_name_collision(name)
            if name is None:                        # user cancelled
                return
        name = self._store_current_analysis(name=name)
        if name is None:
            show_warning(self, "No results", "Run or open a fit first.")
            return
        self._populate_saved_fits()
        self._log(f"Saved analysis '{name}' to project (save the project to persist).")
        show_info(self, "Saved to project",
                  f"Analysis '{name}' added. Save the project (File ▸ Save) to write it to disk.")

    def _prompt_name_collision(self, base):
        """A Kd analysis named `base` already exists. Ask whether to replace it or save
        under a different name. Returns the name to store under (the base to replace, or
        a new name), or None to cancel."""
        box = QMessageBox(self)
        box.setWindowTitle("Analysis already exists")
        box.setText(f"A saved Kd analysis named '{base}' already exists in this project.")
        box.setInformativeText("Replace it, or save under a different name?")
        replace_btn = box.addButton("Replace", QMessageBox.AcceptRole)
        rename_btn = box.addButton("Save as…", QMessageBox.ActionRole)
        box.addButton(QMessageBox.Cancel)
        box.exec()
        clicked = box.clickedButton()
        if clicked is replace_btn:
            return base
        if clicked is rename_btn:
            return self._prompt_new_analysis_name(base)
        return None

    def _prompt_new_analysis_name(self, base):
        """Ask for a new analysis name, pre-filled with a collision-free suggestion.
        Re-prompts on a blank name and confirms before overwriting an existing one.
        Returns the chosen name, or None if the user cancels."""
        existing = getattr(self._lunanmr_main_window(), 'kd_analyses', None) or {}
        suggestion = self._unique_analysis_name(base, existing)
        while True:
            name, ok = QInputDialog.getText(
                self, "Save analysis as", "Name for this Kd analysis:", text=suggestion)
            if not ok:
                return None
            name = name.strip()
            if not name:
                continue
            if name in existing and QMessageBox.question(
                    self, "Name in use", f"'{name}' also exists. Replace it?") != QMessageBox.Yes:
                continue
            return name

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

    def _refit_selected_saved_fit(self):
        item = self.saved_fits_list.currentItem()
        if item is None:
            show_warning(self, "No fit selected", "Select a saved Kd fit to load for refit.")
            return
        mw = self._lunanmr_main_window()
        entry = (getattr(mw, 'kd_analyses', {}) or {}).get(item.text())
        if entry is not None:
            self.load_saved_for_refit(entry, name=item.text())

    def _entry_fit_data(self, entry):
        """The fit JSON dict for a saved analysis, from its bundled file or the
        in-memory copy (saved this session, project not yet written)."""
        import json
        path = entry.get('fit_data_path')
        if path and os.path.exists(path):
            try:
                return json.loads(Path(path).read_text())
            except Exception:
                return None
        return entry.get('fit_data')

    def _reconstruct_tidy_csv(self, fit_data):
        """Rebuild a minimal series_analysis_tidy.csv from a fit JSON's embedded per-point
        series, so a saved fit can be refit when the original CSV is gone. Heights/volumes
        are the stored (already-scaled) values — the caller resets intensity_scales to 1.0
        to avoid double-scaling. Returns the temp CSV path, or None if there is no series."""
        import csv as _csv
        import tempfile
        meta = fit_data.get('metadata', {}) or {}
        points = meta.get('points') or meta.get('concentrations') or []
        if not points or not any(f.get('series') for f in fit_data.get('fits', [])):
            return None
        rows = []
        for f in fit_data.get('fits', []):
            res = f.get('residue')
            s = f.get('series') or {}

            def val(key, i):
                v = s.get(key)
                return v[i] if isinstance(v, list) and i < len(v) and v[i] is not None else ''
            for i in range(len(points)):
                rows.append([str(points[i]), res, val('ppm_x', i), val('ppm_y', i),
                             val('height', i), val('volume', i)])
        tmp = tempfile.NamedTemporaryFile(
            'w', suffix='_series_analysis_tidy.csv', delete=False, newline='')
        w = _csv.writer(tmp)
        w.writerow(['spectrum_name', 'assignment', 'ppm_x', 'ppm_y', 'height', 'volume'])
        w.writerows(rows)
        tmp.close()
        return tmp.name

    def load_saved_for_refit(self, entry, name=None):
        """Load a saved analysis's binding parameters + input data back into the panel so
        it can be edited and re-run. Prefers the original CSV; if it's gone, rebuilds the
        input from the fit's embedded series (resetting intensity scales to 1.0). Does not
        open the viewer — leaves focus on the editable panel."""
        if not isinstance(entry, dict):
            return
        self._current_analysis_name = name
        meta = entry.get('meta') or {}
        if meta:
            self.restore_session_state(meta)
        if not self.input_file or not os.path.exists(self.input_file):
            fit_data = self._entry_fit_data(entry)
            rebuilt = self._reconstruct_tidy_csv(fit_data) if fit_data else None
            if not rebuilt:
                show_warning(self, "Cannot refit",
                             "The original input CSV was not found and the saved fit has no "
                             "embedded data to rebuild from.")
                return
            self.input_file = rebuilt
            self.file_drop.setText(os.path.basename(rebuilt) + "  (rebuilt from saved data)")
            for s in self.scale_spins:
                s.setValue(1.0)
            self._log("Original CSV not found — rebuilt input from the saved series; "
                      "intensity scales reset to 1.0 (stored data is already scaled).")
        self._log(f"Loaded '{name}' for refit — edit parameters, then click Run Kd fit.")
        show_info(self, "Loaded for refit",
                  f"'{name}' loaded. Adjust the binding parameters and click Run Kd fit.")

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
