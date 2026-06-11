# ABOUTME: Methyl bi-exp T2 fit viewer with outlier rejection (mirror of FitViewer for methyl_T2 data).
# ABOUTME: Single panel per residue (bi-exp curve + data); toggleable per-residue eta = T2a/T2b bar chart.

import json
import os
import re
import sys
from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt, QEvent
from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QFrame,
    QLabel, QPushButton, QCheckBox, QFileDialog, QSizePolicy,
    QListWidget, QListWidgetItem, QAbstractItemView, QScrollArea,
    QRadioButton, QButtonGroup,
)

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

sys.path.append(str(Path(__file__).parent.parent))
from constants import (
    BG_COLOR, PANEL_BG_COLOR, FRAME_BG_COLOR,
    PRIMARY_TEXT, SECONDARY_TEXT,
    PRIMARY_BUTTON_BG, PRIMARY_BUTTON_HOVER,
    SPACING_XS, SPACING_SM, SPACING_MD,
    FONT_SIZE_BODY, FONT_SIZE_SMALL,
)
from gui_components import (
    create_primary_button, create_secondary_button, create_label,
    get_font, show_info, show_error, show_warning,
)

# Refit helpers shared with the mono-exp viewer (sidecar IO + JSON updater).
sys.path.insert(0, str(Path(__file__).parent.parent / "dynamiXs_T1_T2"))
from refit import (
    load_outliers, save_outliers,
    refit_residue_methyl, update_json_fit_entry, update_tsv_row_methyl,
)


JSON_GLOB = "*_methylT2_fit_data.json"


class MethylT2FitViewer(QMainWindow):
    """Interactive viewer for methyl bi-exp T2 fits.

    Plots the fitted bi-exp curve over the raw points for one residue at a time,
    with the same outlier-rejection workflow as the mono-exp FitViewer
    (Edit Outliers checkbox + Re-fit / Reset buttons + sidecar persistence).

    A toggleable secondary panel shows eta = T2a / T2b across all loaded residues.
    """

    SUMMARY_QUANTITIES = [
        # (key, label, err_key)
        ("T2_avg", "T2_avg",            "T2_avg_err"),
        ("dT2",    "ΔT2 = T2a − T2b",   "dT2_err"),
        ("t2_a",   "T2a (slow)",        "t2_a_err"),
        ("t2_b",   "T2b (fast)",        "t2_b_err"),
        ("eta",    "η = T2a / T2b",     "eta_err"),
    ]

    def __init__(self, parent=None, json_folder=None, summary_only=False):
        super().__init__(parent)
        self.summary_only = bool(summary_only)
        self.setWindowTitle(
            "Methyl T2 — Across-Residues" if self.summary_only
            else "Methyl T2 Fit Viewer"
        )
        self.setMinimumSize(1300, 850)
        # Set both bg and color: partial stylesheet (bg only) lets descendants
        # fall through to system palette for color (white on macOS Dark Mode).
        self.setStyleSheet(f"background-color: {BG_COLOR}; color: {PRIMARY_TEXT};")

        # State
        self.json_folder = json_folder
        self.data = {}             # {file_key: parsed json}
        self.json_paths = {}       # {file_key: Path}
        self.exclusions = {}       # {file_key: {residue: [indices]}}
        self.available_residues = []
        self.edit_outliers = False
        self.show_t2_avg = False
        self.show_dT2 = False
        self.show_t2a = False
        self.show_t2b = False
        self.show_eta = False
        self._pick_registry = {}
        # In summary-only mode the residue list becomes a checklist for
        # inclusion/exclusion in the bar chart; track which residues are kept.
        self.summary_quantity = "T2_avg"
        self.included_residues = set()
        # Snapshot of loaded JSON file mtimes; used to auto-refresh on focus
        # when a sibling viewer (e.g. the per-residue editor) refits a residue
        # and writes new values to disk.
        self.json_mtimes = {}

        self._create_ui()

        if json_folder and os.path.isdir(json_folder):
            self._load_json_folder(json_folder)

    # ---------- UI ----------

    def _create_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        main = QHBoxLayout(central)
        main.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)
        main.setSpacing(SPACING_SM)

        # Plot panel (left)
        left = QFrame()
        left.setStyleSheet(f"QFrame {{ background-color: {PANEL_BG_COLOR}; border-radius: 8px; }}")
        self._create_plot_panel(left)
        main.addWidget(left, 3)

        # Navigator (right)
        right = QFrame()
        right.setStyleSheet(f"QFrame {{ background-color: {PANEL_BG_COLOR}; border-radius: 8px; }}")
        right.setFixedWidth(330)
        self._create_navigator_panel(right)
        main.addWidget(right, 0)

    def _create_plot_panel(self, parent):
        layout = QVBoxLayout(parent)
        layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)
        layout.setSpacing(SPACING_SM)

        title = create_label(
            "Methyl T2 — across-residues" if self.summary_only
            else "Methyl T2 — bi-exponential fit"
        )
        title.setFont(get_font(18, bold=True))
        layout.addWidget(title)

        plot_frame = QFrame()
        plot_frame.setStyleSheet(f"background-color: {FRAME_BG_COLOR}; border-radius: 8px;")
        pl = QVBoxLayout(plot_frame)
        pl.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)

        self.figure = Figure(figsize=(9, 9), facecolor=PANEL_BG_COLOR)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.canvas.mpl_connect("pick_event", self._on_pick)
        pl.addWidget(self.canvas)

        self.toolbar = NavigationToolbar2QT(self.canvas, plot_frame)
        pl.addWidget(self.toolbar)

        layout.addWidget(plot_frame, 1)
        self._show_blank_state()

    def _create_navigator_panel(self, parent):
        layout = QVBoxLayout(parent)
        layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)
        layout.setSpacing(SPACING_SM)

        title = create_label("Navigator")
        title.setFont(get_font(16, bold=True))
        layout.addWidget(title)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)
        scroll.setStyleSheet("background-color: transparent;")
        body = QWidget()
        bl = QVBoxLayout(body)
        bl.setContentsMargins(0, 0, 0, 0)
        bl.setSpacing(SPACING_MD)

        # Folder section
        folder_section = self._section_frame("JSON Data Folder:")
        fb = folder_section.layout()
        self.folder_label = QLabel(self.json_folder or "No folder selected")
        self.folder_label.setStyleSheet(
            f"QLabel {{ background-color: {BG_COLOR}; color: {SECONDARY_TEXT}; "
            f"padding: 6px; border-radius: 4px; }}"
        )
        self.folder_label.setWordWrap(True)
        fb.addWidget(self.folder_label)
        fb.addWidget(create_primary_button("Browse Folder", clicked=self._browse_folder, width=280))
        bl.addWidget(folder_section)

        # Residue list — single-select navigator in normal mode, checklist
        # for inclusion/exclusion in summary-only mode.
        res_label = "Include residues:" if self.summary_only else "Residue:"
        res_section = self._section_frame(res_label)
        rb = res_section.layout()
        self.residue_list = QListWidget()
        self.residue_list.setMinimumHeight(280)
        self.residue_list.setSelectionMode(QAbstractItemView.SingleSelection)
        self.residue_list.setStyleSheet(
            f"QListWidget {{ background-color: {BG_COLOR}; border: none; border-radius: 4px; }}"
            f"QListWidget::item {{ padding: 6px 8px; border-radius: 4px; }}"
            f"QListWidget::item:selected {{ background-color: {PRIMARY_BUTTON_BG}; color: white; }}"
            f"QListWidget::item:hover:!selected {{ background-color: {PRIMARY_BUTTON_HOVER}; }}"
        )
        if self.summary_only:
            self.residue_list.itemChanged.connect(self._on_inclusion_toggled)
        else:
            self.residue_list.itemClicked.connect(self._on_residue_clicked)
            self.residue_list.currentRowChanged.connect(self._on_residue_row_changed)
        rb.addWidget(self.residue_list)
        if self.summary_only:
            qb_row = QHBoxLayout()
            qb_row.addWidget(create_secondary_button(
                "Include all", clicked=self._include_all, width=110))
            qb_row.addWidget(create_secondary_button(
                "Exclude all", clicked=self._exclude_all, width=110))
            rb.addLayout(qb_row)
            rb.addWidget(create_secondary_button(
                "Exclude flagged (mono-exp regime)",
                clicked=self._exclude_flagged, width=280))
        bl.addWidget(res_section)

        if self.summary_only:
            # Single quantity picker (radio group).
            quantity_section = self._section_frame("Show on plot:")
            qs = quantity_section.layout()
            self._quantity_button_group = QButtonGroup(self)
            self._quantity_radios = {}
            for key, label, _ in self.SUMMARY_QUANTITIES:
                rb_btn = QRadioButton(label)
                if key == self.summary_quantity:
                    rb_btn.setChecked(True)
                rb_btn.toggled.connect(self._on_quantity_changed)
                self._quantity_button_group.addButton(rb_btn)
                self._quantity_radios[key] = rb_btn
                qs.addWidget(rb_btn)
            bl.addWidget(quantity_section)
        else:
            # Across-residues summary panel toggles (multi-select stacked panels).
            summary_section = self._section_frame("Across-residues panels:")
            sb = summary_section.layout()
            self.t2_avg_checkbox = QCheckBox("Show T2_avg per residue")
            self.t2_avg_checkbox.toggled.connect(self._toggle_summary)
            sb.addWidget(self.t2_avg_checkbox)
            self.dT2_checkbox = QCheckBox("Show ΔT2 = T2a − T2b per residue")
            self.dT2_checkbox.toggled.connect(self._toggle_summary)
            sb.addWidget(self.dT2_checkbox)
            self.t2a_checkbox = QCheckBox("Show T2a (slow) per residue")
            self.t2a_checkbox.toggled.connect(self._toggle_summary)
            sb.addWidget(self.t2a_checkbox)
            self.t2b_checkbox = QCheckBox("Show T2b (fast) per residue")
            self.t2b_checkbox.toggled.connect(self._toggle_summary)
            sb.addWidget(self.t2b_checkbox)
            self.eta_checkbox = QCheckBox("Show η = T2a / T2b per residue")
            self.eta_checkbox.toggled.connect(self._toggle_summary)
            sb.addWidget(self.eta_checkbox)
            bl.addWidget(summary_section)

            # Outlier editing — only in fit-viewer mode (no per-residue plot
            # in summary-only mode means there's nothing to click on).
            out_section = self._section_frame("Outlier Editing:")
            ob = out_section.layout()
            self.edit_outliers_checkbox = QCheckBox("Edit Outliers (click points to toggle)")
            self.edit_outliers_checkbox.toggled.connect(self._toggle_edit_outliers)
            ob.addWidget(self.edit_outliers_checkbox)

            self.refit_btn = create_primary_button("Re-fit Residue", clicked=self._on_refit_clicked, width=280)
            self.refit_btn.setToolTip(
                "Re-fit the selected residue with currently-rejected points excluded.\n"
                "Errors are estimated analytically (covariance matrix) for speed.\n"
                "Updates JSON, TSV (if present), and in-memory data."
            )
            ob.addWidget(self.refit_btn)
            self.reset_btn = create_secondary_button("Reset Residue", clicked=self._on_reset_clicked, width=280)
            self.reset_btn.setToolTip(
                "Clear all rejected points for this residue and re-fit on full data\n"
                "(analytical errors)."
            )
            ob.addWidget(self.reset_btn)
            bl.addWidget(out_section)

        bl.addStretch()
        scroll.setWidget(body)
        layout.addWidget(scroll, 1)

    def _section_frame(self, title: str) -> QFrame:
        frame = QFrame()
        frame.setStyleSheet(f"QFrame {{ background-color: {FRAME_BG_COLOR}; border-radius: 8px; }}")
        v = QVBoxLayout(frame)
        v.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)
        v.setSpacing(SPACING_XS)
        lbl = create_label(title)
        lbl.setFont(get_font(FONT_SIZE_BODY, bold=True))
        v.addWidget(lbl)
        return frame

    # ---------- Data loading ----------

    def _browse_folder(self):
        initial = self.json_folder or os.getcwd()
        folder = QFileDialog.getExistingDirectory(self, "Select JSON folder", initial)
        if folder:
            self._load_json_folder(folder)

    def _load_json_folder(self, folder_path: str, silent: bool = False):
        self.json_folder = folder_path
        self.folder_label.setText(folder_path)
        self.data.clear()
        self.json_paths.clear()
        self.exclusions.clear()
        self.available_residues = []

        residues_set = set()
        try:
            for json_file in Path(folder_path).glob(JSON_GLOB):
                try:
                    payload = json.loads(json_file.read_text())
                except Exception as e:
                    print(f"Skipping {json_file}: {e}")
                    continue
                key = json_file.stem.replace("_fit_data", "")
                self.data[key] = payload
                self.json_paths[key] = json_file
                self.exclusions[key] = load_outliers(json_file)
                for fit in payload.get("fits", []):
                    residues_set.add(fit["residue"])

            self._snapshot_mtimes()

            if not self.data:
                if not silent:
                    show_warning(self, "No data",
                                 f"No '{JSON_GLOB}' files found in {folder_path}.")
                self._populate_residue_list([])
                return

            self.available_residues = sorted(residues_set, key=self._residue_sort_key)
            self._populate_residue_list(self.available_residues)
            if not silent:
                show_info(self, "Loaded", f"{len(self.data)} dataset(s), "
                          f"{len(self.available_residues)} residue(s).")
        except Exception as e:
            if not silent:
                show_error(self, "Load error", str(e))

    # ---------- auto-refresh on window focus ----------

    def event(self, qevent):
        """Reload from disk when the window regains focus and a tracked JSON
        file has been modified externally (e.g. a refit happened in a sibling
        viewer window). UI state — inclusion list, selected residue, quantity
        radio — is preserved across the reload.
        """
        if (qevent.type() == QEvent.WindowActivate
                and self.json_folder
                and self._json_files_changed_since_snapshot()):
            self._refresh_preserving_state()
        return super().event(qevent)

    def _snapshot_mtimes(self):
        """Record the modification time of every loaded JSON file."""
        self.json_mtimes = {}
        for path in self.json_paths.values():
            try:
                self.json_mtimes[str(path)] = path.stat().st_mtime_ns
            except OSError:
                pass

    def _json_files_changed_since_snapshot(self) -> bool:
        if not self.json_mtimes:
            return False
        for path in self.json_paths.values():
            try:
                if path.stat().st_mtime_ns != self.json_mtimes.get(str(path)):
                    return True
            except OSError:
                continue
        return False

    def _refresh_preserving_state(self):
        saved_included = (set(self.included_residues)
                          if self.summary_only else None)
        saved_current = self.get_current_residue()

        self._load_json_folder(self.json_folder, silent=True)

        if self.summary_only and saved_included is not None:
            self.residue_list.blockSignals(True)
            new_included = set()
            for i in range(self.residue_list.count()):
                it = self.residue_list.item(i)
                r = it.data(Qt.UserRole)
                if r in saved_included:
                    it.setCheckState(Qt.Checked)
                    new_included.add(r)
                else:
                    it.setCheckState(Qt.Unchecked)
            self.residue_list.blockSignals(False)
            self.included_residues = new_included

        if saved_current:
            for i in range(self.residue_list.count()):
                if self.residue_list.item(i).data(Qt.UserRole) == saved_current:
                    self.residue_list.setCurrentRow(i)
                    break

        self._update_plot()

    def _residue_sort_key(self, name):
        s = str(name)
        m = re.match(r"^\D*(\d+)", s)
        return (int(m.group(1)), s) if m else (10**9, s)

    def _populate_residue_list(self, residues):
        # Block itemChanged signals during repopulation to avoid spurious renders.
        self.residue_list.blockSignals(True)
        self.residue_list.clear()
        self.included_residues = set(str(r) for r in residues)
        for r in residues:
            item = QListWidgetItem(str(r))
            item.setData(Qt.UserRole, str(r))
            if self.summary_only:
                item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
                item.setCheckState(Qt.Checked)
            self.residue_list.addItem(item)
        self.residue_list.blockSignals(False)
        if residues:
            self.residue_list.setCurrentRow(0)
            if self.summary_only:
                self._update_plot()

    # ---------- summary-only inclusion controls ----------

    def _on_inclusion_toggled(self, item):
        residue = item.data(Qt.UserRole)
        if not residue:
            return
        if item.checkState() == Qt.Checked:
            self.included_residues.add(residue)
        else:
            self.included_residues.discard(residue)
        self._update_plot()

    def _set_all_check_states(self, state, predicate=None):
        self.residue_list.blockSignals(True)
        for i in range(self.residue_list.count()):
            it = self.residue_list.item(i)
            residue = it.data(Qt.UserRole)
            if predicate is not None and not predicate(residue):
                continue
            it.setCheckState(state)
            if state == Qt.Checked:
                self.included_residues.add(residue)
            else:
                self.included_residues.discard(residue)
        self.residue_list.blockSignals(False)
        self._update_plot()

    def _include_all(self):
        self._set_all_check_states(Qt.Checked)

    def _exclude_all(self):
        self._set_all_check_states(Qt.Unchecked)

    def _exclude_flagged(self):
        def is_flagged(residue):
            _, fit = self._find_fit_for(residue)
            return bool(fit and fit.get("bi_exp_unidentifiable", False))
        self._set_all_check_states(Qt.Unchecked, predicate=is_flagged)

    def _on_quantity_changed(self, checked):
        if not checked:
            return  # only react to the newly-active radio
        for key, btn in self._quantity_radios.items():
            if btn.isChecked():
                self.summary_quantity = key
                break
        self._update_plot()

    # ---------- Selection / plot dispatch ----------

    def _on_residue_clicked(self, item):
        residue = item.data(Qt.UserRole)
        if residue:
            self._update_plot(residue)

    def _on_residue_row_changed(self, row):
        if row < 0:
            return
        item = self.residue_list.item(row)
        if item:
            self._update_plot(item.data(Qt.UserRole))

    def get_current_residue(self) -> str:
        item = self.residue_list.currentItem()
        return item.data(Qt.UserRole) if item else ""

    def _update_plot(self, residue=None):
        if self.summary_only:
            self._update_plot_summary_only()
            return

        if residue is None:
            residue = self.get_current_residue()
        if not residue:
            self._show_blank_state()
            return

        self.figure.clear()
        self._pick_registry = {}

        # Summary panels stack below the main fit plot.
        summary_specs = []
        if self.show_t2_avg:
            summary_specs.append(("T2_avg", "T2_avg per residue", "T2_avg_err"))
        if self.show_dT2:
            summary_specs.append(("dT2", "ΔT2 = T2a − T2b per residue", "dT2_err"))
        if self.show_t2a:
            summary_specs.append(("t2_a", "T2a (slow) per residue", "t2_a_err"))
        if self.show_t2b:
            summary_specs.append(("t2_b", "T2b (fast) per residue", "t2_b_err"))
        if self.show_eta:
            summary_specs.append(("eta", "η = T2a / T2b per residue", "eta_err"))

        n_panels = 1 + len(summary_specs)
        ax_fit = self.figure.add_subplot(n_panels, 1, 1)
        self._plot_residue_fit(ax_fit, residue)

        for i, (key, title, err_key) in enumerate(summary_specs, start=2):
            ax = self.figure.add_subplot(n_panels, 1, i)
            self._plot_summary_panel(ax, key, err_key, title,
                                     highlight_residue=residue)

        self.figure.tight_layout()
        self.canvas.draw()

    def _update_plot_summary_only(self):
        self.figure.clear()
        self._pick_registry = {}

        if not self.data:
            self._show_blank_state()
            return

        # Find the spec for the currently-chosen quantity.
        spec = next(
            (s for s in self.SUMMARY_QUANTITIES if s[0] == self.summary_quantity),
            self.SUMMARY_QUANTITIES[0],
        )
        key, label, err_key = spec
        ax = self.figure.add_subplot(1, 1, 1)
        self._plot_summary_panel(
            ax, key, err_key, f"{label} per residue",
            included=self.included_residues,
        )
        self.figure.tight_layout()
        self.canvas.draw()

    def _show_blank_state(self):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.text(0.5, 0.5, "Load a folder of *_methylT2_fit_data.json to begin",
                ha="center", va="center", fontsize=14, color=SECONDARY_TEXT,
                transform=ax.transAxes)
        ax.axis("off")
        self.canvas.draw()

    # ---------- Plotting ----------

    def _plot_residue_fit(self, ax, residue):
        ax.set_facecolor(FRAME_BG_COLOR)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        data_key, fit = self._find_fit_for(residue)
        if fit is None:
            ax.text(0.5, 0.5, f"No fit for residue {residue}",
                    ha="center", va="center", fontsize=12, color=SECONDARY_TEXT,
                    transform=ax.transAxes)
            ax.axis("off")
            return

        metadata = self.data[data_key]["metadata"]
        time_points = metadata["time_points"]
        intensities = fit["intensities"]
        fit_time = fit["fit_curve"]["time"]
        fit_intensity = fit["fit_curve"]["intensity"]

        # Data points (with included/excluded styling and pickability).
        self._plot_residue_points(ax, data_key, residue,
                                  color="#1976D2", label="Data",
                                  time_points=time_points, intensities=intensities)
        # Bi-exp curve.
        ax.plot(fit_time, fit_intensity, "-", color="#e63946", linewidth=2,
                label="Bi-exp fit", zorder=2)

        # Optional per-component breakdown (each component carries 0.5*A in the model).
        if all(k in fit for k in ("A", "t2_a", "t2_b")):
            t_dense = np.asarray(fit_time)
            half_A = 0.5 * fit["A"]
            comp_a = half_A * np.exp(-t_dense / fit["t2_a"])
            comp_b = half_A * np.exp(-t_dense / fit["t2_b"])
            ax.plot(t_dense, comp_a, "--", color="#888888", linewidth=1, alpha=0.7,
                    label="slow component (½A)", zorder=1)
            ax.plot(t_dense, comp_b, ":", color="#888888", linewidth=1, alpha=0.7,
                    label="fast component (½A)", zorder=1)

        units = metadata.get("time_units", "ms")
        T2_avg, dT2, T2_avg_err, dT2_err, unident = self._reparam_view(fit)

        if unident:
            # When ΔT2 is at its 0 lower bound, lmfit's stderr is meaningless
            # (regularized infinity). Skip the ± part and label it as bound-active.
            if not np.isfinite(dT2_err):
                dT2_line = f"ΔT2    = {dT2:.2f} {units}   (at bound, not constrained)"
            else:
                dT2_line = f"ΔT2    = {dT2:.2f} ± {dT2_err:.2f} {units}"
            textstr = (
                f"A      = {fit.get('A', float('nan')):.3e} ± {fit.get('A_err', float('nan')):.2e}\n"
                f"T2_avg = {T2_avg:.2f} ± {T2_avg_err:.2f} {units}\n"
                f"{dT2_line}\n"
                f"(mono-exp regime — bi-exp character\nunsupported by data)"
            )
            box_color = "#fff3cd"  # warm yellow for the flagged regime
        else:
            textstr = (
                f"A      = {fit.get('A', float('nan')):.3e} ± {fit.get('A_err', float('nan')):.2e}\n"
                f"T2_avg = {T2_avg:.2f} ± {T2_avg_err:.2f} {units}\n"
                f"ΔT2    = {dT2:.2f} ± {dT2_err:.2f} {units}\n"
                f"T2a    = {fit['t2_a']:.2f} ± {fit['t2_a_err']:.2f} {units}\n"
                f"T2b    = {fit['t2_b']:.2f} ± {fit['t2_b_err']:.2f} {units}\n"
                f"η      = {fit['eta']:.2f} ± {fit['eta_err']:.2f}"
            )
            box_color = "white"
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
                fontsize=10, va="top",
                bbox=dict(boxstyle="round", facecolor=box_color, alpha=0.85))

        ax.set_xlabel(f"Time ({units})", fontsize=11, fontweight="bold")
        ax.set_ylabel("Signal Intensity", fontsize=11, fontweight="bold")
        ax.set_title(f"Residue {residue} — methyl T2 (bi-exp)", fontsize=12, fontweight="bold")
        ax.legend(loc="best", fontsize=9)
        ax.grid(True, alpha=0.3, linestyle="--", zorder=1)

    def _plot_summary_panel(self, ax, value_key, err_key, title,
                            highlight_residue=None, included=None):
        """Bar chart of one fit parameter (e.g. 'T2_avg', 't2_a', 'eta') across
        all loaded residues. The currently-selected residue is highlighted in red.
        Residues flagged bi_exp_unidentifiable are drawn in gray to signal that
        their T2a/T2b/η are not statistically distinguishable from a mono-exp.

        If `included` is provided (set of residue ids), only those residues are
        plotted; this is how the summary-only mode supports user exclusion.
        Bars whose value is non-finite are skipped silently; non-finite errors
        are coerced to 0 so the bar still appears (the value carries info even
        without an error bar).
        """
        ax.set_facecolor(FRAME_BG_COLOR)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        residues, values, errors, flagged = [], [], [], []
        for residue in self.available_residues:
            if included is not None and str(residue) not in included:
                continue
            _, fit = self._find_fit_for(residue)
            if fit is None:
                continue
            v, e = self._panel_value_and_err(fit, value_key, err_key)
            if not np.isfinite(v):
                continue
            if not np.isfinite(e):
                e = 0.0
            residues.append(str(residue))
            values.append(float(v))
            errors.append(float(e))
            flagged.append(bool(fit.get("bi_exp_unidentifiable", False)))

        if not residues:
            ax.text(0.5, 0.5, f"No {value_key} data", ha="center", va="center",
                    transform=ax.transAxes, color=SECONDARY_TEXT)
            ax.axis("off")
            return

        x = np.arange(len(residues))
        # Gray for flagged residues; blue for clean bi-exp; red for selection.
        colors = ["#9e9e9e" if f else "#1976D2" for f in flagged]
        if highlight_residue is not None and str(highlight_residue) in residues:
            colors[residues.index(str(highlight_residue))] = "#e63946"

        ax.bar(x, values, yerr=errors, color=colors, capsize=2)
        ax.set_xticks(x)
        ax.set_xticklabels(residues, rotation=70, fontsize=8)
        if value_key == "eta":
            ax.axhline(1.0, color="grey", linewidth=0.8, linestyle=":")
            ax.set_ylabel("η = T2a / T2b", fontsize=11, fontweight="bold")
        elif value_key == "dT2":
            ax.axhline(0.0, color="grey", linewidth=0.8, linestyle=":")
            ax.set_ylabel("ΔT2 (ms)", fontsize=11, fontweight="bold")
        elif value_key == "T2_avg":
            ax.set_ylabel("T2_avg (ms)", fontsize=11, fontweight="bold")
        else:
            ax.set_ylabel(f"{value_key} (ms)", fontsize=11, fontweight="bold")
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.grid(True, axis="y", alpha=0.3, linestyle="--")

    @staticmethod
    def _panel_value_and_err(fit, value_key, err_key):
        """Resolve (value, err) for a panel, deriving reparam keys from the
        legacy (t2_a, t2_b) entries when older JSON files lack T2_avg/dT2.
        """
        if value_key in ("T2_avg", "dT2"):
            T2_avg, dT2, T2_avg_err, dT2_err, _ = MethylT2FitViewer._reparam_view(fit)
            if value_key == "T2_avg":
                return T2_avg, T2_avg_err
            return dT2, dT2_err
        v = fit.get(value_key, np.nan)
        e = fit.get(err_key, 0.0) or 0.0
        return v, e

    # ---------- Outlier editing ----------

    def _plot_residue_points(self, ax, data_key, residue, color, label,
                             time_points, intensities):
        excluded = set(self.exclusions.get(data_key, {}).get(str(residue), []))
        n = len(time_points)
        inc_idx = [i for i in range(n) if i not in excluded]
        exc_idx = [i for i in range(n) if i in excluded]
        tp = np.asarray(time_points)
        ints = np.asarray(intensities)

        picker_radius = 8 if self.edit_outliers else None

        if inc_idx:
            inc = ax.scatter(tp[inc_idx], ints[inc_idx],
                             color=color, s=64, marker="o",
                             label=label, zorder=3, picker=picker_radius)
            self._pick_registry[id(inc)] = (data_key, str(residue), inc_idx)
        if exc_idx:
            exc = ax.scatter(tp[exc_idx], ints[exc_idx],
                             color="#888888", s=64, marker="x", linewidths=2,
                             label="Excluded", zorder=3, picker=picker_radius)
            self._pick_registry[id(exc)] = (data_key, str(residue), exc_idx)

    def _on_pick(self, event):
        if not self.edit_outliers:
            return
        info = self._pick_registry.get(id(event.artist))
        if info is None:
            return
        data_key, residue, idx_map = info
        excl_for_key = self.exclusions.setdefault(data_key, {})
        excl = excl_for_key.setdefault(residue, [])
        for i in event.ind:
            if 0 <= i < len(idx_map):
                orig = idx_map[i]
                if orig in excl:
                    excl.remove(orig)
                else:
                    excl.append(orig)
        if not excl:
            del excl_for_key[residue]
        json_path = self.json_paths.get(data_key)
        if json_path is not None:
            save_outliers(json_path, excl_for_key)
        self._update_plot()

    def _toggle_edit_outliers(self, checked):
        self.edit_outliers = bool(checked)
        self._update_plot()

    def _toggle_summary(self, _checked=None):
        """Read all summary checkboxes and re-render."""
        self.show_t2_avg = self.t2_avg_checkbox.isChecked()
        self.show_dT2 = self.dT2_checkbox.isChecked()
        self.show_t2a = self.t2a_checkbox.isChecked()
        self.show_t2b = self.t2b_checkbox.isChecked()
        self.show_eta = self.eta_checkbox.isChecked()
        self._update_plot()

    @staticmethod
    def _reparam_view(fit):
        """Return (T2_avg, dT2, T2_avg_err, dT2_err, bi_exp_unidentifiable).

        Handles JSON entries written by older versions of the fitter (which
        lacked the reparam keys) by deriving T2_avg, dT2 from t2_a/t2_b. In
        that case errors fall back to NaN and the unidentifiable flag to False.
        """
        t2_a = fit.get("t2_a", float("nan"))
        t2_b = fit.get("t2_b", float("nan"))
        T2_avg = fit.get("T2_avg")
        if T2_avg is None or not np.isfinite(T2_avg):
            T2_avg = (t2_a + t2_b) / 2.0
        dT2 = fit.get("dT2")
        if dT2 is None or not np.isfinite(dT2):
            dT2 = abs(t2_a - t2_b)
        T2_avg_err = fit.get("T2_avg_err", float("nan")) or float("nan")
        dT2_err = fit.get("dT2_err", float("nan")) or float("nan")
        unident = bool(fit.get("bi_exp_unidentifiable", False))
        return T2_avg, dT2, T2_avg_err, dT2_err, unident

    def _find_fit_for(self, residue):
        """Return (data_key, fit_entry) for the first dataset containing residue.

        Iterates over data_keys in sorted order so the result is deterministic
        when multiple JSON files (different fields) are loaded; otherwise the
        choice depends on Path.glob's filesystem order.
        """
        for data_key in sorted(self.data.keys()):
            payload = self.data[data_key]
            for fit in payload.get("fits", []):
                if str(fit.get("residue")) == str(residue):
                    return data_key, fit
        return None, None

    def _find_tsv_for(self, data_key):
        """Best-effort: locate the *_methylT2_fit_results.txt sibling of the json folder."""
        json_path = self.json_paths.get(data_key)
        if json_path is None:
            return None
        candidate = Path(json_path).parent.parent / f"{data_key}_fit_results.txt"
        return candidate if candidate.exists() else None

    def _refit_one(self, data_key, residue):
        excl = self.exclusions.get(data_key, {}).get(str(residue), [])
        payload = self.data.get(data_key)
        if not payload:
            return None
        metadata = payload["metadata"]
        for fit in payload["fits"]:
            if str(fit["residue"]) == str(residue):
                new_entry = refit_residue_methyl(fit, metadata, excl)
                fit.update(new_entry)
                json_path = self.json_paths.get(data_key)
                if json_path is not None:
                    update_json_fit_entry(json_path, str(residue), new_entry)
                tsv_path = self._find_tsv_for(data_key)
                if tsv_path is not None:
                    try:
                        update_tsv_row_methyl(tsv_path, str(residue), new_entry)
                    except KeyError:
                        pass  # row not present in TSV; skip
                # Refresh our own mtime snapshot so the focus-event handler
                # doesn't see this write as an external change and reload.
                self._snapshot_mtimes()
                return new_entry
        return None

    def _on_refit_clicked(self):
        residue = self.get_current_residue()
        if not residue:
            show_warning(self, "No residue selected", "Select a residue to re-fit.")
            return
        from PySide6.QtWidgets import QApplication
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            refitted = []
            for data_key in list(self.data.keys()):
                if self._refit_one(data_key, residue) is not None:
                    refitted.append(data_key)
            if refitted:
                self._update_plot()
        finally:
            QApplication.restoreOverrideCursor()
        if refitted:
            show_info(self, "Re-fit complete",
                      f"Residue {residue} re-fitted in: {', '.join(refitted)}")

    def _on_reset_clicked(self):
        residue = self.get_current_residue()
        if not residue:
            show_warning(self, "No residue selected", "Select a residue to reset.")
            return
        from PySide6.QtWidgets import QApplication
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            reset = []
            for data_key in list(self.data.keys()):
                excl_map = self.exclusions.get(data_key, {})
                if str(residue) in excl_map:
                    del excl_map[str(residue)]
                    json_path = self.json_paths.get(data_key)
                    if json_path is not None:
                        save_outliers(json_path, excl_map)
                if self._refit_one(data_key, residue) is not None:
                    reset.append(data_key)
            if reset:
                self._update_plot()
        finally:
            QApplication.restoreOverrideCursor()
        if reset:
            show_info(self, "Reset complete",
                      f"Exclusions cleared and residue {residue} re-fitted.")


def open_methyl_t2_fit_viewer(parent=None, json_folder=None):
    viewer = MethylT2FitViewer(parent, json_folder=json_folder)
    viewer.show()
    return viewer


def open_methyl_t2_summary_viewer(parent=None, json_folder=None):
    """Open the across-residues summary plot (Plot T2 vs Residue): a single
    bar chart of T2_avg / ΔT2 / T2a / T2b / η, with a residue checklist for
    excluding contaminating fits.
    """
    viewer = MethylT2FitViewer(parent, json_folder=json_folder, summary_only=True)
    viewer.show()
    return viewer
