# ABOUTME: Viewer for Kd titration fits - per-residue binding curves + summary bars.
# ABOUTME: Reads a *_kd_fit_data.json; click-to-exclude points and refit (dynamiXs pattern).

import json
import os
import re
import sys
from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QListWidget, QListWidgetItem,
    QComboBox, QLabel, QPushButton, QFileDialog, QApplication,
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

# Make the Kd model functions importable (bare-name imports).
_KD_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "dynamiXs_Kd")
if _KD_DIR not in sys.path:
    sys.path.insert(0, _KD_DIR)
from kd_models import csp_model, intensity_decay, compute_csp  # noqa: E402
from kd_fit import fit_residue_csp, fit_residue_intensity, json_safe  # noqa: E402


def _residue_sort_key(name):
    """Order residues by sequence number, not amino-acid letter: 'K14' before
    'A17'. Names without a number sort last, alphabetically among themselves."""
    m = re.search(r'\d+', str(name))
    return (int(m.group()) if m else float('inf'), str(name))


def _peak_present(series, k, value='height'):
    """True if a residue's peak is genuinely present at point k: a detected position
    (finite, non-zero sentinel) AND a real intensity (finite, > 0). Used so CSP and
    intensity treat the same residues as missing (grey)."""
    px = np.asarray(series.get('ppm_x', []), float)
    py = np.asarray(series.get('ppm_y', []), float)
    v = np.asarray(series.get(value, []), float)
    if k < 0 or k >= min(len(px), len(py), len(v)):
        return False
    return (np.isfinite(px[k]) and px[k] != 0.0
            and np.isfinite(py[k]) and py[k] != 0.0
            and np.isfinite(v[k]) and v[k] > 0.0)


def _pair_observable(series, i, j, obs, alpha=0.14, value='height'):
    """Observable for one residue between reference point i and point j.

    obs='csp'  -> CSP (ppm) from the position shift: sqrt(ΔδH² + (alpha·ΔδN)²).
    obs='intensity' -> ratio v[j]/v[i] of the chosen value ('height'/'volume').
    Returns NaN if an endpoint is missing — a 0.0 position is the undetected
    sentinel, a non-positive reference intensity makes the ratio undefined
    (same conventions as kd_input.csp_series / intensity_ratio_series).
    """
    if obs == 'csp':
        px = np.asarray(series.get('ppm_x', []), float)
        py = np.asarray(series.get('ppm_y', []), float)
        if max(i, j) >= len(px) or max(i, j) >= len(py):
            return float('nan')
        xi, xj, yi, yj = px[i], px[j], py[i], py[j]
        if 0.0 in (xi, xj, yi, yj) or not np.all(np.isfinite([xi, xj, yi, yj])):
            return float('nan')
        return float(compute_csp(xj - xi, yj - yi, alpha=alpha))
    v = np.asarray(series.get(value, []), float)
    if max(i, j) >= len(v):
        return float('nan')
    ref = v[i]
    if not np.isfinite(ref) or ref <= 0.0:
        return float('nan')
    return float(v[j] / ref)


class KdTitrationFitViewer(QMainWindow):
    """Inspect per-residue Kd binding curves, exclude points and refit on the go."""

    def __init__(self, parent=None, json_file=None):
        super().__init__(parent)
        self.setWindowTitle("Kd / Titration Fit Viewer")
        self.resize(1100, 720)
        self.data = {}
        self.fits = []
        self.P0 = 50.0
        self.json_file = None
        self.edit_mode = False
        # excluded[obs][residue] = set of point indices (into that fit's L/obs)
        self.excluded = {"csp": {}, "intensity": {}}
        # residues excluded wholesale (problematic values) — click a bar in edit mode.
        # Per-observable: excluding a residue in CSP must not exclude it in Intensity.
        self.excluded_residues = {"csp": set(), "intensity": set()}
        self._pick_registry = {}

        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)

        # Plot
        left = QVBoxLayout()
        self.figure = Figure(figsize=(8, 6))
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.mpl_connect("pick_event", self._on_pick)
        left.addWidget(NavigationToolbar2QT(self.canvas, self))
        left.addWidget(self.canvas, stretch=1)
        root.addLayout(left, stretch=4)

        # Controls
        right = QVBoxLayout()
        right.addWidget(QLabel("Observable"))
        self.obs_combo = QComboBox()
        self.obs_combo.addItems(["CSP", "Intensity"])
        self.obs_combo.currentIndexChanged.connect(self._refresh)
        right.addWidget(self.obs_combo)

        right.addWidget(QLabel("View"))
        self.view_combo = QComboBox()
        self.view_combo.addItems(["Per-residue curve", "Kd vs residue",
                                  "Δδmax (CSP) / I0 (int) vs residue",
                                  "Ref vs point (bars)"])
        self.view_combo.currentIndexChanged.connect(self._refresh)
        right.addWidget(self.view_combo)

        # Ref/compare point selectors — only shown for the "Ref vs point" view.
        self.cmp_box = QWidget()
        cmp_layout = QVBoxLayout(self.cmp_box)
        cmp_layout.setContentsMargins(0, 0, 0, 0)
        cmp_layout.addWidget(QLabel("Reference point"))
        self.ref_point_combo = QComboBox()
        self.ref_point_combo.currentIndexChanged.connect(self._refresh)
        cmp_layout.addWidget(self.ref_point_combo)
        cmp_layout.addWidget(QLabel("Compare point"))
        self.cmp_point_combo = QComboBox()
        self.cmp_point_combo.currentIndexChanged.connect(self._refresh)
        cmp_layout.addWidget(self.cmp_point_combo)
        right.addWidget(self.cmp_box)
        self.cmp_box.setVisible(False)

        # Edit / refit controls (dynamiXs-style outlier rejection)
        self.edit_btn = QPushButton("Exclude points (click)")
        self.edit_btn.setCheckable(True)
        self.edit_btn.toggled.connect(self._toggle_edit)
        right.addWidget(self.edit_btn)
        self.refit_btn = QPushButton("Refit (excl. points)")
        self.refit_btn.clicked.connect(self._on_refit)
        right.addWidget(self.refit_btn)
        self.reset_btn = QPushButton("Reset points")
        self.reset_btn.clicked.connect(self._on_reset)
        right.addWidget(self.reset_btn)

        self.global_label = QLabel("")
        self.global_label.setWordWrap(True)
        right.addWidget(self.global_label)

        right.addWidget(QLabel("Residues"))
        self.residue_list = QListWidget()
        self.residue_list.currentRowChanged.connect(self._refresh)
        right.addWidget(self.residue_list, stretch=1)

        export_btn = QPushButton("Export figure…")
        export_btn.clicked.connect(self._export)
        right.addWidget(export_btn)
        root.addLayout(right, stretch=1)

        if json_file:
            self.load(json_file)

    # ---------- data ----------

    def load(self, json_file):
        self.json_file = json_file
        self.data = json.loads(Path(json_file).read_text())
        self.fits = sorted(self.data.get("fits", []),
                           key=lambda f: _residue_sort_key(f.get("residue", "")))
        self.P0 = float(self.data.get("metadata", {}).get("protein_conc", 50.0))
        # restore any saved exclusions
        self.excluded = {"csp": {}, "intensity": {}}
        for f in self.fits:
            for obs in ("csp", "intensity"):
                excl = f.get(obs, {}).get("excluded")
                if excl:
                    self.excluded[obs][f["residue"]] = set(excl)
        er = self.data.get("excluded_residues", {})
        if isinstance(er, dict):
            self.excluded_residues = {"csp": set(er.get("csp", [])),
                                      "intensity": set(er.get("intensity", []))}
        else:  # legacy flat list → applied to both observables
            self.excluded_residues = {"csp": set(er or []), "intensity": set(er or [])}
        self._update_global_label()
        self._populate_point_combos()
        self.residue_list.clear()
        for f in self.fits:
            self.residue_list.addItem(QListWidgetItem(f["residue"]))
        if self.fits:
            self.residue_list.setCurrentRow(0)
        self._refresh()

    def _point_labels(self):
        """Per-point labels for the comparison selectors — ligand concentrations
        when available, else the raw point labels."""
        meta = self.data.get("metadata", {})
        pts = meta.get("concentrations") or meta.get("points") or []
        return [f"{p:g}" if isinstance(p, (int, float)) else str(p) for p in pts]

    def _populate_point_combos(self):
        """Fill Ref/Compare selectors; default ref = point 0, compare = first
        addition (the next point, lowest non-zero [L] since points are sorted)."""
        labels = self._point_labels()
        for combo in (self.ref_point_combo, self.cmp_point_combo):
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(labels)
        if labels:
            self.ref_point_combo.setCurrentIndex(0)
            self.cmp_point_combo.setCurrentIndex(1 if len(labels) > 1 else 0)
        for combo in (self.ref_point_combo, self.cmp_point_combo):
            combo.blockSignals(False)

    def _update_global_label(self):
        g = self.data.get("global", {}).get("csp", {})
        self.global_label.setText(
            f"Global shared Kd (CSP): {g['Kd']:.4g}  (n={g.get('n_residues','?')}, "
            f"from initial fit — not updated by per-residue refits)"
            if g.get("success") else "Global Kd: n/a")

    def _obs_key(self):
        return "csp" if self.obs_combo.currentIndex() == 0 else "intensity"

    def _current_residue(self):
        row = self.residue_list.currentRow()
        return self.fits[row] if 0 <= row < len(self.fits) else None

    # ---------- plotting ----------

    def _refresh(self, *_):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        view = self.view_combo.currentIndex()
        self.cmp_box.setVisible(view == 3)
        if view == 0:
            self._plot_curve(ax)
        elif view == 1:
            self._plot_summary(ax, "Kd", "Kd_err", "Kd vs residue", "Kd")
        elif view == 3:
            self._plot_comparison(ax)
        else:
            obs = self._obs_key()
            key = "dd_max" if obs == "csp" else "I0"
            self._plot_summary(ax, key, key + "_err", f"{key} vs residue", key)
        self.canvas.draw()

    def _plot_comparison(self, ax):
        """Observable (CSP or I/I₀) per residue between a chosen reference point
        and a comparison point — computed from each residue's raw series."""
        obs = self._obs_key()
        i = self.ref_point_combo.currentIndex()
        j = self.cmp_point_combo.currentIndex()
        labels = self._point_labels()
        if i < 0 or j < 0 or not labels:
            ax.text(0.5, 0.5, "No titration points", ha="center", va="center",
                    transform=ax.transAxes)
            return
        if not any(f.get("series") for f in self.fits):
            ax.text(0.5, 0.5, "no embedded series — re-run the fit to enable this view",
                    ha="center", va="center", transform=ax.transAxes)
            return
        meta = self.data.get("metadata", {})
        alpha = float(meta.get("alpha", 0.14))
        value = meta.get("intensity_value", "height")
        names, vals = [], []
        for f in self.fits:
            series = f.get("series")
            # "Missing" is the SAME for both observables: the peak must be present
            # (detected position + a real intensity) at both the ref and compare point.
            # Otherwise CSP greyed on bad positions while intensity greyed on bad
            # intensity, so the same residue looked missing in one mode but not the other.
            if (series and _peak_present(series, i, value)
                    and _peak_present(series, j, value)):
                v = _pair_observable(series, i, j, obs, alpha=alpha, value=value)
            else:
                v = float("nan")
            names.append(f["residue"])
            vals.append(v)
        bar_color = "seagreen" if obs == "csp" else "indianred"
        # I/I₀ is a ratio → fix the axis at 0–1; CSP is in ppm → auto-scale.
        self._plot_residue_bars(ax, names, vals, bar_color,
                                ymax=None if obs == "csp" else 1.0)
        ax.set_ylabel("CSP (ppm)" if obs == "csp" else "Intensity ratio I/I₀")
        sym = "CSP" if obs == "csp" else "I/I₀"
        ax.set_title(f"{sym}:  {labels[i]} → {labels[j]}  (ref → point)")

    def _plot_residue_bars(self, ax, names, vals, color, errs=None, ymax=None):
        """Bar chart of value-per-residue. Residues with no value (NaN — missing
        assignment or failed fit) draw as a full-height grey bar so the gap is visible;
        residues the user excluded (problematic values) draw hatched. In edit mode the
        bars are pickable so a click toggles that residue's exclusion. `ymax` fixes the
        y-axis top (e.g. 1.0 for an I/I₀ ratio); None auto-scales to the data."""
        self._pick_registry = {}
        excl = self.excluded_residues[self._obs_key()]   # per-observable exclusions
        finite = [v for v, n in zip(vals, names)
                  if np.isfinite(v) and n not in excl]
        data_max = max(finite) if finite else (max([v for v in vals if np.isfinite(v)],
                                                    default=1.0))
        # Honor a fixed top (ratio axis), expanding only if the data exceeds it.
        if ymax is not None:
            top = ymax if data_max <= ymax else data_max * 1.05
        else:
            top = data_max * 1.1
        full = top                         # grey/missing bars span to the axis top
        x = np.arange(len(names))
        heights = [v if np.isfinite(v) else full for v in vals]
        facecolors = []
        for n, v in zip(names, vals):
            if n in excl:
                facecolors.append("#e08585")          # excluded → muted red (hatched below)
            elif np.isfinite(v):
                facecolors.append(color)
            else:
                facecolors.append("#d0d0d0")           # missing → grey full-height
        bar_errs = ([e if np.isfinite(v) and n not in excl else 0.0
                     for v, e, n in zip(vals, errs, names)] if errs is not None else None)
        bars = ax.bar(x, heights, color=facecolors, yerr=bar_errs, capsize=2)
        for n, patch in zip(names, bars):
            if n in excl:
                patch.set_hatch("//")
            if self.edit_mode:
                patch.set_picker(True)
                self._pick_registry[id(patch)] = (n, "RESIDUE", None)
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=90, fontsize=6)
        ax.set_ylim(0, top)
        return full

    def _plot_curve(self, ax):
        self._pick_registry = {}
        f = self._current_residue()
        if f is None:
            return
        obs = self._obs_key()
        fit = f.get(obs, {})
        if not fit.get("success"):
            ax.text(0.5, 0.5, f"No {obs} fit for {f['residue']}", ha="center", va="center",
                    transform=ax.transAxes)
            return
        L = np.asarray(fit["L"], float)
        y = np.asarray(fit["obs"], float)
        excl = self.excluded[obs].get(f["residue"], set())
        n = len(L)
        inc_idx = [i for i in range(n) if i not in excl]
        exc_idx = [i for i in range(n) if i in excl]
        picker = 8 if self.edit_mode else None

        if inc_idx:
            art = ax.scatter(L[inc_idx], y[inc_idx], color="black", s=64, zorder=3,
                             label="data", picker=picker)
            self._pick_registry[id(art)] = (f["residue"], obs, inc_idx)
        if exc_idx:
            art = ax.scatter(L[exc_idx], y[exc_idx], color="#888888", marker="x", s=64,
                             linewidths=2, zorder=3, label="excluded", picker=picker)
            self._pick_registry[id(art)] = (f["residue"], obs, exc_idx)

        Lg = np.linspace(0, L.max() * 1.05, 200)
        if obs == "csp":
            yg = csp_model(Lg, fit["dd_max"], fit["Kd"], self.P0)
            ylab = "CSP (ppm)"
        else:
            yg = intensity_decay(Lg, fit["I0"], fit["I_inf"], fit["Kd"])
            ylab = "Intensity ratio"
        ax.plot(Lg, yg, "r-", lw=2, label="fit")
        kd, kde = fit["Kd"], fit.get("Kd_err", float("nan"))
        ax.set_title(f"{f['residue']}   Kd = {kd:.3g} ± {kde:.2g}   R² = {fit.get('r_squared', 0):.3f}")
        ax.set_xlabel("[Ligand]")
        ax.set_ylabel(ylab)
        ax.legend()

    def _plot_summary(self, ax, key, err_key, title, ylab):
        obs = self._obs_key()
        names, vals, errs = [], [], []
        for f in self.fits:
            fit = f.get(obs, {})
            ok = fit.get("success") and isinstance(fit.get(key), (int, float))
            names.append(f["residue"])
            vals.append(float(fit[key]) if ok else float("nan"))
            e = fit.get(err_key)
            errs.append(float(e) if ok and isinstance(e, (int, float)) else 0.0)
        if not any(np.isfinite(v) for v in vals):
            ax.text(0.5, 0.5, f"No {obs} fits", ha="center", va="center", transform=ax.transAxes)
            return
        # Missing/failed residues show as full-height grey bars (see _plot_residue_bars).
        self._plot_residue_bars(ax, names, vals, "steelblue", errs=errs)
        if key == "Kd":
            g = self.data.get("global", {}).get("csp", {})
            if obs == "csp" and g.get("success"):
                ax.axhline(g["Kd"], color="red", ls="--", label=f"global Kd={g['Kd']:.3g}")
                ax.legend()
        ax.set_ylabel(ylab)
        ax.set_title(f"{title}  ({obs})")

    # ---------- exclude / refit (dynamiXs pattern) ----------

    def _toggle_edit(self, checked):
        self.edit_mode = bool(checked)
        self._refresh()

    def _on_pick(self, event):
        if not self.edit_mode:
            return
        info = self._pick_registry.get(id(event.artist))
        if info is None:
            return
        residue, obs, idx_map = info
        if obs == "RESIDUE":
            # Toggle whole-residue exclusion for the CURRENT observable only — CSP and
            # Intensity exclusions are independent.
            cur = self._obs_key()
            s = self.excluded_residues[cur]
            s.discard(residue) if residue in s else s.add(residue)
            self.data["excluded_residues"] = {
                k: sorted(v) for k, v in self.excluded_residues.items()}
            self._save_json()
            self._refresh()
            return
        excl = self.excluded[obs].setdefault(residue, set())
        for i in event.ind:
            if 0 <= i < len(idx_map):
                orig = idx_map[i]
                excl.discard(orig) if orig in excl else excl.add(orig)
        if not excl:
            self.excluded[obs].pop(residue, None)
        self._refresh()

    def _refit_residue(self, f, obs):
        """Refit one residue/observable on the included points; keep full L/obs for
        display, update only the parameter fields. Returns True on success."""
        fit = f.get(obs, {})
        if not fit.get("L"):
            return False
        L = np.asarray(fit["L"], float)
        y = np.asarray(fit["obs"], float)
        excl = self.excluded[obs].get(f["residue"], set())
        inc = [i for i in range(len(L)) if i not in excl]
        L_inc, y_inc = L[inc], y[inc]
        n_boot = int(self.data.get("metadata", {}).get("n_bootstrap", 0))
        if obs == "csp":
            new = fit_residue_csp(L_inc, y_inc, self.P0, n_bootstrap=n_boot)
        else:
            new = fit_residue_intensity(L_inc, y_inc, n_bootstrap=n_boot)
        if not new.get("success"):
            return False
        # keep the full display points + record the exclusions; update params only
        for k, v in new.items():
            if k not in ("L", "obs"):
                fit[k] = v
        fit["excluded"] = sorted(excl)
        return True

    def _on_refit(self):
        f = self._current_residue()
        if f is None:
            return
        obs = self._obs_key()
        ok = False
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            ok = self._refit_residue(f, obs)
            if ok:
                self._save_json()
                self._refresh()
        finally:
            QApplication.restoreOverrideCursor()
        if not ok:
            from PySide6.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Refit failed",
                                "Could not refit — too few points left after exclusion "
                                "(need at least 3), or the fit did not converge.")

    def _on_reset(self):
        f = self._current_residue()
        if f is None:
            return
        obs = self._obs_key()
        self.excluded[obs].pop(f["residue"], None)
        ok = False
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            ok = self._refit_residue(f, obs)
            self._save_json()
            self._refresh()
        finally:
            QApplication.restoreOverrideCursor()
        if not ok:
            from PySide6.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Reset incomplete",
                                "Exclusions were cleared but the all-points refit did "
                                "not converge — parameters left unchanged.")

    def _save_json(self):
        if self.json_file:
            Path(self.json_file).write_text(json.dumps(json_safe(self.data), indent=2))

    def _export(self):
        path, selected = QFileDialog.getSaveFileName(
            self, "Export figure", "",
            "SVG — editable vector (*.svg);;PDF — vector (*.pdf);;PNG — raster (*.png)")
        if not path:
            return
        # Ensure the filename has an extension matching the chosen filter.
        ext = os.path.splitext(path)[1].lower()
        if ext not in (".svg", ".pdf", ".png"):
            ext = ".svg" if "svg" in selected else ".pdf" if "pdf" in selected else ".png"
            path += ext
        # Keep text editable in vector output so it can be edited in Illustrator:
        # svg.fonttype='none' keeps labels as <text> (not outlined paths); pdf.fonttype=42
        # embeds TrueType so glyphs stay selectable. rc_context avoids touching globals.
        import matplotlib as mpl
        with mpl.rc_context({"svg.fonttype": "none", "pdf.fonttype": 42}):
            self.figure.savefig(path, bbox_inches="tight", dpi=300)


def open_kd_titration_viewer(parent=None, json_file=None):
    viewer = KdTitrationFitViewer(parent=parent, json_file=json_file)
    viewer.show()
    return viewer
