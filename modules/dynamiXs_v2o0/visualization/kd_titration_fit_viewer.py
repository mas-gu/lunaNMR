# ABOUTME: Viewer for Kd titration fits - per-residue binding curves + summary bars.
# ABOUTME: Reads a *_kd_fit_data.json; click-to-exclude points and refit (dynamiXs pattern).

import json
import os
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
from kd_models import csp_model, intensity_decay  # noqa: E402
from kd_fit import fit_residue_csp, fit_residue_intensity  # noqa: E402


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
        self.view_combo.addItems(["Per-residue curve", "Kd vs residue", "amplitude vs residue"])
        self.view_combo.currentIndexChanged.connect(self._refresh)
        right.addWidget(self.view_combo)

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
        self.fits = self.data.get("fits", [])
        self.P0 = float(self.data.get("metadata", {}).get("protein_conc", 50.0))
        # restore any saved exclusions
        self.excluded = {"csp": {}, "intensity": {}}
        for f in self.fits:
            for obs in ("csp", "intensity"):
                excl = f.get(obs, {}).get("excluded")
                if excl:
                    self.excluded[obs][f["residue"]] = set(excl)
        self._update_global_label()
        self.residue_list.clear()
        for f in self.fits:
            self.residue_list.addItem(QListWidgetItem(f["residue"]))
        if self.fits:
            self.residue_list.setCurrentRow(0)
        self._refresh()

    def _update_global_label(self):
        g = self.data.get("global", {}).get("csp", {})
        self.global_label.setText(
            f"Global shared Kd (CSP): {g['Kd']:.4g}  (n={g.get('n_residues','?')})"
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
        if view == 0:
            self._plot_curve(ax)
        elif view == 1:
            self._plot_summary(ax, "Kd", "Kd_err", "Kd vs residue", "Kd")
        else:
            obs = self._obs_key()
            key = "dd_max" if obs == "csp" else "I0"
            self._plot_summary(ax, key, key + "_err", f"{key} vs residue", key)
        self.canvas.draw()

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
            if fit.get("success") and key in fit:
                names.append(f["residue"])
                vals.append(fit[key])
                errs.append(fit.get(err_key) if isinstance(fit.get(err_key), (int, float)) else 0.0)
        if not names:
            ax.text(0.5, 0.5, f"No {obs} fits", ha="center", va="center", transform=ax.transAxes)
            return
        x = np.arange(len(names))
        ax.bar(x, vals, yerr=errs, color="steelblue", capsize=2)
        if key == "Kd":
            g = self.data.get("global", {}).get("csp", {})
            if obs == "csp" and g.get("success"):
                ax.axhline(g["Kd"], color="red", ls="--", label=f"global Kd={g['Kd']:.3g}")
                ax.legend()
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=90, fontsize=6)
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
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            self._refit_residue(f, obs)
            self._save_json()
            self._refresh()
        finally:
            QApplication.restoreOverrideCursor()

    def _save_json(self):
        if self.json_file:
            Path(self.json_file).write_text(json.dumps(self.data, indent=2))

    def _export(self):
        path, _ = QFileDialog.getSaveFileName(self, "Export figure", "", "PDF (*.pdf);;PNG (*.png)")
        if path:
            self.figure.savefig(path, bbox_inches="tight", dpi=300)


def open_kd_titration_viewer(parent=None, json_file=None):
    viewer = KdTitrationFitViewer(parent=parent, json_file=json_file)
    viewer.show()
    return viewer
