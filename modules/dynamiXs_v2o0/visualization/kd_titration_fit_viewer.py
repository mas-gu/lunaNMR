# ABOUTME: Viewer for Kd titration fits - per-residue binding curves + summary bars.
# ABOUTME: Reads a *_kd_fit_data.json produced by dynamiXs_Kd.kd_fit.

import json
import os
import sys
from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QListWidget, QListWidgetItem,
    QComboBox, QLabel, QPushButton,
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

# Make the Kd model functions importable (bare-name imports).
_KD_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "dynamiXs_Kd")
if _KD_DIR not in sys.path:
    sys.path.insert(0, _KD_DIR)
from kd_models import csp_model, intensity_model  # noqa: E402


class KdTitrationFitViewer(QMainWindow):
    """Inspect per-residue Kd binding curves and across-residue summaries."""

    def __init__(self, parent=None, json_file=None):
        super().__init__(parent)
        self.setWindowTitle("Kd / Titration Fit Viewer")
        self.resize(1100, 720)
        self.data = {}
        self.fits = []
        self.P0 = 50.0

        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)

        # Plot
        left = QVBoxLayout()
        self.figure = Figure(figsize=(8, 6))
        self.canvas = FigureCanvasQTAgg(self.figure)
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
        self.view_combo.addItems(["Per-residue curve", "Kd vs residue", "Δδmax / amp vs residue"])
        self.view_combo.currentIndexChanged.connect(self._refresh)
        right.addWidget(self.view_combo)

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
        self.data = json.loads(Path(json_file).read_text())
        self.fits = self.data.get("fits", [])
        self.P0 = float(self.data.get("metadata", {}).get("protein_conc", 50.0))
        g = self.data.get("global", {}).get("csp", {})
        self.global_label.setText(
            f"Global shared Kd (CSP): {g['Kd']:.4g}  (n={g.get('n_residues','?')})"
            if g.get("success") else "Global Kd: n/a")
        self.residue_list.clear()
        for f in self.fits:
            item = QListWidgetItem(f["residue"])
            self.residue_list.addItem(item)
        if self.fits:
            self.residue_list.setCurrentRow(0)
        self._refresh()

    def _obs_key(self):
        return "csp" if self.obs_combo.currentIndex() == 0 else "intensity"

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
            key = "dd_max" if obs == "csp" else "amp"
            self._plot_summary(ax, key, key + "_err", f"{key} vs residue", key)
        self.canvas.draw()

    def _plot_curve(self, ax):
        row = self.residue_list.currentRow()
        if row < 0 or row >= len(self.fits):
            return
        f = self.fits[row]
        obs = self._obs_key()
        fit = f.get(obs, {})
        if not fit.get("success"):
            ax.text(0.5, 0.5, f"No {obs} fit for {f['residue']}", ha="center", va="center",
                    transform=ax.transAxes)
            return
        L = np.asarray(fit["L"], float)
        y = np.asarray(fit["obs"], float)
        ax.plot(L, y, "ko", ms=7, label="data")
        Lg = np.linspace(0, L.max() * 1.05, 200)
        if obs == "csp":
            yg = csp_model(Lg, fit["dd_max"], fit["Kd"], self.P0)
            ylab = "CSP (ppm)"
        else:
            yg = intensity_model(Lg, fit["baseline"], fit["amp"], fit["Kd"], self.P0)
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

    def _export(self):
        from PySide6.QtWidgets import QFileDialog
        path, _ = QFileDialog.getSaveFileName(self, "Export figure", "", "PDF (*.pdf);;PNG (*.png)")
        if path:
            self.figure.savefig(path, bbox_inches="tight", dpi=300)


def open_kd_titration_viewer(parent=None, json_file=None):
    viewer = KdTitrationFitViewer(parent=parent, json_file=json_file)
    viewer.show()
    return viewer
