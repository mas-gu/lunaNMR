# ABOUTME: Viewer for Kd titration fits - per-residue binding curves + summary bars.
# ABOUTME: Reads a *_kd_fit_data.json; click-to-exclude points and refit (dynamiXs pattern).

import json
import os
import sys
from contextlib import contextmanager
from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QListWidget, QListWidgetItem,
    QComboBox, QLabel, QPushButton, QApplication,
    QDialog, QCheckBox, QDialogButtonBox, QMessageBox,
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure

# Make the Kd model functions importable (bare-name imports).
_KD_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "dynamiXs_Kd")
if _KD_DIR not in sys.path:
    sys.path.insert(0, _KD_DIR)
from kd_models import (csp_model, intensity_decay,  # noqa: E402
                       pair_observable as _pair_observable,  # re-exported for tests
                       ref_point_values, ref_vs_point_table,
                       residue_sort_key as _residue_sort_key)
from kd_fit import (fit_residue_csp, fit_residue_intensity, json_safe,  # noqa: E402
                    _GLOBAL_R2_MIN)


_KD_OUTLIER_K = 3.5   # Hampel: drop a Kd more than k robust-SD (1.4826·MAD) from the median

# Matches lunaNMR/cli.py's _INTENSITY_YLIM exactly: every intensity curve panel (raw-scale
# or already-normalized) shares this 0-1 I/I(0) axis, so CLI- and GUI-exported figures are
# directly comparable/overlay-able regardless of which surface produced them.
_INTENSITY_YLIM = (-0.05, 1.05)


def _fmt_num(value, spec):
    """Format a fit value, or 'n/a' when it is missing or non-finite. A degenerate
    fit (e.g. too few titration points) can yield an unbounded covariance error;
    json_safe writes such non-finite floats as null, so value may be None or inf."""
    if value is None or not np.isfinite(value):
        return "n/a"
    return format(value, spec)


def _hampel_keep(vals, k=_KD_OUTLIER_K):
    """Drop values more than k robust-SDs (1.4826·MAD) ABOVE the median — the weak /
    undetermined-Kd residues (Kd off by multiples) that inflate the mean/std even when
    their own fit is good. One-sided on purpose: low Kd is genuine tight binding, not an
    artifact, so it is kept. No-op for <3 values or a zero MAD."""
    if len(vals) < 3:
        return list(vals)
    med = float(np.median(vals))
    mad = float(np.median([abs(x - med) for x in vals]))
    if mad <= 0:
        return list(vals)
    hi = med + k * 1.4826 * mad
    return [x for x in vals if x <= hi]


def _kd_with_err(kd, err, kd_spec=".4g", err_spec=".3g"):
    """Render 'Kd ± err', dropping the '± err' when the error is missing/non-finite
    (e.g. an older fit with no stored Kd_err, or a singular covariance)."""
    s = _fmt_num(kd, kd_spec)
    if err is not None and np.isfinite(err):
        s += f" ± {format(err, err_spec)}"
    return s


class _ExportChoiceDialog(QDialog):
    """Pick which fit figures to export as editable-vector PDF. All checked by
    default so 'export figures → OK' writes the full set."""

    def __init__(self, parent, dest):
        super().__init__(parent)
        self.setWindowTitle("Export figures")
        lay = QVBoxLayout(self)
        lay.addWidget(QLabel("Export editable-vector PDF (ref→point also writes CSV + JSON):"))
        self.cb_intensity_fit = QCheckBox("1  Intensity vs titration (per-residue fits)")
        self.cb_intensity_ref = QCheckBox("2  Intensity I/I₀: reference → point (+ CSV/JSON)")
        self.cb_csp_ref = QCheckBox("3  CSP: reference → point (+ CSV/JSON)")
        self.cb_csp_kd = QCheckBox("4  CSP: Kd vs residue (+ global Kd line)")
        self.cb_intensity_kd = QCheckBox("5  Intensity: Kd vs residue (+ global apparent Kd line)")
        self.cb_csp_global = QCheckBox("6  CSP: global shared-Kd fit over data (per residue)")
        self.cb_intensity_global = QCheckBox("7  Intensity: global shared-Kd fit over data (per residue)")
        for cb in (self.cb_intensity_fit, self.cb_intensity_ref, self.cb_csp_ref,
                   self.cb_csp_kd, self.cb_intensity_kd, self.cb_csp_global,
                   self.cb_intensity_global):
            cb.setChecked(True)
            lay.addWidget(cb)
        lay.addWidget(QLabel(f"Destination:\n{dest}"))
        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(self.accept)
        btns.rejected.connect(self.reject)
        lay.addWidget(btns)

    def choices(self):
        return {"intensity_fit": self.cb_intensity_fit.isChecked(),
                "intensity_ref": self.cb_intensity_ref.isChecked(),
                "csp_ref": self.cb_csp_ref.isChecked(),
                "csp_kd": self.cb_csp_kd.isChecked(),
                "intensity_kd": self.cb_intensity_kd.isChecked(),
                "csp_global": self.cb_csp_global.isChecked(),
                "intensity_global": self.cb_intensity_global.isChecked()}


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

        export_btn = QPushButton("Export figures…")
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

    def _clean_kds(self, obs):
        """Per-residue Kd for `obs` used in the summary: successful, not user-excluded,
        well-fit (R² ≥ _GLOBAL_R2_MIN), then a Hampel reject to drop residues whose Kd is
        off by multiples of the median even though their own fit is good."""
        excl = self.excluded_residues.get(obs, set())
        vals = []
        for f in self.fits:
            fit = f.get(obs, {})
            v = fit.get("Kd")
            r2 = fit.get("r_squared")
            if (fit.get("success") and f.get("residue") not in excl
                    and isinstance(v, (int, float)) and np.isfinite(v)
                    and isinstance(r2, (int, float)) and r2 >= _GLOBAL_R2_MIN):
                vals.append(float(v))
        return _hampel_keep(vals)

    def _mean_kd(self, obs):
        """Mean of the cleaned per-residue Kd (see _clean_kds). Returns (mean, n) or
        (None, 0) if there are none."""
        vals = self._clean_kds(obs)
        return (float(np.mean(vals)), len(vals)) if vals else (None, 0)

    def _update_global_label(self):
        """Show the global shared Kd AND the per-residue average for the CURRENTLY
        selected observable, so the label tracks the Observable selector (intensity mode
        shows the intensity Kd, not CSP). The global fit is from the initial run —
        per-residue refits don't update it."""
        obs = self._obs_key()
        g = self.data.get("global", {}).get(obs, {})
        mean_kd, n_mean = self._mean_kd(obs)
        avg = f"   |   per-residue avg: {mean_kd:.4g} (n={n_mean})" if mean_kd is not None else ""
        if g.get("success"):
            lbl = ("Global shared Kd (CSP)" if obs == "csp"
                   else "Global shared apparent Kd (intensity)")
            # A Kd the titration could not resolve looks identical to one it could,
            # so the verdict belongs next to the number, not only in the JSON.
            verdict = "" if g.get("reliable") else "   ⚠ NOT reliable"
            self.global_label.setText(
                f"{lbl}: {_kd_with_err(g.get('Kd'), g.get('Kd_err'))} "
                f"(n={g.get('n_residues','?')})" + verdict + avg +
                "   — global from initial fit, not updated by per-residue refits")
        else:
            which = "CSP" if obs == "csp" else "intensity"
            self.global_label.setText(
                f"Global shared Kd ({which}): n/a  (re-run the fit to compute it)" + avg)

    def _obs_key(self):
        return "csp" if self.obs_combo.currentIndex() == 0 else "intensity"

    def _current_residue(self):
        row = self.residue_list.currentRow()
        return self.fits[row] if 0 <= row < len(self.fits) else None

    # ---------- plotting ----------

    def _refresh(self, *_):
        self._update_global_label()          # track the Observable selector
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
        """Interactive comparison view: draw the observable between the reference
        and compare points chosen in the selectors."""
        self._draw_comparison(ax, self.ref_point_combo.currentIndex(),
                              self.cmp_point_combo.currentIndex(), self._obs_key())

    def _draw_comparison(self, ax, i, j, obs):
        """Observable (CSP or I/I₀) per residue between reference point i and
        comparison point j — computed from each residue's raw series."""
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
        # "Missing" is the SAME for both observables: the peak must be present (detected
        # position + a real intensity) at both endpoints, so a residue isn't shown as
        # present in one mode and grey in the other. (ref_point_values applies that rule.)
        names, vals = ref_point_values(self.fits, i, j, obs, alpha=alpha, value=value)
        bar_color = "seagreen" if obs == "csp" else "indianred"
        # I/I₀ is a ratio → fix the axis at 0–1; CSP is in ppm → auto-scale.
        self._plot_residue_bars(ax, names, vals, bar_color,
                                ymax=None if obs == "csp" else 1.0, obs=obs)
        ax.set_ylabel("CSP (ppm)" if obs == "csp" else "Intensity ratio I/I₀")
        sym = "CSP" if obs == "csp" else "I/I₀"
        ax.set_title(f"{sym}:  {labels[i]} → {labels[j]}  (ref → point)")

    def _plot_residue_bars(self, ax, names, vals, color, errs=None, ymax=None, obs=None):
        """Bar chart of value-per-residue. Residues with no value (NaN — missing
        assignment or failed fit) draw as a full-height grey bar so the gap is visible;
        residues the user excluded (problematic values) draw hatched. In edit mode the
        bars are pickable so a click toggles that residue's exclusion. `obs` selects
        which observable's exclusion set to honor (defaults to the current view). `ymax` fixes the
        y-axis top (e.g. 1.0 for an I/I₀ ratio); None auto-scales to the data."""
        self._pick_registry = {}
        excl = self.excluded_residues[obs or self._obs_key()]   # per-observable exclusions
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
        if obs == "intensity":
            y = y / fit["I0"]
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
            # Normalized to I/I(0), matching the exported figures, so the interactive
            # view and the PDF output share the same 0-1 axis.
            yg = intensity_decay(Lg, fit["I0"], fit["I_inf"], fit["Kd"]) / fit["I0"]
            ylab = "I / I(0)"
        ax.plot(Lg, yg, "r-", lw=2, label="fit")
        kd, kde, r2 = fit.get("Kd"), fit.get("Kd_err"), fit.get("r_squared")
        ax.set_title(f"{f['residue']}   Kd = {_fmt_num(kd, '.3g')} ± {_fmt_num(kde, '.2g')}"
                     f"   R² = {_fmt_num(r2, '.3f')}")
        ax.set_xlabel("[Ligand]")
        ax.set_ylabel(ylab)
        if obs == "intensity":
            ax.set_ylim(*_INTENSITY_YLIM)
        ax.legend()

    def _plot_summary(self, ax, key, err_key, title, ylab, obs=None):
        obs = obs or self._obs_key()
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
        self._plot_residue_bars(ax, names, vals, "steelblue", errs=errs, obs=obs)
        if key == "Kd":
            g = self.data.get("global", {}).get(obs, {})
            gkd = g.get("Kd")
            drew = False
            if g.get("success") and gkd is not None and np.isfinite(gkd):
                lbl = "global Kd" if obs == "csp" else "global apparent Kd"
                ax.axhline(gkd, color="red", ls="--",
                           label=f"{lbl}={_kd_with_err(gkd, g.get('Kd_err'), '.3g', '.2g')}")
                drew = True
            mean_kd, n_mean = self._mean_kd(obs)
            if mean_kd is not None:
                ax.axhline(mean_kd, color="green", ls=":",
                           label=f"per-residue avg={mean_kd:.3g} (n={n_mean})")
                drew = True
            if drew:
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

    # ---------- figure export ----------

    def _export_base(self):
        """Series-name prefix for exported files, taken from the fit's filename
        (`<series>_kd_fit_data.json` → `<series>`). For a generic filename (e.g. a
        project-bundled `fit_data.json`) it falls back to the name stored in metadata,
        then to 'kd'."""
        stem = os.path.basename(self.json_file or "")
        suffix = "_kd_fit_data.json"
        if stem.endswith(suffix):
            base = stem[:-len(suffix)]
            if base:
                return base
        name = (self.data.get("metadata", {}) or {}).get("name")
        if name:
            return str(name)
        return "kd"

    def _export(self):
        """Export the fit figures as editable-vector PDFs, written into the folder
        where the fit results live (next to the *_kd_fit_data.json)."""
        if not self.json_file:
            QMessageBox.warning(self, "No data", "Load a fit first.")
            return
        dest = os.path.dirname(os.path.abspath(self.json_file))
        dlg = _ExportChoiceDialog(self, dest)
        if not dlg.exec():
            return
        ch = dlg.choices()
        if not any(ch.values()):
            return
        base = self._export_base()
        jobs = [
            (ch["intensity_fit"], f"{base}_intensity_titration_fits.pdf",
             lambda p: self._export_intensity_fits(p)),
            (ch["intensity_ref"], f"{base}_intensity_ref_vs_point.pdf",
             lambda p: self._export_ref_vs_point(p, "intensity")),
            (ch["csp_ref"], f"{base}_csp_ref_vs_point.pdf",
             lambda p: self._export_ref_vs_point(p, "csp")),
            (ch["csp_kd"], f"{base}_csp_kd_vs_residue.pdf",
             lambda p: self._export_kd_vs_residue(p, "csp")),
            (ch["intensity_kd"], f"{base}_intensity_kd_vs_residue.pdf",
             lambda p: self._export_kd_vs_residue(p, "intensity")),
            (ch["csp_global"], f"{base}_csp_global_fit.pdf",
             lambda p: self._export_global_fit(p, "csp")),
            (ch["intensity_global"], f"{base}_intensity_global_fit.pdf",
             lambda p: self._export_global_fit(p, "intensity")),
        ]
        written = []
        try:
            for wanted, fname, fn in jobs:
                if wanted and fn(os.path.join(dest, fname)):
                    written.append(fname)
            # The ref→point observable data as Excel-ready CSV + JSON, alongside the PDFs.
            for wanted, obs in [(ch["intensity_ref"], "intensity"), (ch["csp_ref"], "csp")]:
                if wanted:
                    for p in self._write_ref_vs_point_data(dest, base, obs):
                        written.append(os.path.basename(p))
            # The global shared-Kd per-residue params as CSV + JSON, alongside the PDFs.
            for wanted, obs in [(ch["intensity_global"], "intensity"), (ch["csp_global"], "csp")]:
                if wanted:
                    for p in self._write_global_fit_data(dest, base, obs):
                        written.append(os.path.basename(p))
        except Exception as e:
            QMessageBox.critical(self, "Export failed", str(e))
            return
        if written:
            QMessageBox.information(self, "Export complete",
                                    "Wrote to:\n" + dest + "\n\n" + "\n".join(written))
        else:
            QMessageBox.warning(self, "Nothing exported",
                                "The selected figures had no data to plot "
                                "(need embedded series and/or successful fits).")

    def _export_intensity_fits(self, path):
        """Per-residue intensity-decay curves (I vs [L] + fit), 20 per page like the
        T1/T2 export. Returns the number of PDF pages written (0 if no fits)."""
        fits = [f for f in self.fits
                if f.get("intensity", {}).get("success") and f.get("intensity", {}).get("L")]
        if not fits:
            return 0
        per_page = 20                          # 5 rows × 4 cols, matching the CLI export
        pages = 0
        with self._vector_pdf(path) as pdf:
            for start in range(0, len(fits), per_page):
                chunk = fits[start:start + per_page]
                fig = Figure(figsize=(16, 15))    # matches lunaNMR/cli.py's (cols*4, rows*3)
                axes = np.asarray(fig.subplots(5, 4)).flatten()
                for ax, f in zip(axes, chunk):
                    self._draw_intensity_fit(ax, f)
                for ax in axes[len(chunk):]:
                    ax.set_visible(False)
                fig.tight_layout()
                pdf.savefig(fig)
                pages += 1
        return pages

    def _export_kd_vs_residue(self, path, obs):
        """Per-residue Kd bar chart (+ the global-Kd line) for one observable, one
        page. Returns 1 if written, 0 if that observable has no successful fits."""
        if not any(f.get(obs, {}).get("success") for f in self.fits):
            return 0
        with self._vector_pdf(path) as pdf:
            fig = Figure(figsize=(11, 6))
            ax = fig.add_subplot(111)
            self._plot_summary(ax, "Kd", "Kd_err", "Kd vs residue", "Kd", obs=obs)
            fig.tight_layout()
            pdf.savefig(fig)
        return 1

    def _export_global_fit(self, path, obs):
        """Per-residue observed data overlaid with the GLOBAL (shared-Kd) model curve,
        so one can judge how well a single Kd fits every residue. 20 panels/page. The
        per-panel R² is the observed vs the shared-Kd curve (a per-residue goodness under
        the global model). Returns pages written (0 if no global fit for this obs)."""
        g = self.data.get("global", {}).get(obs, {})
        amp = g.get("dd_max") if obs == "csp" else g.get("I0")
        if not g.get("success") or not amp:
            return 0
        panels = [f for f in self.fits
                  if f.get("residue") in amp and f.get(obs, {}).get("L")]
        if not panels:
            return 0
        per_page = 20
        pages = 0
        with self._vector_pdf(path) as pdf:
            for start in range(0, len(panels), per_page):
                chunk = panels[start:start + per_page]
                fig = Figure(figsize=(16, 15))    # matches lunaNMR/cli.py's (cols*4, rows*3)
                axes = np.asarray(fig.subplots(5, 4)).flatten()
                for ax, f in zip(axes, chunk):
                    self._draw_global_fit_panel(ax, f, obs, g)
                for ax in axes[len(chunk):]:
                    ax.set_visible(False)
                fig.suptitle(f"Global shared-Kd fit ({obs}): one Kd = "
                             f"{_fmt_num(g.get('Kd'), '.4g')} for all residues "
                             "(per-residue amplitudes)", fontsize=15)
                fig.tight_layout(rect=[0, 0, 1, 0.99])
                pdf.savefig(fig)
                pages += 1
        return pages

    def _draw_global_fit_panel(self, ax, f, obs, g):
        """One residue: observed points + the shared-Kd global-model curve. Intensity is
        normalized to I/I(0) (dividing by the global fit's per-residue I0) so it shares
        the same 0-1 axis as every other intensity panel, matching lunaNMR/cli.py."""
        res = f["residue"]
        fit = f[obs]
        L = np.asarray(fit["L"], float)
        y = np.asarray(fit["obs"], float)
        Lg = np.linspace(0, L.max() * 1.05, 200)
        kd = g["Kd"]
        if obs == "csp":
            yg = csp_model(Lg, g["dd_max"][res], kd, self.P0)
            yhat = csp_model(L, g["dd_max"][res], kd, self.P0)
            ylab = "CSP (ppm)"
        else:
            yg = intensity_decay(Lg, g["I0"][res], g["I_inf"][res], kd)
            yhat = intensity_decay(L, g["I0"][res], g["I_inf"][res], kd)
            ylab = "I / I(0)"
        ss_res = float(np.sum((y - yhat) ** 2))
        ss_tot = float(np.sum((y - np.mean(y)) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        if obs == "intensity":
            y = y / g["I0"][res]
            yg = yg / g["I0"][res]
            ax.set_ylim(*_INTENSITY_YLIM)
        ax.scatter(L, y, color="black", s=36, zorder=3)
        ax.plot(Lg, yg, "-", color="#1f77b4", lw=1.5)   # blue = global (vs red per-residue)
        ax.set_title(f"{res}   R²(global)={_fmt_num(r2, '.2f')}", fontsize=9)
        ax.set_xlabel("[Ligand]", fontsize=8)
        ax.set_ylabel(ylab, fontsize=8)

    def _export_ref_vs_point(self, path, obs):
        """Reference (point 0) → point bar charts of the observable per residue, one
        page per later titration point. Returns the number of pages written."""
        labels = self._point_labels()
        if len(labels) < 2 or not any(f.get("series") for f in self.fits):
            return 0
        pages = 0
        with self._vector_pdf(path) as pdf:
            for j in range(1, len(labels)):
                fig = Figure(figsize=(11, 6))
                ax = fig.add_subplot(111)
                self._draw_comparison(ax, 0, j, obs)
                fig.tight_layout()
                pdf.savefig(fig)
                pages += 1
        return pages

    def _write_ref_vs_point_data(self, dest, base, obs):
        """Write the wide CSV + JSON for one observable's reference→point data, so the
        bars are reproducible in Excel. Returns the written paths ([] if no data)."""
        labels = self._point_labels()
        if len(labels) < 2 or not any(f.get("series") for f in self.fits):
            return []
        meta = self.data.get("metadata", {})
        alpha = float(meta.get("alpha", 0.14))
        value = meta.get("intensity_value", "height")
        from kd_export import export_ref_vs_point
        return export_ref_vs_point(os.path.join(dest, f"{base}_{obs}_ref_vs_point"),
                                   self.fits, labels, obs, alpha=alpha, value=value)

    def _write_global_fit_data(self, dest, base, obs):
        """Write the CSV + JSON of the per-residue global shared-Kd params for one
        observable, so the global-fit figure is reproducible as data. Returns the written
        paths ([] if there is no successful global fit for this observable)."""
        from kd_export import export_global_fit
        return export_global_fit(os.path.join(dest, f"{base}_{obs}_global_fit"),
                                 self.fits, self.data.get("global", {}) or {},
                                 obs, self.P0)

    @contextmanager
    def _vector_pdf(self, path):
        """A PdfPages context that keeps text editable in the output (pdf.fonttype=42
        embeds TrueType so glyphs stay selectable in Illustrator)."""
        import matplotlib as mpl
        from matplotlib.backends.backend_pdf import PdfPages
        with mpl.rc_context({"svg.fonttype": "none", "pdf.fonttype": 42}):
            with PdfPages(path) as pdf:
                yield pdf

    def _draw_intensity_fit(self, ax, f):
        """Draw one residue's intensity data + fitted decay curve onto ax (compact,
        for the multi-plot export page). Normalized to I/I(0) (dividing by the fitted
        I0 amplitude) and pinned to the shared 0-1 axis, matching lunaNMR/cli.py, so
        raw-scale residues (I0 in the hundreds) don't auto-scale to invisible flat
        lines and every panel/document is directly comparable."""
        fit = f["intensity"]
        L = np.asarray(fit["L"], float)
        y = np.asarray(fit["obs"], float) / fit["I0"]
        excl = self.excluded["intensity"].get(f["residue"], set())
        inc = [i for i in range(len(L)) if i not in excl]
        exc = [i for i in range(len(L)) if i in excl]
        if inc:
            ax.scatter(L[inc], y[inc], color="black", s=36, zorder=3)
        if exc:
            ax.scatter(L[exc], y[exc], color="#888888", marker="x", s=36,
                       linewidths=2, zorder=3)
        Lg = np.linspace(0, L.max() * 1.05, 200)
        yg = intensity_decay(Lg, fit["I0"], fit["I_inf"], fit["Kd"]) / fit["I0"]
        ax.plot(Lg, yg, "r-", lw=1.5)
        kd, kde, r2 = fit.get("Kd"), fit.get("Kd_err"), fit.get("r_squared")
        ax.set_title(f"{f['residue']}   Kd={_fmt_num(kd, '.3g')}±{_fmt_num(kde, '.2g')}"
                     f"  R²={_fmt_num(r2, '.2f')}", fontsize=9)
        ax.set_xlabel("[Ligand]", fontsize=8)
        ax.set_ylabel("I / I(0)", fontsize=8)
        ax.set_ylim(*_INTENSITY_YLIM)


def open_kd_titration_viewer(parent=None, json_file=None):
    viewer = KdTitrationFitViewer(parent=parent, json_file=json_file)
    viewer.show()
    return viewer
