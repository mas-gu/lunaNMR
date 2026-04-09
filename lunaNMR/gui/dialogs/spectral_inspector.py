# ABOUTME: Standalone spectral inspector window for visual comparison of multiple NMR spectra
# ABOUTME: Supports overlay display, peak list management, and visual peak repositioning

from __future__ import annotations

import os
import uuid
import logging
from dataclasses import dataclass, field
from typing import List, Optional

from PySide6.QtWidgets import (
    QWidget, QSplitter, QLabel, QVBoxLayout, QHBoxLayout,
    QPushButton, QCheckBox, QToolBar, QTreeWidget, QTreeWidgetItem,
    QFileDialog, QInputDialog, QMessageBox, QDoubleSpinBox, QSpinBox,
    QGroupBox, QColorDialog, QMenu, QAbstractItemView
)
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor

from lunaNMR.gui.base.base_window import BaseWindow

logger = logging.getLogger(__name__)

# Default color palette — mirrors MultiSpectrumViewerDialog
DEFAULT_COLORS = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
]


# ---------------------------------------------------------------------------
# Data Model
# ---------------------------------------------------------------------------

@dataclass
class InspectorPeak:
    """One peak in the inspector's local copy.

    Fields mirror the fitted_peaks dict format used throughout the codebase
    so that existing plot helpers can be reused without modification.
    """
    peak_id: str = field(default_factory=lambda: str(uuid.uuid4()))
    ppm_x: float = 0.0
    ppm_y: float = 0.0
    assignment: str = ''
    quality: str = 'Unknown'

    def to_dict(self) -> dict:
        """Return dict compatible with SpectrumPlotter peak drawing methods."""
        return {
            'peak_id': self.peak_id,
            'ppm_x': self.ppm_x,
            'ppm_y': self.ppm_y,
            'center_x': self.ppm_x,   # alias expected by SpectrumPlotter
            'center_y': self.ppm_y,   # alias expected by SpectrumPlotter
            'assignment': self.assignment,
            'quality': self.quality,
        }

    @classmethod
    def from_dict(cls, d: dict) -> InspectorPeak:
        """Construct from any fitted_peaks / peak list dict."""
        def _first_not_none(*keys):
            for k in keys:
                v = d.get(k)
                if v is not None:
                    return v
            return None

        x = _first_not_none('ppm_x', 'center_x', 'peak_x')
        y = _first_not_none('ppm_y', 'center_y', 'peak_y')
        assignment = _first_not_none('assignment', 'Assignment')
        quality = (_first_not_none('quality', 'Quality', 'fitting_quality') or 'Unknown')
        return cls(
            ppm_x=float(x) if x is not None else 0.0,
            ppm_y=float(y) if y is not None else 0.0,
            assignment=str(assignment) if assignment is not None else '',
            quality=str(quality),
        )


@dataclass
class InspectorSpectrum:
    """All inspector-local state for one loaded spectrum."""
    spectrum_id: str = field(default_factory=lambda: str(uuid.uuid4()))
    display_name: str = ''
    file_path: str = ''
    color: str = field(default_factory=lambda: DEFAULT_COLORS[0])
    visible: bool = True

    # Per-spectrum contour display parameters
    contour_levels: int = 10
    contour_min_pct: float = 5.0
    contour_increment: float = 1.3

    # Loaded via EnhancedVoigtIntegrator (same pattern as MultiSpectrumViewerDialog)
    integrator: Optional[object] = None
    loaded: bool = False

    # Local peak list — inspector-only copy, never pushed back to main window
    peaks: List[InspectorPeak] = field(default_factory=list)


@dataclass
class InspectorGroup:
    """A named folder in the spectrum library tree."""
    group_id: str = field(default_factory=lambda: str(uuid.uuid4()))
    name: str = 'Group'
    spectra: List[InspectorSpectrum] = field(default_factory=list)


# Tree item type constants
ITEM_GROUP = 0
ITEM_SPECTRUM = 1
ITEM_PEAKLIST = 2

# QTreeWidgetItem data roles
_ROLE_ID = Qt.UserRole        # stores group_id / spectrum_id
_ROLE_TYPE = Qt.UserRole + 1  # stores ITEM_GROUP / ITEM_SPECTRUM / ITEM_PEAKLIST


# ---------------------------------------------------------------------------
# Spectrum row widget (used as a tree item widget)
# ---------------------------------------------------------------------------

class SpectrumRowWidget(QWidget):
    """Compact row for one spectrum: visibility checkbox, color swatch, name, peak-list button."""

    visibility_changed = Signal(str, bool)   # spectrum_id, visible
    color_changed = Signal(str, str)          # spectrum_id, color_hex
    peak_list_requested = Signal(str)         # spectrum_id

    def __init__(self, spectrum: InspectorSpectrum, parent=None):
        super().__init__(parent)
        self._spectrum_id = spectrum.spectrum_id
        self._current_color = spectrum.color

        layout = QHBoxLayout()
        layout.setContentsMargins(4, 2, 4, 2)
        layout.setSpacing(4)

        self._checkbox = QCheckBox()
        self._checkbox.setChecked(spectrum.visible)
        self._checkbox.setToolTip('Toggle visibility')
        self._checkbox.toggled.connect(self._on_visibility)
        layout.addWidget(self._checkbox)

        self._color_btn = QPushButton()
        self._color_btn.setFixedSize(18, 18)
        self._color_btn.setStyleSheet(
            f"background-color: {spectrum.color}; border: 1px solid gray; border-radius: 2px;"
        )
        self._color_btn.setToolTip('Change color')
        self._color_btn.clicked.connect(self._on_color)
        layout.addWidget(self._color_btn)

        self._label = QLabel(spectrum.display_name or os.path.basename(spectrum.file_path))
        layout.addWidget(self._label, stretch=1)

        self._peak_btn = QPushButton('+')
        self._peak_btn.setFixedSize(20, 20)
        self._peak_btn.setToolTip('Attach peak list from file')
        self._peak_btn.clicked.connect(lambda: self.peak_list_requested.emit(self._spectrum_id))
        layout.addWidget(self._peak_btn)

        self.setLayout(layout)

    def _on_visibility(self, checked: bool):
        self.visibility_changed.emit(self._spectrum_id, checked)

    def _on_color(self):
        color = QColorDialog.getColor(QColor(self._current_color), self, 'Choose Color')
        if color.isValid():
            self._current_color = color.name()
            self._color_btn.setStyleSheet(
                f"background-color: {self._current_color}; border: 1px solid gray; border-radius: 2px;"
            )
            self.color_changed.emit(self._spectrum_id, self._current_color)


# ---------------------------------------------------------------------------
# Drag-aware tree widget
# ---------------------------------------------------------------------------

class _LibraryTree(QTreeWidget):
    """QTreeWidget that emits spectrum_moved when a spectrum is dropped onto a different group."""

    spectrum_moved = Signal(str, str)  # spectrum_id, target_group_id

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setDragEnabled(True)
        self.setAcceptDrops(True)
        self.setDropIndicatorShown(True)
        self.setDragDropMode(QAbstractItemView.InternalMove)
        self.setDefaultDropAction(Qt.MoveAction)

    def dropEvent(self, event):
        dragged = self.currentItem()
        if dragged is None or dragged.data(0, _ROLE_TYPE) != ITEM_SPECTRUM:
            event.ignore()
            return

        spectrum_id = dragged.data(0, _ROLE_ID)
        target_item = self.itemAt(event.position().toPoint())

        # Resolve drop target to a group item
        target_group_item = None
        if target_item is not None:
            if target_item.data(0, _ROLE_TYPE) == ITEM_GROUP:
                target_group_item = target_item
            elif target_item.data(0, _ROLE_TYPE) == ITEM_SPECTRUM:
                target_group_item = target_item.parent()

        if target_group_item is None:
            event.ignore()
            return

        source_group_item = dragged.parent()
        if source_group_item is target_group_item:
            event.ignore()
            return

        target_group_id = target_group_item.data(0, _ROLE_ID)
        self.spectrum_moved.emit(spectrum_id, target_group_id)
        event.accept()


# ---------------------------------------------------------------------------
# Spectrum Library Panel
# ---------------------------------------------------------------------------

class SpectrumLibraryPanel(QWidget):
    """Left panel: groups/spectra tree with load controls and contour settings."""

    new_group_requested = Signal(str)            # group name
    load_spectrum_requested = Signal(str, str)  # file_path, group_id
    spectrum_visibility_changed = Signal(str, bool)
    spectrum_color_changed = Signal(str, str)
    peak_list_attach_requested = Signal(str)    # spectrum_id
    contour_update_requested = Signal(int, float, float)  # levels, min_pct, increment
    active_spectrum_changed = Signal(str)       # spectrum_id
    group_delete_requested = Signal(str)        # group_id
    spectrum_delete_requested = Signal(str)     # spectrum_id
    group_renamed = Signal(str, str)            # group_id, new_name
    spectrum_renamed = Signal(str, str)         # spectrum_id, new_name
    group_properties_requested = Signal(str)    # group_id
    spectrum_contour_requested = Signal(str)    # spectrum_id
    spectrum_moved = Signal(str, str)           # spectrum_id, target_group_id

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumWidth(220)
        self.setMaximumWidth(400)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)

        # Top buttons
        btn_row = QHBoxLayout()
        self._new_group_btn = QPushButton('New Group')
        self._new_group_btn.clicked.connect(self._on_new_group)
        btn_row.addWidget(self._new_group_btn)

        self._load_btn = QPushButton('Load Spectrum')
        self._load_btn.clicked.connect(self._on_load_spectrum)
        btn_row.addWidget(self._load_btn)
        layout.addLayout(btn_row)

        # Tree
        self._tree = _LibraryTree()
        self._tree.setHeaderHidden(True)
        self._tree.setIndentation(12)
        self._tree.setAnimated(True)
        self._tree.setContextMenuPolicy(Qt.CustomContextMenu)
        self._tree.customContextMenuRequested.connect(self._on_context_menu)
        self._tree.currentItemChanged.connect(self._on_selection_changed)
        self._tree.spectrum_moved.connect(self.spectrum_moved)
        layout.addWidget(self._tree, stretch=1)

        # Contour controls
        contour_box = QGroupBox('Contour')
        cl = QVBoxLayout(contour_box)
        cl.setSpacing(2)

        row1 = QHBoxLayout()
        row1.addWidget(QLabel('Levels:'))
        self._levels_spin = QSpinBox()
        self._levels_spin.setRange(1, 30)
        self._levels_spin.setValue(10)
        row1.addWidget(self._levels_spin)
        cl.addLayout(row1)

        row2 = QHBoxLayout()
        row2.addWidget(QLabel('Min %:'))
        self._min_spin = QDoubleSpinBox()
        self._min_spin.setRange(0.01, 100.0)
        self._min_spin.setDecimals(2)
        self._min_spin.setValue(5.0)
        row2.addWidget(self._min_spin)
        cl.addLayout(row2)

        row3 = QHBoxLayout()
        row3.addWidget(QLabel('Step:'))
        self._incr_spin = QDoubleSpinBox()
        self._incr_spin.setRange(1.01, 5.0)
        self._incr_spin.setDecimals(2)
        self._incr_spin.setValue(1.3)
        row3.addWidget(self._incr_spin)
        cl.addLayout(row3)

        self._update_btn = QPushButton('Apply to All')
        self._update_btn.clicked.connect(self._on_contour_update)
        cl.addWidget(self._update_btn)

        layout.addWidget(contour_box)

        # Internal map: group_id → QTreeWidgetItem
        self._group_items: dict = {}

    # ------------------------------------------------------------------
    # Public API — called by SpectralInspector
    # ------------------------------------------------------------------

    def add_group(self, group: InspectorGroup) -> QTreeWidgetItem:
        """Add a group node to the tree and return it."""
        item = QTreeWidgetItem(self._tree, [group.name])
        item.setData(0, _ROLE_ID, group.group_id)
        item.setData(0, _ROLE_TYPE, ITEM_GROUP)
        item.setExpanded(True)
        self._group_items[group.group_id] = item
        return item

    def add_spectrum(self, spectrum: InspectorSpectrum, group_id: str):
        """Add a spectrum row under its group item."""
        group_item = self._group_items.get(group_id)
        if group_item is None:
            return

        child = QTreeWidgetItem(group_item)
        child.setData(0, _ROLE_ID, spectrum.spectrum_id)
        child.setData(0, _ROLE_TYPE, ITEM_SPECTRUM)

        row_widget = SpectrumRowWidget(spectrum)
        row_widget.visibility_changed.connect(self.spectrum_visibility_changed)
        row_widget.color_changed.connect(self.spectrum_color_changed)
        row_widget.peak_list_requested.connect(self.peak_list_attach_requested)
        self._tree.setItemWidget(child, 0, row_widget)

        group_item.setExpanded(True)

    def _selected_group_id(self) -> Optional[str]:
        """Return group_id of the currently selected group, or first group if any."""
        for item in self._tree.selectedItems():
            type_ = item.data(0, _ROLE_TYPE)
            if type_ == ITEM_GROUP:
                return item.data(0, _ROLE_ID)
            if type_ == ITEM_SPECTRUM:
                parent = item.parent()
                if parent and parent.data(0, _ROLE_TYPE) == ITEM_GROUP:
                    return parent.data(0, _ROLE_ID)
        # Fall back to first group if nothing is selected
        root = self._tree.invisibleRootItem()
        for i in range(root.childCount()):
            item = root.child(i)
            if item.data(0, _ROLE_TYPE) == ITEM_GROUP:
                return item.data(0, _ROLE_ID)
        return None

    def remove_group(self, group_id: str):
        """Remove a group tree item."""
        item = self._group_items.pop(group_id, None)
        if item is not None:
            idx = self._tree.indexOfTopLevelItem(item)
            if idx >= 0:
                self._tree.takeTopLevelItem(idx)

    def remove_spectrum(self, spectrum_id: str):
        """Remove a spectrum tree item from its parent group."""
        root = self._tree.invisibleRootItem()
        for i in range(root.childCount()):
            group_item = root.child(i)
            for j in range(group_item.childCount()):
                child = group_item.child(j)
                if child.data(0, _ROLE_ID) == spectrum_id:
                    group_item.removeChild(child)
                    return

    # ------------------------------------------------------------------
    # Slot handlers
    # ------------------------------------------------------------------

    def _on_selection_changed(self, current, previous):
        if current is None:
            return
        type_ = current.data(0, _ROLE_TYPE)
        if type_ == ITEM_SPECTRUM:
            self.active_spectrum_changed.emit(current.data(0, _ROLE_ID))
        elif type_ == ITEM_GROUP and current.childCount() > 0:
            first = current.child(0)
            if first.data(0, _ROLE_TYPE) == ITEM_SPECTRUM:
                self.active_spectrum_changed.emit(first.data(0, _ROLE_ID))

    def _on_context_menu(self, pos):
        item = self._tree.itemAt(pos)
        if item is None:
            return
        type_ = item.data(0, _ROLE_TYPE)
        menu = QMenu(self)
        if type_ == ITEM_GROUP:
            group_id = item.data(0, _ROLE_ID)
            props_act = menu.addAction('Group Properties…')
            menu.addSeparator()
            rename_act = menu.addAction('Rename')
            delete_act = menu.addAction('Delete Group')
            action = menu.exec(self._tree.viewport().mapToGlobal(pos))
            if action == props_act:
                self.group_properties_requested.emit(group_id)
            elif action == rename_act:
                name, ok = QInputDialog.getText(self, 'Rename Group', 'New name:', text=item.text(0))
                if ok and name.strip():
                    item.setText(0, name.strip())
                    self.group_renamed.emit(group_id, name.strip())
            elif action == delete_act:
                self.group_delete_requested.emit(group_id)
        elif type_ == ITEM_SPECTRUM:
            spectrum_id = item.data(0, _ROLE_ID)
            contour_act = menu.addAction('Contour Settings…')
            menu.addSeparator()
            rename_act = menu.addAction('Rename')
            delete_act = menu.addAction('Remove Spectrum')
            action = menu.exec(self._tree.viewport().mapToGlobal(pos))
            if action == contour_act:
                self.spectrum_contour_requested.emit(spectrum_id)
            elif action == rename_act:
                widget = self._tree.itemWidget(item, 0)
                current_name = widget._label.text() if widget else ''
                name, ok = QInputDialog.getText(self, 'Rename Spectrum', 'New name:', text=current_name)
                if ok and name.strip():
                    if widget:
                        widget._label.setText(name.strip())
                    self.spectrum_renamed.emit(spectrum_id, name.strip())
            elif action == delete_act:
                self.spectrum_delete_requested.emit(spectrum_id)

    def _on_new_group(self):
        name, ok = QInputDialog.getText(self, 'New Group', 'Group name:')
        if ok and name.strip():
            self.new_group_requested.emit(name.strip())

    def _on_load_spectrum(self):
        # Auto-create a default group if none exists yet
        group_id = self._selected_group_id()
        if group_id is None:
            self.new_group_requested.emit('Default')
            group_id = self._selected_group_id()
        if group_id is None:
            return  # still none after auto-create attempt (shouldn't happen)

        paths, _ = QFileDialog.getOpenFileNames(
            self, 'Load NMR Spectra', '',
            'NMR files (*.ft2 *.ft1 *.ft *.pipe *.ucsf *.fid);;All files (*)'
        )
        for path in paths:
            self.load_spectrum_requested.emit(path, group_id)

    def _on_contour_update(self):
        self.contour_update_requested.emit(
            self._levels_spin.value(),
            self._min_spin.value(),
            self._incr_spin.value(),
        )


# ---------------------------------------------------------------------------
# Peak list file parsing
# ---------------------------------------------------------------------------

def _parse_peak_list_file(file_path: str) -> List[dict]:
    """Parse a CSV or JSON peak list file and return a list of dicts.

    CSV files must have a header row. JSON files must contain a list of objects.
    """
    ext = os.path.splitext(file_path)[1].lower()
    if ext == '.json':
        import json
        with open(file_path, 'r') as fh:
            data = json.load(fh)
        return data if isinstance(data, list) else []
    else:
        import csv
        with open(file_path, 'r', newline='') as fh:
            return list(csv.DictReader(fh))


# ---------------------------------------------------------------------------
# Inspector Canvas
# ---------------------------------------------------------------------------

class InspectorCanvas(QWidget):
    """Overlay canvas: draws N spectra as contour plots with optional peak markers."""

    peak_add_requested = Signal(str, float, float)  # spectrum_id, ppm_x, ppm_y
    peak_delete_requested = Signal(str, str)         # spectrum_id, peak_id

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        from lunaNMR.gui.components.matplotlib_widget import MatplotlibWidget
        self._plot = MatplotlibWidget(figsize=(8, 6), toolbar=False)
        layout.addWidget(self._plot)

        # Peak interaction state
        self._current_spectra: List[InspectorSpectrum] = []
        self._active_spectrum_id: Optional[str] = None
        self._selected_peak: Optional[InspectorPeak] = None
        self._selected_spectrum_id: Optional[str] = None
        self._highlight = None  # matplotlib Ellipse artist

        # Navigation handler (pan, zoom, keyboard, peak callbacks)
        from lunaNMR.gui.components.nmr_navigation_handler import NMRNavigationHandler
        self._nav_handler = NMRNavigationHandler()
        self._nav_handler.attach(self._plot)
        self._nav_handler.on_peak_select = self._on_nav_peak_select
        self._nav_handler.on_peak_edit = self._on_nav_peak_edit
        self._nav_handler.on_delete = self._on_nav_delete
        self._nav_handler.on_escape = self._on_nav_escape
        self._nav_handler.on_reset_zoom = self._reset_zoom

        self._show_empty_state()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def set_active_spectrum(self, spectrum_id: Optional[str]):
        """Set which spectrum receives new peaks from Shift+click."""
        self._active_spectrum_id = spectrum_id

    def update_plot(self, spectra: List[InspectorSpectrum]):
        """Redraw all visible loaded spectra using each spectrum's own contour settings."""
        import numpy as np
        import matplotlib.pyplot as plt

        self._current_spectra = spectra
        ax = self._plot.axes

        # Preserve zoom state across redraws
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        has_zoom = not (xlim[0] == 0.0 and xlim[1] == 1.0
                        and ylim[0] == 0.0 and ylim[1] == 1.0)

        self._plot.clear()

        visible = [s for s in spectra if s.visible and s.loaded and s.integrator is not None]

        if not visible:
            self._show_empty_state()
            self._plot.refresh()
            return

        legend_handles = []
        legend_labels = []

        for spec in visible:
            integrator = spec.integrator
            try:
                data = integrator.nmr_data
                ppm_x = integrator.ppm_x_axis
                ppm_y = integrator.ppm_y_axis
                if data is None or ppm_x is None or ppm_y is None:
                    continue

                max_int = np.max(np.abs(data))
                if max_int == 0:
                    continue

                levels = self._calculate_contour_levels(
                    max_int, spec.contour_levels, spec.contour_min_pct, spec.contour_increment
                )
                if not levels:
                    continue

                ax.contour(ppm_x, ppm_y, data,
                           levels=levels, colors=spec.color,
                           linewidths=1.0, alpha=0.8)

                if spec.peaks:
                    xs = [p.ppm_x for p in spec.peaks]
                    ys = [p.ppm_y for p in spec.peaks]
                    ax.scatter(xs, ys, marker='+', s=40, c=spec.color,
                               linewidths=1.2, zorder=5)

                legend_handles.append(
                    plt.Line2D([0], [0], color=spec.color, linewidth=2)
                )
                legend_labels.append(spec.display_name)

            except Exception as e:
                logger.warning(f"Failed to draw {spec.display_name}: {e}")
                continue

        ax.set_xlabel('¹H (ppm)', fontsize=10)
        ax.set_ylabel('¹⁵N / ¹³C (ppm)', fontsize=10)

        if has_zoom:
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
        else:
            ax.invert_xaxis()
            ax.invert_yaxis()

        if legend_handles:
            ax.legend(legend_handles, legend_labels, loc='best', fontsize=8)

        # Restore selection highlight after axes are cleared by redraw
        self._draw_selection_highlight()
        self._plot.refresh()

    @staticmethod
    def _calculate_contour_levels(max_intensity: float, num_levels: int,
                                   min_pct: float, increment: float) -> list:
        """Return multiplicatively-spaced contour levels up to max_intensity."""
        min_level = (min_pct / 100.0) * max_intensity
        levels = []
        current = min_level
        for _ in range(num_levels):
            if current > max_intensity:
                break
            levels.append(current)
            current *= increment
        return levels

    # ------------------------------------------------------------------
    # Peak interaction callbacks (wired to NMRNavigationHandler)
    # ------------------------------------------------------------------

    def _on_nav_peak_select(self, x: float, y: float):
        """Middle-click: select the nearest peak within threshold."""
        ax = self._plot.axes
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        xlim_range = abs(xlim[1] - xlim[0])
        ylim_range = abs(ylim[1] - ylim[0])
        result = self._find_nearest_peak_in_spectra(
            x, y, self._current_spectra, xlim_range, ylim_range
        )
        if result:
            self._selected_peak, self._selected_spectrum_id = result
        else:
            self._selected_peak = None
            self._selected_spectrum_id = None
        self._draw_selection_highlight()
        self._plot.canvas.draw_idle()

    def _on_nav_peak_edit(self, x: float, y: float, modifiers):
        """Shift+click: add a new peak to the active spectrum."""
        if self._active_spectrum_id:
            self.peak_add_requested.emit(self._active_spectrum_id, x, y)

    def _on_nav_delete(self):
        """Delete key: remove the selected peak."""
        if self._selected_peak and self._selected_spectrum_id:
            self.peak_delete_requested.emit(
                self._selected_spectrum_id, self._selected_peak.peak_id
            )
            self._selected_peak = None
            self._selected_spectrum_id = None
            self._highlight = None

    def _on_nav_escape(self):
        """Escape: deselect current peak."""
        self._selected_peak = None
        self._selected_spectrum_id = None
        self._clear_selection_highlight()
        self._plot.canvas.draw_idle()

    def _reset_zoom(self):
        self._plot.axes.autoscale()
        self._plot.refresh()

    # ------------------------------------------------------------------
    # Selection highlight
    # ------------------------------------------------------------------

    def _draw_selection_highlight(self):
        """Draw a dashed ellipse around the selected peak."""
        self._clear_selection_highlight()
        if self._selected_peak is None:
            return
        from matplotlib.patches import Ellipse
        ax = self._plot.axes
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        width = abs(xlim[1] - xlim[0]) * 0.02
        height = abs(ylim[1] - ylim[0]) * 0.02
        ell = Ellipse(
            (self._selected_peak.ppm_x, self._selected_peak.ppm_y),
            width=width, height=height,
            fill=False, edgecolor='yellow', linewidth=2.0,
            linestyle='--', zorder=10,
        )
        ax.add_patch(ell)
        self._highlight = ell

    def _clear_selection_highlight(self):
        if self._highlight is not None:
            try:
                self._highlight.remove()
            except Exception:
                pass
            self._highlight = None

    # ------------------------------------------------------------------
    # Peak lookup
    # ------------------------------------------------------------------

    @staticmethod
    def _find_nearest_peak_in_spectra(
        x: float, y: float,
        spectra: List[InspectorSpectrum],
        xlim_range: float, ylim_range: float,
        threshold: float = 0.05,
    ) -> Optional[tuple]:
        """Return (InspectorPeak, spectrum_id) for the nearest visible peak within threshold.

        Distance is normalized by the current axis ranges so that the threshold
        is a fraction of the displayed view (default 5%).
        """
        if xlim_range <= 0 or ylim_range <= 0:
            return None
        best_dist = threshold
        best = None
        for spec in spectra:
            if not spec.visible:
                continue
            for peak in spec.peaks:
                dx = (peak.ppm_x - x) / xlim_range
                dy = (peak.ppm_y - y) / ylim_range
                dist = (dx * dx + dy * dy) ** 0.5
                if dist < best_dist:
                    best_dist = dist
                    best = (peak, spec.spectrum_id)
        return best

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _show_empty_state(self):
        self._plot.axes.text(
            0.5, 0.5,
            'Load spectra to begin\n(New Group  →  Load Spectrum)',
            ha='center', va='center',
            transform=self._plot.axes.transAxes,
            fontsize=12, color='gray', style='italic',
        )


# ---------------------------------------------------------------------------
# Spectral Inspector Window
# ---------------------------------------------------------------------------

class SpectralInspector(BaseWindow):
    """Standalone window for visual inspection and comparison of NMR spectra.

    Displays multiple spectra as overlaid contour plots, each in a distinct
    color. Supports local peak list management and visual peak repositioning.
    All edits are local — nothing is propagated back to the main window.
    """

    def __init__(self, parent=None):
        super().__init__(
            title='Spectral Inspector',
            default_size=(1400, 900),
            min_size=(900, 600),
            enable_status_bar=True,
            parent=parent,
        )

        # Data model
        self.groups: List[InspectorGroup] = []
        self._color_index = 0  # Cycles through DEFAULT_COLORS

        # Contour parameters (match SpectrumLibraryPanel defaults)
        self._contour_levels = 10
        self._contour_min_pct = 5.0
        self._contour_increment = 1.3

        # Toolbar action references (spectrum_id → QAction) — populated in Increment 5
        self._toolbar_actions: dict = {}
        self._toolbar_buttons: dict = {}

        self._build_ui()
        self.update_status('Spectral Inspector ready — load spectra to begin')

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_ui(self):
        """Assemble the top-level layout: toolbar + left panel + canvas."""
        self._create_top_toolbar()

        # Central splitter: left panel | canvas
        splitter = QSplitter(Qt.Horizontal)
        splitter.setChildrenCollapsible(False)

        self._library_panel = SpectrumLibraryPanel()
        self._library_panel.new_group_requested.connect(self._on_new_group_requested)
        self._library_panel.load_spectrum_requested.connect(self._on_load_spectrum_requested)
        self._library_panel.spectrum_visibility_changed.connect(self._on_spectrum_visibility)
        self._library_panel.spectrum_color_changed.connect(self._on_spectrum_color)
        self._library_panel.peak_list_attach_requested.connect(self._on_peak_list_attach_requested)
        self._library_panel.contour_update_requested.connect(self._on_contour_update_requested)
        self._library_panel.active_spectrum_changed.connect(self._on_active_spectrum_changed)
        self._library_panel.group_delete_requested.connect(self._on_group_delete_requested)
        self._library_panel.spectrum_delete_requested.connect(self._on_spectrum_delete_requested)
        self._library_panel.group_renamed.connect(self._on_group_renamed)
        self._library_panel.spectrum_renamed.connect(self._on_spectrum_renamed)
        self._library_panel.group_properties_requested.connect(self._on_group_properties_requested)
        self._library_panel.spectrum_contour_requested.connect(self._on_spectrum_contour_requested)
        self._library_panel.spectrum_moved.connect(self._on_spectrum_moved)
        splitter.addWidget(self._library_panel)

        self._canvas = InspectorCanvas()
        self._canvas.peak_add_requested.connect(self._on_canvas_peak_add_requested)
        self._canvas.peak_delete_requested.connect(self._on_canvas_peak_delete_requested)
        splitter.addWidget(self._canvas)

        splitter.setSizes([280, 1120])
        self.setCentralWidget(splitter)

    def _create_top_toolbar(self):
        """Create the top toolbar for per-spectrum visibility toggles."""
        tb = QToolBar('Spectra', self)
        tb.setMovable(False)
        tb.setFloatable(False)

        label = QLabel('  Visible spectra: ')
        tb.addWidget(label)
        self._top_toolbar = tb
        self.addToolBar(Qt.TopToolBarArea, tb)

    # ------------------------------------------------------------------
    # Data model helpers
    # ------------------------------------------------------------------

    def _next_color(self) -> str:
        """Return the next color from the palette, cycling as needed."""
        color = DEFAULT_COLORS[self._color_index % len(DEFAULT_COLORS)]
        self._color_index += 1
        return color

    def _find_spectrum(self, spectrum_id: str) -> Optional[InspectorSpectrum]:
        """Return spectrum by id, searching all groups."""
        for group in self.groups:
            for spec in group.spectra:
                if spec.spectrum_id == spectrum_id:
                    return spec
        return None

    def _all_spectra(self) -> List[InspectorSpectrum]:
        """Return flat list of all spectra across all groups."""
        return [s for g in self.groups for s in g.spectra]

    def _add_group(self, name: str) -> InspectorGroup:
        """Create a new named group, append it, and return it."""
        group = InspectorGroup(name=name)
        self.groups.append(group)
        return group

    def _remove_group(self, group_id: str):
        """Remove a group by id (no-op if not found)."""
        self.groups = [g for g in self.groups if g.group_id != group_id]

    def _remove_spectrum(self, spectrum_id: str):
        """Remove a spectrum from whichever group contains it (no-op if not found)."""
        for group in self.groups:
            group.spectra = [s for s in group.spectra if s.spectrum_id != spectrum_id]

    def _add_peak(self, spectrum_id: str, ppm_x: float, ppm_y: float) -> Optional[InspectorPeak]:
        """Add a new peak to a spectrum's local list and return it."""
        spec = self._find_spectrum(spectrum_id)
        if spec is None:
            return None
        peak = InspectorPeak(ppm_x=ppm_x, ppm_y=ppm_y)
        spec.peaks.append(peak)
        return peak

    def _delete_peak(self, spectrum_id: str, peak_id: str):
        """Remove a peak by id from its spectrum (no-op if not found)."""
        spec = self._find_spectrum(spectrum_id)
        if spec is None:
            return
        spec.peaks = [p for p in spec.peaks if p.peak_id != peak_id]

    def _load_peak_list_from_file(self, file_path: str, spectrum_id: str) -> int:
        """Parse a CSV/JSON peak list file and append parsed peaks to the spectrum.

        Returns the number of peaks added.
        """
        spec = self._find_spectrum(spectrum_id)
        if spec is None or not os.path.exists(file_path):
            return 0
        try:
            raw = _parse_peak_list_file(file_path)
            peaks = [InspectorPeak.from_dict(d) for d in raw]
            spec.peaks.extend(peaks)
            return len(peaks)
        except Exception as e:
            logger.error(f"Error loading peak list from {file_path}: {e}")
            return 0

    def _load_spectrum_file(self, file_path: str, group: InspectorGroup) -> Optional[InspectorSpectrum]:
        """Load an NMR file, create an InspectorSpectrum, and add it to group.

        Returns the new InspectorSpectrum on success, None on failure.
        Mirrors MultiSpectrumViewerDialog._load_spectrum_data loading pattern.
        """
        if not file_path or not os.path.exists(file_path):
            logger.warning(f"File not found: {file_path}")
            return None

        try:
            from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator

            integrator = EnhancedVoigtIntegrator()
            success = integrator.load_nmr_file(file_path)
            if not success:
                logger.warning(f"Failed to load NMR file: {file_path}")
                return None

            if integrator.nmr_data is None:
                logger.warning(f"NMR data is None after loading: {file_path}")
                return None

            if not hasattr(integrator, 'ppm_x_axis') or not hasattr(integrator, 'ppm_y_axis'):
                logger.warning(f"Missing PPM axes after loading: {file_path}")
                return None

            display_name = os.path.splitext(os.path.basename(file_path))[0]
            spec = InspectorSpectrum(
                display_name=display_name,
                file_path=file_path,
                color=self._next_color(),
                contour_levels=self._contour_levels,
                contour_min_pct=self._contour_min_pct,
                contour_increment=self._contour_increment,
                integrator=integrator,
                loaded=True,
            )
            group.spectra.append(spec)
            logger.info(f"Loaded spectrum: {display_name}")
            return spec

        except Exception as e:
            logger.error(f"Error loading {file_path}: {e}")
            return None

    # ------------------------------------------------------------------
    # Slot handlers (connected to SpectrumLibraryPanel signals)
    # ------------------------------------------------------------------

    def _refresh_canvas(self):
        """Redraw the overlay canvas."""
        self._canvas.update_plot(self._all_spectra())

    def _update_spectrum_contour(self, spectrum_id: str,
                                  levels: int, min_pct: float, increment: float):
        """Update contour display params for a single spectrum."""
        spec = self._find_spectrum(spectrum_id)
        if spec is None:
            return
        spec.contour_levels = levels
        spec.contour_min_pct = min_pct
        spec.contour_increment = increment

    def _propagate_settings(self, source_id: str, target_ids: List[str], fields: List[str]):
        """Copy field values from source spectrum to each target spectrum."""
        source = self._find_spectrum(source_id)
        if source is None:
            return
        for target_id in target_ids:
            target = self._find_spectrum(target_id)
            if target is None:
                continue
            for field_name in fields:
                if hasattr(source, field_name) and hasattr(target, field_name):
                    setattr(target, field_name, getattr(source, field_name))

    def _on_new_group_requested(self, name: str):
        group = self._add_group(name)
        self._library_panel.add_group(group)
        self.update_status(f"Group '{name}' created")

    def _on_load_spectrum_requested(self, file_path: str, group_id: str):
        group = next((g for g in self.groups if g.group_id == group_id), None)
        if group is None:
            return
        spec = self._load_spectrum_file(file_path, group)
        if spec is not None:
            self._library_panel.add_spectrum(spec, group_id)
            self._add_toolbar_spectrum_button(spec)
            self._refresh_canvas()
            self.update_status(f"Loaded: {spec.display_name}")
        else:
            self.update_status(f"Failed to load: {os.path.basename(file_path)}")

    def _on_spectrum_visibility(self, spectrum_id: str, visible: bool):
        spec = self._find_spectrum(spectrum_id)
        if spec:
            spec.visible = visible
            btn = self._toolbar_buttons.get(spectrum_id)
            if btn and btn.isChecked() != visible:
                btn.blockSignals(True)
                btn.setChecked(visible)
                btn.blockSignals(False)
            self._refresh_canvas()

    def _on_spectrum_color(self, spectrum_id: str, color: str):
        spec = self._find_spectrum(spectrum_id)
        if spec:
            spec.color = color
            self._refresh_canvas()

    def _on_peak_list_attach_requested(self, spectrum_id: str):
        spec = self._find_spectrum(spectrum_id)
        if spec is None:
            return
        path, _ = QFileDialog.getOpenFileName(
            self, 'Load Peak List', '',
            'Peak list files (*.csv *.json);;All files (*)'
        )
        if not path:
            return
        count = self._load_peak_list_from_file(path, spectrum_id)
        if count > 0:
            self._refresh_canvas()
            self.update_status(f"Loaded {count} peaks for {spec.display_name}")
        else:
            self.update_status(f"No peaks loaded from {os.path.basename(path)}")

    def _on_contour_update_requested(self, levels: int, min_pct: float, increment: float):
        # Store as defaults for new spectra
        self._contour_levels = levels
        self._contour_min_pct = min_pct
        self._contour_increment = increment
        # Apply to all existing spectra
        for spec in self._all_spectra():
            spec.contour_levels = levels
            spec.contour_min_pct = min_pct
            spec.contour_increment = increment
        self._refresh_canvas()
        self.update_status(f'Contour: {levels} levels, min {min_pct:.1f}%, step ×{increment:.2f} — applied to all')

    def _on_active_spectrum_changed(self, spectrum_id: str):
        self._canvas.set_active_spectrum(spectrum_id)
        spec = self._find_spectrum(spectrum_id)
        if spec:
            self.update_status(f"Active: {spec.display_name}  |  Shift+click to add peaks")

    def _on_group_delete_requested(self, group_id: str):
        # Remove toolbar buttons for all spectra in this group
        group = next((g for g in self.groups if g.group_id == group_id), None)
        if group:
            for spec in group.spectra:
                self._remove_toolbar_button(spec.spectrum_id)
        self._remove_group(group_id)
        self._library_panel.remove_group(group_id)
        self._refresh_canvas()

    def _on_spectrum_delete_requested(self, spectrum_id: str):
        self._remove_toolbar_button(spectrum_id)
        self._remove_spectrum(spectrum_id)
        self._library_panel.remove_spectrum(spectrum_id)
        self._refresh_canvas()

    def _on_spectrum_moved(self, spectrum_id: str, target_group_id: str):
        """Move a spectrum to a different group (data model + tree)."""
        target_group = next((g for g in self.groups if g.group_id == target_group_id), None)
        if target_group is None:
            return

        spectrum = None
        source_group = None
        for g in self.groups:
            for s in g.spectra:
                if s.spectrum_id == spectrum_id:
                    spectrum = s
                    source_group = g
                    break
            if spectrum is not None:
                break

        if spectrum is None or source_group is target_group:
            return

        source_group.spectra.remove(spectrum)
        target_group.spectra.append(spectrum)

        self._library_panel.remove_spectrum(spectrum_id)
        self._library_panel.add_spectrum(spectrum, target_group_id)
        self._refresh_canvas()

    def _on_group_renamed(self, group_id: str, new_name: str):
        group = next((g for g in self.groups if g.group_id == group_id), None)
        if group:
            group.name = new_name

    def _on_spectrum_renamed(self, spectrum_id: str, new_name: str):
        spec = self._find_spectrum(spectrum_id)
        if spec:
            spec.display_name = new_name
            btn = self._toolbar_buttons.get(spectrum_id)
            if btn:
                btn.setText(new_name)

    def _on_canvas_peak_add_requested(self, spectrum_id: str, ppm_x: float, ppm_y: float):
        peak = self._add_peak(spectrum_id, ppm_x, ppm_y)
        if peak:
            self._refresh_canvas()
            self.update_status(f"Peak added at ({ppm_x:.3f}, {ppm_y:.3f})")

    def _on_canvas_peak_delete_requested(self, spectrum_id: str, peak_id: str):
        self._delete_peak(spectrum_id, peak_id)
        self._refresh_canvas()
        self.update_status("Peak deleted")

    # ------------------------------------------------------------------
    # Toolbar button management
    # ------------------------------------------------------------------

    def _add_toolbar_spectrum_button(self, spec: InspectorSpectrum):
        """Add a visibility toggle button for this spectrum to the top toolbar."""
        btn = QPushButton(spec.display_name)
        btn.setCheckable(True)
        btn.setChecked(spec.visible)
        btn.setStyleSheet(
            f"QPushButton {{ background-color: {spec.color}; color: white; "
            f"border: none; padding: 2px 8px; border-radius: 3px; font-weight: bold; }}"
            f"QPushButton:!checked {{ background-color: #888; color: #ccc; }}"
        )
        btn.toggled.connect(
            lambda checked, sid=spec.spectrum_id: self._on_spectrum_visibility(sid, checked)
        )
        action = self._top_toolbar.addWidget(btn)
        self._toolbar_buttons[spec.spectrum_id] = btn
        self._toolbar_actions[spec.spectrum_id] = action

    def _remove_toolbar_button(self, spectrum_id: str):
        """Remove the toolbar toggle button for a spectrum."""
        action = self._toolbar_actions.pop(spectrum_id, None)
        if action:
            self._top_toolbar.removeAction(action)
        self._toolbar_buttons.pop(spectrum_id, None)

    def _on_spectrum_contour_requested(self, spectrum_id: str):
        spec = self._find_spectrum(spectrum_id)
        if spec is None:
            return
        dlg = SpectrumContourDialog(spec, parent=self)
        if dlg.exec():
            self._refresh_canvas()

    def _on_group_properties_requested(self, group_id: str):
        group = next((g for g in self.groups if g.group_id == group_id), None)
        if group is None or not group.spectra:
            self.update_status('No spectra in this group')
            return
        dlg = GroupPropertiesDialog(group, parent=self)
        dlg.exec()
        self._refresh_canvas()

# ---------------------------------------------------------------------------
# Per-spectrum contour settings dialog
# ---------------------------------------------------------------------------

class SpectrumContourDialog(QWidget):
    """Small dialog to edit contour display settings for one spectrum."""

    def __init__(self, spectrum: InspectorSpectrum, parent=None):
        from PySide6.QtWidgets import QDialog, QDialogButtonBox
        # Use QDialog via dynamic super — keep class independent of import order
        self._dialog = QDialog(parent)
        self._dialog.setWindowTitle(f'Contour — {spectrum.display_name}')
        self._dialog.setMinimumWidth(280)
        self._spec = spectrum

        layout = QVBoxLayout(self._dialog)
        layout.setSpacing(6)

        def _row(label, widget):
            row = QHBoxLayout()
            row.addWidget(QLabel(label))
            row.addWidget(widget)
            return row

        self._levels_spin = QSpinBox()
        self._levels_spin.setRange(1, 30)
        self._levels_spin.setValue(spectrum.contour_levels)
        layout.addLayout(_row('Levels:', self._levels_spin))

        self._min_spin = QDoubleSpinBox()
        self._min_spin.setRange(0.01, 100.0)
        self._min_spin.setDecimals(2)
        self._min_spin.setValue(spectrum.contour_min_pct)
        layout.addLayout(_row('Min %:', self._min_spin))

        self._incr_spin = QDoubleSpinBox()
        self._incr_spin.setRange(1.01, 5.0)
        self._incr_spin.setDecimals(2)
        self._incr_spin.setValue(spectrum.contour_increment)
        layout.addLayout(_row('Step ×:', self._incr_spin))

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(self._apply)
        btns.accepted.connect(self._dialog.accept)
        btns.rejected.connect(self._dialog.reject)
        layout.addWidget(btns)

    def exec(self) -> bool:
        return self._dialog.exec() == 1

    def _apply(self):
        self._spec.contour_levels = self._levels_spin.value()
        self._spec.contour_min_pct = self._min_spin.value()
        self._spec.contour_increment = self._incr_spin.value()


# ---------------------------------------------------------------------------
# Group properties dialog
# ---------------------------------------------------------------------------

class GroupPropertiesDialog:
    """Dialog showing all spectra in a group with editable properties and propagation."""

    def __init__(self, group: InspectorGroup, parent=None):
        from PySide6.QtWidgets import (
            QDialog, QDialogButtonBox, QTableWidget, QTableWidgetItem,
            QComboBox, QHeaderView, QAbstractItemView
        )
        self._group = group
        self._dialog = QDialog(parent)
        self._dialog.setWindowTitle(f'Group Properties — {group.name}')
        self._dialog.setMinimumSize(680, 360)

        layout = QVBoxLayout(self._dialog)
        layout.setSpacing(6)

        # Table: one row per spectrum
        cols = ['Name', 'Visible', 'Color', 'Levels', 'Min %', 'Step ×']
        self._table = QTableWidget(len(group.spectra), len(cols), self._dialog)
        self._table.setHorizontalHeaderLabels(cols)
        self._table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self._table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self._table.verticalHeader().setVisible(False)

        self._row_widgets = []  # list of (checkbox, color_btn, levels_spin, min_spin, incr_spin)
        for row, spec in enumerate(group.spectra):
            # Name (read-only)
            name_item = QTableWidgetItem(spec.display_name)
            name_item.setFlags(name_item.flags() & ~Qt.ItemIsEditable)
            self._table.setItem(row, 0, name_item)

            # Visible checkbox
            vis_cb = QCheckBox()
            vis_cb.setChecked(spec.visible)
            cell = QWidget()
            cl = QHBoxLayout(cell)
            cl.addWidget(vis_cb)
            cl.setAlignment(Qt.AlignCenter)
            cl.setContentsMargins(0, 0, 0, 0)
            self._table.setCellWidget(row, 1, cell)

            # Color button
            color_btn = QPushButton()
            color_btn.setFixedSize(24, 24)
            color_btn.setStyleSheet(
                f'background-color: {spec.color}; border: 1px solid gray; border-radius: 2px;'
            )
            color_btn.clicked.connect(lambda _, b=color_btn, s=spec: self._pick_color(b, s))
            self._table.setCellWidget(row, 2, color_btn)

            # Levels
            lev = QSpinBox()
            lev.setRange(1, 30)
            lev.setValue(spec.contour_levels)
            self._table.setCellWidget(row, 3, lev)

            # Min %
            min_sp = QDoubleSpinBox()
            min_sp.setRange(0.01, 100.0)
            min_sp.setDecimals(2)
            min_sp.setValue(spec.contour_min_pct)
            self._table.setCellWidget(row, 4, min_sp)

            # Step
            incr_sp = QDoubleSpinBox()
            incr_sp.setRange(1.01, 5.0)
            incr_sp.setDecimals(2)
            incr_sp.setValue(spec.contour_increment)
            self._table.setCellWidget(row, 5, incr_sp)

            self._row_widgets.append((vis_cb, color_btn, lev, min_sp, incr_sp))

        layout.addWidget(self._table)

        # Propagation row
        prop_row = QHBoxLayout()
        prop_row.addWidget(QLabel('Propagate from:'))
        self._source_combo = QComboBox()
        for spec in group.spectra:
            self._source_combo.addItem(spec.display_name, spec.spectrum_id)
        prop_row.addWidget(self._source_combo)
        prop_row.addWidget(QLabel('→ all other spectra'))
        prop_btn = QPushButton('Propagate')
        prop_btn.clicked.connect(self._on_propagate)
        prop_row.addWidget(prop_btn)
        prop_row.addStretch()
        layout.addLayout(prop_row)

        # Dialog buttons
        btns = QDialogButtonBox(QDialogButtonBox.Apply | QDialogButtonBox.Close)
        btns.button(QDialogButtonBox.Apply).clicked.connect(self._apply)
        btns.button(QDialogButtonBox.Close).clicked.connect(self._dialog.accept)
        layout.addWidget(btns)

        # Stored colors (mutated by _pick_color)
        self._colors = [spec.color for spec in group.spectra]

    def exec(self):
        return self._dialog.exec()

    def _pick_color(self, btn: QPushButton, spec: InspectorSpectrum):
        idx = self._group.spectra.index(spec)
        color = QColorDialog.getColor(QColor(self._colors[idx]), self._dialog, 'Choose Color')
        if color.isValid():
            self._colors[idx] = color.name()
            btn.setStyleSheet(
                f'background-color: {color.name()}; border: 1px solid gray; border-radius: 2px;'
            )

    def _apply(self):
        for row, (spec, (vis_cb, color_btn, lev, min_sp, incr_sp)) in enumerate(
            zip(self._group.spectra, self._row_widgets)
        ):
            spec.visible = vis_cb.isChecked()
            spec.color = self._colors[row]
            spec.contour_levels = lev.value()
            spec.contour_min_pct = min_sp.value()
            spec.contour_increment = incr_sp.value()

    def _on_propagate(self):
        source_id = self._source_combo.currentData()
        source_idx = next(
            (i for i, s in enumerate(self._group.spectra) if s.spectrum_id == source_id), None
        )
        if source_idx is None:
            return
        src_widgets = self._row_widgets[source_idx]
        _, _, src_lev, src_min, src_incr = src_widgets
        src_color = self._colors[source_idx]

        for row, (vis_cb, color_btn, lev, min_sp, incr_sp) in enumerate(self._row_widgets):
            if row == source_idx:
                continue
            lev.setValue(src_lev.value())
            min_sp.setValue(src_min.value())
            incr_sp.setValue(src_incr.value())
            self._colors[row] = src_color
            color_btn.setStyleSheet(
                f'background-color: {src_color}; border: 1px solid gray; border-radius: 2px;'
            )
