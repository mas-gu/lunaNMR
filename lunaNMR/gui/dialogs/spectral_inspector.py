# ABOUTME: Standalone spectral inspector window for visual comparison of multiple NMR spectra
# ABOUTME: Supports overlay display, peak list management, peak add/delete, and coordinate shifting

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
from PySide6.QtWidgets import QSizePolicy, QFormLayout, QSlider
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor

from lunaNMR.gui.base.base_window import BaseWindow
from lunaNMR.gui.styles.design_system import (
    PRIMARY_TEXT, FRAME_BG_COLOR,
    SPACING_XS, SPACING_SM, SPACING_MD,
    FONT_SIZE_BODY, FONT_SIZE_SECTION_LABEL,
    PRIMARY_BUTTON_BG, PRIMARY_BUTTON_HOVER, PRIMARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BG, SECONDARY_BUTTON_HOVER, SECONDARY_BUTTON_TEXT, SECONDARY_BUTTON_BORDER,
    BUTTON_CORNER_RADIUS,
)

logger = logging.getLogger(__name__)

# Largest 1D-trace amplitude, in percent of the visible axis range. Weak cross-sections
# in a zoomed region need well over 1× the axis range to be readable.
TRACE_SCALE_MAX_PERCENT = 2000

# Default color palette — mirrors MultiSpectrumViewerDialog
DEFAULT_COLORS = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
]

# Shared widget styles — mirror the main window's Spectrum Display Controls / Shift sections
_GROUPBOX_STYLE = f"""
    QGroupBox {{
        font-size: {FONT_SIZE_SECTION_LABEL}pt;
        font-weight: bold;
        color: {PRIMARY_TEXT};
        border: 1px solid {SECONDARY_BUTTON_BORDER};
        border-radius: {BUTTON_CORNER_RADIUS}px;
        margin-top: {SPACING_SM}px;
        padding-top: {SPACING_MD}px;
    }}
    QGroupBox::title {{
        subcontrol-origin: margin;
        left: {SPACING_SM}px;
        padding: 0 {SPACING_XS}px;
    }}
"""

_SECONDARY_BTN_STYLE = f"""
    QPushButton {{
        background-color: {SECONDARY_BUTTON_BG};
        color: {SECONDARY_BUTTON_TEXT};
        border: 1px solid {SECONDARY_BUTTON_BORDER};
        border-radius: {BUTTON_CORNER_RADIUS}px;
        padding: 4px 6px;
        font-size: {FONT_SIZE_BODY}pt;
    }}
    QPushButton:hover {{ background-color: {SECONDARY_BUTTON_HOVER}; }}
    QPushButton:pressed {{ background-color: {SECONDARY_BUTTON_BORDER}; }}
"""

_PRIMARY_BTN_STYLE = f"""
    QPushButton {{
        background-color: {PRIMARY_BUTTON_BG};
        color: {PRIMARY_BUTTON_TEXT};
        border: none;
        border-radius: {BUTTON_CORNER_RADIUS}px;
        padding: 4px 6px;
        font-size: {FONT_SIZE_BODY}pt;
        font-weight: bold;
    }}
    QPushButton:hover {{ background-color: {PRIMARY_BUTTON_HOVER}; }}
    QPushButton:disabled {{ background-color: #cccccc; color: #666666; }}
"""

_SPIN_STYLE = f"""
    QSpinBox, QDoubleSpinBox {{
        background-color: {FRAME_BG_COLOR};
        color: {PRIMARY_TEXT};
        border: 1px solid {SECONDARY_BUTTON_BORDER};
        border-radius: {BUTTON_CORNER_RADIUS // 2}px;
        padding: {SPACING_XS}px;
        font-size: {FONT_SIZE_BODY}pt;
    }}
"""


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

    # Net coordinate shift applied to this peak (for Shift Peak List undo).
    # Tracked per-peak so a Reset only moves peaks that were actually shifted.
    applied_shift_h: float = 0.0
    applied_shift_n: float = 0.0

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
                # Treat NaN (blank cell from pandas) as missing, not a value
                if v is not None and not (isinstance(v, float) and v != v):
                    return v
            return None

        x = _first_not_none('ppm_x', 'center_x', 'peak_x', 'Position_X', 'Position_x')
        y = _first_not_none('ppm_y', 'center_y', 'peak_y', 'Position_Y', 'Position_y')
        assignment = _first_not_none('assignment', 'Assignment')
        quality = (_first_not_none('quality', 'Quality', 'fitting_quality') or 'Unknown')
        return cls(
            ppm_x=float(x) if x is not None else 0.0,
            ppm_y=float(y) if y is not None else 0.0,
            assignment=str(assignment) if assignment is not None else '',
            quality=str(quality),
        )

    def to_state(self) -> dict:
        """Full round-trippable state for project persistence."""
        return {
            'peak_id': self.peak_id,
            'ppm_x': self.ppm_x, 'ppm_y': self.ppm_y,
            'assignment': self.assignment, 'quality': self.quality,
            'applied_shift_h': self.applied_shift_h,
            'applied_shift_n': self.applied_shift_n,
        }

    @classmethod
    def from_state(cls, d: dict) -> InspectorPeak:
        """Reconstruct a peak from :meth:`to_state` output."""
        return cls(
            peak_id=d.get('peak_id') or str(uuid.uuid4()),
            ppm_x=float(d.get('ppm_x', 0.0)),
            ppm_y=float(d.get('ppm_y', 0.0)),
            assignment=str(d.get('assignment', '')),
            quality=str(d.get('quality', 'Unknown')),
            applied_shift_h=float(d.get('applied_shift_h', 0.0)),
            applied_shift_n=float(d.get('applied_shift_n', 0.0)),
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
    # contour_min_level is a fraction of the spectrum's max intensity (0.0001–1.0)
    contour_levels: int = 10
    contour_min_level: float = 0.05
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
    contour_update_requested = Signal(int, float, float)  # levels, min_level, increment
    min_scale_requested = Signal(float)         # scale each spectrum's own min level
    propagate_requested = Signal()
    horizontal_trace_toggled = Signal(bool)
    vertical_trace_toggled = Signal(bool)
    trace_scale_changed = Signal(float)         # trace amplitude fraction (0–1)
    reset_zoom_requested = Signal()
    peak_markers_toggled = Signal(bool)         # show markers
    shift_apply_requested = Signal(float, float)  # h_shift, n_shift
    shift_reset_requested = Signal()
    active_spectrum_changed = Signal(str)       # spectrum_id
    group_delete_requested = Signal(str)        # group_id
    spectrum_delete_requested = Signal(str)     # spectrum_id
    group_renamed = Signal(str, str)            # group_id, new_name
    spectrum_renamed = Signal(str, str)         # spectrum_id, new_name
    spectrum_reorder_requested = Signal(str, int)  # spectrum_id, delta (-1 up / +1 down)
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

        # Peak marker visibility toggle
        self._peaks_checkbox = QCheckBox('Show peak markers')
        self._peaks_checkbox.setChecked(True)
        self._peaks_checkbox.setToolTip('Turn peak-list visualization on/off for all spectra')
        self._peaks_checkbox.toggled.connect(self.peak_markers_toggled)
        layout.addWidget(self._peaks_checkbox)

        layout.addWidget(self._build_display_controls())
        layout.addWidget(self._build_traces_section())
        layout.addWidget(self._build_shift_section())

        # Internal map: group_id → QTreeWidgetItem
        self._group_items: dict = {}

    def _build_traces_section(self) -> QGroupBox:
        """1D Traces: toggle horizontal/vertical cross-sections + amplitude scale slider."""
        group = QGroupBox('1D Traces')
        group.setStyleSheet(_GROUPBOX_STYLE)
        layout = QVBoxLayout(group)
        layout.setContentsMargins(SPACING_SM, SPACING_MD, SPACING_SM, SPACING_SM)
        layout.setSpacing(SPACING_XS)

        cb_row = QHBoxLayout()
        cb_row.setSpacing(SPACING_XS)
        self._h_trace_cb = QCheckBox('Horizontal')
        self._h_trace_cb.setToolTip('Show a horizontal cross-section following the cursor')
        self._h_trace_cb.toggled.connect(self.horizontal_trace_toggled)
        cb_row.addWidget(self._h_trace_cb)
        self._v_trace_cb = QCheckBox('Vertical')
        self._v_trace_cb.setToolTip('Show a vertical cross-section following the cursor')
        self._v_trace_cb.toggled.connect(self.vertical_trace_toggled)
        cb_row.addWidget(self._v_trace_cb)
        layout.addLayout(cb_row)

        scale_row = QHBoxLayout()
        scale_row.setSpacing(SPACING_XS)
        scale_label = QLabel('Scale:')
        scale_label.setStyleSheet(f'font-size: {FONT_SIZE_BODY}pt; color: {PRIMARY_TEXT};')
        scale_row.addWidget(scale_label)
        self._trace_scale_slider = QSlider(Qt.Horizontal)
        self._trace_scale_slider.setRange(1, TRACE_SCALE_MAX_PERCENT)  # → amplitude 0.01–20×
        self._trace_scale_slider.setValue(15)
        self._trace_scale_slider.setToolTip(
            'Trace height as a multiple of the visible axis range (1.00 = full range)')
        self._trace_scale_slider.valueChanged.connect(self._on_trace_scale)
        scale_row.addWidget(self._trace_scale_slider, stretch=1)
        self._trace_scale_spin = QDoubleSpinBox()
        self._trace_scale_spin.setRange(0.01, TRACE_SCALE_MAX_PERCENT / 100.0)
        self._trace_scale_spin.setDecimals(2)
        self._trace_scale_spin.setSingleStep(0.25)
        self._trace_scale_spin.setValue(0.15)
        self._trace_scale_spin.setFixedWidth(70)
        self._trace_scale_spin.setStyleSheet(_SPIN_STYLE)
        self._trace_scale_spin.setToolTip('Type an exact trace amplitude (axis-range multiples)')
        self._trace_scale_spin.valueChanged.connect(self._on_trace_scale_typed)
        scale_row.addWidget(self._trace_scale_spin)
        layout.addLayout(scale_row)
        return group

    def _build_display_controls(self) -> QGroupBox:
        """Spectrum Display Controls: collapsible contour params + Reset Zoom.

        Contour changes apply live to all spectra (no explicit Apply button).
        """
        group = QGroupBox('Spectrum Display Controls')
        group.setStyleSheet(_GROUPBOX_STYLE)
        layout = QVBoxLayout(group)
        layout.setContentsMargins(SPACING_SM, SPACING_MD, SPACING_SM, SPACING_SM)
        layout.setSpacing(SPACING_SM)

        # Button row: Contour Settings toggle + Reset Zoom
        button_row = QHBoxLayout()
        button_row.setSpacing(SPACING_XS)

        self._contour_toggle_btn = QPushButton('Contour Settings ▼')
        self._contour_toggle_btn.setToolTip('Show/hide contour display settings')
        self._contour_toggle_btn.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        self._contour_toggle_btn.setStyleSheet(_SECONDARY_BTN_STYLE)
        self._contour_toggle_btn.clicked.connect(self._toggle_contour_settings)
        button_row.addWidget(self._contour_toggle_btn)

        reset_zoom_btn = QPushButton('Reset Zoom')
        reset_zoom_btn.setToolTip('Reset view to full spectrum')
        reset_zoom_btn.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        reset_zoom_btn.setStyleSheet(_SECONDARY_BTN_STYLE)
        reset_zoom_btn.clicked.connect(self.reset_zoom_requested)
        button_row.addWidget(reset_zoom_btn)
        layout.addLayout(button_row)

        # Collapsible contour parameters
        self._contour_params_widget = QWidget()
        params_layout = QFormLayout(self._contour_params_widget)
        params_layout.setContentsMargins(0, SPACING_SM, 0, 0)
        params_layout.setSpacing(SPACING_XS)
        params_layout.setLabelAlignment(Qt.AlignRight)

        self._levels_spin = QSpinBox()
        self._levels_spin.setRange(5, 50)
        self._levels_spin.setValue(10)
        self._levels_spin.setStyleSheet(_SPIN_STYLE)
        self._levels_spin.valueChanged.connect(self._on_contour_changed)
        params_layout.addRow('Levels:', self._levels_spin)

        # Min Level row: spinbox + ×2 / ÷2 quick-scale buttons
        min_row = QWidget()
        min_layout = QHBoxLayout(min_row)
        min_layout.setContentsMargins(0, 0, 0, 0)
        min_layout.setSpacing(SPACING_XS)
        self._min_spin = QDoubleSpinBox()
        self._min_spin.setRange(0.0001, 1.0)
        self._min_spin.setDecimals(5)
        self._min_spin.setSingleStep(0.002)
        self._min_spin.setValue(0.05)
        self._min_spin.setStyleSheet(_SPIN_STYLE)
        self._min_spin.valueChanged.connect(self._on_contour_changed)
        min_layout.addWidget(self._min_spin, stretch=1)
        x2_btn = QPushButton('×2')
        x2_btn.setToolTip("Double each spectrum's own min level")
        x2_btn.setFixedWidth(34)
        x2_btn.setStyleSheet(_SECONDARY_BTN_STYLE)
        x2_btn.clicked.connect(self._on_min_times2)
        min_layout.addWidget(x2_btn)
        div2_btn = QPushButton('÷2')
        div2_btn.setToolTip("Halve each spectrum's own min level")
        div2_btn.setFixedWidth(34)
        div2_btn.setStyleSheet(_SECONDARY_BTN_STYLE)
        div2_btn.clicked.connect(self._on_min_div2)
        min_layout.addWidget(div2_btn)
        params_layout.addRow('Min Level:', min_row)

        self._incr_spin = QDoubleSpinBox()
        # Increment must stay > 1 or contour levels stop increasing (empty/broken render)
        self._incr_spin.setRange(1.01, 10.0)
        self._incr_spin.setDecimals(2)
        self._incr_spin.setSingleStep(0.1)
        self._incr_spin.setValue(1.3)
        self._incr_spin.setStyleSheet(_SPIN_STYLE)
        self._incr_spin.valueChanged.connect(self._on_contour_changed)
        params_layout.addRow('Increment:', self._incr_spin)

        # Propagate one spectrum's contour settings to chosen others
        propagate_btn = QPushButton('Propagate…')
        propagate_btn.setToolTip('Copy Levels / Min Level / Increment from one spectrum to others')
        propagate_btn.setStyleSheet(_SECONDARY_BTN_STYLE)
        propagate_btn.clicked.connect(self.propagate_requested)
        params_layout.addRow(propagate_btn)

        self._contour_params_widget.setVisible(False)
        layout.addWidget(self._contour_params_widget)
        return group

    def _build_shift_section(self) -> QGroupBox:
        """Shift Peak List: offset the active spectrum's peaks in 1H and 15N/13C."""
        group = QGroupBox('Shift Peak List')
        group.setStyleSheet(_GROUPBOX_STYLE)
        layout = QVBoxLayout(group)
        layout.setContentsMargins(SPACING_SM, SPACING_MD, SPACING_SM, SPACING_SM)
        layout.setSpacing(SPACING_XS)

        # Row 1: offset spinboxes
        spin_row = QHBoxLayout()
        spin_row.setSpacing(SPACING_XS)
        h_label = QLabel('1H:')
        h_label.setStyleSheet(f'font-size: {FONT_SIZE_BODY}pt; color: {PRIMARY_TEXT};')
        spin_row.addWidget(h_label)
        self._shift_h_spin = QDoubleSpinBox()
        self._shift_h_spin.setRange(-10.0, 10.0)
        self._shift_h_spin.setSingleStep(0.001)
        self._shift_h_spin.setDecimals(3)
        self._shift_h_spin.setSuffix(' ppm')
        self._shift_h_spin.setStyleSheet(_SPIN_STYLE)
        spin_row.addWidget(self._shift_h_spin)

        n_label = QLabel('15N/13C:')
        n_label.setStyleSheet(f'font-size: {FONT_SIZE_BODY}pt; color: {PRIMARY_TEXT};')
        spin_row.addWidget(n_label)
        self._shift_n_spin = QDoubleSpinBox()
        self._shift_n_spin.setRange(-50.0, 50.0)
        self._shift_n_spin.setSingleStep(0.01)
        self._shift_n_spin.setDecimals(2)
        self._shift_n_spin.setSuffix(' ppm')
        self._shift_n_spin.setStyleSheet(_SPIN_STYLE)
        spin_row.addWidget(self._shift_n_spin)
        layout.addLayout(spin_row)

        # Row 2: Apply / Reset
        btn_row = QHBoxLayout()
        btn_row.setSpacing(SPACING_XS)
        apply_btn = QPushButton('Apply')
        apply_btn.setToolTip('Shift the active spectrum\'s peaks by these offsets')
        apply_btn.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        apply_btn.setStyleSheet(_PRIMARY_BTN_STYLE)
        apply_btn.clicked.connect(self._on_shift_apply)
        btn_row.addWidget(apply_btn)

        reset_btn = QPushButton('Reset')
        reset_btn.setToolTip('Undo the applied shift and zero the offsets')
        reset_btn.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        reset_btn.setStyleSheet(_SECONDARY_BTN_STYLE)
        reset_btn.clicked.connect(self._on_shift_reset)
        btn_row.addWidget(reset_btn)
        btn_row.addStretch()
        layout.addLayout(btn_row)
        return group

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

    def rebuild_group_spectra(self, group: InspectorGroup):
        """Rebuild a group's spectrum rows to match the model order (after a reorder)."""
        group_item = self._group_items.get(group.group_id)
        if group_item is None:
            return
        while group_item.childCount():
            group_item.removeChild(group_item.child(0))
        for spec in group.spectra:
            self.add_spectrum(spec, group.group_id)

    def select_spectrum(self, spectrum_id: str):
        """Select the tree row for a spectrum (keeps focus on a just-moved item)."""
        root = self._tree.invisibleRootItem()
        for i in range(root.childCount()):
            group_item = root.child(i)
            for j in range(group_item.childCount()):
                child = group_item.child(j)
                if child.data(0, _ROLE_ID) == spectrum_id:
                    self._tree.setCurrentItem(child)
                    return

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
            up_act = menu.addAction('Move Up (behind)')
            down_act = menu.addAction('Move Down (in front)')
            menu.addSeparator()
            rename_act = menu.addAction('Rename')
            delete_act = menu.addAction('Remove Spectrum')
            action = menu.exec(self._tree.viewport().mapToGlobal(pos))
            if action == contour_act:
                self.spectrum_contour_requested.emit(spectrum_id)
            elif action == up_act:
                self.spectrum_reorder_requested.emit(spectrum_id, -1)
            elif action == down_act:
                self.spectrum_reorder_requested.emit(spectrum_id, +1)
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

    def _on_contour_changed(self, _value=None):
        """Emit contour params live whenever a spinbox changes."""
        self.contour_update_requested.emit(
            self._levels_spin.value(),
            self._min_spin.value(),
            self._incr_spin.value(),
        )

    def _on_trace_scale(self, value: int):
        """Slider step → trace amplitude in axis-range multiples (1 step = 0.01)."""
        scale = value / 100.0
        self._set_spin_silently(self._trace_scale_spin, scale)
        self.trace_scale_changed.emit(scale)

    def _on_trace_scale_typed(self, scale: float):
        """Typed amplitude → move the slider (silently) and emit."""
        slider = self._trace_scale_slider
        slider.blockSignals(True)
        slider.setValue(int(round(scale * 100)))
        slider.blockSignals(False)
        self.trace_scale_changed.emit(scale)

    def set_trace_scale_silently(self, scale: float):
        """Show ``scale`` on both the slider and the spin box without emitting."""
        slider = self._trace_scale_slider
        slider.blockSignals(True)
        slider.setValue(int(round(scale * 100)))
        slider.blockSignals(False)
        self._set_spin_silently(self._trace_scale_spin, scale)

    @staticmethod
    def _set_spin_silently(spin, value: float):
        spin.blockSignals(True)
        spin.setValue(value)
        spin.blockSignals(False)

    def _on_min_times2(self):
        """Scale each spectrum's own min level up by 2× (preserves per-spectrum differences)."""
        self.min_scale_requested.emit(2.0)

    def _on_min_div2(self):
        """Scale each spectrum's own min level down by 2×."""
        self.min_scale_requested.emit(0.5)

    def _toggle_contour_settings(self):
        """Show/hide the collapsible contour parameters."""
        visible = self._contour_params_widget.isVisible()
        self._contour_params_widget.setVisible(not visible)
        arrow = '▼' if visible else '▶'
        self._contour_toggle_btn.setText(f'Contour Settings {arrow}')

    def _on_shift_apply(self):
        self.shift_apply_requested.emit(
            self._shift_h_spin.value(), self._shift_n_spin.value()
        )

    def _on_shift_reset(self):
        self.shift_reset_requested.emit()
        self._shift_h_spin.setValue(0.0)
        self._shift_n_spin.setValue(0.0)


# ---------------------------------------------------------------------------
# Peak list file parsing
# ---------------------------------------------------------------------------

def _parse_peak_list_file(file_path: str) -> List[dict]:
    """Parse a peak list file and return a list of dicts.

    JSON files must contain a list of objects. Everything else (.txt, .csv,
    .peaks) is read as the comma-separated, header-row format the main window
    accepts (``Assignment, Position_X, Position_Y[, Height]``), using the same
    pandas parsing so acceptance matches exactly.
    """
    ext = os.path.splitext(file_path)[1].lower()
    if ext == '.json':
        import json
        with open(file_path, 'r') as fh:
            data = json.load(fh)
        return data if isinstance(data, list) else []
    else:
        import pandas as pd
        df = pd.read_csv(file_path, skipinitialspace=True)
        df.columns = df.columns.str.strip()
        return df.to_dict('records')


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
        self._rubberband = None  # matplotlib Rectangle artist for the zoom box
        self._peaks_visible = True  # global peak-marker visibility toggle
        self._has_plotted = False   # True once real data has been drawn (zoom preservation)

        # Navigation handler (pan, zoom, keyboard, peak callbacks)
        from lunaNMR.gui.components.nmr_navigation_handler import NMRNavigationHandler
        self._nav_handler = NMRNavigationHandler()
        self._nav_handler.attach(self._plot)
        self._nav_handler.on_peak_select = self._on_nav_peak_select
        self._nav_handler.on_peak_edit = self._on_nav_peak_edit
        self._nav_handler.on_delete = self._on_nav_delete
        self._nav_handler.on_escape = self._on_nav_escape
        self._nav_handler.on_reset_zoom = self._reset_zoom
        # Middle-drag draws a zoom box (independent x/y → reshape the aspect);
        # middle-click (no drag) still selects a peak.
        self._nav_handler.on_rect_drag = self._on_zoom_rubberband
        self._nav_handler.on_area_select = self._on_zoom_box

        # 1D trace (cross-section) state — hover-following, blitted for smoothness
        self._h_trace_on = False
        self._v_trace_on = False
        self._trace_scale = 0.15    # trace amplitude as a fraction of the visible axis range
        self._trace_bg = None       # cached clean background for blitting (no trace lines)
        self._trace_lines = []      # reusable animated Line2D pool
        self._last_cursor = None    # (cx, cy) so a scale change can redraw without mouse motion
        self._trace_max_cache = {}  # spectrum_id → max|intensity| (avoids O(N) per mouse move)
        self._plot.canvas.mpl_connect('motion_notify_event', self._on_cursor_move)
        self._plot.canvas.mpl_connect('draw_event', self._on_trace_draw_event)

        self._show_empty_state()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def set_active_spectrum(self, spectrum_id: Optional[str]):
        """Set which spectrum receives new peaks from Shift+click."""
        self._active_spectrum_id = spectrum_id

    def set_peaks_visible(self, visible: bool):
        """Toggle whether peak markers are drawn, then redraw."""
        self._peaks_visible = visible
        self.update_plot(self._current_spectra)

    def reset_zoom(self):
        """Restore the full-spectrum view on the next redraw."""
        self._has_plotted = False
        self.update_plot(self._current_spectra)

    # ------------------------------------------------------------------
    # Zoom box (middle-drag) — reshape the view to any x/y aspect
    # ------------------------------------------------------------------

    def _clear_rubberband(self):
        if self._rubberband is not None:
            try:
                self._rubberband.remove()
            except Exception:
                pass
            self._rubberband = None

    def _on_zoom_rubberband(self, x0, y0, x1, y1):
        """Draw/update the dashed zoom-box rectangle while the middle button drags."""
        if None in (x0, y0, x1, y1):
            return
        from matplotlib.patches import Rectangle
        self._clear_rubberband()
        self._rubberband = Rectangle(
            (min(x0, x1), min(y0, y1)), abs(x1 - x0), abs(y1 - y0),
            fill=False, edgecolor='dimgray', linewidth=1.0, linestyle='--', zorder=20,
        )
        self._plot.axes.add_patch(self._rubberband)
        self._plot.canvas.draw_idle()

    def _on_zoom_box(self, x0, y0, x1, y1):
        """Middle-drag release: zoom the view to the drawn box (independent x and y).

        This sets the ¹H and ¹⁵N/¹³C ranges independently, so the displayed x/y
        aspect ratio follows the box. Reset Zoom restores the full spectrum.
        """
        self._clear_rubberband()
        ax = self._plot.axes
        if None in (x0, y0, x1, y1) or x0 == x1 or y0 == y1:
            self._plot.canvas.draw_idle()
            return
        # Keep the NMR-inverted orientation (high ppm first) regardless of drag direction
        ax.set_xlim(max(x0, x1), min(x0, x1))
        ax.set_ylim(max(y0, y1), min(y0, y1))
        self._has_plotted = True
        self._plot.canvas.draw_idle()

    # ------------------------------------------------------------------
    # 1D traces (cross-sections) — follow the cursor, drawn in each spectrum's color
    # ------------------------------------------------------------------

    def set_traces(self, horizontal=None, vertical=None):
        """Turn the horizontal and/or vertical cross-section traces on or off."""
        if horizontal is not None:
            self._h_trace_on = horizontal
        if vertical is not None:
            self._v_trace_on = vertical
        if not (self._h_trace_on or self._v_trace_on):
            self._clear_traces()
        elif self._last_cursor is not None:
            self._draw_traces_at(*self._last_cursor)
        else:
            self._plot.canvas.draw_idle()

    def set_trace_scale(self, fraction: float):
        """Set trace amplitude as a fraction of the visible axis range and redraw."""
        self._trace_scale = fraction
        if (self._h_trace_on or self._v_trace_on) and self._last_cursor is not None:
            self._draw_traces_at(*self._last_cursor)

    @staticmethod
    def _horizontal_trace(data, ppm_x, ppm_y, cy):
        """Return (ppm_x, intensity_row) for the data row nearest to ppm_y == cy."""
        import numpy as np
        iy = int(np.argmin(np.abs(np.asarray(ppm_y) - cy)))
        return np.asarray(ppm_x), np.asarray(data)[iy, :]

    @staticmethod
    def _vertical_trace(data, ppm_x, ppm_y, cx):
        """Return (ppm_y, intensity_col) for the data column nearest to ppm_x == cx."""
        import numpy as np
        ix = int(np.argmin(np.abs(np.asarray(ppm_x) - cx)))
        return np.asarray(ppm_y), np.asarray(data)[:, ix]

    @staticmethod
    def _offset_intensity(intensity, max_int, scale, axis_range, baseline):
        """Map raw intensity to plotted coordinates: baseline + scaled normalized amplitude."""
        import numpy as np
        denom = max_int if max_int else 1.0
        return baseline + (scale * axis_range) * (np.asarray(intensity) / denom)

    def _on_trace_draw_event(self, _event):
        """Capture a clean background after every full draw (animated traces excluded)."""
        try:
            self._trace_bg = self._plot.canvas.copy_from_bbox(self._plot.axes.bbox)
        except Exception:
            self._trace_bg = None

    def _on_cursor_move(self, event):
        if not (self._h_trace_on or self._v_trace_on):
            return
        if event.inaxes != self._plot.axes or event.xdata is None or event.ydata is None:
            return
        self._last_cursor = (event.xdata, event.ydata)
        self._draw_traces_at(event.xdata, event.ydata)

    def _trace_line(self, i):
        """Reuse (or lazily create) the i-th animated trace line."""
        while i >= len(self._trace_lines):
            (line,) = self._plot.axes.plot([], [], linewidth=0.8, animated=True, zorder=15)
            self._trace_lines.append(line)
        return self._trace_lines[i]

    def _draw_traces_at(self, cx, cy):
        """Blit the cross-section traces for every visible spectrum at cursor (cx, cy)."""
        import numpy as np
        canvas = self._plot.canvas
        ax = self._plot.axes
        if self._trace_bg is None:
            canvas.draw()          # triggers draw_event → captures a clean background
            if self._trace_bg is None:
                return
        canvas.restore_region(self._trace_bg)

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        # Signed spans so positive intensity bumps toward the top (horizontal) / right
        # (vertical) regardless of the axes being inverted (NMR high-ppm-first convention).
        x_range = xlim[1] - xlim[0]
        y_range = ylim[1] - ylim[0]

        idx = 0
        for spec in self._current_spectra:
            if not (spec.visible and spec.loaded and spec.integrator is not None):
                continue
            integrator = spec.integrator
            data = getattr(integrator, 'nmr_data', None)
            px = getattr(integrator, 'ppm_x_axis', None)
            py = getattr(integrator, 'ppm_y_axis', None)
            if data is None or px is None or py is None:
                continue
            max_int = self._trace_max_cache.get(spec.spectrum_id)
            if not max_int:
                max_int = float(np.max(np.abs(data))) or 1.0
                self._trace_max_cache[spec.spectrum_id] = max_int

            if self._h_trace_on:
                xarr, inten = self._horizontal_trace(data, px, py, cy)
                yplot = self._offset_intensity(inten, max_int, self._trace_scale, y_range, cy)
                line = self._trace_line(idx); idx += 1
                line.set_data(xarr, yplot)
                line.set_color(spec.color); line.set_visible(True)
                ax.draw_artist(line)

            if self._v_trace_on:
                yarr, inten = self._vertical_trace(data, px, py, cx)
                xplot = self._offset_intensity(inten, max_int, self._trace_scale, x_range, cx)
                line = self._trace_line(idx); idx += 1
                line.set_data(xplot, yarr)
                line.set_color(spec.color); line.set_visible(True)
                ax.draw_artist(line)

        for j in range(idx, len(self._trace_lines)):
            self._trace_lines[j].set_visible(False)
        canvas.blit(ax.bbox)

    def _clear_traces(self):
        for line in self._trace_lines:
            line.set_visible(False)
        self._plot.canvas.draw_idle()

    # ------------------------------------------------------------------
    # PDF export
    # ------------------------------------------------------------------

    def export_pdf(self, file_path: str):
        """Save the current view to a vector PDF with editable (TrueType) text.

        Sets ``pdf.fonttype = 42`` so glyphs are embedded as TrueType and stay
        selectable/editable in Illustrator — the matplotlib default (Type 3)
        writes each glyph as an un-editable path. Contours and markers are
        already vector geometry. The transient rubber-band and selection ring
        are excluded from the output.
        """
        import matplotlib as mpl
        self._clear_rubberband()
        self._clear_selection_highlight()
        saved = {k: mpl.rcParams[k] for k in ('pdf.fonttype', 'ps.fonttype')}
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42
        try:
            self._plot.figure.savefig(
                file_path, format='pdf', bbox_inches='tight',
                facecolor='white', edgecolor='none',
            )
        finally:
            mpl.rcParams.update(saved)
            # restore the on-screen selection ring even if the save failed
            self._draw_selection_highlight()
            self._plot.canvas.draw_idle()

    def update_plot(self, spectra: List[InspectorSpectrum]):
        """Redraw all visible loaded spectra using each spectrum's own contour settings."""
        import numpy as np
        import matplotlib.pyplot as plt

        self._current_spectra = spectra
        ax = self._plot.axes

        # Preserve the user's zoom across redraws once real data has been drawn once
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        has_zoom = self._has_plotted

        self._plot.clear()
        self._rubberband = None   # removed by axes clear; drop the stale reference
        self._trace_lines = []    # animated trace artists were removed by the clear too
        self._trace_bg = None     # background is stale until the next draw
        self._trace_max_cache = {}  # recomputed below for the spectra actually drawn

        visible = [s for s in spectra if s.visible and s.loaded and s.integrator is not None]

        if not visible:
            # Empty state leaves the axes at matplotlib's default (0,1) range;
            # forget the prior view so the next real draw re-fits to full data
            # instead of preserving that (0,1) as if it were a user zoom.
            self._has_plotted = False
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
                self._trace_max_cache[spec.spectrum_id] = float(max_int)

                levels = self._calculate_contour_levels(
                    max_int, spec.contour_levels, spec.contour_min_level, spec.contour_increment
                )
                if not levels:
                    continue

                ax.contour(ppm_x, ppm_y, data,
                           levels=levels, colors=spec.color,
                           linewidths=1.0, alpha=0.8)

                if self._peaks_visible and spec.peaks:
                    xs = [p.ppm_x for p in spec.peaks]
                    ys = [p.ppm_y for p in spec.peaks]
                    ax.scatter(xs, ys, marker='+', s=40, c=spec.color,
                               linewidths=1.2, zorder=5)
                    # Annotate each peak with its assignment (falling back to its number)
                    for i, peak in enumerate(spec.peaks):
                        label = peak.assignment or str(i + 1)
                        ax.annotate(
                            label, (peak.ppm_x, peak.ppm_y),
                            xytext=(5, 5), textcoords='offset points',
                            fontsize='small', color=spec.color, fontweight='bold',
                            zorder=6,
                            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                      alpha=0.8, edgecolor=spec.color),
                        )

                legend_handles.append(
                    plt.Line2D([0], [0], color=spec.color, linewidth=2)
                )
                legend_labels.append(spec.display_name)

            except Exception as e:
                logger.warning(f"Failed to draw {spec.display_name}: {e}")
                continue

        if not legend_handles:
            # No spectrum drew successfully (e.g. every contour failed) — show the
            # empty hint and DON'T lock this blank (0,1) view as a preserved zoom.
            self._has_plotted = False
            self._plot.clear()
            self._show_empty_state()
            self._plot.refresh()
            return

        ax.set_xlabel('¹H (ppm)', fontsize=10)
        ax.set_ylabel('¹⁵N / ¹³C (ppm)', fontsize=10)

        if has_zoom:
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
        else:
            ax.invert_xaxis()
            ax.invert_yaxis()
        self._has_plotted = True

        if legend_handles:
            ax.legend(legend_handles, legend_labels, loc='best', fontsize=8)

        # Restore selection highlight after axes are cleared by redraw
        self._draw_selection_highlight()
        self._plot.refresh()

    @staticmethod
    def _calculate_contour_levels(max_intensity: float, num_levels: int,
                                   min_level_fraction: float, increment: float) -> list:
        """Return multiplicatively-spaced contour levels up to max_intensity.

        min_level_fraction is the first level as a fraction of max_intensity
        (e.g. 0.05 → first contour at 5% of the peak).
        """
        min_level = min_level_fraction * max_intensity
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
    color. Supports local peak list management (load, add, delete) and
    coordinate shifting of a spectrum's peaks. All edits are local — nothing
    is propagated back to the main window.
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
        self._active_spectrum_id: Optional[str] = None

        # Contour parameters (match SpectrumLibraryPanel defaults)
        # _contour_min_level is a fraction of max intensity (0.0001–1.0)
        self._contour_levels = 10
        self._contour_min_level = 0.05
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
        self._library_panel.min_scale_requested.connect(self._on_min_scale_requested)
        self._library_panel.propagate_requested.connect(self._on_propagate_requested)
        self._library_panel.horizontal_trace_toggled.connect(self._on_horizontal_trace_toggled)
        self._library_panel.vertical_trace_toggled.connect(self._on_vertical_trace_toggled)
        self._library_panel.trace_scale_changed.connect(self._on_trace_scale_changed)
        self._library_panel.reset_zoom_requested.connect(self._on_reset_zoom_requested)
        self._library_panel.peak_markers_toggled.connect(self._on_peak_markers_toggled)
        self._library_panel.shift_apply_requested.connect(self._on_shift_apply_requested)
        self._library_panel.shift_reset_requested.connect(self._on_shift_reset_requested)
        self._library_panel.active_spectrum_changed.connect(self._on_active_spectrum_changed)
        self._library_panel.group_delete_requested.connect(self._on_group_delete_requested)
        self._library_panel.spectrum_delete_requested.connect(self._on_spectrum_delete_requested)
        self._library_panel.group_renamed.connect(self._on_group_renamed)
        self._library_panel.spectrum_renamed.connect(self._on_spectrum_renamed)
        self._library_panel.spectrum_reorder_requested.connect(self._on_spectrum_reorder_requested)
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
        """Create the top toolbar: PDF export + per-spectrum visibility toggles."""
        tb = QToolBar('Spectra', self)
        tb.setMovable(False)
        tb.setFloatable(False)

        export_btn = QPushButton('Export PDF…')
        export_btn.setToolTip('Export the current view to a vector PDF '
                              '(editable text, ready for Illustrator)')
        export_btn.setStyleSheet(_SECONDARY_BTN_STYLE)
        export_btn.clicked.connect(self._on_export_pdf)
        tb.addWidget(export_btn)
        tb.addSeparator()

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

    # ------------------------------------------------------------------
    # State serialization (project persistence + session close/reopen)
    # ------------------------------------------------------------------

    STATE_SCHEMA = 1

    def get_state(self) -> dict:
        """Serialize all inspector content to a JSON-safe dict (no integrators)."""
        return {
            'schema': self.STATE_SCHEMA,
            'color_index': self._color_index,
            'active_spectrum_id': self._active_spectrum_id,
            'contour_defaults': {
                'levels': self._contour_levels,
                'min_level': self._contour_min_level,
                'increment': self._contour_increment,
            },
            'peaks_visible': self._canvas._peaks_visible,
            'traces': {
                'horizontal': self._canvas._h_trace_on,
                'vertical': self._canvas._v_trace_on,
                'scale': self._canvas._trace_scale,
            },
            'groups': [
                {
                    'group_id': g.group_id,
                    'name': g.name,
                    'spectra': [self._spectrum_to_state(s) for s in g.spectra],
                }
                for g in self.groups
            ],
        }

    @staticmethod
    def _spectrum_to_state(s: InspectorSpectrum) -> dict:
        return {
            'spectrum_id': s.spectrum_id,
            'display_name': s.display_name,
            'file_path': s.file_path,
            'color': s.color,
            'visible': s.visible,
            'contour_levels': s.contour_levels,
            'contour_min_level': s.contour_min_level,
            'contour_increment': s.contour_increment,
            'peaks': [p.to_state() for p in s.peaks],
        }

    def load_state(self, state: dict):
        """Replace all inspector content from a :meth:`get_state` dict.

        Spectra are reloaded from their stored file paths; a spectrum whose file
        is missing is kept in the tree but stays unloaded (degrades gracefully).
        """
        self._clear_all()
        self._color_index = state.get('color_index', 0)
        cd = state.get('contour_defaults', {})
        self._contour_levels = cd.get('levels', self._contour_levels)
        self._contour_min_level = cd.get('min_level', self._contour_min_level)
        self._contour_increment = cd.get('increment', self._contour_increment)

        for gd in state.get('groups', []):
            group = InspectorGroup(
                group_id=gd.get('group_id') or str(uuid.uuid4()),
                name=gd.get('name', 'Group'),
            )
            self.groups.append(group)
            self._library_panel.add_group(group)
            for sd in gd.get('spectra', []):
                spec = self._spectrum_from_state(sd)
                group.spectra.append(spec)
                self._library_panel.add_spectrum(spec, group.group_id)
                self._add_toolbar_spectrum_button(spec)

        self._active_spectrum_id = state.get('active_spectrum_id')
        self._canvas.set_active_spectrum(self._active_spectrum_id)

        # Restore peak-marker + trace toggles (block signals so we don't double-apply)
        pv = bool(state.get('peaks_visible', True))
        self._canvas._peaks_visible = pv
        self._set_checkbox_silently(self._library_panel._peaks_checkbox, pv)

        tr = state.get('traces', {})
        self._canvas._h_trace_on = bool(tr.get('horizontal', False))
        self._canvas._v_trace_on = bool(tr.get('vertical', False))
        self._canvas._trace_scale = float(tr.get('scale', self._canvas._trace_scale))
        self._set_checkbox_silently(self._library_panel._h_trace_cb, self._canvas._h_trace_on)
        self._set_checkbox_silently(self._library_panel._v_trace_cb, self._canvas._v_trace_on)
        self._library_panel.set_trace_scale_silently(self._canvas._trace_scale)

        self._refresh_canvas()

    @staticmethod
    def _set_checkbox_silently(checkbox, checked: bool):
        checkbox.blockSignals(True)
        checkbox.setChecked(checked)
        checkbox.blockSignals(False)

    def _spectrum_from_state(self, sd: dict) -> InspectorSpectrum:
        spec = InspectorSpectrum(
            spectrum_id=sd.get('spectrum_id') or str(uuid.uuid4()),
            display_name=sd.get('display_name', ''),
            file_path=sd.get('file_path', ''),
            color=sd.get('color', DEFAULT_COLORS[0]),
            visible=bool(sd.get('visible', True)),
            contour_levels=sd.get('contour_levels', self._contour_levels),
            contour_min_level=sd.get('contour_min_level', self._contour_min_level),
            contour_increment=sd.get('contour_increment', self._contour_increment),
        )
        spec.peaks = [InspectorPeak.from_state(p) for p in sd.get('peaks', [])]
        if spec.file_path and os.path.exists(spec.file_path):
            try:
                from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
                integrator = EnhancedVoigtIntegrator()
                if integrator.load_nmr_file(spec.file_path) and integrator.nmr_data is not None:
                    spec.integrator = integrator
                    spec.loaded = True
            except Exception as e:
                logger.warning(f"Could not reload spectrum {spec.file_path}: {e}")
        return spec

    def _clear_all(self):
        """Remove every group/spectrum from the model, tree, and toolbar."""
        for group in list(self.groups):
            for spec in list(group.spectra):
                self._remove_toolbar_button(spec.spectrum_id)
            self._library_panel.remove_group(group.group_id)
        self.groups = []

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

    def _shift_spectrum_peaks(self, spectrum_id: str, h_shift: float, n_shift: float) -> int:
        """Add coordinate offsets to every peak of a spectrum, tracking each peak's net shift.

        Returns the number of peaks shifted (0 if the spectrum is missing or empty).
        """
        spec = self._find_spectrum(spectrum_id)
        if spec is None:
            return 0
        for peak in spec.peaks:
            peak.ppm_x += h_shift
            peak.ppm_y += n_shift
            peak.applied_shift_h += h_shift
            peak.applied_shift_n += n_shift
        return len(spec.peaks)

    def _reset_spectrum_shift(self, spectrum_id: str) -> int:
        """Undo each peak's applied shift, restoring pre-shift positions.

        Peaks added or loaded after a shift carry no applied shift, so Reset
        leaves them where the user put them. Returns the number of peaks in the
        spectrum (0 if the spectrum is missing).
        """
        spec = self._find_spectrum(spectrum_id)
        if spec is None:
            return 0
        for peak in spec.peaks:
            peak.ppm_x -= peak.applied_shift_h
            peak.ppm_y -= peak.applied_shift_n
            peak.applied_shift_h = 0.0
            peak.applied_shift_n = 0.0
        return len(spec.peaks)

    def _load_peak_list_from_file(self, file_path: str, spectrum_id: str) -> int:
        """Parse a peak list file (.txt/.csv/.peaks/.json) and append its peaks.

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
                contour_min_level=self._contour_min_level,
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
                                  levels: int, min_level: float, increment: float):
        """Update contour display params for a single spectrum."""
        spec = self._find_spectrum(spectrum_id)
        if spec is None:
            return
        spec.contour_levels = levels
        spec.contour_min_level = min_level
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
            self._style_toolbar_button(spectrum_id)
            self._refresh_canvas()

    def _on_peak_list_attach_requested(self, spectrum_id: str):
        spec = self._find_spectrum(spectrum_id)
        if spec is None:
            return
        path, _ = QFileDialog.getOpenFileName(
            self, 'Load Peak List', '',
            'Peak list files (*.txt *.csv *.peaks *.json);;All files (*)'
        )
        if not path:
            return
        count = self._load_peak_list_from_file(path, spectrum_id)
        if count > 0:
            self._refresh_canvas()
            self.update_status(f"Loaded {count} peaks for {spec.display_name}")
        else:
            self.update_status(f"No peaks loaded from {os.path.basename(path)}")

    def _on_contour_update_requested(self, levels: int, min_level: float, increment: float):
        # Store as defaults for new spectra
        self._contour_levels = levels
        self._contour_min_level = min_level
        self._contour_increment = increment
        # Apply to all existing spectra
        for spec in self._all_spectra():
            spec.contour_levels = levels
            spec.contour_min_level = min_level
            spec.contour_increment = increment
        self._refresh_canvas()
        self.update_status(f'Contour: {levels} levels, min {min_level:.5f}, increment ×{increment:.2f} — applied to all')

    def _on_min_scale_requested(self, factor: float):
        """Multiply each spectrum's OWN min level by factor, preserving differences.

        Unlike the Min Level spinbox (which sets one absolute value on every
        spectrum), ×2/÷2 scales each spectrum relative to its current min level.
        """
        specs = self._all_spectra()
        if not specs:
            return
        for spec in specs:
            # clamp to the Min Level spinbox range (0.0001–1.0)
            spec.contour_min_level = min(max(spec.contour_min_level * factor, 0.0001), 1.0)
        self._refresh_canvas()
        self.update_status(f'Scaled each spectrum\'s min level ×{factor:g} ({len(specs)} spectra)')

    def _on_active_spectrum_changed(self, spectrum_id: str):
        self._active_spectrum_id = spectrum_id
        self._canvas.set_active_spectrum(spectrum_id)
        spec = self._find_spectrum(spectrum_id)
        if spec:
            self.update_status(
                f"Active: {spec.display_name}  |  Shift+click add peak · "
                f"middle-drag zoom box · Reset Zoom to restore"
            )

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

    def _reorder_spectrum(self, spectrum_id: str, delta: int) -> bool:
        """Move a spectrum by delta within its group. Later in the list draws on top.

        Returns True if the order changed, False at the ends / if not found.
        """
        for group in self.groups:
            ids = [s.spectrum_id for s in group.spectra]
            if spectrum_id not in ids:
                continue
            i = ids.index(spectrum_id)
            j = i + delta
            if j < 0 or j >= len(group.spectra):
                return False
            group.spectra[i], group.spectra[j] = group.spectra[j], group.spectra[i]
            return True
        return False

    def _rebuild_toolbar(self):
        """Rebuild the top-toolbar spectrum buttons to match the model (draw) order."""
        for sid in list(self._toolbar_actions.keys()):
            self._remove_toolbar_button(sid)
        for spec in self._all_spectra():
            self._add_toolbar_spectrum_button(spec)

    def _on_spectrum_reorder_requested(self, spectrum_id: str, delta: int):
        """Reorder a spectrum's draw stacking and update the tree, toolbar, and canvas."""
        if not self._reorder_spectrum(spectrum_id, delta):
            return
        group = next((g for g in self.groups
                      if any(s.spectrum_id == spectrum_id for s in g.spectra)), None)
        if group is not None:
            self._library_panel.rebuild_group_spectra(group)
            self._library_panel.select_spectrum(spectrum_id)
        self._rebuild_toolbar()
        self._refresh_canvas()
        spec = self._find_spectrum(spectrum_id)
        if spec:
            self.update_status(
                f"Reordered {spec.display_name} "
                f"({'toward front' if delta > 0 else 'toward back'})"
            )

    def _on_canvas_peak_add_requested(self, spectrum_id: str, ppm_x: float, ppm_y: float):
        peak = self._add_peak(spectrum_id, ppm_x, ppm_y)
        if peak:
            self._refresh_canvas()
            self.update_status(f"Peak added at ({ppm_x:.3f}, {ppm_y:.3f})")

    def _on_canvas_peak_delete_requested(self, spectrum_id: str, peak_id: str):
        self._delete_peak(spectrum_id, peak_id)
        self._refresh_canvas()
        self.update_status("Peak deleted")

    def _on_peak_markers_toggled(self, visible: bool):
        self._canvas.set_peaks_visible(visible)
        self.update_status('Peak markers ' + ('shown' if visible else 'hidden'))

    def _on_reset_zoom_requested(self):
        self._canvas.reset_zoom()
        self.update_status('View reset to full spectrum')

    def _on_horizontal_trace_toggled(self, on: bool):
        self._canvas.set_traces(horizontal=on)
        self.update_status('Horizontal 1D trace ' + ('on — hover the spectrum' if on else 'off'))

    def _on_vertical_trace_toggled(self, on: bool):
        self._canvas.set_traces(vertical=on)
        self.update_status('Vertical 1D trace ' + ('on — hover the spectrum' if on else 'off'))

    def _on_trace_scale_changed(self, fraction: float):
        self._canvas.set_trace_scale(fraction)

    def _on_export_pdf(self):
        """Export the current canvas (visible spectra + peaks) to a vector PDF."""
        has_content = any(
            s.visible and s.loaded and s.integrator is not None
            for s in self._all_spectra()
        )
        if not has_content:
            self.update_status('Nothing to export — load and show a spectrum first')
            return
        path, _ = QFileDialog.getSaveFileName(
            self, 'Export PDF', '', 'PDF files (*.pdf);;All files (*)'
        )
        if not path:
            return
        if not path.lower().endswith('.pdf'):
            path += '.pdf'
        try:
            self._canvas.export_pdf(path)
            self.update_status(f'Exported PDF: {os.path.basename(path)}')
        except Exception as e:
            logger.error(f'PDF export failed: {e}')
            self.update_status(f'PDF export failed: {e}')

    def _on_shift_apply_requested(self, h_shift: float, n_shift: float):
        """Shift the active spectrum's peak list by the given offsets."""
        if self._active_spectrum_id is None:
            self.update_status('Select a spectrum first to shift its peaks')
            return
        spec = self._find_spectrum(self._active_spectrum_id)
        if spec is None:
            self.update_status('Select a spectrum first to shift its peaks')
            return
        if not spec.peaks:
            self.update_status(f'{spec.display_name} has no peaks to shift')
            return
        count = self._shift_spectrum_peaks(self._active_spectrum_id, h_shift, n_shift)
        self._refresh_canvas()
        self.update_status(
            f'Shifted {count} peaks of {spec.display_name} by '
            f'{h_shift:+.3f} ppm (1H), {n_shift:+.2f} ppm (15N/13C)'
        )

    def _on_shift_reset_requested(self):
        """Undo the active spectrum's applied shift and restore original positions."""
        if self._active_spectrum_id is None:
            self.update_status('Select a spectrum first to reset its shift')
            return
        spec = self._find_spectrum(self._active_spectrum_id)
        if spec is None:
            self.update_status('Select a spectrum first to reset its shift')
            return
        count = self._reset_spectrum_shift(self._active_spectrum_id)
        self._refresh_canvas()
        self.update_status(f'Reset shift for {spec.display_name} ({count} peaks restored)')

    # ------------------------------------------------------------------
    # Toolbar button management
    # ------------------------------------------------------------------

    def _add_toolbar_spectrum_button(self, spec: InspectorSpectrum):
        """Add a visibility toggle button for this spectrum to the top toolbar."""
        btn = QPushButton(spec.display_name)
        btn.setCheckable(True)
        btn.setChecked(spec.visible)
        self._toolbar_buttons[spec.spectrum_id] = btn
        self._style_toolbar_button(spec.spectrum_id)
        btn.toggled.connect(
            lambda checked, sid=spec.spectrum_id: self._on_spectrum_visibility(sid, checked)
        )
        action = self._top_toolbar.addWidget(btn)
        self._toolbar_actions[spec.spectrum_id] = action

    def _style_toolbar_button(self, spectrum_id: str):
        """(Re)apply the spectrum's current color to its toolbar button."""
        btn = self._toolbar_buttons.get(spectrum_id)
        spec = self._find_spectrum(spectrum_id)
        if btn is None or spec is None:
            return
        btn.setStyleSheet(
            f"QPushButton {{ background-color: {spec.color}; color: white; "
            f"border: none; padding: 2px 8px; border-radius: 3px; font-weight: bold; }}"
            f"QPushButton:!checked {{ background-color: #888; color: #ccc; }}"
        )

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

    def _on_propagate_requested(self):
        """Open the picker to copy one spectrum's contour settings to chosen others."""
        spectra = self._all_spectra()
        if len(spectra) < 2:
            self.update_status('Need at least two spectra to propagate')
            return
        dlg = ContourPropagateDialog(
            spectra, active_id=self._active_spectrum_id, parent=self,
            on_apply=self._apply_propagation,
        )
        dlg.exec()

    def _apply_propagation(self, source_id: str, target_ids: List[str], fields: List[str]):
        """Copy the chosen contour fields from source to targets, then redraw."""
        if not target_ids or not fields:
            self.update_status('Propagate: nothing selected')
            return
        self._propagate_settings(source_id, target_ids, fields)
        self._refresh_canvas()
        self.update_status(
            f'Propagated {len(fields)} setting(s) to {len(target_ids)} spectrum(s)'
        )

    def _on_group_properties_requested(self, group_id: str):
        group = next((g for g in self.groups if g.group_id == group_id), None)
        if group is None or not group.spectra:
            self.update_status('No spectra in this group')
            return
        dlg = GroupPropertiesDialog(group, parent=self, on_apply=self._refresh_canvas)
        dlg.exec()
        self._refresh_canvas()

# ---------------------------------------------------------------------------
# Per-spectrum contour settings dialog
# ---------------------------------------------------------------------------

class SpectrumContourDialog:
    """Small dialog to edit contour display settings for one spectrum."""

    def __init__(self, spectrum: InspectorSpectrum, parent=None):
        from PySide6.QtWidgets import QDialog, QDialogButtonBox
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
        self._levels_spin.setRange(1, 50)
        self._levels_spin.setValue(spectrum.contour_levels)
        layout.addLayout(_row('Levels:', self._levels_spin))

        self._min_spin = QDoubleSpinBox()
        self._min_spin.setRange(0.0001, 1.0)
        self._min_spin.setDecimals(5)
        self._min_spin.setSingleStep(0.002)
        self._min_spin.setValue(spectrum.contour_min_level)
        layout.addLayout(_row('Min Level:', self._min_spin))

        self._incr_spin = QDoubleSpinBox()
        self._incr_spin.setRange(1.01, 10.0)
        self._incr_spin.setDecimals(2)
        self._incr_spin.setValue(spectrum.contour_increment)
        layout.addLayout(_row('Increment:', self._incr_spin))

        btns = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btns.accepted.connect(self._apply)
        btns.accepted.connect(self._dialog.accept)
        btns.rejected.connect(self._dialog.reject)
        layout.addWidget(btns)

    def exec(self) -> bool:
        return self._dialog.exec() == 1

    def _apply(self):
        self._spec.contour_levels = self._levels_spin.value()
        self._spec.contour_min_level = self._min_spin.value()
        self._spec.contour_increment = self._incr_spin.value()


# ---------------------------------------------------------------------------
# Group properties dialog
# ---------------------------------------------------------------------------

class GroupPropertiesDialog:
    """Dialog showing all spectra in a group with editable properties and propagation."""

    def __init__(self, group: InspectorGroup, parent=None, on_apply=None):
        from PySide6.QtWidgets import (
            QDialog, QDialogButtonBox, QTableWidget, QTableWidgetItem,
            QComboBox, QHeaderView, QAbstractItemView
        )
        self._group = group
        self._on_apply = on_apply
        self._dialog = QDialog(parent)
        self._dialog.setWindowTitle(f'Group Properties — {group.name}')
        self._dialog.setMinimumSize(680, 360)

        layout = QVBoxLayout(self._dialog)
        layout.setSpacing(6)

        # Table: one row per spectrum
        cols = ['Name', 'Visible', 'Color', 'Levels', 'Min Level', 'Increment']
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
            lev.setRange(1, 50)
            lev.setValue(spec.contour_levels)
            self._table.setCellWidget(row, 3, lev)

            # Min Level (fraction of max intensity)
            min_sp = QDoubleSpinBox()
            min_sp.setRange(0.0001, 1.0)
            min_sp.setDecimals(5)
            min_sp.setSingleStep(0.002)
            min_sp.setValue(spec.contour_min_level)
            self._table.setCellWidget(row, 4, min_sp)

            # Increment
            incr_sp = QDoubleSpinBox()
            incr_sp.setRange(1.01, 10.0)
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
            spec.contour_min_level = min_sp.value()
            spec.contour_increment = incr_sp.value()
        if self._on_apply is not None:
            self._on_apply()

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


# ---------------------------------------------------------------------------
# Contour settings propagation dialog
# ---------------------------------------------------------------------------

# Contour fields that can be propagated, in display order
_PROPAGATE_FIELDS = (
    ('Levels', 'contour_levels'),
    ('Min Level', 'contour_min_level'),
    ('Increment', 'contour_increment'),
)


class ContourPropagateDialog:
    """Copy contour settings from one spectrum to a user-chosen set of others.

    Calls ``on_apply(source_id, target_ids, fields)`` where ``fields`` are
    InspectorSpectrum attribute names (e.g. 'contour_min_level').
    """

    def __init__(self, spectra, active_id=None, parent=None, on_apply=None):
        from PySide6.QtWidgets import (
            QDialog, QDialogButtonBox, QComboBox, QListWidget, QListWidgetItem, QGroupBox
        )
        self._spectra = list(spectra)
        self._on_apply = on_apply
        self._dialog = QDialog(parent)
        self._dialog.setWindowTitle('Propagate Contour Settings')
        self._dialog.setMinimumWidth(320)

        layout = QVBoxLayout(self._dialog)
        layout.setSpacing(SPACING_SM)

        # Source spectrum
        src_row = QHBoxLayout()
        src_row.addWidget(QLabel('From:'))
        self._source_combo = QComboBox()
        for s in self._spectra:
            self._source_combo.addItem(s.display_name, s.spectrum_id)
        if active_id is not None:
            idx = next((i for i, s in enumerate(self._spectra)
                        if s.spectrum_id == active_id), -1)
            if idx >= 0:
                self._source_combo.setCurrentIndex(idx)
        self._source_combo.currentIndexChanged.connect(self._refresh_targets)
        src_row.addWidget(self._source_combo, stretch=1)
        layout.addLayout(src_row)

        # Which settings to copy
        param_box = QGroupBox('Settings to copy')
        pl = QVBoxLayout(param_box)
        self._param_checks = {}
        for label, attr in _PROPAGATE_FIELDS:
            cb = QCheckBox(label)
            cb.setChecked(True)
            pl.addWidget(cb)
            self._param_checks[attr] = cb
        layout.addWidget(param_box)

        # Target spectra
        layout.addWidget(QLabel('To:'))
        self._target_list = QListWidget()
        for s in self._spectra:
            item = QListWidgetItem(s.display_name)
            item.setData(Qt.UserRole, s.spectrum_id)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Checked)
            self._target_list.addItem(item)
        layout.addWidget(self._target_list)
        self._refresh_targets()

        btns = QDialogButtonBox(QDialogButtonBox.Apply | QDialogButtonBox.Close)
        btns.button(QDialogButtonBox.Apply).clicked.connect(self._apply)
        btns.button(QDialogButtonBox.Close).clicked.connect(self._dialog.accept)
        layout.addWidget(btns)

    def _refresh_targets(self):
        """Disable the source spectrum as a target; default all others to checked.

        Re-checking non-source rows on every source change keeps the sensible
        'propagate to all other spectra' default, including a row that has just
        stopped being the source.
        """
        source_id = self._source_combo.currentData()
        for i in range(self._target_list.count()):
            item = self._target_list.item(i)
            if item.data(Qt.UserRole) == source_id:
                item.setCheckState(Qt.Unchecked)
                item.setFlags(item.flags() & ~Qt.ItemIsEnabled)
            else:
                item.setFlags(item.flags() | Qt.ItemIsEnabled)
                item.setCheckState(Qt.Checked)

    def exec(self):
        return self._dialog.exec()

    def selected_targets(self) -> list:
        source_id = self._source_combo.currentData()
        return [
            self._target_list.item(i).data(Qt.UserRole)
            for i in range(self._target_list.count())
            if self._target_list.item(i).checkState() == Qt.Checked
            and self._target_list.item(i).data(Qt.UserRole) != source_id
        ]

    def selected_fields(self) -> list:
        return [attr for attr, cb in self._param_checks.items() if cb.isChecked()]

    def _apply(self):
        if self._on_apply is not None:
            self._on_apply(
                self._source_combo.currentData(),
                self.selected_targets(),
                self.selected_fields(),
            )
