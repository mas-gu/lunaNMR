# ABOUTME: Qt-based multi-spectrum overlay viewer dialog for comparing multiple NMR spectra
# ABOUTME: Port of Tkinter multi_spectrum_viewer.py to PySide6 for v1.0 Qt interface

import os
import logging
from typing import Optional, Dict, List, Any, TYPE_CHECKING

from PySide6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QCheckBox,
    QGroupBox, QScrollArea, QWidget, QColorDialog, QMessageBox,
    QSpinBox, QDoubleSpinBox, QSplitter, QListWidget, QListWidgetItem,
    QTabWidget, QSlider, QComboBox, QTextEdit
)
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor, QFont

from lunaNMR.gui.base.base_dialog import BaseDialog
from lunaNMR.gui.components.peak_navigator_table import PeakNavigatorTable
from lunaNMR.gui.styles.design_system import (
    SPACING_SM, SPACING_MD, SPACING_LG,
    FONT_SIZE_BODY, FONT_SIZE_SECTION_LABEL, FONT_SIZE_SMALL,
    PRIMARY_TEXT, SECONDARY_TEXT,
    PRIMARY_BUTTON_BG, PRIMARY_BUTTON_HOVER, PRIMARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BG, SECONDARY_BUTTON_HOVER, SECONDARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BORDER, BUTTON_CORNER_RADIUS, BUTTON_HEIGHT_DIALOG,
    FRAME_BG_COLOR
)

logger = logging.getLogger(__name__)


class SpectrumListItem(QWidget):
    """Custom widget for spectrum item with checkbox and color picker."""

    visibility_changed = Signal(str, bool)  # spectrum_name, visible
    color_changed = Signal(str, str)  # spectrum_name, color_hex

    def __init__(self, spectrum_name: str, color: str = "#1f77b4",
                 visible: bool = False, parent=None):
        super().__init__(parent)
        self.spectrum_name = spectrum_name
        self.color = color

        layout = QHBoxLayout()
        layout.setContentsMargins(4, 2, 4, 2)
        layout.setSpacing(SPACING_SM)

        # Visibility checkbox
        self.checkbox = QCheckBox()
        self.checkbox.setChecked(visible)
        self.checkbox.toggled.connect(self._on_visibility_changed)
        layout.addWidget(self.checkbox)

        # Color button
        self.color_btn = QPushButton()
        self.color_btn.setFixedSize(24, 24)
        self.color_btn.setStyleSheet(f"background-color: {color}; border: 1px solid gray;")
        self.color_btn.clicked.connect(self._on_color_clicked)
        layout.addWidget(self.color_btn)

        # Name label
        name_label = QLabel(spectrum_name)
        name_label.setStyleSheet(f"font-size: {FONT_SIZE_BODY}px;")
        layout.addWidget(name_label, stretch=1)

        self.setLayout(layout)

    def _on_visibility_changed(self, checked: bool):
        self.visibility_changed.emit(self.spectrum_name, checked)

    def _on_color_clicked(self):
        color = QColorDialog.getColor(QColor(self.color), self, "Choose Color")
        if color.isValid():
            self.color = color.name()
            self.color_btn.setStyleSheet(f"background-color: {self.color}; border: 1px solid gray;")
            self.color_changed.emit(self.spectrum_name, self.color)


class MultiSpectrumViewerDialog(BaseDialog):
    """Dialog for viewing multiple spectra overlaid for comparison.

    This is a Qt port of the Tkinter MultiSpectrumViewer from multi_spectrum_viewer.py.

    Features:
        - List of spectra with visibility toggles
        - Color customization per spectrum
        - Contour control settings
        - Peak markers from selected spectrum

    Based on v0.9 MultiSpectrumViewer (multi_spectrum_viewer.py:137-1306)

    Example:
        ```python
        dialog = MultiSpectrumViewerDialog(parent, all_results, file_manager)
        dialog.show()
        ```
    """

    def __init__(self, parent=None, all_results: List[Dict] = None,
                 file_manager=None):
        """Initialize the multi-spectrum viewer dialog.

        Args:
            parent: Parent widget
            all_results: List of result dictionaries for all spectra
            file_manager: NMRFileManager instance for loading spectra
        """
        super().__init__(
            parent=parent,
            title="Multi-Spectrum Overlay Viewer",
            default_size=(1200, 800),
            modal=False  # Non-modal to allow interaction
        )

        self.all_results = all_results or []
        self.file_manager = file_manager

        # Spectrum data management
        self.spectra = []  # List of spectrum info dicts
        self.visible_spectra = set()  # Names of visible spectra

        # Default color palette
        self.default_colors = [
            '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
            '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
        ]

        # Voigt analysis state (v0.9 multi_spectrum_viewer.py:198-200)
        self.selected_spectrum_index = 0
        self.selected_peak_index = None
        self.voigt_2d_plotter = None
        self.voigt_3d_plotter = None

        # Initialize spectra data
        self._initialize_spectra()

        # Build UI
        self.setup_ui()

        # Center on parent
        if parent:
            self.center_on_parent()
        else:
            self.center_on_screen()

        logger.debug(f"MultiSpectrumViewerDialog initialized with {len(self.spectra)} spectra")

    def _initialize_spectra(self):
        """Initialize spectrum data from results."""
        for idx, result in enumerate(self.all_results):
            spectrum_name = result.get('spectrum_name', f'spectrum_{idx+1:03d}')
            file_path = result.get('spectrum_file', '')
            color = self.default_colors[idx % len(self.default_colors)]

            # Extract fitted peaks from result_data
            # v0.9 multi_spectrum_viewer.py:117-125 - check integration_results first
            fitted_peaks = result.get('integration_results', [])
            if not fitted_peaks:
                fitted_peaks = result.get('fitted_peaks', [])
            if fitted_peaks is None:
                fitted_peaks = []

            self.spectra.append({
                'name': spectrum_name,
                'file_path': file_path,
                'color': color,
                'visible': idx == 0,  # First spectrum visible by default
                'result_data': result,
                'fitted_peaks': fitted_peaks,  # Store extracted peaks
                'loaded': False,
                'data': None
            })

            if idx == 0:
                self.visible_spectra.add(spectrum_name)

            logger.debug(f"Initialized spectrum {spectrum_name} with {len(fitted_peaks)} fitted peaks")

    def _get_fitted_peaks(self, spec_idx: int) -> List[Dict]:
        """Get fitted peaks for a spectrum.

        Follows v0.9 pattern: check integration_results first, then fitted_peaks.
        Based on v0.9 multi_spectrum_viewer.py:117-125.
        """
        if spec_idx < 0 or spec_idx >= len(self.spectra):
            return []

        spec = self.spectra[spec_idx]

        # First check cached fitted_peaks
        if spec.get('fitted_peaks'):
            return spec['fitted_peaks']

        # Fall back to result_data extraction
        result_data = spec.get('result_data', {})
        fitted_peaks = result_data.get('integration_results', [])
        if not fitted_peaks:
            fitted_peaks = result_data.get('fitted_peaks', [])
        if fitted_peaks is None:
            fitted_peaks = []

        return fitted_peaks

    def setup_ui(self):
        """Setup the dialog user interface."""
        layout = QVBoxLayout()
        layout.setSpacing(SPACING_SM)
        layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)

        # Main content splitter (no title - maximizes space)
        splitter = QSplitter(Qt.Horizontal)

        # Left panel: Spectrum list + Contour controls (compact)
        left_panel = self.create_spectrum_list_panel()
        splitter.addWidget(left_panel)

        # Right panel: Plot area
        right_panel = self.create_plot_panel()
        splitter.addWidget(right_panel)

        # Set splitter proportions (left column 25% smaller: 280 * 0.75 = 210)
        splitter.setSizes([210, 990])

        layout.addWidget(splitter, stretch=1)

        self.setLayout(layout)

    def create_spectrum_list_panel(self) -> QWidget:
        """Create the spectrum list panel with contour controls.

        Returns:
            QWidget containing spectrum list and contour controls
        """
        # Main container widget
        container = QWidget()
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(SPACING_SM)

        # Spectra group
        spectra_group = QGroupBox("Spectra")
        spectra_group.setStyleSheet(self._get_group_style())

        spectra_layout = QVBoxLayout()
        spectra_layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)
        spectra_layout.setSpacing(2)

        # Info label
        info_label = QLabel(f"{len(self.spectra)} spectra available")
        info_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}px; color: {SECONDARY_TEXT};")
        spectra_layout.addWidget(info_label)

        # Scroll area for spectrum list
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setStyleSheet("QScrollArea { border: none; }")

        scroll_content = QWidget()
        scroll_layout = QVBoxLayout()
        scroll_layout.setSpacing(2)
        scroll_layout.setContentsMargins(0, 0, 0, 0)

        # Add spectrum items
        self.spectrum_items = {}
        for spec in self.spectra:
            item = SpectrumListItem(
                spectrum_name=spec['name'],
                color=spec['color'],
                visible=spec['visible']
            )
            item.visibility_changed.connect(self._on_visibility_changed)
            item.color_changed.connect(self._on_color_changed)
            scroll_layout.addWidget(item)
            self.spectrum_items[spec['name']] = item

        scroll_layout.addStretch()
        scroll_content.setLayout(scroll_layout)
        scroll_area.setWidget(scroll_content)

        spectra_layout.addWidget(scroll_area, stretch=1)

        # Quick selection buttons
        btn_layout = QHBoxLayout()
        btn_layout.setSpacing(SPACING_SM)

        select_all_btn = QPushButton("Select All")
        select_all_btn.clicked.connect(self._select_all)
        btn_layout.addWidget(select_all_btn)

        select_none_btn = QPushButton("Select None")
        select_none_btn.clicked.connect(self._select_none)
        btn_layout.addWidget(select_none_btn)

        spectra_layout.addLayout(btn_layout)
        spectra_group.setLayout(spectra_layout)
        main_layout.addWidget(spectra_group, stretch=1)

        # Contour controls group (compact vertical layout)
        contour_group = QGroupBox("Contour")
        contour_group.setStyleSheet(self._get_group_style())

        contour_layout = QVBoxLayout()
        contour_layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)
        contour_layout.setSpacing(2)

        # Row 1: Levels
        levels_row = QHBoxLayout()
        levels_row.addWidget(QLabel("Levels:"))
        self.levels_spin = QSpinBox()
        self.levels_spin.setRange(5, 100)
        self.levels_spin.setValue(10)
        levels_row.addWidget(self.levels_spin)
        contour_layout.addLayout(levels_row)

        # Row 2: Min Level
        min_row = QHBoxLayout()
        min_row.addWidget(QLabel("Min:"))
        self.min_spin = QDoubleSpinBox()
        self.min_spin.setRange(0.01, 10.0)
        self.min_spin.setValue(0.2)
        self.min_spin.setSingleStep(0.01)
        min_row.addWidget(self.min_spin)
        contour_layout.addLayout(min_row)

        # Row 3: Increment
        inc_row = QHBoxLayout()
        inc_row.addWidget(QLabel("Inc:"))
        self.increment_spin = QDoubleSpinBox()
        self.increment_spin.setRange(1.01, 2.0)
        self.increment_spin.setValue(2.0)
        self.increment_spin.setSingleStep(0.01)
        inc_row.addWidget(self.increment_spin)
        contour_layout.addLayout(inc_row)

        # Update button
        update_btn = QPushButton("Update")
        update_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {PRIMARY_BUTTON_BG};
                color: {PRIMARY_BUTTON_TEXT};
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px;
                font-size: {FONT_SIZE_BODY}px;
            }}
            QPushButton:hover {{
                background-color: {PRIMARY_BUTTON_HOVER};
            }}
        """)
        update_btn.clicked.connect(self._update_plot)
        contour_layout.addWidget(update_btn)

        contour_group.setLayout(contour_layout)
        main_layout.addWidget(contour_group)

        container.setLayout(main_layout)
        return container

    def create_plot_panel(self) -> QWidget:
        """Create the plot panel with tabbed views: Overlay, Voigt Analysis, 3D Voigt.

        Based on v0.9 multi_spectrum_viewer.py tab structure (lines 304-312).

        Returns:
            QWidget containing tabbed plot area
        """
        panel = QWidget()
        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        # Create tab widget
        self.plot_tabs = QTabWidget()

        # Tab 1: Overlay View
        overlay_tab = self._create_overlay_tab()
        self.plot_tabs.addTab(overlay_tab, "Overlay View")

        # Tab 2: Voigt Analysis
        voigt_tab = self._create_voigt_analysis_tab()
        self.plot_tabs.addTab(voigt_tab, "Voigt Analysis")

        # Tab 3: 3D Voigt Analysis
        voigt_3d_tab = self._create_voigt_3d_analysis_tab()
        self.plot_tabs.addTab(voigt_3d_tab, "3D Voigt Analysis")

        # Tab 4: Peak Parameters
        peak_params_tab = self._create_peak_parameters_tab()
        self.plot_tabs.addTab(peak_params_tab, "Peak Parameters")

        layout.addWidget(self.plot_tabs)
        panel.setLayout(layout)
        return panel

    def _create_overlay_tab(self) -> QWidget:
        """Create the overlay view tab with contour plot and peak navigator.

        Based on v0.9 multi_spectrum_viewer.py overlay tab structure.
        Includes Select Spectrum and Peak Navigator panels on the right.
        """
        from lunaNMR.gui.components.matplotlib_widget import MatplotlibWidget

        tab = QWidget()
        layout = QHBoxLayout()
        layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)
        layout.setSpacing(SPACING_SM)

        # Left panel (75%): Contour plot
        left_panel = QWidget()
        left_layout = QVBoxLayout()
        left_layout.setContentsMargins(0, 0, 0, 0)

        # Create matplotlib widget for contour plot
        self.plot_widget = MatplotlibWidget(
            parent=left_panel,
            toolbar=True,
            figsize=(8, 6),
            dpi=100
        )

        # Initialize plot with placeholder text
        self.plot_widget.axes.set_title("Select spectra to display overlay")
        self.plot_widget.axes.text(
            0.5, 0.5,
            "Select spectra from the list\nto display contour overlay",
            ha='center', va='center',
            transform=self.plot_widget.axes.transAxes,
            fontsize=12, color='gray'
        )
        self.plot_widget.refresh()

        left_layout.addWidget(self.plot_widget, stretch=1)
        left_panel.setLayout(left_layout)
        layout.addWidget(left_panel, stretch=7)

        # Right panel (30%): Spectrum selector and Peak Navigator (20% larger than 25%)
        right_panel = QWidget()
        right_layout = QVBoxLayout()
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.setSpacing(SPACING_SM)

        # Spectrum selector for overlay
        spec_group = QGroupBox("Select Spectrum")
        spec_group.setStyleSheet(self._get_group_style())
        spec_layout = QVBoxLayout()

        self.overlay_spectrum_list = QListWidget()
        self.overlay_spectrum_list.setMaximumHeight(75)
        for i, spec in enumerate(self.spectra):
            item = QListWidgetItem(spec['name'])
            item.setData(Qt.UserRole, i)
            self.overlay_spectrum_list.addItem(item)
        self.overlay_spectrum_list.itemClicked.connect(self._on_overlay_spectrum_selected)
        spec_layout.addWidget(self.overlay_spectrum_list)
        spec_group.setLayout(spec_layout)
        right_layout.addWidget(spec_group)

        # Peak navigator for overlay (table-based)
        peak_group = QGroupBox("Peak Navigator")
        peak_group.setStyleSheet(self._get_group_style())
        peak_layout = QVBoxLayout()

        # Use table-based PeakNavigatorTable for consistent UI with main GUI
        self.overlay_peak_table = PeakNavigatorTable()
        self.overlay_peak_table.peak_selected.connect(self._on_overlay_peak_selected_by_index)
        peak_layout.addWidget(self.overlay_peak_table)

        # Show Peaks checkbox - toggle peak markers on overlay plot
        self.show_peaks_checkbox = QCheckBox("Show Peak Markers")
        self.show_peaks_checkbox.setChecked(True)
        self.show_peaks_checkbox.toggled.connect(self._on_toggle_peak_markers)
        peak_layout.addWidget(self.show_peaks_checkbox)

        peak_group.setLayout(peak_layout)
        right_layout.addWidget(peak_group, stretch=1)

        right_panel.setLayout(right_layout)
        layout.addWidget(right_panel, stretch=3)

        # Track selected spectrum for overlay peak display
        self.overlay_selected_spectrum_idx = 0
        self.overlay_show_peaks = True

        tab.setLayout(layout)
        return tab

    def _create_voigt_analysis_tab(self) -> QWidget:
        """Create the Voigt analysis tab with spectrum/peak selection and VoigtAnalysisPlotter.

        Based on v0.9 setup_voigt_analysis_tab() (multi_spectrum_viewer.py:819-902).
        """
        from lunaNMR.gui.components.voigt_analysis_plotter import VoigtAnalysisPlotter

        tab = QWidget()
        layout = QHBoxLayout()
        layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)
        layout.setSpacing(SPACING_SM)

        # Left panel (70%): Voigt analysis plot
        left_panel = QWidget()
        left_layout = QVBoxLayout()
        left_layout.setContentsMargins(0, 0, 0, 0)

        plot_group = QGroupBox("Voigt Analysis")
        plot_group.setStyleSheet(self._get_group_style())
        plot_layout = QVBoxLayout()

        # Create VoigtAnalysisPlotter
        self.voigt_2d_plotter = VoigtAnalysisPlotter(parent=plot_group, toolbar=True)
        plot_layout.addWidget(self.voigt_2d_plotter)
        plot_group.setLayout(plot_layout)
        left_layout.addWidget(plot_group)

        left_panel.setLayout(left_layout)
        layout.addWidget(left_panel, stretch=7)

        # Right panel (30%): Spectrum and peak selection (20% larger)
        right_panel = QWidget()
        right_layout = QVBoxLayout()
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.setSpacing(SPACING_SM)

        # Spectrum selector
        spec_group = QGroupBox("Select Spectrum")
        spec_group.setStyleSheet(self._get_group_style())
        spec_layout = QVBoxLayout()

        self.voigt_spectrum_list = QListWidget()
        self.voigt_spectrum_list.setMaximumHeight(75)
        for i, spec in enumerate(self.spectra):
            item = QListWidgetItem(spec['name'])
            item.setData(Qt.UserRole, i)
            self.voigt_spectrum_list.addItem(item)
        self.voigt_spectrum_list.itemClicked.connect(self._on_voigt_spectrum_selected)
        spec_layout.addWidget(self.voigt_spectrum_list)
        spec_group.setLayout(spec_layout)
        right_layout.addWidget(spec_group)

        # Peak navigator (table-based)
        peak_group = QGroupBox("Peak Navigator")
        peak_group.setStyleSheet(self._get_group_style())
        peak_layout = QVBoxLayout()

        # Use table-based PeakNavigatorTable for consistent UI
        self.voigt_peak_table = PeakNavigatorTable()
        self.voigt_peak_table.peak_selected.connect(self._on_voigt_peak_selected_by_index)
        peak_layout.addWidget(self.voigt_peak_table)

        peak_group.setLayout(peak_layout)
        right_layout.addWidget(peak_group, stretch=1)

        right_panel.setLayout(right_layout)
        layout.addWidget(right_panel, stretch=3)

        # Initialize with placeholder
        self.voigt_2d_plotter.show_placeholder("Select spectrum and peak to view Voigt analysis")

        tab.setLayout(layout)
        return tab

    def _create_voigt_3d_analysis_tab(self) -> QWidget:
        """Create the 3D Voigt analysis tab with VoigtAnalysisPlotter.

        Based on v0.9 main_gui.py Tab 3: 3D Voigt Analysis (lines 1829-1944).
        Full controls matching v0.9 implementation.
        """
        from lunaNMR.gui.components.voigt_analysis_plotter import VoigtAnalysisPlotter
        from PySide6.QtWidgets import QRadioButton, QButtonGroup

        tab = QWidget()
        layout = QHBoxLayout()
        layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)
        layout.setSpacing(SPACING_SM)

        # Left panel: Controls + 3D plot
        left_panel = QWidget()
        left_layout = QVBoxLayout()
        left_layout.setContentsMargins(0, 0, 0, 0)

        # === Row 1: Layer Visibility Controls (compact, bold text) ===
        row1_layout = QHBoxLayout()
        row1_layout.setContentsMargins(0, 0, 0, 0)
        row1_layout.setSpacing(SPACING_SM)

        layer_group = QGroupBox("Layer Visibility")
        layer_group.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                font-size: 13px;
                margin-top: 8px;
                padding-top: 8px;
                padding-bottom: 2px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 6px;
                padding: 0 3px;
            }
            QCheckBox { font-weight: bold; font-size: 13px; }
        """)
        layer_layout = QHBoxLayout()
        layer_layout.setContentsMargins(6, 2, 6, 2)
        layer_layout.setSpacing(10)

        # Initialize control state variables
        self.show_exp_3d = True
        self.show_fit_3d = True
        self.show_individual_3d = False
        self.show_peak_labels_3d = False
        self.show_resid_3d = False
        self.limit_peak_display_3d = True

        self.exp_3d_checkbox = QCheckBox("Exp")
        self.exp_3d_checkbox.setChecked(True)
        self.exp_3d_checkbox.toggled.connect(self._on_toggle_exp_3d)
        layer_layout.addWidget(self.exp_3d_checkbox)

        self.fit_3d_checkbox = QCheckBox("Fit")
        self.fit_3d_checkbox.setChecked(True)
        self.fit_3d_checkbox.toggled.connect(self._on_toggle_fit_3d)
        layer_layout.addWidget(self.fit_3d_checkbox)

        self.individual_3d_checkbox = QCheckBox("Individual")
        self.individual_3d_checkbox.setChecked(False)
        self.individual_3d_checkbox.toggled.connect(self._on_toggle_individual_3d)
        layer_layout.addWidget(self.individual_3d_checkbox)

        self.resid_3d_checkbox = QCheckBox("Residuals")
        self.resid_3d_checkbox.setChecked(False)
        self.resid_3d_checkbox.toggled.connect(self._on_toggle_resid_3d)
        layer_layout.addWidget(self.resid_3d_checkbox)

        self.limit_3d_checkbox = QCheckBox("Limit Peak")
        self.limit_3d_checkbox.setChecked(True)
        self.limit_3d_checkbox.toggled.connect(self._on_toggle_limit_3d)
        layer_layout.addWidget(self.limit_3d_checkbox)

        layer_group.setLayout(layer_layout)
        row1_layout.addWidget(layer_group)

        # Color scheme dropdown (next to Layer Visibility, compact)
        color_group = QGroupBox("Color Scheme")
        color_group.setStyleSheet("""
            QGroupBox {
                font-weight: bold;
                font-size: 13px;
                margin-top: 8px;
                padding-top: 8px;
                padding-bottom: 2px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 6px;
                padding: 0 3px;
            }
        """)
        color_layout = QHBoxLayout()
        color_layout.setContentsMargins(6, 2, 6, 2)

        self.color_scheme_3d = 'Clean'
        self.colormap_3d_combo = QComboBox()
        self.colormap_3d_combo.addItems(['Classic', 'Clean', 'Dark', 'Warm'])
        self.colormap_3d_combo.setCurrentText('Clean')
        self.colormap_3d_combo.setMinimumWidth(100)
        self.colormap_3d_combo.setStyleSheet("QComboBox { font-size: 12px; padding: 2px 6px; }")
        self.colormap_3d_combo.currentTextChanged.connect(self._on_3d_colormap_changed)
        color_layout.addWidget(self.colormap_3d_combo)

        color_group.setLayout(color_layout)
        row1_layout.addWidget(color_group)
        row1_layout.addStretch()
        left_layout.addLayout(row1_layout)

        # Keep residual mode and intensity scale state variables functional (hidden controls)
        self.residual_mode_3d = 'separate'
        self.intensity_scale_3d = 100

        # === Plot area ===
        plot_group = QGroupBox("3D Voigt Surface Analysis")
        plot_group.setStyleSheet(self._get_group_style())
        plot_layout = QVBoxLayout()

        # Create VoigtAnalysisPlotter for 3D
        self.voigt_3d_plotter = VoigtAnalysisPlotter(parent=plot_group, toolbar=True)
        plot_layout.addWidget(self.voigt_3d_plotter)

        plot_group.setLayout(plot_layout)
        left_layout.addWidget(plot_group, stretch=1)

        left_panel.setLayout(left_layout)
        layout.addWidget(left_panel, stretch=7)

        # Right panel (30%): Spectrum and peak selection (20% larger)
        right_panel = QWidget()
        right_layout = QVBoxLayout()
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.setSpacing(SPACING_SM)

        # Spectrum selector for 3D
        spec_group = QGroupBox("Select Spectrum")
        spec_group.setStyleSheet(self._get_group_style())
        spec_layout = QVBoxLayout()

        self.voigt_3d_spectrum_list = QListWidget()
        self.voigt_3d_spectrum_list.setMaximumHeight(75)
        for i, spec in enumerate(self.spectra):
            item = QListWidgetItem(spec['name'])
            item.setData(Qt.UserRole, i)
            self.voigt_3d_spectrum_list.addItem(item)
        self.voigt_3d_spectrum_list.itemClicked.connect(self._on_voigt_3d_spectrum_selected)
        spec_layout.addWidget(self.voigt_3d_spectrum_list)
        spec_group.setLayout(spec_layout)
        right_layout.addWidget(spec_group)

        # Peak navigator for 3D (table-based)
        peak_group = QGroupBox("Peak Navigator")
        peak_group.setStyleSheet(self._get_group_style())
        peak_layout = QVBoxLayout()

        # Use table-based PeakNavigatorTable for consistent UI
        self.voigt_3d_peak_table = PeakNavigatorTable()
        self.voigt_3d_peak_table.peak_selected.connect(self._on_voigt_3d_peak_selected_by_index)
        peak_layout.addWidget(self.voigt_3d_peak_table)

        peak_group.setLayout(peak_layout)
        right_layout.addWidget(peak_group, stretch=1)

        right_panel.setLayout(right_layout)
        layout.addWidget(right_panel, stretch=3)

        # Initialize with placeholder
        self.voigt_3d_plotter.show_placeholder("Select spectrum and peak to view 3D Voigt analysis")

        tab.setLayout(layout)
        return tab

    def _create_peak_parameters_tab(self) -> QWidget:
        """Create the Peak Parameters tab with QTextEdit for parameter display.

        Based on main_window.py _create_peak_parameters_tab().

        Returns:
            QWidget containing a read-only QTextEdit showing peak fit results
        """
        tab_widget = QWidget()
        layout = QVBoxLayout(tab_widget)
        layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)
        layout.setSpacing(0)

        # Create text editor for peak parameters
        self.peak_params_text = QTextEdit()
        self.peak_params_text.setReadOnly(True)
        self.peak_params_text.setFont(QFont("Courier", 9))
        self.peak_params_text.setPlainText(
            "Peak Parameters\n"
            "===============\n\n"
            "Peak parameters will appear here after selecting a peak.\n\n"
            "Expected output format:\n"
            "Peak #    X (ppm)    Y (ppm)    Width_X    Width_Y    Intensity    R²\n"
            "----  -----------  -----------  ---------  ---------  -----------  ------\n"
        )
        self.peak_params_text.setStyleSheet(f"""
            QTextEdit {{
                background-color: {FRAME_BG_COLOR};
                color: {PRIMARY_TEXT};
                border: none;
                padding: {SPACING_MD}px;
                font-family: Courier;
                font-size: 9px;
            }}
        """)
        layout.addWidget(self.peak_params_text)

        logger.debug("Peak Parameters tab created with QTextEdit")
        return tab_widget

    def update_peak_parameters(self, voigt_result: dict):
        """Update Peak Parameters tab with Voigt fitting results.

        Based on main_window.py update_peak_parameters().

        Args:
            voigt_result: Dictionary containing Voigt fitting results
        """
        if not hasattr(self, 'peak_params_text') or self.peak_params_text is None:
            return

        # Extract result data
        assignment = voigt_result.get('assignment', 'Unknown')
        quality = voigt_result.get('fitting_quality', 'Unknown')
        method = voigt_result.get('method', '')
        all_peaks = voigt_result.get('all_peaks', [])
        r_squared = voigt_result.get('avg_r_squared', voigt_result.get('r_squared', 0.0))

        # Color palette matching visualization
        colors = ['red', 'orange', 'purple', 'brown', 'pink', 'olive', 'cyan', 'magenta']

        # Build text content
        lines = []

        # Header
        lines.append(f"Peak Parameters - {assignment}")
        lines.append("=" * 60)
        lines.append("")

        # Overall quality info
        lines.append(f"Overall Quality: {quality}")
        lines.append(f"R² = {r_squared:.4f}")
        lines.append(f"Fitting Method: {method}")
        lines.append(f"Number of Peaks: {len(all_peaks)}")
        lines.append("")

        # Individual peak parameters
        if len(all_peaks) > 0:
            lines.append("Individual Peak Details:")
            lines.append("=" * 60)
            lines.append("")

            for i, peak in enumerate(all_peaks):
                color = colors[i % len(colors)]
                peak_assignment = peak.get('assignment', f'Peak {i+1}')

                # Extract parameters
                pos_f2 = peak.get('pos_f2', peak.get('center_x', 0.0))
                pos_f1 = peak.get('pos_f1', peak.get('center_y', 0.0))
                volume = peak.get('volume', peak.get('intensity', 0.0))
                height = peak.get('height', peak.get('amplitude', 0.0))

                # Linewidths
                lw_gau_f2 = peak.get('lw_gau_f2', peak.get('sigma_x', 0.0))
                lw_lor_f2 = peak.get('lw_lor_f2', peak.get('gamma_x', 0.0))
                lw_gau_f1 = peak.get('lw_gau_f1', peak.get('sigma_y', 0.0))
                lw_lor_f1 = peak.get('lw_lor_f1', peak.get('gamma_y', 0.0))

                lw_f2 = lw_gau_f2 + lw_lor_f2
                lw_f1 = lw_gau_f1 + lw_lor_f1

                # Lorentz/Gauss ratio
                lg_ratio_f2 = lw_lor_f2 / lw_f2 if lw_f2 > 0 else 0.0
                lg_ratio_f1 = lw_lor_f1 / lw_f1 if lw_f1 > 0 else 0.0

                # Peak header
                lines.append(f"{peak_assignment} ({color}):")
                lines.append("-" * 50)

                # Position
                lines.append(f"  Position (F2, F1):     ({pos_f2:.4f}, {pos_f1:.2f}) ppm")

                # Intensity
                lines.append(f"  Height:                {height:.3e}")
                lines.append(f"  Volume:                {volume:.3e}")
                lines.append("")

                # Linewidths F2 (1H)
                lines.append(f"  Linewidth F2 (1H):     {lw_f2:.4f} ppm")
                lines.append(f"    - Gaussian (sigma):  {lw_gau_f2:.4f} ppm")
                lines.append(f"    - Lorentzian (gamma):{lw_lor_f2:.4f} ppm")
                lines.append(f"    - L/G Ratio:         {lg_ratio_f2:.2%}")
                lines.append("")

                # Linewidths F1 (15N/13C)
                lines.append(f"  Linewidth F1 (15N):    {lw_f1:.3f} ppm")
                lines.append(f"    - Gaussian (sigma):  {lw_gau_f1:.3f} ppm")
                lines.append(f"    - Lorentzian (gamma):{lw_lor_f1:.3f} ppm")
                lines.append(f"    - L/G Ratio:         {lg_ratio_f1:.2%}")

                if i < len(all_peaks) - 1:
                    lines.append("")
        else:
            lines.append("No peak details available.")

        # Update text widget
        self.peak_params_text.setPlainText("\n".join(lines))
        logger.debug(f"Updated Peak Parameters tab for {assignment}")

    def _get_group_style(self) -> str:
        """Get standard group box style."""
        return f"""
            QGroupBox {{
                font-size: {FONT_SIZE_BODY}px;
                font-weight: bold;
                color: {PRIMARY_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                margin-top: {SPACING_SM}px;
                padding-top: {SPACING_MD}px;
                background-color: {FRAME_BG_COLOR};
            }}
            QGroupBox::title {{
                subcontrol-origin: margin;
                left: {SPACING_SM}px;
                padding: 0 {SPACING_SM}px;
            }}
        """

    def create_contour_controls(self) -> QGroupBox:
        """Create contour control settings.

        Returns:
            QGroupBox containing contour controls
        """
        group = QGroupBox("Contour Controls (Applied to All Spectra)")
        group.setStyleSheet(f"""
            QGroupBox {{
                font-size: {FONT_SIZE_BODY}px;
                font-weight: bold;
                color: {PRIMARY_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                margin-top: {SPACING_SM}px;
                padding-top: {SPACING_MD}px;
                background-color: {FRAME_BG_COLOR};
            }}
            QGroupBox::title {{
                subcontrol-origin: margin;
                left: {SPACING_SM}px;
                padding: 0 {SPACING_SM}px;
            }}
        """)

        layout = QHBoxLayout()
        layout.setSpacing(SPACING_MD)
        layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)

        # Levels
        layout.addWidget(QLabel("Levels:"))
        self.levels_spin = QSpinBox()
        self.levels_spin.setRange(5, 100)
        self.levels_spin.setValue(10)
        self.levels_spin.setMaximumWidth(80)
        layout.addWidget(self.levels_spin)

        # Min level
        layout.addWidget(QLabel("Min Level:"))
        self.min_spin = QDoubleSpinBox()
        self.min_spin.setRange(0.01, 10.0)
        self.min_spin.setValue(0.2)
        self.min_spin.setSingleStep(0.01)
        self.min_spin.setMaximumWidth(80)
        layout.addWidget(self.min_spin)

        # Increment
        layout.addWidget(QLabel("Increment:"))
        self.increment_spin = QDoubleSpinBox()
        self.increment_spin.setRange(1.01, 2.0)
        self.increment_spin.setValue(2.0)
        self.increment_spin.setSingleStep(0.01)
        self.increment_spin.setMaximumWidth(80)
        layout.addWidget(self.increment_spin)

        layout.addStretch()

        # Update plot button
        update_btn = QPushButton("Update Plot")
        update_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {PRIMARY_BUTTON_BG};
                color: {PRIMARY_BUTTON_TEXT};
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px {SPACING_MD}px;
                font-size: {FONT_SIZE_BODY}px;
            }}
            QPushButton:hover {{
                background-color: {PRIMARY_BUTTON_HOVER};
            }}
        """)
        update_btn.clicked.connect(self._update_plot)
        layout.addWidget(update_btn)

        # Export button
        export_btn = QPushButton("Export PNG")
        export_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {SECONDARY_BUTTON_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px {SPACING_MD}px;
                font-size: {FONT_SIZE_BODY}px;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
        """)
        export_btn.clicked.connect(self._export_plot)
        layout.addWidget(export_btn)

        group.setLayout(layout)
        return group

    def create_button_row(self) -> QHBoxLayout:
        """Create the bottom button row.

        Returns:
            QHBoxLayout containing buttons
        """
        layout = QHBoxLayout()
        layout.addStretch()

        # Close button
        close_btn = QPushButton("Close")
        close_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        close_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {SECONDARY_BUTTON_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px {SPACING_MD}px;
                font-size: {FONT_SIZE_BODY}px;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
        """)
        close_btn.clicked.connect(self.close)
        layout.addWidget(close_btn)

        return layout

    def _on_visibility_changed(self, spectrum_name: str, visible: bool):
        """Handle spectrum visibility toggle."""
        if visible:
            self.visible_spectra.add(spectrum_name)
        else:
            self.visible_spectra.discard(spectrum_name)

        # Update spectrum data
        for spec in self.spectra:
            if spec['name'] == spectrum_name:
                spec['visible'] = visible
                break

        logger.debug(f"Spectrum '{spectrum_name}' visibility: {visible}")
        self._update_plot()

    def _on_color_changed(self, spectrum_name: str, color: str):
        """Handle spectrum color change."""
        for spec in self.spectra:
            if spec['name'] == spectrum_name:
                spec['color'] = color
                break

        logger.debug(f"Spectrum '{spectrum_name}' color: {color}")
        self._update_plot()

    def _select_all(self):
        """Select all spectra."""
        for name, item in self.spectrum_items.items():
            item.checkbox.setChecked(True)

    def _select_none(self):
        """Deselect all spectra."""
        for name, item in self.spectrum_items.items():
            item.checkbox.setChecked(False)

    def _update_plot(self):
        """Update the overlay plot with contours for visible spectra.

        Based on v0.9 update_overlay_plot() (multi_spectrum_viewer.py:647-747)
        """
        import numpy as np

        visible_count = len(self.visible_spectra)
        logger.debug(f"Updating plot with {visible_count} visible spectra")

        # Save current zoom state before clearing
        xlim = self.plot_widget.axes.get_xlim()
        ylim = self.plot_widget.axes.get_ylim()
        has_valid_zoom = (xlim[0] != 0.0 or xlim[1] != 1.0) and (ylim[0] != 0.0 or ylim[1] != 1.0)

        # Clear current plot
        self.plot_widget.clear()

        if visible_count == 0:
            # Show placeholder message
            self.plot_widget.axes.set_title("Select spectra to display overlay")
            self.plot_widget.axes.text(
                0.5, 0.5,
                "No spectra selected\nCheck spectrum boxes to display",
                ha='center', va='center',
                transform=self.plot_widget.axes.transAxes,
                fontsize=12, color='gray'
            )
            self.plot_widget.refresh()
            return

        # Get contour parameters
        num_levels = self.levels_spin.value()
        min_level = self.min_spin.value()

        # Track if any spectra were successfully plotted
        plotted_count = 0
        legend_labels = []
        legend_colors = []

        # Plot each visible spectrum
        for spec in self.spectra:
            if spec['name'] not in self.visible_spectra:
                continue

            # Load spectrum data if not already loaded
            if not spec['loaded']:
                if not self._load_spectrum_data(spec):
                    logger.warning(f"Failed to load spectrum: {spec['name']}")
                    continue

            # Get data and axes from integrator
            integrator = spec.get('integrator')
            if integrator is None:
                continue

            try:
                data = integrator.nmr_data
                ppm_x = integrator.ppm_x_axis
                ppm_y = integrator.ppm_y_axis

                if data is None or ppm_x is None or ppm_y is None:
                    logger.warning(f"Missing data for {spec['name']}")
                    continue

                # Calculate contour levels
                data_abs = np.abs(data)
                max_intensity = np.max(data_abs)
                min_intensity = min_level * max_intensity

                # Use geomspace for better contour distribution
                level_values = np.geomspace(min_intensity, max_intensity * 0.8, num_levels)

                # Plot contours with spectrum color
                self.plot_widget.axes.contour(
                    ppm_x, ppm_y, data,
                    levels=level_values,
                    colors=spec['color'],
                    linewidths=1.0,
                    alpha=0.8
                )

                plotted_count += 1
                legend_labels.append(spec['name'])
                legend_colors.append(spec['color'])

                logger.debug(f"Plotted {spec['name']} (color: {spec['color']})")

            except Exception as e:
                logger.error(f"Error plotting {spec['name']}: {e}")
                continue

        if plotted_count > 0:
            # Configure axes for NMR convention
            self.plot_widget.axes.set_xlabel('¹H (ppm)', fontsize=10)
            self.plot_widget.axes.set_ylabel('¹⁵N/¹³C (ppm)', fontsize=10)

            # Restore zoom if valid, otherwise use autoscale with NMR axis inversion
            if has_valid_zoom:
                self.plot_widget.axes.set_xlim(xlim)
                self.plot_widget.axes.set_ylim(ylim)
            else:
                self.plot_widget.axes.invert_xaxis()
                self.plot_widget.axes.invert_yaxis()

            self.plot_widget.axes.set_title(
                f'Multi-Spectrum Overlay ({plotted_count} spectra)',
                fontsize=11
            )

            # Add legend
            if legend_labels:
                import matplotlib.pyplot as plt
                legend_handles = [
                    plt.Line2D([0], [0], color=color, linewidth=2)
                    for color in legend_colors
                ]
                self.plot_widget.axes.legend(
                    legend_handles, legend_labels,
                    loc='best', fontsize=8
                )

            # Plot peak markers for selected spectrum (v0.9 multi_spectrum_viewer.py:764-817)
            if getattr(self, 'overlay_show_peaks', True):
                self._plot_overlay_peak_markers()

        self.plot_widget.refresh()
        logger.info(f"Plot updated with {plotted_count} spectra")

    def _load_spectrum_data(self, spec: Dict) -> bool:
        """Load spectrum data for a spectrum entry.

        Based on v0.9 SpectrumData.load_data() (multi_spectrum_viewer.py:84-134)

        Args:
            spec: Spectrum info dictionary

        Returns:
            True if loaded successfully, False otherwise
        """
        import os

        file_path = spec.get('file_path', '')
        if not file_path or not os.path.exists(file_path):
            logger.warning(f"File not found: {file_path}")
            return False

        try:
            # Create integrator and load data
            from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator

            integrator = EnhancedVoigtIntegrator()
            success = integrator.load_nmr_file(file_path)

            if not success:
                logger.warning(f"Failed to load NMR file: {file_path}")
                return False

            # Verify data loaded correctly
            if integrator.nmr_data is None:
                logger.warning(f"NMR data is None after loading: {file_path}")
                return False

            if not hasattr(integrator, 'ppm_x_axis') or not hasattr(integrator, 'ppm_y_axis'):
                logger.warning(f"Missing PPM axes after loading: {file_path}")
                return False

            # Store in spectrum info
            spec['integrator'] = integrator
            spec['data'] = integrator.nmr_data
            spec['loaded'] = True

            logger.info(f"Loaded spectrum: {spec['name']}")
            return True

        except Exception as e:
            logger.error(f"Error loading {spec['name']}: {e}")
            return False

    def _export_plot(self):
        """Export the overlay plot as PNG."""
        from PySide6.QtWidgets import QFileDialog

        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Export Overlay Plot",
            "",
            "PNG files (*.png);;All files (*.*)"
        )

        if filename:
            try:
                self.plot_widget.save_figure(filename, dpi=300)
                QMessageBox.information(
                    self,
                    "Export Successful",
                    f"Plot exported to:\n{filename}"
                )
                logger.info(f"Exported plot to {filename}")
            except Exception as e:
                QMessageBox.critical(
                    self,
                    "Export Error",
                    f"Failed to export plot:\n{str(e)}"
                )

    # ===== Voigt Analysis Methods (v0.9 multi_spectrum_viewer.py:1113-1305) =====

    def _on_voigt_spectrum_selected(self, item: QListWidgetItem):
        """Handle spectrum selection in Voigt analysis tab."""
        spec_idx = item.data(Qt.UserRole)
        self.selected_spectrum_index = spec_idx

        # Load peaks into the table
        fitted_peaks = self._get_fitted_peaks(spec_idx)
        self.voigt_peak_table.load_peaks(fitted_peaks)

        logger.debug(f"Selected spectrum {spec_idx} for Voigt analysis")

    def _on_voigt_peak_selected_by_index(self, peak_idx: int):
        """Handle peak selection from Voigt PeakNavigatorTable.

        Plots the Voigt analysis for the selected peak.
        """
        self.selected_peak_index = peak_idx
        self._plot_voigt_analysis()

    def _on_voigt_3d_spectrum_selected(self, item: QListWidgetItem):
        """Handle spectrum selection in 3D Voigt analysis tab."""
        spec_idx = item.data(Qt.UserRole)
        self.selected_spectrum_index = spec_idx

        # Load peaks into the table
        fitted_peaks = self._get_fitted_peaks(spec_idx)
        self.voigt_3d_peak_table.load_peaks(fitted_peaks)

        logger.debug(f"Selected spectrum {spec_idx} for 3D Voigt analysis")

    def _on_voigt_3d_peak_selected_by_index(self, peak_idx: int):
        """Handle peak selection from 3D Voigt PeakNavigatorTable.

        Plots the 3D Voigt analysis for the selected peak.
        """
        self.selected_peak_index = peak_idx
        self._plot_voigt_analysis_3d()

    def _plot_voigt_analysis(self):
        """Plot Voigt analysis for currently selected peak.

        Based on v0.9 plot_voigt_analysis() (multi_spectrum_viewer.py:1222-1274).
        Uses _get_fitted_peaks() which checks integration_results first.
        """
        if self.selected_peak_index is None:
            logger.debug("No peak selected")
            return

        if self.selected_spectrum_index < 0 or self.selected_spectrum_index >= len(self.spectra):
            logger.debug("Invalid spectrum index")
            return

        spec = self.spectra[self.selected_spectrum_index]
        fitted_peaks = self._get_fitted_peaks(self.selected_spectrum_index)

        if self.selected_peak_index >= len(fitted_peaks):
            logger.debug("Invalid peak index")
            return

        peak = fitted_peaks[self.selected_peak_index]
        assignment = peak.get('assignment', peak.get('Assignment', 'Unknown'))

        # Find stored Voigt result - the peak itself may contain the full result
        stored_result = self._find_stored_voigt_result(spec, assignment, peak)

        try:
            if stored_result and self.voigt_2d_plotter:
                self.voigt_2d_plotter.plot_voigt_analysis(stored_result)
                # Also update Peak Parameters tab
                self.update_peak_parameters(stored_result)
                logger.debug(f"Plotted Voigt analysis for {assignment}")
            else:
                if self.voigt_2d_plotter:
                    self.voigt_2d_plotter.show_placeholder(
                        f"No Voigt fitting results available for peak {assignment}"
                    )
        except Exception as e:
            logger.error(f"Error plotting Voigt analysis: {e}")
            if self.voigt_2d_plotter:
                self.voigt_2d_plotter.show_placeholder(f"Error: {str(e)[:50]}")

    def _plot_voigt_analysis_3d(self):
        """Plot 3D Voigt analysis for currently selected peak.

        Uses _get_fitted_peaks() which checks integration_results first.
        """
        if self.selected_peak_index is None:
            return

        if self.selected_spectrum_index < 0 or self.selected_spectrum_index >= len(self.spectra):
            return

        spec = self.spectra[self.selected_spectrum_index]
        fitted_peaks = self._get_fitted_peaks(self.selected_spectrum_index)

        if self.selected_peak_index >= len(fitted_peaks):
            return

        peak = fitted_peaks[self.selected_peak_index]
        assignment = peak.get('assignment', peak.get('Assignment', 'Unknown'))

        # Find stored Voigt result - the peak itself may contain the full result
        stored_result = self._find_stored_voigt_result(spec, assignment, peak)

        try:
            if stored_result and self.voigt_3d_plotter:
                self.voigt_3d_plotter.plot_voigt_analysis_3d(stored_result)
                # Also update Peak Parameters tab
                self.update_peak_parameters(stored_result)
                logger.debug(f"Plotted 3D Voigt analysis for {assignment}")
            else:
                if self.voigt_3d_plotter:
                    self.voigt_3d_plotter.show_placeholder(
                        f"No Voigt fitting results available for peak {assignment}"
                    )
        except Exception as e:
            logger.error(f"Error plotting 3D Voigt analysis: {e}")
            if self.voigt_3d_plotter:
                self.voigt_3d_plotter.show_placeholder(f"Error: {str(e)[:50]}")

    def _find_stored_voigt_result(self, spec: Dict, assignment: str, peak: Dict = None) -> Optional[Dict]:
        """Find stored Voigt fitting results for a peak by assignment.

        Based on v0.9 find_stored_voigt_result() (multi_spectrum_viewer.py:1276-1305).
        Updated to check integration_results first and accept peak directly.

        Args:
            spec: Spectrum info dictionary
            assignment: Peak assignment string
            peak: Optional peak dict that may contain the full result

        Returns:
            Complete peak result dict with region_2d and fit data, or None
        """
        if not assignment:
            return None

        # First check if peak itself contains region_2d (it's the full result)
        # This is the case when integration_results contains the full Voigt data
        if peak and peak.get('region_2d') is not None:
            logger.debug(f"Using peak directly as Voigt result for {assignment}")
            return peak

        result_data = spec.get('result_data', {})

        # Check integration_results first (v0.9 pattern)
        integration_results = result_data.get('integration_results', [])
        if integration_results:
            for result in integration_results:
                result_assignment = result.get('assignment') or result.get('Assignment')
                if result_assignment == assignment:
                    r_squared = (result.get('r_squared') or
                                result.get('R_Squared') or
                                result.get('avg_r_squared') or 0)
                    logger.debug(f"Found integration_results Voigt result for {assignment}: R²={r_squared:.3f}")
                    return result

        # Check fitted_results (detailed Voigt analysis)
        fitted_results = result_data.get('fitted_results', [])
        if fitted_results:
            for result in fitted_results:
                result_assignment = result.get('assignment') or result.get('Assignment')
                if result_assignment == assignment:
                    r_squared = (result.get('r_squared') or
                                result.get('R_Squared') or
                                result.get('avg_r_squared') or 0)
                    logger.debug(f"Found fitted_results Voigt result for {assignment}: R²={r_squared:.3f}")
                    return result

        # Fallback: check integrator fitted_peaks
        integrator = spec.get('integrator')
        if integrator and hasattr(integrator, 'fitted_peaks'):
            if integrator.fitted_peaks:
                for result in integrator.fitted_peaks:
                    if result.get('assignment') == assignment:
                        logger.debug(f"Found integrator Voigt result for {assignment}")
                        return result

        # Also check fitted_peaks in result_data
        fitted_peaks = result_data.get('fitted_peaks', [])
        for result in fitted_peaks:
            result_assignment = result.get('assignment') or result.get('Assignment')
            if result_assignment == assignment:
                return result

        logger.debug(f"No stored Voigt result found for {assignment}")
        return None

    def _on_3d_intensity_changed(self, value: int):
        """Handle intensity scale slider change for 3D Voigt plot.

        Based on v0.9 main_gui.py - uses set_intensity_scale() method.
        """
        # Update the label
        if hasattr(self, 'intensity_3d_label'):
            self.intensity_3d_label.setText(f"{value}%")

        # Use set_intensity_scale() method on the plotter
        if self.voigt_3d_plotter:
            self.voigt_3d_plotter.set_intensity_scale(float(value))

    def _on_3d_colormap_changed(self, scheme: str):
        """Handle color scheme change for 3D Voigt plot.

        Based on v0.9 main_gui.py - uses set_color_scheme() method.
        """
        self.color_scheme_3d = scheme
        if self.voigt_3d_plotter:
            self.voigt_3d_plotter.set_color_scheme(scheme)

    # === 3D Voigt Control Handlers (v0.9 main_gui.py:816-829) ===
    def _on_toggle_exp_3d(self, checked):
        """Toggle experimental layer visibility."""
        self.show_exp_3d = checked
        if self.voigt_3d_plotter:
            self.voigt_3d_plotter.toggle_experimental(checked)

    def _on_toggle_fit_3d(self, checked):
        """Toggle fitted layer visibility."""
        self.show_fit_3d = checked
        if self.voigt_3d_plotter:
            self.voigt_3d_plotter.toggle_fitted(checked)

    def _on_toggle_individual_3d(self, checked):
        """Toggle individual peaks visibility."""
        self.show_individual_3d = checked
        if self.voigt_3d_plotter:
            self.voigt_3d_plotter.toggle_individual_peaks(checked)

    def _on_toggle_labels_3d(self, checked):
        """Toggle peak labels visibility."""
        self.show_peak_labels_3d = checked
        if self.voigt_3d_plotter:
            self.voigt_3d_plotter.toggle_peak_labels(checked)

    def _on_toggle_resid_3d(self, checked):
        """Toggle residuals visibility."""
        self.show_resid_3d = checked
        if self.voigt_3d_plotter:
            self.voigt_3d_plotter.toggle_residuals(checked)

    def _on_toggle_limit_3d(self, checked):
        """Toggle peak extent limiting."""
        self.limit_peak_display_3d = checked
        if self.voigt_3d_plotter:
            self.voigt_3d_plotter.toggle_peak_clipping(checked)

    def _on_residual_mode_changed(self, mode):
        """Handle residual mode change."""
        self.residual_mode_3d = mode
        if self.voigt_3d_plotter:
            self.voigt_3d_plotter.set_residual_mode(mode)

    # ===== Overlay Tab Methods =====

    def _on_overlay_spectrum_selected(self, item: QListWidgetItem):
        """Handle spectrum selection in overlay tab.

        Loads peaks into the PeakNavigatorTable for the selected spectrum.
        """
        spec_idx = item.data(Qt.UserRole)
        self.overlay_selected_spectrum_idx = spec_idx

        # Load peaks into the table
        fitted_peaks = self._get_fitted_peaks(spec_idx)
        self.overlay_peak_table.load_peaks(fitted_peaks)

        logger.debug(f"Selected spectrum {spec_idx} for overlay peak display")
        # Refresh the plot to show peaks for newly selected spectrum
        self._update_plot()

    def _on_overlay_peak_selected_by_index(self, peak_idx: int):
        """Handle peak selection from overlay PeakNavigatorTable.

        Zooms to the selected peak on the overlay plot.
        """
        if self.overlay_selected_spectrum_idx < 0 or self.overlay_selected_spectrum_idx >= len(self.spectra):
            return

        fitted_peaks = self._get_fitted_peaks(self.overlay_selected_spectrum_idx)
        if peak_idx >= len(fitted_peaks):
            return

        peak = fitted_peaks[peak_idx]

        # Get peak position
        pos_x = (peak.get('center_x') or peak.get('peak_x') or peak.get('ppm_x') or
                peak.get('pos_f2') or peak.get('position_x') or 0)
        pos_y = (peak.get('center_y') or peak.get('peak_y') or peak.get('ppm_y') or
                peak.get('pos_f1') or peak.get('position_y') or 0)

        # Zoom to peak region (typical zoom window: ±0.2 ppm 1H, ±2 ppm 15N)
        zoom_margin_x = 0.2  # 1H dimension
        zoom_margin_y = 2.0  # 15N dimension

        try:
            ax = self.plot_widget.axes
            ax.set_xlim(pos_x + zoom_margin_x, pos_x - zoom_margin_x)  # Inverted for NMR
            ax.set_ylim(pos_y + zoom_margin_y, pos_y - zoom_margin_y)  # Inverted for NMR
            self.plot_widget.refresh()
            logger.debug(f"Zoomed to peak at ({pos_x:.3f}, {pos_y:.2f})")
        except Exception as e:
            logger.error(f"Error zooming to peak: {e}")

    def _on_toggle_peak_markers(self, checked: bool):
        """Toggle display of peak markers on overlay plot."""
        self.overlay_show_peaks = checked
        self._update_plot()

    def _plot_overlay_peak_markers(self):
        """Plot peak markers on overlay spectrum.

        Based on v0.9 multi_spectrum_viewer.py:764-817 (plot_peak_markers).
        Shows peaks from the selected spectrum with quality-based coloring.
        """
        spec_idx = getattr(self, 'overlay_selected_spectrum_idx', 0)

        if spec_idx < 0 or spec_idx >= len(self.spectra):
            return

        fitted_peaks = self._get_fitted_peaks(spec_idx)
        if not fitted_peaks:
            return

        try:
            x_coords = []
            y_coords = []
            labels = []
            colors = []

            # Extract peak coordinates and quality info
            for peak in fitted_peaks:
                # Try multiple possible keys for coordinates
                x_ppm = (peak.get('center_x') or peak.get('peak_x') or peak.get('ppm_x') or
                        peak.get('pos_f2') or peak.get('Position_X') or 0)
                y_ppm = (peak.get('center_y') or peak.get('peak_y') or peak.get('ppm_y') or
                        peak.get('pos_f1') or peak.get('Position_Y') or 0)
                assignment = peak.get('assignment', peak.get('Assignment', ''))
                quality = peak.get('quality', peak.get('fitting_quality', peak.get('Quality', 'Unknown')))

                if x_ppm and y_ppm and x_ppm != 0 and y_ppm != 0:
                    x_coords.append(float(x_ppm))
                    y_coords.append(float(y_ppm))
                    labels.append(assignment)

                    # Color code by quality (same as v0.9)
                    if quality in ['Excellent', 'excellent']:
                        colors.append('lime')
                    elif quality in ['Good', 'good']:
                        colors.append('green')
                    elif quality in ['Fair', 'fair']:
                        colors.append('orange')
                    elif quality in ['Poor', 'poor']:
                        colors.append('red')
                    else:
                        colors.append('yellow')  # Default for unknown

            # Plot peak markers
            if x_coords and y_coords:
                self.plot_widget.axes.scatter(
                    x_coords, y_coords, c=colors, s=40, marker='o',
                    alpha=0.8, edgecolors='black', linewidths=1,
                    label='Peaks', zorder=10
                )

                # Add peak labels
                for x, y, label in zip(x_coords, y_coords, labels):
                    if label and label != 'Unknown':
                        self.plot_widget.axes.annotate(
                            label, (x, y), xytext=(5, 5),
                            textcoords='offset points', fontsize=7,
                            color='black', weight='bold',
                            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7)
                        )

                logger.debug(f"Plotted {len(x_coords)} peak markers for spectrum {spec_idx}")

        except Exception as e:
            logger.error(f"Error plotting peak markers: {e}")
