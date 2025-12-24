# ABOUTME: Qt-based detailed spectrum viewer dialog for individual spectrum analysis
# ABOUTME: Port of Tkinter SpectrumViewer class to PySide6 with contour plot, peak nav, Voigt analysis

import os
import logging
from typing import Optional, Dict, List, Any

import numpy as np

from PySide6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QCheckBox,
    QGroupBox, QScrollArea, QWidget, QSplitter, QSpinBox,
    QDoubleSpinBox, QComboBox, QMessageBox, QListWidget, QListWidgetItem,
    QTabWidget, QSlider
)
from PySide6.QtCore import Qt, Signal

from lunaNMR.gui.base.base_dialog import BaseDialog
from lunaNMR.gui.components.matplotlib_widget import MatplotlibWidget, MatplotlibMultiAxesWidget
from lunaNMR.gui.components.nmr_navigation_handler import NMRNavigationHandler
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


class SpectrumViewerDialog(BaseDialog):
    """Detailed spectrum viewer dialog with contour plot, peak navigator, and Voigt analysis.

    This is a Qt port of the Tkinter SpectrumViewer from spectrum_browser.py.

    Features:
        - 2D contour plot with navigation toolbar
        - Peak navigator list
        - 1D cross-section plots (F1 and F2)
        - Voigt analysis display
        - Quality metrics
        - Contour level controls

    Based on v0.9 SpectrumViewer (spectrum_browser.py:691-1306)

    Example:
        ```python
        dialog = SpectrumViewerDialog(
            parent,
            spectrum_name="spectrum_001",
            spectrum_file_path="/path/to/spectrum.ft2",
            result_data={'peaks': [...], 'avg_r_squared': 0.95}
        )
        dialog.show()
        ```
    """

    peak_selected = Signal(int)  # Emits peak index when selected

    def __init__(self, parent=None, spectrum_name: str = "",
                 spectrum_file_path: str = "", result_data: Dict = None,
                 series_processor=None):
        """Initialize the spectrum viewer dialog.

        Args:
            parent: Parent widget
            spectrum_name: Name of spectrum being viewed
            spectrum_file_path: Path to NMR spectrum file
            result_data: Result dictionary containing peaks, metrics, etc.
            series_processor: Optional series processor for additional context
        """
        super().__init__(
            parent=parent,
            title=f"Spectrum Viewer - {spectrum_name}",
            default_size=(1400, 900),
            modal=False  # Non-modal to allow multiple viewers
        )

        self.spectrum_name = spectrum_name
        self.spectrum_file_path = spectrum_file_path
        self.result_data = result_data or {}
        self.series_processor = series_processor

        # State variables
        self.integrator = None
        self.selected_peak_idx = None
        self.fitted_peaks = []
        self.fitted_results = self.result_data.get('fitted_results', [])

        # Contour control values
        self.contour_levels = 20
        self.contour_min = 0.15
        self.contour_increment = 1.2

        # Voigt analysis plotters (initialized in setup_ui)
        self.voigt_2d_plotter = None
        self.voigt_3d_plotter = None

        # Build UI
        self.setup_ui()

        # Load spectrum data
        self._load_spectrum()

        # Update displays
        self._update_contour_plot()
        self._populate_peak_list()

        # Center on parent
        if parent:
            self.center_on_parent()
        else:
            self.center_on_screen()

        logger.debug(f"SpectrumViewerDialog initialized for {spectrum_name}")

    def setup_ui(self):
        """Setup the dialog user interface."""
        layout = QVBoxLayout()
        layout.setSpacing(SPACING_SM)
        layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)

        # Title bar with info
        title_layout = self._create_title_bar()
        layout.addLayout(title_layout)

        # Main splitter for plots and peak navigator
        main_splitter = QSplitter(Qt.Horizontal)

        # Left side: Plots
        plots_widget = self._create_plots_panel()
        main_splitter.addWidget(plots_widget)

        # Right side: Peak navigator and controls
        right_panel = self._create_right_panel()
        main_splitter.addWidget(right_panel)

        # Set splitter proportions (plots get more space)
        main_splitter.setSizes([1000, 400])

        layout.addWidget(main_splitter, stretch=1)

        # Contour controls bar
        contour_bar = self._create_contour_controls()
        layout.addWidget(contour_bar)

        # Button row
        button_layout = self._create_button_row()
        layout.addLayout(button_layout)

        self.setLayout(layout)

    def _create_title_bar(self) -> QHBoxLayout:
        """Create the title bar with spectrum info."""
        layout = QHBoxLayout()

        # Spectrum name
        title_label = QLabel(f"Spectrum: {self.spectrum_name}")
        title_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_SECTION_LABEL}px;
                font-weight: bold;
                color: {PRIMARY_TEXT};
            }}
        """)
        layout.addWidget(title_label)

        layout.addStretch()

        # Quality info
        if self.result_data:
            status = self.result_data.get('status', 'Unknown')
            detection_rate = self.result_data.get('detection_rate', 0)
            avg_r2 = self.result_data.get('avg_r_squared', 0)

            info_text = f"Status: {status} | Detection: {detection_rate:.1f}% | Avg R²: {avg_r2:.4f}"
            info_label = QLabel(info_text)
            info_label.setStyleSheet(f"""
                QLabel {{
                    font-size: {FONT_SIZE_BODY}px;
                    color: {'#34C759' if status == 'success' else '#FF3B30'};
                }}
            """)
            layout.addWidget(info_label)

        return layout

    def _create_plots_panel(self) -> QWidget:
        """Create the plots panel with tabbed views: Main, Voigt Analysis, 3D Voigt.

        Based on v0.9 spectrum_browser.py tab structure (lines 762-775).
        """
        panel = QWidget()
        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        # Create tab widget
        self.plot_tabs = QTabWidget()

        # Tab 1: Main Spectrum
        main_tab = self._create_main_spectrum_tab()
        self.plot_tabs.addTab(main_tab, "Main Spectrum")

        # Tab 2: Voigt Analysis
        voigt_tab = self._create_voigt_analysis_tab()
        self.plot_tabs.addTab(voigt_tab, "Voigt Analysis")

        # Tab 3: 3D Voigt Analysis
        voigt_3d_tab = self._create_voigt_3d_analysis_tab()
        self.plot_tabs.addTab(voigt_3d_tab, "3D Voigt Analysis")

        layout.addWidget(self.plot_tabs)
        panel.setLayout(layout)
        return panel

    def _create_main_spectrum_tab(self) -> QWidget:
        """Create the main spectrum tab with contour and cross-section plots."""
        tab = QWidget()
        layout = QVBoxLayout()
        layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)
        layout.setSpacing(SPACING_SM)

        # Main contour plot
        contour_group = QGroupBox("2D Contour Plot")
        contour_group.setStyleSheet(self._get_group_style())
        contour_layout = QVBoxLayout()

        self.contour_widget = MatplotlibWidget(
            parent=contour_group,
            toolbar=False,
            figsize=(8, 6),
            dpi=100
        )
        contour_layout.addWidget(self.contour_widget)

        # Attach navigation handler for pan/zoom
        self._nav_handler = NMRNavigationHandler()
        self._nav_handler.attach(self.contour_widget)
        contour_group.setLayout(contour_layout)

        layout.addWidget(contour_group, stretch=2)

        # Cross-section plots (horizontal splitter)
        cross_section_splitter = QSplitter(Qt.Horizontal)

        # F2 (1H) cross-section
        f2_group = QGroupBox("F2 Cross-Section (¹H)")
        f2_group.setStyleSheet(self._get_group_style())
        f2_layout = QVBoxLayout()
        self.f2_widget = MatplotlibWidget(
            parent=f2_group,
            toolbar=False,
            figsize=(4, 2),
            dpi=80
        )
        f2_layout.addWidget(self.f2_widget)
        f2_group.setLayout(f2_layout)
        cross_section_splitter.addWidget(f2_group)

        # F1 (15N/13C) cross-section
        f1_group = QGroupBox("F1 Cross-Section (¹⁵N/¹³C)")
        f1_group.setStyleSheet(self._get_group_style())
        f1_layout = QVBoxLayout()
        self.f1_widget = MatplotlibWidget(
            parent=f1_group,
            toolbar=False,
            figsize=(4, 2),
            dpi=80
        )
        f1_layout.addWidget(self.f1_widget)
        f1_group.setLayout(f1_layout)
        cross_section_splitter.addWidget(f1_group)

        layout.addWidget(cross_section_splitter, stretch=1)

        tab.setLayout(layout)
        return tab

    def _create_voigt_analysis_tab(self) -> QWidget:
        """Create the Voigt analysis tab with VoigtAnalysisPlotter.

        Based on v0.9 spectrum_browser.py setup_voigt_analysis_tab() (lines 873-908).
        """
        from lunaNMR.gui.components.voigt_analysis_plotter import VoigtAnalysisPlotter

        tab = QWidget()
        layout = QVBoxLayout()
        layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)
        layout.setSpacing(SPACING_SM)

        # Voigt analysis plot container
        plot_group = QGroupBox("Voigt Analysis")
        plot_group.setStyleSheet(self._get_group_style())
        plot_layout = QVBoxLayout()

        # Create VoigtAnalysisPlotter with 2×2 grid
        self.voigt_2d_plotter = VoigtAnalysisPlotter(parent=plot_group, toolbar=True)
        plot_layout.addWidget(self.voigt_2d_plotter)

        plot_group.setLayout(plot_layout)
        layout.addWidget(plot_group)

        # Initialize with placeholder message
        self.voigt_2d_plotter.show_placeholder("Select a peak to view Voigt analysis")

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
        layout = QVBoxLayout()
        layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)
        layout.setSpacing(SPACING_SM)

        # === Row 1: Layer Visibility Controls (bold text, 25% larger font) ===
        row1_layout = QHBoxLayout()
        layer_group = QGroupBox("Layer Visibility")
        layer_group.setStyleSheet("""
            QGroupBox { font-weight: bold; font-size: 15px; }
            QCheckBox { font-weight: bold; font-size: 15px; }
        """)
        layer_layout = QHBoxLayout()
        layer_layout.setContentsMargins(8, 3, 8, 3)
        layer_layout.setSpacing(15)

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

        # Color scheme dropdown (next to Layer Visibility)
        color_group = QGroupBox("Color Scheme")
        color_group.setStyleSheet(self._get_group_style())
        color_layout = QHBoxLayout()

        self.color_scheme_3d = 'Clean'
        self.colormap_3d_combo = QComboBox()
        self.colormap_3d_combo.addItems(['Classic', 'Clean', 'Dark', 'Warm'])
        self.colormap_3d_combo.setCurrentText('Clean')
        # Make dropdown larger with larger text
        self.colormap_3d_combo.setMinimumWidth(120)
        self.colormap_3d_combo.setStyleSheet("QComboBox { font-size: 13px; padding: 4px 8px; }")
        self.colormap_3d_combo.currentTextChanged.connect(self._on_3d_colormap_changed)
        color_layout.addWidget(self.colormap_3d_combo)

        color_group.setLayout(color_layout)
        row1_layout.addWidget(color_group)
        row1_layout.addStretch()
        layout.addLayout(row1_layout)

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
        layout.addWidget(plot_group, stretch=1)

        # Initialize with placeholder message
        self.voigt_3d_plotter.show_placeholder("Select a peak to view 3D Voigt analysis")

        tab.setLayout(layout)
        return tab

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

    def _create_right_panel(self) -> QWidget:
        """Create the right panel with peak navigator and info."""
        panel = QWidget()
        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(SPACING_SM)

        # Peak Navigator
        peak_group = QGroupBox("Peak Navigator")
        peak_group.setStyleSheet(self._get_group_style())
        peak_layout = QVBoxLayout()

        # Peak list
        self.peak_list = QListWidget()
        self.peak_list.setStyleSheet("""
            QListWidget {
                background-color: white;
                border: 1px solid #C7C7CC;
                border-radius: 6px;
            }
            QListWidget::item {
                padding: 4px;
            }
            QListWidget::item:selected {
                background-color: #007AFF;
                color: white;
            }
        """)
        self.peak_list.itemClicked.connect(self._on_peak_selected)
        peak_layout.addWidget(self.peak_list)

        # Navigation buttons
        nav_layout = QHBoxLayout()
        prev_btn = QPushButton("◀ Previous")
        prev_btn.clicked.connect(self._prev_peak)
        nav_layout.addWidget(prev_btn)

        next_btn = QPushButton("Next ▶")
        next_btn.clicked.connect(self._next_peak)
        nav_layout.addWidget(next_btn)

        peak_layout.addLayout(nav_layout)
        peak_group.setLayout(peak_layout)

        layout.addWidget(peak_group, stretch=1)

        # Peak Info panel
        info_group = QGroupBox("Selected Peak Info")
        info_group.setStyleSheet(self._get_group_style())
        info_layout = QVBoxLayout()

        self.peak_info_label = QLabel("Select a peak to view details")
        self.peak_info_label.setWordWrap(True)
        self.peak_info_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_SMALL}px;
                color: {SECONDARY_TEXT};
                padding: {SPACING_SM}px;
                background-color: white;
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
            }}
        """)
        info_layout.addWidget(self.peak_info_label)
        info_group.setLayout(info_layout)

        layout.addWidget(info_group)

        panel.setLayout(layout)
        return panel

    def _create_contour_controls(self) -> QGroupBox:
        """Create contour control bar."""
        group = QGroupBox("Contour Settings")
        group.setStyleSheet(self._get_group_style())

        layout = QHBoxLayout()
        layout.setSpacing(SPACING_MD)

        # Levels
        layout.addWidget(QLabel("Levels:"))
        self.levels_spin = QSpinBox()
        self.levels_spin.setRange(5, 100)
        self.levels_spin.setValue(self.contour_levels)
        self.levels_spin.setMaximumWidth(80)
        self.levels_spin.valueChanged.connect(self._on_contour_changed)
        layout.addWidget(self.levels_spin)

        # Min level
        layout.addWidget(QLabel("Min Level:"))
        self.min_spin = QDoubleSpinBox()
        self.min_spin.setRange(0.01, 10.0)
        self.min_spin.setValue(self.contour_min)
        self.min_spin.setSingleStep(0.01)
        self.min_spin.setMaximumWidth(80)
        self.min_spin.valueChanged.connect(self._on_contour_changed)
        layout.addWidget(self.min_spin)

        # Increment
        layout.addWidget(QLabel("Increment:"))
        self.increment_spin = QDoubleSpinBox()
        self.increment_spin.setRange(1.01, 2.0)
        self.increment_spin.setValue(self.contour_increment)
        self.increment_spin.setSingleStep(0.01)
        self.increment_spin.setMaximumWidth(80)
        self.increment_spin.valueChanged.connect(self._on_contour_changed)
        layout.addWidget(self.increment_spin)

        layout.addStretch()

        # Colormap
        layout.addWidget(QLabel("Colormap:"))
        self.colormap_combo = QComboBox()
        self.colormap_combo.addItems(['viridis', 'plasma', 'inferno', 'magma', 'coolwarm', 'Blues', 'Reds'])
        self.colormap_combo.setMaximumWidth(100)
        self.colormap_combo.currentTextChanged.connect(self._on_contour_changed)
        layout.addWidget(self.colormap_combo)

        # Refresh button
        refresh_btn = QPushButton("Refresh Plot")
        refresh_btn.clicked.connect(self._update_contour_plot)
        layout.addWidget(refresh_btn)

        group.setLayout(layout)
        return group

    def _create_button_row(self) -> QHBoxLayout:
        """Create the bottom button row."""
        layout = QHBoxLayout()

        # Export button
        export_btn = QPushButton("Export PNG")
        export_btn.setStyleSheet(self._get_secondary_button_style())
        export_btn.clicked.connect(self._export_plot)
        layout.addWidget(export_btn)

        layout.addStretch()

        # Close button
        close_btn = QPushButton("Close")
        close_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        close_btn.setStyleSheet(self._get_secondary_button_style())
        close_btn.clicked.connect(self.close)
        layout.addWidget(close_btn)

        return layout

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

    def _get_secondary_button_style(self) -> str:
        """Get standard secondary button style."""
        return f"""
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
        """

    def _load_spectrum(self):
        """Load the spectrum data."""
        if not self.spectrum_file_path or not os.path.exists(self.spectrum_file_path):
            logger.warning(f"Spectrum file not found: {self.spectrum_file_path}")
            return

        try:
            from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator

            self.integrator = EnhancedVoigtIntegrator()
            success = self.integrator.load_nmr_file(self.spectrum_file_path)

            if not success:
                logger.error(f"Failed to load spectrum: {self.spectrum_file_path}")
                return

            # Load fitted peaks from result data
            # v0.9 pattern: check integration_results first, then fitted_peaks
            self.fitted_peaks = self.result_data.get('integration_results', [])
            if not self.fitted_peaks:
                self.fitted_peaks = self.result_data.get('fitted_peaks', [])
            if self.fitted_peaks is None:
                self.fitted_peaks = []

            logger.info(f"Loaded spectrum with {len(self.fitted_peaks)} peaks")

        except Exception as e:
            logger.error(f"Error loading spectrum: {e}")

    def _update_contour_plot(self):
        """Update the main contour plot."""
        if self.integrator is None or self.integrator.nmr_data is None:
            self.contour_widget.axes.clear()
            self.contour_widget.axes.text(
                0.5, 0.5, "No spectrum data loaded",
                ha='center', va='center',
                transform=self.contour_widget.axes.transAxes,
                fontsize=12, color='gray'
            )
            self.contour_widget.refresh()
            return

        try:
            data = self.integrator.nmr_data
            ppm_x = self.integrator.ppm_x_axis
            ppm_y = self.integrator.ppm_y_axis

            # Clear and redraw
            self.contour_widget.clear()

            # Get contour parameters
            num_levels = self.levels_spin.value()
            min_level = self.min_spin.value()
            cmap = self.colormap_combo.currentText()

            # Calculate contour levels
            data_abs = np.abs(data)
            max_intensity = np.max(data_abs)
            min_intensity = min_level * max_intensity

            level_values = np.geomspace(min_intensity, max_intensity * 0.8, num_levels)

            # Plot contours
            self.contour_widget.axes.contour(
                ppm_x, ppm_y, data,
                levels=level_values,
                cmap=cmap,
                linewidths=0.5
            )

            # Mark peaks
            if self.fitted_peaks:
                for i, peak in enumerate(self.fitted_peaks):
                    pos_x = peak.get('pos_f2', peak.get('position_x', 0))
                    pos_y = peak.get('pos_f1', peak.get('position_y', 0))

                    color = 'red' if i == self.selected_peak_idx else 'blue'
                    size = 50 if i == self.selected_peak_idx else 30

                    self.contour_widget.axes.scatter(
                        [pos_x], [pos_y],
                        c=color, s=size, marker='x', linewidths=2
                    )

            # Configure axes
            self.contour_widget.axes.set_xlabel('¹H (ppm)', fontsize=10)
            self.contour_widget.axes.set_ylabel('¹⁵N/¹³C (ppm)', fontsize=10)
            self.contour_widget.axes.invert_xaxis()
            self.contour_widget.axes.invert_yaxis()
            self.contour_widget.axes.set_title(self.spectrum_name, fontsize=11)

            self.contour_widget.refresh()

        except Exception as e:
            logger.error(f"Error updating contour plot: {e}")

    def _populate_peak_list(self):
        """Populate the peak navigator list."""
        self.peak_list.clear()

        for i, peak in enumerate(self.fitted_peaks):
            pos_x = peak.get('pos_f2', peak.get('position_x', 0))
            pos_y = peak.get('pos_f1', peak.get('position_y', 0))
            r2 = peak.get('r_squared', peak.get('fit_quality', 0))

            label = f"Peak {i+1}: ({pos_x:.3f}, {pos_y:.2f}) R²={r2:.3f}"
            item = QListWidgetItem(label)
            item.setData(Qt.UserRole, i)
            self.peak_list.addItem(item)

    def _on_peak_selected(self, item: QListWidgetItem):
        """Handle peak selection from list."""
        peak_idx = item.data(Qt.UserRole)
        self.selected_peak_idx = peak_idx
        self._show_peak_info(peak_idx)
        self._update_cross_sections(peak_idx)
        self._update_contour_plot()  # Refresh to highlight selected peak
        self._update_voigt_analysis(peak_idx)  # Update Voigt analysis plots
        self.peak_selected.emit(peak_idx)

    def _show_peak_info(self, peak_idx: int):
        """Display info for selected peak."""
        if peak_idx < 0 or peak_idx >= len(self.fitted_peaks):
            self.peak_info_label.setText("Invalid peak index")
            return

        peak = self.fitted_peaks[peak_idx]

        info = f"Peak {peak_idx + 1}\n"
        info += "-" * 20 + "\n"
        info += f"Position F2 (¹H): {peak.get('pos_f2', 'N/A'):.4f} ppm\n"
        info += f"Position F1: {peak.get('pos_f1', 'N/A'):.2f} ppm\n"
        info += f"Intensity: {peak.get('intensity', 'N/A'):.2e}\n"
        info += f"Volume: {peak.get('volume', 'N/A'):.2e}\n"
        info += f"R²: {peak.get('r_squared', 'N/A'):.4f}\n"
        info += f"SNR: {peak.get('snr', 'N/A'):.2f}\n"

        if 'sigma_x' in peak:
            info += f"σx: {peak['sigma_x']:.4f}\n"
        if 'sigma_y' in peak:
            info += f"σy: {peak['sigma_y']:.4f}\n"

        self.peak_info_label.setText(info)

    def _update_cross_sections(self, peak_idx: int):
        """Update cross-section plots for selected peak."""
        if self.integrator is None or peak_idx < 0 or peak_idx >= len(self.fitted_peaks):
            return

        peak = self.fitted_peaks[peak_idx]
        pos_x = peak.get('pos_f2', peak.get('position_x', 0))
        pos_y = peak.get('pos_f1', peak.get('position_y', 0))

        try:
            data = self.integrator.nmr_data
            ppm_x = self.integrator.ppm_x_axis
            ppm_y = self.integrator.ppm_y_axis

            # Find indices for cross-sections
            idx_x = np.argmin(np.abs(ppm_x - pos_x))
            idx_y = np.argmin(np.abs(ppm_y - pos_y))

            # F2 cross-section (horizontal, at selected F1)
            self.f2_widget.clear()
            f2_slice = data[idx_y, :]
            self.f2_widget.axes.plot(ppm_x, f2_slice, 'b-', linewidth=1)
            self.f2_widget.axes.axvline(pos_x, color='red', linestyle='--', alpha=0.7)
            self.f2_widget.axes.set_xlabel('¹H (ppm)')
            self.f2_widget.axes.invert_xaxis()
            self.f2_widget.refresh()

            # F1 cross-section (vertical, at selected F2)
            self.f1_widget.clear()
            f1_slice = data[:, idx_x]
            self.f1_widget.axes.plot(ppm_y, f1_slice, 'b-', linewidth=1)
            self.f1_widget.axes.axvline(pos_y, color='red', linestyle='--', alpha=0.7)
            self.f1_widget.axes.set_xlabel('¹⁵N/¹³C (ppm)')
            self.f1_widget.axes.invert_xaxis()
            self.f1_widget.refresh()

        except Exception as e:
            logger.error(f"Error updating cross-sections: {e}")

    def _prev_peak(self):
        """Navigate to previous peak."""
        if not self.fitted_peaks:
            return

        if self.selected_peak_idx is None:
            self.selected_peak_idx = 0
        else:
            self.selected_peak_idx = (self.selected_peak_idx - 1) % len(self.fitted_peaks)

        self.peak_list.setCurrentRow(self.selected_peak_idx)
        self._show_peak_info(self.selected_peak_idx)
        self._update_cross_sections(self.selected_peak_idx)
        self._update_contour_plot()
        self._update_voigt_analysis(self.selected_peak_idx)

    def _next_peak(self):
        """Navigate to next peak."""
        if not self.fitted_peaks:
            return

        if self.selected_peak_idx is None:
            self.selected_peak_idx = 0
        else:
            self.selected_peak_idx = (self.selected_peak_idx + 1) % len(self.fitted_peaks)

        self.peak_list.setCurrentRow(self.selected_peak_idx)
        self._show_peak_info(self.selected_peak_idx)
        self._update_cross_sections(self.selected_peak_idx)
        self._update_contour_plot()
        self._update_voigt_analysis(self.selected_peak_idx)

    def _on_contour_changed(self):
        """Handle contour parameter changes."""
        self.contour_levels = self.levels_spin.value()
        self.contour_min = self.min_spin.value()
        self.contour_increment = self.increment_spin.value()
        # Don't auto-update to avoid performance issues during dragging
        # User can click Refresh Plot

    def _export_plot(self):
        """Export the contour plot."""
        from PySide6.QtWidgets import QFileDialog

        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Export Spectrum Plot",
            f"{self.spectrum_name}.png",
            "PNG files (*.png);;All files (*.*)"
        )

        if filename:
            try:
                self.contour_widget.save_figure(filename, dpi=300)
                QMessageBox.information(
                    self,
                    "Export Successful",
                    f"Plot exported to:\n{filename}"
                )
            except Exception as e:
                QMessageBox.critical(
                    self,
                    "Export Error",
                    f"Failed to export plot:\n{str(e)}"
                )

    # ===== Voigt Analysis Methods =====

    def _update_voigt_analysis(self, peak_idx: int):
        """Update Voigt analysis plots for selected peak.

        Based on v0.9 plot_voigt_analysis() (spectrum_browser.py:1222-1274).
        """
        if peak_idx < 0 or peak_idx >= len(self.fitted_peaks):
            return

        peak = self.fitted_peaks[peak_idx]
        assignment = peak.get('assignment', peak.get('Assignment', f'Peak_{peak_idx+1}'))

        # Find stored Voigt result with full analysis data
        stored_result = self._find_stored_voigt_result(assignment)

        try:
            if stored_result and self.voigt_2d_plotter:
                # Plot 2D Voigt analysis
                self.voigt_2d_plotter.plot_voigt_analysis(stored_result)

                # Plot 3D Voigt analysis
                if self.voigt_3d_plotter:
                    self.voigt_3d_plotter.plot_voigt_analysis_3d(stored_result)

                logger.debug(f"Plotted Voigt analysis for peak {assignment}")
            else:
                # No stored result - show placeholder
                if self.voigt_2d_plotter:
                    self.voigt_2d_plotter.show_placeholder(
                        f"No Voigt fitting results available for peak {assignment}"
                    )
                if self.voigt_3d_plotter:
                    self.voigt_3d_plotter.show_placeholder(
                        f"No Voigt fitting results available for peak {assignment}"
                    )
                logger.debug(f"No stored Voigt result found for peak {assignment}")

        except Exception as e:
            logger.error(f"Error plotting Voigt analysis: {e}")
            if self.voigt_2d_plotter:
                self.voigt_2d_plotter.show_placeholder(f"Error: {str(e)[:50]}")

    def _find_stored_voigt_result(self, assignment: str, peak: Dict = None) -> Optional[Dict]:
        """Find stored Voigt fitting results for a peak by assignment.

        Based on v0.9 find_stored_voigt_result() (multi_spectrum_viewer.py:1276-1305).
        Updated to check integration_results first and accept peak directly.

        Args:
            assignment: Peak assignment string
            peak: Optional peak dict that may contain the full result

        Returns:
            Complete peak result dict with region_2d and fit data, or None
        """
        if not assignment:
            return None

        # First check if peak itself contains region_2d (it's the full result)
        if peak and peak.get('region_2d') is not None:
            logger.debug(f"Using peak directly as Voigt result for {assignment}")
            return peak

        # Check integration_results first (v0.9 pattern)
        integration_results = self.result_data.get('integration_results', [])
        if integration_results:
            for result in integration_results:
                result_assignment = result.get('assignment') or result.get('Assignment')
                if result_assignment == assignment:
                    r_squared = (result.get('r_squared') or
                                result.get('R_Squared') or
                                result.get('avg_r_squared') or 0)
                    logger.debug(f"Found integration_results Voigt result for {assignment}: R²={r_squared:.3f}")
                    return result

        # Check fitted_results from result_data (detailed Voigt analysis)
        if self.fitted_results:
            for result in self.fitted_results:
                result_assignment = result.get('assignment') or result.get('Assignment')
                if result_assignment == assignment:
                    r_squared = (result.get('r_squared') or
                                result.get('R_Squared') or
                                result.get('avg_r_squared') or 0)
                    logger.debug(f"Found fitted_results Voigt result for {assignment}: R²={r_squared:.3f}")
                    return result

        # Fallback: check integrator fitted_peaks
        if self.integrator and hasattr(self.integrator, 'fitted_peaks'):
            if self.integrator.fitted_peaks:
                for result in self.integrator.fitted_peaks:
                    if result.get('assignment') == assignment:
                        logger.debug(f"Found integrator Voigt result for {assignment}")
                        return result

        # Also check fitted_peaks in result_data
        fitted_peaks = self.result_data.get('fitted_peaks', [])
        for result in fitted_peaks:
            result_assignment = result.get('assignment') or result.get('Assignment')
            if result_assignment == assignment:
                return result

        logger.debug(f"No stored Voigt result found for {assignment}")
        return None

    def _on_3d_colormap_changed(self, scheme: str):
        """Handle color scheme change for 3D Voigt plot.

        Based on v0.9 _on_color_scheme_change_3d() (main_gui.py:829-837).
        Maps scheme names (Classic, Clean, Dark, Warm) to matplotlib colormaps.
        """
        # Map v0.9 color scheme names to matplotlib colormaps
        scheme_map = {
            'Classic': 'viridis',
            'Clean': 'RdBu_r',
            'Dark': 'inferno',
            'Warm': 'plasma'
        }
        colormap = scheme_map.get(scheme, 'viridis')

        self.color_scheme_3d = scheme
        if self.voigt_3d_plotter and hasattr(self.voigt_3d_plotter, 'current_result'):
            if self.voigt_3d_plotter.current_result:
                self.voigt_3d_plotter.plot_voigt_analysis_3d(
                    self.voigt_3d_plotter.current_result,
                    colormap=colormap
                )
