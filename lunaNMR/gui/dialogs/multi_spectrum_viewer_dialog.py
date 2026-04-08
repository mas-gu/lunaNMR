# ABOUTME: Qt-based multi-spectrum overlay viewer dialog for comparing multiple NMR spectra
# ABOUTME: Port of Tkinter multi_spectrum_viewer.py to PySide6 for v1.0 Qt interface

import logging
from typing import Optional, Dict, List

from PySide6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QCheckBox,
    QGroupBox, QWidget, QColorDialog, QMessageBox,
    QSpinBox, QDoubleSpinBox, QSplitter, QListWidget, QListWidgetItem,
    QTabWidget, QComboBox, QTextEdit, QStackedWidget
)
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor, QFont, QIcon

from lunaNMR.gui.base.base_dialog import BaseDialog
from lunaNMR.gui.components.peak_navigator_table import PeakNavigatorTable
from lunaNMR.gui.styles.design_system import (
    SPACING_SM, SPACING_MD, FONT_SIZE_BODY, FONT_SIZE_SMALL,
    PRIMARY_TEXT, SECONDARY_TEXT,
    PRIMARY_BUTTON_BG, PRIMARY_BUTTON_HOVER, PRIMARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BG, SECONDARY_BUTTON_HOVER, SECONDARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BORDER, BUTTON_CORNER_RADIUS, BUTTON_HEIGHT_DIALOG,
    FRAME_BG_COLOR
)

logger = logging.getLogger(__name__)


def _to_scalar(value):
    """Convert a value to scalar, handling lists from JSON deserialization.

    JSON deserialization may convert numpy scalars to single-element lists.
    This helper extracts the scalar value.

    Args:
        value: A scalar, list, or None

    Returns:
        Scalar value (first element if list, 0 if None/empty)
    """
    if value is None:
        return 0
    if isinstance(value, (list, tuple)):
        return value[0] if len(value) > 0 else 0
    return value


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
                 file_manager=None, initial_assignment: str = None,
                 exponential_fit_data: dict = None):
        """Initialize the multi-spectrum viewer dialog.

        Args:
            parent: Parent widget
            all_results: List of result dictionaries for all spectra
            file_manager: NMRFileManager instance for loading spectra
            initial_assignment: Optional assignment to auto-select in Peak Mode
                               (e.g., "142" or "A.142.ALA.H"). If provided, the
                               dialog will start in Peak Mode with this peak selected.
            exponential_fit_data: Dict mapping residue ID to exponential fit data from DynamiXs.
                E.g., {'142': {fit_data}, '143': {fit_data}, ...}
                Each fit_data contains: time_points, intensities, fit_curve, t_value, t_error, etc.
        """
        super().__init__(
            parent=parent,
            title="Multi-Spectrum Overlay Viewer",
            default_size=(1200, 800),
            modal=False  # Non-modal to allow interaction
        )

        self.all_results = all_results or []
        self.file_manager = file_manager
        self._initial_assignment = initial_assignment  # Store for post-init use
        self._exponential_fit_data = exponential_fit_data  # T1/T2 fit from DynamiXs

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

        # Peak Mode navigation state
        self.navigation_mode = 'spectrum'  # 'spectrum' or 'peak'
        self.master_peak_list = []  # All unique peaks by assignment
        self.peak_spectrum_map = {}  # assignment -> list of spectrum indices
        self.selected_master_peak_idx = None  # Currently selected peak in Peak Mode

        # Peak Mode axis locking (for comparing intensity/position across spectra)
        self.lock_z_scale_to_reference = True  # ON by default
        self.reference_z_min = None  # Captured from first spectrum
        self.reference_z_max = None
        self.lock_xy_to_reference = True  # ON by default
        self.reference_x_min = None  # F2 (1H) ppm range
        self.reference_x_max = None
        self.reference_y_min = None  # F1 (15N) ppm range
        self.reference_y_max = None

        # Initialize spectra data
        self._initialize_spectra()

        # Build peak-spectrum mapping for Peak Mode
        self._build_peak_spectrum_map()

        # Build UI
        self.setup_ui()

        # Center on parent
        if parent:
            self.center_on_parent()
        else:
            self.center_on_screen()

        # Handle initial assignment (for DynamiXs → lunaNMR drill-down)
        if self._initial_assignment:
            self._select_initial_assignment(self._initial_assignment)

        logger.debug(f"MultiSpectrumViewerDialog initialized with {len(self.spectra)} spectra")

    def _initialize_spectra(self):
        """Initialize spectrum data from results, sorted by delay value."""
        from lunaNMR.gui.components.intensity_decay_widget import extract_delay_from_spectrum_name

        # First pass: collect all spectrum data
        unsorted_spectra = []
        for idx, result in enumerate(self.all_results):
            spectrum_name = result.get('spectrum_name', f'spectrum_{idx+1:03d}')
            file_path = result.get('spectrum_file', '')

            # Extract fitted peaks from result_data
            # v0.9 multi_spectrum_viewer.py:117-125 - check integration_results first
            fitted_peaks = result.get('integration_results', [])
            if not fitted_peaks:
                fitted_peaks = result.get('fitted_peaks', [])
            if fitted_peaks is None:
                fitted_peaks = []

            # Extract delay for sorting
            delay_ms = extract_delay_from_spectrum_name(spectrum_name)

            unsorted_spectra.append({
                'name': spectrum_name,
                'file_path': file_path,
                'result_data': result,
                'fitted_peaks': fitted_peaks,
                'loaded': False,
                'data': None,
                '_delay_ms': delay_ms  # Temporary key for sorting
            })

        # Sort by delay (None values go to end)
        def sort_key(spec):
            delay = spec.get('_delay_ms')
            if delay is None:
                return (1, spec['name'])  # Sort unparseable names alphabetically at end
            return (0, delay)

        unsorted_spectra.sort(key=sort_key)

        # Second pass: assign colors and visibility after sorting
        for idx, spec in enumerate(unsorted_spectra):
            spec['color'] = self.default_colors[idx % len(self.default_colors)]
            spec['visible'] = (idx == 0)  # First spectrum visible by default
            del spec['_delay_ms']  # Remove temporary key
            self.spectra.append(spec)

            if idx == 0:
                self.visible_spectra.add(spec['name'])

            logger.debug(f"Initialized spectrum {spec['name']} with {len(spec['fitted_peaks'])} fitted peaks")

    def _build_peak_spectrum_map(self):
        """Build mapping from peak assignments to spectrum indices.

        Creates:
            - master_peak_list: List of unique peaks by assignment (from first spectrum)
            - peak_spectrum_map: Dict mapping assignment -> list of spectrum indices where peak exists

        Coordinates for master_peak_list come from first spectrum's peaks.
        R² is left blank in master list (shown as "—").
        """
        self.peak_spectrum_map = {}
        self.master_peak_list = []
        reference_peaks = {}  # assignment -> peak data from first spectrum

        # First pass: collect all assignments and their spectrum indices
        for spec_idx, spec in enumerate(self.spectra):
            fitted_peaks = spec.get('fitted_peaks', [])
            for peak in fitted_peaks:
                assignment = peak.get('assignment', peak.get('Assignment', ''))
                if not assignment:
                    continue

                # Track which spectra have this peak
                if assignment not in self.peak_spectrum_map:
                    self.peak_spectrum_map[assignment] = []
                self.peak_spectrum_map[assignment].append(spec_idx)

                # Use first spectrum's coordinates as reference
                if assignment not in reference_peaks:
                    reference_peaks[assignment] = peak

        # Build master_peak_list from reference peaks
        # Use first spectrum's order if possible
        if self.spectra and self.spectra[0].get('fitted_peaks'):
            first_spectrum_peaks = self.spectra[0]['fitted_peaks']
            seen_assignments = set()

            for peak in first_spectrum_peaks:
                assignment = peak.get('assignment', peak.get('Assignment', ''))
                if assignment and assignment not in seen_assignments:
                    seen_assignments.add(assignment)
                    # Get coordinates using _to_scalar for JSON deserialization safety
                    cx = _to_scalar(peak.get('center_x') or peak.get('peak_x') or peak.get('ppm_x') or peak.get('pos_f2') or 0)
                    cy = _to_scalar(peak.get('center_y') or peak.get('peak_y') or peak.get('ppm_y') or peak.get('pos_f1') or 0)
                    # Create master peak entry (R² left blank)
                    master_peak = {
                        'assignment': assignment,
                        'center_x': cx,
                        'center_y': cy,
                        'r_squared': None,  # Leave blank for master list
                        'quality': ''
                    }
                    self.master_peak_list.append(master_peak)

            # Add any peaks from other spectra not in first spectrum
            for assignment, peak in reference_peaks.items():
                if assignment not in seen_assignments:
                    cx = _to_scalar(peak.get('center_x') or peak.get('peak_x') or peak.get('ppm_x') or peak.get('pos_f2') or 0)
                    cy = _to_scalar(peak.get('center_y') or peak.get('peak_y') or peak.get('ppm_y') or peak.get('pos_f1') or 0)
                    master_peak = {
                        'assignment': assignment,
                        'center_x': cx,
                        'center_y': cy,
                        'r_squared': None,
                        'quality': ''
                    }
                    self.master_peak_list.append(master_peak)

        logger.debug(f"Built peak_spectrum_map with {len(self.master_peak_list)} unique peaks")

    def _select_initial_assignment(self, assignment: str):
        """Select a peak by assignment and switch to Peak Mode in 3D Voigt tab.

        Used for DynamiXs → lunaNMR drill-down: when user clicks "Inspect Peak"
        in FitViewer, this opens the MultiSpectrumViewer with the corresponding
        peak already selected.

        Args:
            assignment: Peak assignment (e.g., "142", "A.8.LEU.H") to select.
                       Matches against assignments in master_peak_list using:
                       1. Exact match
                       2. Substring match
                       3. Residue number extraction and match
        """
        import re

        if not self.master_peak_list:
            logger.warning(f"Cannot select assignment '{assignment}': no peaks loaded")
            return

        # Find matching peak index
        target_idx = None
        assignment_clean = str(assignment).strip()

        # Extract residue number from assignment (e.g., "A.8.LEU.H" -> "8")
        def extract_residue_number(s):
            """Extract the first number from an assignment string."""
            numbers = re.findall(r'\d+', str(s))
            return numbers[0] if numbers else None

        search_resnum = extract_residue_number(assignment_clean)

        # First try exact match
        for idx, peak in enumerate(self.master_peak_list):
            peak_assignment = peak.get('assignment', '')
            if peak_assignment == assignment_clean:
                target_idx = idx
                logger.debug(f"Exact match: '{assignment_clean}' == '{peak_assignment}'")
                break

        # If no exact match, try partial match (assignment as substring)
        if target_idx is None:
            for idx, peak in enumerate(self.master_peak_list):
                peak_assignment = peak.get('assignment', '')
                # Check if search term is in peak assignment
                if assignment_clean in peak_assignment:
                    target_idx = idx
                    logger.debug(f"Substring match: '{assignment_clean}' in '{peak_assignment}'")
                    break
                # Also check if peak assignment is in the search string
                if peak_assignment in assignment_clean:
                    target_idx = idx
                    logger.debug(f"Reverse substring match: '{peak_assignment}' in '{assignment_clean}'")
                    break

        # If still no match, try matching by residue number
        if target_idx is None and search_resnum:
            for idx, peak in enumerate(self.master_peak_list):
                peak_assignment = peak.get('assignment', '')
                peak_resnum = extract_residue_number(peak_assignment)
                if peak_resnum and peak_resnum == search_resnum:
                    target_idx = idx
                    logger.debug(f"Residue number match: {search_resnum} in '{peak_assignment}'")
                    break

        if target_idx is None:
            logger.warning(f"Assignment '{assignment}' not found in master_peak_list")
            # Log available assignments for debugging
            available = [p.get('assignment', '') for p in self.master_peak_list[:10]]
            logger.debug(f"Available assignments (first 10): {available}")
            return

        # Switch to 3D Voigt Analysis tab (index 2)
        self.plot_tabs.setCurrentIndex(2)

        # Enable Peak Mode
        self.voigt_3d_mode_toggle_btn.setChecked(True)
        self._toggle_voigt_3d_navigation_mode()

        # Select the peak in the peak table
        self.voigt_3d_peak_mode_peak_table.select_peak(target_idx)

        logger.debug(f"Selected initial assignment '{assignment}' (index {target_idx}) in Peak Mode")

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
        """Setup the dialog user interface.

        Tabs contain their own layout - no global left panel.
        Spectrum list is now inside Overlay View tab only.
        """
        layout = QVBoxLayout()
        layout.setSpacing(SPACING_SM)
        layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)

        # Plot area with tabs (full width - no global left panel)
        plot_panel = self.create_plot_panel()
        layout.addWidget(plot_panel, stretch=1)

        self.setLayout(layout)

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
        """Create the overlay view tab with spectrum list, contour plot, and peak navigator.

        Layout: [Spectra + Contour] | [Plot] | [Peak Navigator]

        The spectrum list with visibility toggles and contour controls is now
        integrated into this tab (moved from global left panel).
        """
        from lunaNMR.gui.components.matplotlib_widget import MatplotlibWidget
        from lunaNMR.gui.components.nmr_navigation_handler import NMRNavigationHandler
        from PySide6.QtWidgets import QScrollArea

        tab = QWidget()
        layout = QHBoxLayout()
        layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)
        layout.setSpacing(SPACING_SM)

        # === Left panel: Spectrum list + Contour controls ===
        left_panel = QWidget()
        left_layout = QVBoxLayout()
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.setSpacing(SPACING_SM)

        # Spectra group with visibility toggles
        spectra_group = QGroupBox("Spectra")
        spectra_group.setStyleSheet(self._get_group_style())
        spectra_layout = QVBoxLayout()
        spectra_layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)
        spectra_layout.setSpacing(2)

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

        # Add spectrum items with visibility toggles
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
        left_layout.addWidget(spectra_group, stretch=1)

        # Contour controls group
        contour_group = QGroupBox("Contour")
        contour_group.setStyleSheet(self._get_group_style())
        contour_layout = QVBoxLayout()
        contour_layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)
        contour_layout.setSpacing(2)

        levels_row = QHBoxLayout()
        levels_row.addWidget(QLabel("Levels:"))
        self.levels_spin = QSpinBox()
        self.levels_spin.setRange(5, 100)
        self.levels_spin.setValue(10)
        levels_row.addWidget(self.levels_spin)
        contour_layout.addLayout(levels_row)

        min_row = QHBoxLayout()
        min_row.addWidget(QLabel("Min:"))
        self.min_spin = QDoubleSpinBox()
        self.min_spin.setRange(0.01, 10.0)
        self.min_spin.setValue(0.2)
        self.min_spin.setSingleStep(0.01)
        min_row.addWidget(self.min_spin)
        contour_layout.addLayout(min_row)

        inc_row = QHBoxLayout()
        inc_row.addWidget(QLabel("Inc:"))
        self.increment_spin = QDoubleSpinBox()
        self.increment_spin.setRange(1.01, 2.0)
        self.increment_spin.setValue(2.0)
        self.increment_spin.setSingleStep(0.01)
        inc_row.addWidget(self.increment_spin)
        contour_layout.addLayout(inc_row)

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
        left_layout.addWidget(contour_group)

        left_panel.setLayout(left_layout)
        layout.addWidget(left_panel, stretch=2)

        # === Center panel: Contour plot ===
        center_panel = QWidget()
        center_layout = QVBoxLayout()
        center_layout.setContentsMargins(0, 0, 0, 0)

        self.plot_widget = MatplotlibWidget(
            parent=center_panel,
            toolbar=False,
            figsize=(8, 6),
            dpi=100
        )

        self._nav_handler = NMRNavigationHandler()
        self._nav_handler.attach(self.plot_widget)

        self.plot_widget.axes.set_title("Select spectra to display overlay")
        self.plot_widget.axes.text(
            0.5, 0.5,
            "Select spectra from the list\nto display contour overlay",
            ha='center', va='center',
            transform=self.plot_widget.axes.transAxes,
            fontsize=12, color='gray'
        )
        self.plot_widget.refresh()

        center_layout.addWidget(self.plot_widget, stretch=1)
        center_panel.setLayout(center_layout)
        layout.addWidget(center_panel, stretch=6)

        # === Right panel: Spectrum selector and peak navigator ===
        right_panel = QWidget()
        right_layout = QVBoxLayout()
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.setSpacing(SPACING_SM)

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

        peak_group = QGroupBox("Peak Navigator")
        peak_group.setStyleSheet(self._get_group_style())
        peak_layout = QVBoxLayout()

        self.overlay_peak_table = PeakNavigatorTable()
        self.overlay_peak_table.peak_selected.connect(self._on_overlay_peak_selected_by_index)
        peak_layout.addWidget(self.overlay_peak_table)

        self.show_peaks_checkbox = QCheckBox("Show Peak Markers")
        self.show_peaks_checkbox.setChecked(True)
        self.show_peaks_checkbox.toggled.connect(self._on_toggle_peak_markers)
        peak_layout.addWidget(self.show_peaks_checkbox)

        peak_group.setLayout(peak_layout)
        right_layout.addWidget(peak_group, stretch=1)

        right_panel.setLayout(right_layout)
        layout.addWidget(right_panel, stretch=2)

        # Track selected spectrum for overlay peak display
        self.overlay_selected_spectrum_idx = 0
        self.overlay_show_peaks = True

        tab.setLayout(layout)
        return tab

    def _create_voigt_analysis_tab(self) -> QWidget:
        """Create the Voigt analysis tab with spectrum/peak selection and VoigtAnalysisPlotter.

        Based on v0.9 setup_voigt_analysis_tab() (multi_spectrum_viewer.py:819-902).
        Includes mode toggle for Spectrum Mode and Peak Mode navigation.
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

        # Right panel (30%): Mode toggle + swappable selector panels
        right_panel = QWidget()
        right_layout = QVBoxLayout()
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.setSpacing(SPACING_SM)

        # Mode toggle button at top
        self.voigt_mode_toggle_btn = QPushButton("Spectrum Mode")
        self.voigt_mode_toggle_btn.setCheckable(True)
        self.voigt_mode_toggle_btn.setChecked(False)  # Start in Spectrum Mode
        self.voigt_mode_toggle_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {SECONDARY_BUTTON_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px {SPACING_MD}px;
                font-size: {FONT_SIZE_BODY}px;
                font-weight: bold;
            }}
            QPushButton:checked {{
                background-color: {PRIMARY_BUTTON_BG};
                color: {PRIMARY_BUTTON_TEXT};
                border: none;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
            QPushButton:checked:hover {{
                background-color: {PRIMARY_BUTTON_HOVER};
            }}
        """)
        self.voigt_mode_toggle_btn.clicked.connect(self._toggle_voigt_navigation_mode)
        right_layout.addWidget(self.voigt_mode_toggle_btn)

        # QStackedWidget for swappable panels
        self.voigt_mode_stack = QStackedWidget()

        # === Page 0: Spectrum Mode (Select Spectrum on top, Peak Navigator below) ===
        spectrum_mode_page = QWidget()
        spectrum_mode_layout = QVBoxLayout()
        spectrum_mode_layout.setContentsMargins(0, 0, 0, 0)
        spectrum_mode_layout.setSpacing(SPACING_SM)

        # Spectrum selector for Spectrum Mode
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
        spectrum_mode_layout.addWidget(spec_group)

        # Peak navigator for Spectrum Mode
        peak_group = QGroupBox("Peak Navigator")
        peak_group.setStyleSheet(self._get_group_style())
        peak_layout = QVBoxLayout()

        self.voigt_peak_table = PeakNavigatorTable()
        self.voigt_peak_table.peak_selected.connect(self._on_voigt_peak_selected_by_index)
        peak_layout.addWidget(self.voigt_peak_table)

        peak_group.setLayout(peak_layout)
        spectrum_mode_layout.addWidget(peak_group, stretch=1)

        spectrum_mode_page.setLayout(spectrum_mode_layout)
        self.voigt_mode_stack.addWidget(spectrum_mode_page)

        # === Page 1: Peak Mode (Select Peak on top, Spectrum Navigator below) ===
        peak_mode_page = QWidget()
        peak_mode_layout = QVBoxLayout()
        peak_mode_layout.setContentsMargins(0, 0, 0, 0)
        peak_mode_layout.setSpacing(SPACING_SM)

        # Peak selector for Peak Mode (master peak list)
        peak_select_group = QGroupBox("Select Peak")
        peak_select_group.setStyleSheet(self._get_group_style())
        peak_select_layout = QVBoxLayout()

        self.voigt_peak_mode_peak_table = PeakNavigatorTable()
        self.voigt_peak_mode_peak_table.peak_selected.connect(self._on_voigt_peak_mode_peak_selected)
        self.voigt_peak_mode_peak_table.load_peaks(self.master_peak_list)
        peak_select_layout.addWidget(self.voigt_peak_mode_peak_table)

        peak_select_group.setLayout(peak_select_layout)
        peak_mode_layout.addWidget(peak_select_group, stretch=1)

        # Spectrum navigator for Peak Mode
        spectrum_nav_group = QGroupBox("Spectrum Navigator")
        spectrum_nav_group.setStyleSheet(self._get_group_style())
        spectrum_nav_layout = QVBoxLayout()

        self.voigt_peak_mode_spectrum_list = QListWidget()
        self.voigt_peak_mode_spectrum_list.setMaximumHeight(120)
        # Connect both click and keyboard navigation (arrow keys)
        self.voigt_peak_mode_spectrum_list.itemClicked.connect(self._on_voigt_peak_mode_spectrum_selected)
        self.voigt_peak_mode_spectrum_list.currentRowChanged.connect(self._on_voigt_peak_mode_spectrum_row_changed)
        spectrum_nav_layout.addWidget(self.voigt_peak_mode_spectrum_list)

        # Navigation buttons for Peak Mode spectrum navigation
        nav_btn_layout = QHBoxLayout()
        nav_btn_layout.setSpacing(SPACING_SM)

        self.voigt_peak_mode_prev_btn = QPushButton("◀ Prev")
        self.voigt_peak_mode_prev_btn.setStyleSheet(self._get_nav_button_style())
        self.voigt_peak_mode_prev_btn.clicked.connect(self._on_voigt_peak_mode_prev_spectrum)
        self.voigt_peak_mode_prev_btn.setEnabled(False)
        nav_btn_layout.addWidget(self.voigt_peak_mode_prev_btn)

        self.voigt_peak_mode_next_btn = QPushButton("Next ▶")
        self.voigt_peak_mode_next_btn.setStyleSheet(self._get_nav_button_style())
        self.voigt_peak_mode_next_btn.clicked.connect(self._on_voigt_peak_mode_next_spectrum)
        self.voigt_peak_mode_next_btn.setEnabled(False)
        nav_btn_layout.addWidget(self.voigt_peak_mode_next_btn)

        spectrum_nav_layout.addLayout(nav_btn_layout)
        spectrum_nav_group.setLayout(spectrum_nav_layout)
        peak_mode_layout.addWidget(spectrum_nav_group)

        peak_mode_page.setLayout(peak_mode_layout)
        self.voigt_mode_stack.addWidget(peak_mode_page)

        right_layout.addWidget(self.voigt_mode_stack, stretch=1)

        right_panel.setLayout(right_layout)
        layout.addWidget(right_panel, stretch=3)

        # Initialize with placeholder
        self.voigt_2d_plotter.show_placeholder("Select spectrum and peak to view Voigt analysis")

        tab.setLayout(layout)
        return tab

    def _create_voigt_3d_analysis_tab(self) -> QWidget:
        """Create the 3D Voigt analysis tab with VoigtAnalysisPlotter and IntensityDecayWidget.

        In Peak Mode, shows side-by-side: 3D Voigt surface + Intensity vs Delay plot.
        In Spectrum Mode, shows only the 3D Voigt surface (full width).

        Based on v0.9 main_gui.py Tab 3: 3D Voigt Analysis (lines 1829-1944).
        """
        from lunaNMR.gui.components.voigt_analysis_plotter import VoigtAnalysisPlotter
        from lunaNMR.gui.components.intensity_decay_widget import IntensityDecayWidget

        tab = QWidget()
        layout = QHBoxLayout()
        layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)
        layout.setSpacing(SPACING_SM)

        # Left panel: Controls + 3D plot + IntensityDecay (in Peak Mode)
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

        # === Plot area with splitter for 3D Voigt + Intensity Decay ===
        plot_splitter = QSplitter(Qt.Horizontal)

        # Left side: 3D Voigt Surface Analysis
        plot_group = QGroupBox("3D Voigt Surface Analysis")
        plot_group.setStyleSheet(self._get_group_style())
        plot_layout = QVBoxLayout()

        # Create VoigtAnalysisPlotter for 3D
        self.voigt_3d_plotter = VoigtAnalysisPlotter(parent=plot_group, toolbar=True)
        plot_layout.addWidget(self.voigt_3d_plotter)

        plot_group.setLayout(plot_layout)
        plot_splitter.addWidget(plot_group)

        # Right side: Intensity Decay Widget (shown in Peak Mode only)
        self.intensity_decay_widget = IntensityDecayWidget(parent=tab)
        self.intensity_decay_widget.point_clicked.connect(self._on_intensity_decay_point_clicked)
        self.intensity_decay_widget.setVisible(False)  # Hidden by default (Spectrum Mode)
        plot_splitter.addWidget(self.intensity_decay_widget)

        # Set initial splitter sizes (3D plot takes full width when decay is hidden)
        plot_splitter.setSizes([1000, 0])
        self._plot_splitter_3d = plot_splitter  # Store reference for mode toggle

        left_layout.addWidget(plot_splitter, stretch=1)

        left_panel.setLayout(left_layout)
        layout.addWidget(left_panel, stretch=7)

        # Right panel (30%): Mode toggle + swappable selector panels
        right_panel = QWidget()
        right_layout = QVBoxLayout()
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.setSpacing(SPACING_SM)

        # Mode toggle button at top
        self.voigt_3d_mode_toggle_btn = QPushButton("Spectrum Mode")
        self.voigt_3d_mode_toggle_btn.setCheckable(True)
        self.voigt_3d_mode_toggle_btn.setChecked(False)  # Start in Spectrum Mode
        self.voigt_3d_mode_toggle_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {SECONDARY_BUTTON_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px {SPACING_MD}px;
                font-size: {FONT_SIZE_BODY}px;
                font-weight: bold;
            }}
            QPushButton:checked {{
                background-color: {PRIMARY_BUTTON_BG};
                color: {PRIMARY_BUTTON_TEXT};
                border: none;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
            QPushButton:checked:hover {{
                background-color: {PRIMARY_BUTTON_HOVER};
            }}
        """)
        self.voigt_3d_mode_toggle_btn.clicked.connect(self._toggle_voigt_3d_navigation_mode)
        right_layout.addWidget(self.voigt_3d_mode_toggle_btn)

        # QStackedWidget for swappable panels
        self.voigt_3d_mode_stack = QStackedWidget()

        # === Page 0: Spectrum Mode (Select Spectrum on top, Peak Navigator below) ===
        spectrum_mode_page = QWidget()
        spectrum_mode_layout = QVBoxLayout()
        spectrum_mode_layout.setContentsMargins(0, 0, 0, 0)
        spectrum_mode_layout.setSpacing(SPACING_SM)

        # Spectrum selector for Spectrum Mode
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
        spectrum_mode_layout.addWidget(spec_group)

        # Peak navigator for Spectrum Mode
        peak_group = QGroupBox("Peak Navigator")
        peak_group.setStyleSheet(self._get_group_style())
        peak_layout = QVBoxLayout()

        self.voigt_3d_peak_table = PeakNavigatorTable()
        self.voigt_3d_peak_table.peak_selected.connect(self._on_voigt_3d_peak_selected_by_index)
        peak_layout.addWidget(self.voigt_3d_peak_table)

        peak_group.setLayout(peak_layout)
        spectrum_mode_layout.addWidget(peak_group, stretch=1)

        spectrum_mode_page.setLayout(spectrum_mode_layout)
        self.voigt_3d_mode_stack.addWidget(spectrum_mode_page)

        # === Page 1: Peak Mode (Select Peak on top, Spectrum Navigator below) ===
        peak_mode_page = QWidget()
        peak_mode_layout = QVBoxLayout()
        peak_mode_layout.setContentsMargins(0, 0, 0, 0)
        peak_mode_layout.setSpacing(SPACING_SM)

        # Peak selector for Peak Mode (master peak list)
        peak_select_group = QGroupBox("Select Peak")
        peak_select_group.setStyleSheet(self._get_group_style())
        peak_select_layout = QVBoxLayout()

        self.voigt_3d_peak_mode_peak_table = PeakNavigatorTable()
        self.voigt_3d_peak_mode_peak_table.peak_selected.connect(self._on_voigt_3d_peak_mode_peak_selected)
        self.voigt_3d_peak_mode_peak_table.load_peaks(self.master_peak_list)
        peak_select_layout.addWidget(self.voigt_3d_peak_mode_peak_table)

        peak_select_group.setLayout(peak_select_layout)
        peak_mode_layout.addWidget(peak_select_group, stretch=1)

        # Spectrum navigator for Peak Mode
        spectrum_nav_group = QGroupBox("Spectrum Navigator")
        spectrum_nav_group.setStyleSheet(self._get_group_style())
        spectrum_nav_layout = QVBoxLayout()

        self.voigt_3d_peak_mode_spectrum_list = QListWidget()
        self.voigt_3d_peak_mode_spectrum_list.setMaximumHeight(120)
        # Connect both click and keyboard navigation (arrow keys)
        self.voigt_3d_peak_mode_spectrum_list.itemClicked.connect(self._on_voigt_3d_peak_mode_spectrum_selected)
        self.voigt_3d_peak_mode_spectrum_list.currentRowChanged.connect(self._on_voigt_3d_peak_mode_spectrum_row_changed)
        spectrum_nav_layout.addWidget(self.voigt_3d_peak_mode_spectrum_list)

        # Navigation buttons for Peak Mode spectrum navigation
        nav_btn_layout = QHBoxLayout()
        nav_btn_layout.setSpacing(SPACING_SM)

        self.voigt_3d_peak_mode_prev_btn = QPushButton("◀ Prev")
        self.voigt_3d_peak_mode_prev_btn.setStyleSheet(self._get_nav_button_style())
        self.voigt_3d_peak_mode_prev_btn.clicked.connect(self._on_voigt_3d_peak_mode_prev_spectrum)
        self.voigt_3d_peak_mode_prev_btn.setEnabled(False)
        nav_btn_layout.addWidget(self.voigt_3d_peak_mode_prev_btn)

        self.voigt_3d_peak_mode_next_btn = QPushButton("Next ▶")
        self.voigt_3d_peak_mode_next_btn.setStyleSheet(self._get_nav_button_style())
        self.voigt_3d_peak_mode_next_btn.clicked.connect(self._on_voigt_3d_peak_mode_next_spectrum)
        self.voigt_3d_peak_mode_next_btn.setEnabled(False)
        nav_btn_layout.addWidget(self.voigt_3d_peak_mode_next_btn)

        spectrum_nav_layout.addLayout(nav_btn_layout)

        # Axis lock checkboxes for comparing across spectra
        lock_layout = QHBoxLayout()
        lock_layout.setSpacing(SPACING_SM)

        self.lock_z_scale_checkbox = QCheckBox("Lock Z")
        self.lock_z_scale_checkbox.setChecked(True)
        self.lock_z_scale_checkbox.setToolTip("Lock intensity axis to first spectrum's scale")
        self.lock_z_scale_checkbox.toggled.connect(self._on_lock_z_scale_toggled)
        lock_layout.addWidget(self.lock_z_scale_checkbox)

        self.lock_xy_checkbox = QCheckBox("Lock X/Y")
        self.lock_xy_checkbox.setChecked(True)
        self.lock_xy_checkbox.setToolTip("Lock chemical shift axes to first spectrum's range")
        self.lock_xy_checkbox.toggled.connect(self._on_lock_xy_toggled)
        lock_layout.addWidget(self.lock_xy_checkbox)

        lock_layout.addStretch()
        spectrum_nav_layout.addLayout(lock_layout)

        spectrum_nav_group.setLayout(spectrum_nav_layout)
        peak_mode_layout.addWidget(spectrum_nav_group)

        peak_mode_page.setLayout(peak_mode_layout)
        self.voigt_3d_mode_stack.addWidget(peak_mode_page)

        right_layout.addWidget(self.voigt_3d_mode_stack, stretch=1)

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

                # Linewidths
                lines.append(f"  Linewidth F2 (1H):     {lw_f2:.4f} ppm")
                lines.append(f"  Linewidth F1 (15N):    {lw_f1:.3f} ppm")

                if i < len(all_peaks) - 1:
                    lines.append("")
        else:
            lines.append("No peak details available.")

        # Update text widget
        self.peak_params_text.setPlainText("\n".join(lines))
        logger.debug(f"Updated Peak Parameters tab for {assignment}")

    def _create_color_icon(self, color: str, size: int = 12) -> QIcon:
        """Create a small colored square icon for list items.

        Args:
            color: Color as hex string (e.g., '#1f77b4')
            size: Icon size in pixels

        Returns:
            QIcon with colored square
        """
        from PySide6.QtGui import QPixmap, QPainter, QBrush, QPen

        pixmap = QPixmap(size, size)
        pixmap.fill(QColor('transparent'))

        painter = QPainter(pixmap)
        painter.setRenderHint(QPainter.Antialiasing)
        painter.setBrush(QBrush(QColor(color)))
        painter.setPen(QPen(QColor('#333333'), 1))
        painter.drawRect(1, 1, size - 2, size - 2)
        painter.end()

        return QIcon(pixmap)

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

    def _get_nav_button_style(self) -> str:
        """Get navigation button style for Peak Mode."""
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
            QPushButton:disabled {{
                background-color: #f0f0f0;
                color: #aaaaaa;
                border-color: #cccccc;
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
        self.min_spin.setRange(0.0001, 10.0)
        self.min_spin.setValue(0.2)
        self.min_spin.setSingleStep(0.002)
        self.min_spin.setDecimals(5)
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
                # Try multiple possible keys for coordinates, using _to_scalar for JSON safety
                x_ppm = _to_scalar(peak.get('center_x') or peak.get('peak_x') or peak.get('ppm_x') or
                        peak.get('pos_f2') or peak.get('Position_X') or 0)
                y_ppm = _to_scalar(peak.get('center_y') or peak.get('peak_y') or peak.get('ppm_y') or
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

    # ===== Voigt Tab Peak Mode Methods =====

    def _toggle_voigt_navigation_mode(self):
        """Toggle between Spectrum Mode and Peak Mode for Voigt Analysis tab."""
        if self.voigt_mode_toggle_btn.isChecked():
            # Switch to Peak Mode
            self.navigation_mode = 'peak'
            self.voigt_mode_toggle_btn.setText("Peak Mode")
            self.voigt_mode_stack.setCurrentIndex(1)
            # Reset to first peak
            if self.master_peak_list:
                self.voigt_peak_mode_peak_table.select_peak(0)
            logger.debug("Voigt tab: Switched to Peak Mode")
        else:
            # Switch to Spectrum Mode
            self.navigation_mode = 'spectrum'
            self.voigt_mode_toggle_btn.setText("Spectrum Mode")
            self.voigt_mode_stack.setCurrentIndex(0)
            # Reset to first spectrum
            if self.voigt_spectrum_list.count() > 0:
                self.voigt_spectrum_list.setCurrentRow(0)
                item = self.voigt_spectrum_list.item(0)
                if item:
                    self._on_voigt_spectrum_selected(item)
            logger.debug("Voigt tab: Switched to Spectrum Mode")

    def _on_voigt_peak_mode_peak_selected(self, peak_idx: int):
        """Handle peak selection in Peak Mode for Voigt tab."""
        if peak_idx < 0 or peak_idx >= len(self.master_peak_list):
            return

        self.selected_master_peak_idx = peak_idx
        peak = self.master_peak_list[peak_idx]
        assignment = peak.get('assignment', '')

        # Get list of spectra that have this peak
        spectra_with_peak = self.peak_spectrum_map.get(assignment, [])

        # Update spectrum navigator list with grayed-out items
        self._update_voigt_peak_mode_spectrum_list(spectra_with_peak)

        logger.debug(f"Voigt Peak Mode: Selected peak {assignment}, found in {len(spectra_with_peak)} spectra")

    def _update_voigt_peak_mode_spectrum_list(self, spectra_with_peak: List[int]):
        """Update the Peak Mode spectrum list for Voigt tab."""
        self.voigt_peak_mode_spectrum_list.clear()
        self._voigt_peak_mode_available_spectra = []

        for i, spec in enumerate(self.spectra):
            # Create colored icon for spectrum (matches overlay view colors)
            spectrum_color = spec.get('color', self.default_colors[i % len(self.default_colors)])
            icon = self._create_color_icon(spectrum_color)

            item = QListWidgetItem(icon, spec['name'])
            item.setData(Qt.UserRole, i)

            if i in spectra_with_peak:
                item.setFlags(item.flags() | Qt.ItemIsSelectable | Qt.ItemIsEnabled)
                self._voigt_peak_mode_available_spectra.append(i)
            else:
                item.setFlags(item.flags() & ~Qt.ItemIsSelectable & ~Qt.ItemIsEnabled)
                item.setForeground(QColor('#aaaaaa'))

            self.voigt_peak_mode_spectrum_list.addItem(item)

        # Update navigation button states
        has_multiple = len(self._voigt_peak_mode_available_spectra) > 1
        self.voigt_peak_mode_prev_btn.setEnabled(has_multiple)
        self.voigt_peak_mode_next_btn.setEnabled(has_multiple)

        # Select first available spectrum and plot
        if self._voigt_peak_mode_available_spectra:
            first_idx = self._voigt_peak_mode_available_spectra[0]
            self.voigt_peak_mode_spectrum_list.setCurrentRow(first_idx)
            self._select_voigt_spectrum_for_peak_mode(first_idx)

    def _on_voigt_peak_mode_spectrum_selected(self, item: QListWidgetItem):
        """Handle spectrum selection in Peak Mode for Voigt tab (click)."""
        spec_idx = item.data(Qt.UserRole)

        if not hasattr(self, '_voigt_peak_mode_available_spectra'):
            return
        if spec_idx not in self._voigt_peak_mode_available_spectra:
            return

        self._select_voigt_spectrum_for_peak_mode(spec_idx)

    def _on_voigt_peak_mode_spectrum_row_changed(self, row: int):
        """Handle spectrum selection in Peak Mode for Voigt tab (keyboard navigation)."""
        if row < 0:
            return

        item = self.voigt_peak_mode_spectrum_list.item(row)
        if item is None:
            return

        spec_idx = item.data(Qt.UserRole)

        if not hasattr(self, '_voigt_peak_mode_available_spectra'):
            return
        if spec_idx not in self._voigt_peak_mode_available_spectra:
            return

        self._select_voigt_spectrum_for_peak_mode(spec_idx)

    def _select_voigt_spectrum_for_peak_mode(self, spec_idx: int):
        """Select a spectrum in Peak Mode and plot Voigt analysis."""
        self.selected_spectrum_index = spec_idx

        if self.selected_master_peak_idx is None or self.selected_master_peak_idx >= len(self.master_peak_list):
            return

        peak = self.master_peak_list[self.selected_master_peak_idx]
        assignment = peak.get('assignment', '')

        # Find the peak index in this spectrum's peak list
        fitted_peaks = self._get_fitted_peaks(spec_idx)
        for i, p in enumerate(fitted_peaks):
            p_assignment = p.get('assignment', p.get('Assignment', ''))
            if p_assignment == assignment:
                self.selected_peak_index = i
                break

        # Plot the Voigt analysis
        self._plot_voigt_analysis()
        logger.debug(f"Voigt Peak Mode: Selected spectrum {spec_idx} for peak {assignment}")

    def _on_voigt_peak_mode_prev_spectrum(self):
        """Navigate to previous available spectrum in Voigt Peak Mode."""
        if not hasattr(self, '_voigt_peak_mode_available_spectra') or not self._voigt_peak_mode_available_spectra:
            return

        current_row = self.voigt_peak_mode_spectrum_list.currentRow()
        try:
            current_pos = self._voigt_peak_mode_available_spectra.index(current_row)
            new_pos = (current_pos - 1) % len(self._voigt_peak_mode_available_spectra)
            new_idx = self._voigt_peak_mode_available_spectra[new_pos]
            self.voigt_peak_mode_spectrum_list.setCurrentRow(new_idx)
            self._select_voigt_spectrum_for_peak_mode(new_idx)
        except ValueError:
            if self._voigt_peak_mode_available_spectra:
                new_idx = self._voigt_peak_mode_available_spectra[0]
                self.voigt_peak_mode_spectrum_list.setCurrentRow(new_idx)
                self._select_voigt_spectrum_for_peak_mode(new_idx)

    def _on_voigt_peak_mode_next_spectrum(self):
        """Navigate to next available spectrum in Voigt Peak Mode."""
        if not hasattr(self, '_voigt_peak_mode_available_spectra') or not self._voigt_peak_mode_available_spectra:
            return

        current_row = self.voigt_peak_mode_spectrum_list.currentRow()
        try:
            current_pos = self._voigt_peak_mode_available_spectra.index(current_row)
            new_pos = (current_pos + 1) % len(self._voigt_peak_mode_available_spectra)
            new_idx = self._voigt_peak_mode_available_spectra[new_pos]
            self.voigt_peak_mode_spectrum_list.setCurrentRow(new_idx)
            self._select_voigt_spectrum_for_peak_mode(new_idx)
        except ValueError:
            if self._voigt_peak_mode_available_spectra:
                new_idx = self._voigt_peak_mode_available_spectra[0]
                self.voigt_peak_mode_spectrum_list.setCurrentRow(new_idx)
                self._select_voigt_spectrum_for_peak_mode(new_idx)

    # ===== 3D Voigt Tab Peak Mode Methods =====

    def _toggle_voigt_3d_navigation_mode(self):
        """Toggle between Spectrum Mode and Peak Mode for 3D Voigt Analysis tab."""
        if self.voigt_3d_mode_toggle_btn.isChecked():
            # Switch to Peak Mode
            self.navigation_mode = 'peak'
            self.voigt_3d_mode_toggle_btn.setText("Peak Mode")
            self.voigt_3d_mode_stack.setCurrentIndex(1)
            # Show intensity decay widget and adjust splitter
            if hasattr(self, 'intensity_decay_widget'):
                self.intensity_decay_widget.setVisible(True)
                self._plot_splitter_3d.setSizes([600, 400])  # 60/40 split
            # Reset to first peak
            if self.master_peak_list:
                self.voigt_3d_peak_mode_peak_table.select_peak(0)
            logger.debug("3D Voigt tab: Switched to Peak Mode")
        else:
            # Switch to Spectrum Mode
            self.navigation_mode = 'spectrum'
            self.voigt_3d_mode_toggle_btn.setText("Spectrum Mode")
            self.voigt_3d_mode_stack.setCurrentIndex(0)
            # Hide intensity decay widget
            if hasattr(self, 'intensity_decay_widget'):
                self.intensity_decay_widget.setVisible(False)
                self._plot_splitter_3d.setSizes([1000, 0])  # Full width for 3D plot
            # Clear reference scales (not used in Spectrum Mode)
            self.reference_z_min = None
            self.reference_z_max = None
            self.reference_x_min = None
            self.reference_x_max = None
            self.reference_y_min = None
            self.reference_y_max = None
            self.voigt_3d_plotter.clear_fixed_z_limits()
            self.voigt_3d_plotter.clear_fixed_xy_limits()
            # Reset to first spectrum
            if self.voigt_3d_spectrum_list.count() > 0:
                self.voigt_3d_spectrum_list.setCurrentRow(0)
                item = self.voigt_3d_spectrum_list.item(0)
                if item:
                    self._on_voigt_3d_spectrum_selected(item)
            logger.debug("3D Voigt tab: Switched to Spectrum Mode")

    def _on_voigt_3d_peak_mode_peak_selected(self, peak_idx: int):
        """Handle peak selection in Peak Mode for 3D Voigt tab."""
        if peak_idx < 0 or peak_idx >= len(self.master_peak_list):
            return

        self.selected_master_peak_idx = peak_idx
        peak = self.master_peak_list[peak_idx]
        assignment = peak.get('assignment', '')

        # Reset view orientation when selecting a new peak
        self.voigt_3d_plotter.reset_view_orientation()

        # Get list of spectra that have this peak
        spectra_with_peak = self.peak_spectrum_map.get(assignment, [])

        # Update intensity decay widget with new peak data
        if hasattr(self, 'intensity_decay_widget'):
            self.intensity_decay_widget.set_data(
                self.spectra,
                assignment,
                current_spectrum_index=None  # Will be set when spectrum is selected
            )
            # Pass exponential fit data from DynamiXs if available
            # Look up fit for this assignment from the all_fit_data dict
            # Keys can be full assignments (A.142.ALA.H) or just residue numbers (142)
            if self._exponential_fit_data:
                # Try full assignment first, then extracted residue number
                fit_data = self._exponential_fit_data.get(assignment)
                if fit_data is None:
                    residue_id = self._extract_residue_from_assignment(assignment)
                    fit_data = self._exponential_fit_data.get(residue_id) if residue_id else None

                if fit_data:
                    self.intensity_decay_widget.set_exponential_fit(fit_data)
                else:
                    self.intensity_decay_widget.clear_exponential_fit()

        # Update spectrum navigator list with grayed-out items
        self._update_voigt_3d_peak_mode_spectrum_list(spectra_with_peak)

        logger.debug(f"3D Voigt Peak Mode: Selected peak {assignment}, found in {len(spectra_with_peak)} spectra")

    def _update_voigt_3d_peak_mode_spectrum_list(self, spectra_with_peak: List[int]):
        """Update the Peak Mode spectrum list for 3D Voigt tab."""
        self.voigt_3d_peak_mode_spectrum_list.clear()
        self._voigt_3d_peak_mode_available_spectra = []

        for i, spec in enumerate(self.spectra):
            # Create colored icon for spectrum (matches overlay view colors)
            spectrum_color = spec.get('color', self.default_colors[i % len(self.default_colors)])
            icon = self._create_color_icon(spectrum_color)

            item = QListWidgetItem(icon, spec['name'])
            item.setData(Qt.UserRole, i)

            if i in spectra_with_peak:
                item.setFlags(item.flags() | Qt.ItemIsSelectable | Qt.ItemIsEnabled)
                self._voigt_3d_peak_mode_available_spectra.append(i)
            else:
                item.setFlags(item.flags() & ~Qt.ItemIsSelectable & ~Qt.ItemIsEnabled)
                item.setForeground(QColor('#aaaaaa'))

            self.voigt_3d_peak_mode_spectrum_list.addItem(item)

        # Update navigation button states
        has_multiple = len(self._voigt_3d_peak_mode_available_spectra) > 1
        self.voigt_3d_peak_mode_prev_btn.setEnabled(has_multiple)
        self.voigt_3d_peak_mode_next_btn.setEnabled(has_multiple)

        # Select first available spectrum and plot
        if self._voigt_3d_peak_mode_available_spectra:
            first_idx = self._voigt_3d_peak_mode_available_spectra[0]
            # Capture reference z-scale from first spectrum
            self._capture_reference_z_scale(first_idx)
            self.voigt_3d_peak_mode_spectrum_list.setCurrentRow(first_idx)
            self._select_voigt_3d_spectrum_for_peak_mode(first_idx)

    def _extract_residue_from_assignment(self, assignment: str) -> str:
        """Extract residue number from a peak assignment string.

        Handles various assignment formats:
        - "A.142.ALA.H" -> "142"
        - "142" -> "142"
        - "B.55.GLY.H" -> "55"

        Args:
            assignment: Peak assignment string

        Returns:
            Residue number as string, or None if not extractable
        """
        import re

        if not assignment:
            return None

        # Try to extract number from dotted format (e.g., "A.142.ALA.H")
        parts = assignment.split('.')
        if len(parts) >= 2:
            # Second part should be the residue number
            if parts[1].isdigit():
                return parts[1]

        # Try to find any number in the assignment
        match = re.search(r'\d+', assignment)
        if match:
            return match.group()

        return None

    def _on_voigt_3d_peak_mode_spectrum_selected(self, item: QListWidgetItem):
        """Handle spectrum selection in Peak Mode for 3D Voigt tab (click)."""
        spec_idx = item.data(Qt.UserRole)

        if not hasattr(self, '_voigt_3d_peak_mode_available_spectra'):
            return
        if spec_idx not in self._voigt_3d_peak_mode_available_spectra:
            return

        self._select_voigt_3d_spectrum_for_peak_mode(spec_idx)

    def _on_voigt_3d_peak_mode_spectrum_row_changed(self, row: int):
        """Handle spectrum selection in Peak Mode for 3D Voigt tab (keyboard navigation)."""
        if row < 0:
            return

        item = self.voigt_3d_peak_mode_spectrum_list.item(row)
        if item is None:
            return

        spec_idx = item.data(Qt.UserRole)

        if not hasattr(self, '_voigt_3d_peak_mode_available_spectra'):
            return
        if spec_idx not in self._voigt_3d_peak_mode_available_spectra:
            return

        self._select_voigt_3d_spectrum_for_peak_mode(spec_idx)

    def _select_voigt_3d_spectrum_for_peak_mode(self, spec_idx: int):
        """Select a spectrum in Peak Mode and plot 3D Voigt analysis."""
        self.selected_spectrum_index = spec_idx

        if self.selected_master_peak_idx is None or self.selected_master_peak_idx >= len(self.master_peak_list):
            return

        peak = self.master_peak_list[self.selected_master_peak_idx]
        assignment = peak.get('assignment', '')

        # Find the peak index in this spectrum's peak list
        fitted_peaks = self._get_fitted_peaks(spec_idx)
        for i, p in enumerate(fitted_peaks):
            p_assignment = p.get('assignment', p.get('Assignment', ''))
            if p_assignment == assignment:
                self.selected_peak_index = i
                break

        # Update intensity decay widget highlight
        if hasattr(self, 'intensity_decay_widget'):
            self.intensity_decay_widget.set_highlight(spec_idx)

        # Apply reference z-scale if locked (for intensity comparison across spectra)
        if self.lock_z_scale_to_reference and self.reference_z_min is not None:
            self.voigt_3d_plotter.set_fixed_z_limits(
                self.reference_z_min,
                self.reference_z_max
            )
        else:
            self.voigt_3d_plotter.clear_fixed_z_limits()

        # Apply reference x/y if locked (for position comparison across spectra)
        if self.lock_xy_to_reference and self.reference_x_min is not None:
            self.voigt_3d_plotter.set_fixed_xy_limits(
                self.reference_x_min,
                self.reference_x_max,
                self.reference_y_min,
                self.reference_y_max
            )
        else:
            self.voigt_3d_plotter.clear_fixed_xy_limits()

        # Plot the 3D Voigt analysis
        self._plot_voigt_analysis_3d()
        logger.debug(f"3D Voigt Peak Mode: Selected spectrum {spec_idx} for peak {assignment}")

    def _on_voigt_3d_peak_mode_prev_spectrum(self):
        """Navigate to previous available spectrum in 3D Voigt Peak Mode."""
        if not hasattr(self, '_voigt_3d_peak_mode_available_spectra') or not self._voigt_3d_peak_mode_available_spectra:
            return

        current_row = self.voigt_3d_peak_mode_spectrum_list.currentRow()
        try:
            current_pos = self._voigt_3d_peak_mode_available_spectra.index(current_row)
            new_pos = (current_pos - 1) % len(self._voigt_3d_peak_mode_available_spectra)
            new_idx = self._voigt_3d_peak_mode_available_spectra[new_pos]
            self.voigt_3d_peak_mode_spectrum_list.setCurrentRow(new_idx)
            self._select_voigt_3d_spectrum_for_peak_mode(new_idx)
        except ValueError:
            if self._voigt_3d_peak_mode_available_spectra:
                new_idx = self._voigt_3d_peak_mode_available_spectra[0]
                self.voigt_3d_peak_mode_spectrum_list.setCurrentRow(new_idx)
                self._select_voigt_3d_spectrum_for_peak_mode(new_idx)

    def _on_voigt_3d_peak_mode_next_spectrum(self):
        """Navigate to next available spectrum in 3D Voigt Peak Mode."""
        if not hasattr(self, '_voigt_3d_peak_mode_available_spectra') or not self._voigt_3d_peak_mode_available_spectra:
            return

        current_row = self.voigt_3d_peak_mode_spectrum_list.currentRow()
        try:
            current_pos = self._voigt_3d_peak_mode_available_spectra.index(current_row)
            new_pos = (current_pos + 1) % len(self._voigt_3d_peak_mode_available_spectra)
            new_idx = self._voigt_3d_peak_mode_available_spectra[new_pos]
            self.voigt_3d_peak_mode_spectrum_list.setCurrentRow(new_idx)
            self._select_voigt_3d_spectrum_for_peak_mode(new_idx)
        except ValueError:
            if self._voigt_3d_peak_mode_available_spectra:
                new_idx = self._voigt_3d_peak_mode_available_spectra[0]
                self.voigt_3d_peak_mode_spectrum_list.setCurrentRow(new_idx)
                self._select_voigt_3d_spectrum_for_peak_mode(new_idx)

    def _on_intensity_decay_point_clicked(self, spec_idx: int):
        """Handle click on intensity decay plot point to navigate to that spectrum.

        Args:
            spec_idx: Index of the spectrum corresponding to the clicked point
        """
        if not hasattr(self, '_voigt_3d_peak_mode_available_spectra'):
            return
        if spec_idx not in self._voigt_3d_peak_mode_available_spectra:
            return

        # Update spectrum list selection
        self.voigt_3d_peak_mode_spectrum_list.setCurrentRow(spec_idx)
        # Navigate to spectrum
        self._select_voigt_3d_spectrum_for_peak_mode(spec_idx)
        logger.debug(f"Intensity decay: Clicked point for spectrum {spec_idx}")

    # ===== Z-Scale Reference Methods for Peak Mode =====

    def _capture_reference_z_scale(self, spec_idx: int):
        """Capture z-axis limits from specified spectrum for reference scaling.

        Called when a peak is selected in Peak Mode to set the reference
        z-scale from the first spectrum. This allows visual comparison of
        intensity changes across spectra.

        Args:
            spec_idx: Index of the spectrum to capture z-limits from
        """
        import numpy as np

        if self.selected_master_peak_idx is None:
            return

        peak = self.master_peak_list[self.selected_master_peak_idx]
        assignment = peak.get('assignment', '')

        # Find the peak in this spectrum using same method as _plot_voigt_analysis_3d
        spec = self.spectra[spec_idx]
        fitted_peaks = self._get_fitted_peaks(spec_idx)

        for p in fitted_peaks:
            p_assignment = p.get('assignment', p.get('Assignment', ''))
            if p_assignment == assignment:
                # Use same lookup as _plot_voigt_analysis_3d
                stored_result = self._find_stored_voigt_result(spec, assignment, p)

                if stored_result is not None:
                    region_2d = stored_result.get('region_2d')
                    if region_2d is not None:
                        # Capture z (intensity) limits
                        experimental = region_2d['intensity']
                        self.reference_z_min = float(np.min(experimental))
                        self.reference_z_max = float(np.max(experimental))

                        # Capture x/y (chemical shift) limits
                        f1_ppm = region_2d['f1_ppm']
                        f2_ppm = region_2d['f2_ppm']
                        self.reference_x_min = float(np.min(f2_ppm))
                        self.reference_x_max = float(np.max(f2_ppm))
                        self.reference_y_min = float(np.min(f1_ppm))
                        self.reference_y_max = float(np.max(f1_ppm))

                        logger.debug(
                            f"Captured reference scales from spectrum {spec_idx}: "
                            f"z=[{self.reference_z_min:.2e}, {self.reference_z_max:.2e}], "
                            f"x=[{self.reference_x_min:.3f}, {self.reference_x_max:.3f}], "
                            f"y=[{self.reference_y_min:.2f}, {self.reference_y_max:.2f}]"
                        )
                break

    def _on_lock_z_scale_toggled(self, checked: bool):
        """Handle z-scale lock toggle checkbox.

        Args:
            checked: True if lock is enabled, False otherwise
        """
        self.lock_z_scale_to_reference = checked

        if not checked:
            self.voigt_3d_plotter.clear_fixed_z_limits()

        # Refresh current view if we have a selection
        if self.selected_master_peak_idx is not None:
            self._plot_voigt_analysis_3d()

    def _on_lock_xy_toggled(self, checked: bool):
        """Handle x/y axis lock toggle checkbox.

        Args:
            checked: True if lock is enabled, False otherwise
        """
        self.lock_xy_to_reference = checked

        if not checked:
            self.voigt_3d_plotter.clear_fixed_xy_limits()

        # Refresh current view if we have a selection
        if self.selected_master_peak_idx is not None:
            self._plot_voigt_analysis_3d()
