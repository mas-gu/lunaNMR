#!/usr/bin/env python3
# ABOUTME: Main window for lunaNMR v1.0 Qt interface - coordinates all GUI components
# ABOUTME: Integrates spectrum display, peak navigation, analysis panels, and file management

"""
Main Window for lunaNMR v1.0 (Qt/PySide6)

This is the central coordination point for the entire application, managing:
- Left panel: File selection and processing controls
- Center panel: Spectrum display and visualization
- Right panel: Peak Navigator and analysis results
- Bottom panel: Voigt Analysis tab (when active)
- Top menu bar: File operations, settings, and help

Architecture:
    - Inherits from BaseWindow (Phase 0) for common window functionality
    - Uses design system constants from styles module (Phase 0)
    - Integrates components from Phase 1-2 (dialogs, navigator, plotters)
    - Coordinates core processing engines (integrator, processors)

State Management:
    - Core components: integrator, file_manager, config_manager
    - Processing state: processing_active, integration_active
    - File paths: current_nmr_file, current_peak_file
    - Selection state: selected_peak_info
    - Workflow mode: peak_list vs sn_threshold

Author: Guillaume Mas
Date: 2025
"""

import logging
from pathlib import Path
from typing import Optional

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QScrollArea,
    QTabWidget, QLabel, QPushButton, QMessageBox, QGroupBox,
    QDoubleSpinBox, QSpinBox, QFormLayout, QCheckBox, QRadioButton, QButtonGroup,
    QTextEdit, QDialog, QDialogButtonBox, QComboBox, QSlider, QFrame, QGridLayout,
    QSizePolicy
)
from PySide6.QtCore import Qt, Signal, QThread, Slot
from PySide6.QtGui import QFont

from lunaNMR.gui.base.base_window import BaseWindow
from lunaNMR.gui.components.peak_navigator import PeakNavigator
from lunaNMR.gui.styles.design_system import (
    PANEL_BG_COLOR, FRAME_BG_COLOR, PRIMARY_TEXT, SECONDARY_TEXT,
    SPACING_SM, SPACING_MD, SPACING_LG, SPACING_XS,
    FONT_SIZE_BODY, FONT_SIZE_SECTION_LABEL, FONT_SIZE_SMALL,
    PRIMARY_BUTTON_BG, PRIMARY_BUTTON_HOVER, PRIMARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BG, SECONDARY_BUTTON_HOVER, SECONDARY_BUTTON_TEXT, SECONDARY_BUTTON_BORDER,
    SUCCESS_GREEN, DESTRUCTIVE_BUTTON_BG, DESTRUCTIVE_BUTTON_HOVER, DESTRUCTIVE_BUTTON_TEXT,
    BUTTON_CORNER_RADIUS, BUTTON_HEIGHT_DIALOG
)

# Import core processing components
from lunaNMR.core.core_integrator import EnhancedVoigtIntegrator
from lunaNMR.core.ps2d_config import set_ps2d_config, get_ps2d_config
from lunaNMR.utils.file_manager import NMRFileManager
from lunaNMR.utils.config_manager import ConfigurationManager, UserPreferences, ProcessingParameters
from lunaNMR.utils.parameter_manager import NMRParameterManager
from lunaNMR.utils.project_manager import ProjectManager

# Import processors (will be instantiated when needed)
from lunaNMR.processors.single_spectrum_processor import SingleSpectrumProcessor
from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor

logger = logging.getLogger(__name__)


class PeakDetectionWorker(QThread):
    """Worker thread for peak detection to prevent GUI freezing.

    Signals:
        progress_updated: (int, str) - Progress percentage and status message
        detection_complete: (list) - Detected peaks list
        detection_failed: (str) - Error message if detection fails
    """

    progress_updated = Signal(int, str)  # progress_percent, status_message
    detection_complete = Signal(list)     # detected_peaks
    detection_failed = Signal(str)        # error_message

    def __init__(self, integrator, workflow_mode: str):
        """Initialize the peak detection worker.

        Args:
            integrator: Core integrator instance with nmr_data loaded
            workflow_mode: 'peak_list' or 'sn_threshold'
        """
        super().__init__()
        self.integrator = integrator
        self.workflow_mode = workflow_mode

    def run(self):
        """Run peak detection in background thread."""
        try:
            logger.info(f"Starting peak detection in background (mode={self.workflow_mode})")
            self.progress_updated.emit(10, "Preparing detection...")

            # Set OutputManager callback for this thread so log_* functions route to GUI
            from lunaNMR.utils.output_manager import OutputManager
            def progress_callback(progress, task, log_msg=None, failed=False):
                self.progress_updated.emit(int(progress) if progress >= 0 else 50, log_msg or task)
            OutputManager.set_callback(progress_callback)

            # Check if NMR data is loaded
            if not hasattr(self.integrator, 'nmr_data') or self.integrator.nmr_data is None:
                self.detection_failed.emit("No NMR data loaded. Please load a spectrum first.")
                return

            self.progress_updated.emit(20, "Running peak detection...")

            # Call appropriate detection method based on workflow mode
            if self.workflow_mode == 'sn_threshold':
                # S/N native detection
                detected_peaks = self.integrator.detect_peaks_sn_native()
            else:
                # Peak list mode - use process_peaks
                detected_peaks = self.integrator.process_peaks()

            self.progress_updated.emit(80, "Processing results...")

            # Validate results
            if detected_peaks is None or len(detected_peaks) == 0:
                self.detection_failed.emit("No peaks detected. Check parameters or spectrum quality.")
                return

            self.progress_updated.emit(100, "Detection complete")
            logger.info(f"Peak detection complete: {len(detected_peaks)} peaks detected")

            # Emit results
            self.detection_complete.emit(detected_peaks)

        except Exception as e:
            logger.error(f"Peak detection failed: {e}", exc_info=True)
            self.detection_failed.emit(f"Detection error: {str(e)}")
        finally:
            # Clear OutputManager callback for this thread
            from lunaNMR.utils.output_manager import OutputManager
            OutputManager.set_callback(None)


class VoigtFittingWorker(QThread):
    """Worker thread for Voigt fitting to prevent GUI freezing.

    Uses SingleSpectrumProcessor which handles PS2D multi-peak fitting automatically.
    This is the exact logic from lunaNMR_v0o9.

    Signals:
        progress_updated: (int, str, str) - Progress percentage, task, log message
        fitting_complete: (dict) - Fitting summary with results
        fitting_failed: (str) - Error message if fitting fails
    """

    progress_updated = Signal(int, str, str)  # progress, task, log_msg
    fitting_complete = Signal(dict)            # summary dict
    fitting_failed = Signal(str)               # error_message

    def __init__(self, integrator, param_manager, processing_options):
        """Initialize the Voigt fitting worker.

        Args:
            integrator: Core integrator instance with peak_list and nmr_data
            param_manager: Parameter manager for configuration
            processing_options: Dict with keys:
                - use_parallel: bool
                - use_global_optimization: bool (optional, default False)
                - use_voigt_fitting: bool
        """
        super().__init__()
        self.integrator = integrator
        self.param_manager = param_manager
        self.processing_options = processing_options

    def run(self):
        """Run Voigt fitting in background thread using SingleSpectrumProcessor."""
        try:
            logger.info("Starting Voigt fitting in background")
            self.progress_updated.emit(5, "Initializing processor...", "")

            # Import SingleSpectrumProcessor (v0o9 pattern)
            from lunaNMR.processors.single_spectrum_processor import SingleSpectrumProcessor
            from lunaNMR.utils.output_manager import OutputManager

            # Create processor with integrator and parameter manager
            processor = SingleSpectrumProcessor(self.integrator, self.param_manager)

            # Define progress callback (Qt signals are thread-safe)
            def progress_callback(progress, task, log_msg=None, failed=False):
                log_text = log_msg if log_msg else ""
                self.progress_updated.emit(int(progress), task, log_text)

            # Set OutputManager callback for this thread so log_* functions route to GUI
            OutputManager.set_callback(progress_callback)

            self.progress_updated.emit(10, "Starting peak fitting...", "")

            # Process peak list using SingleSpectrumProcessor (exact v0o9 logic)
            # This internally handles PS2D multi-peak fitting
            # Returns tuple: (fitted_results, learned_statistics)
            fitted_results, learned_statistics = processor.process_peak_list(
                self.integrator.peak_list,
                self.processing_options,
                progress_callback
            )

            self.progress_updated.emit(90, "Getting summary...", "")

            # Get comprehensive summary (v0o9 pattern)
            summary = processor.get_processing_summary(
                fitted_results,
                len(self.integrator.peak_list)
            )

            # Add fitted_results to summary for GUI update
            summary['fitted_results'] = fitted_results
            summary['learned_statistics'] = learned_statistics

            self.progress_updated.emit(100, "Fitting complete", "")
            logger.info(f"Voigt fitting complete: {summary['successful_peaks']}/{summary['total_peaks']} successful")

            # Emit results
            self.fitting_complete.emit(summary)

        except Exception as e:
            logger.error(f"Voigt fitting failed: {e}", exc_info=True)
            self.fitting_failed.emit(f"Fitting error: {str(e)}")
        finally:
            # Clear OutputManager callback for this thread
            from lunaNMR.utils.output_manager import OutputManager
            OutputManager.set_callback(None)


class LunaNMRMainWindow(BaseWindow):
    """Main application window for lunaNMR v1.0.

    This is the central coordination point for the entire application, managing:
    - Left panel: File selection and processing controls
    - Center panel: Spectrum display and visualization
    - Right panel: Peak Navigator and analysis results
    - Bottom panel: Voigt Analysis tab (when active)
    - Top menu bar: File operations, settings, and help

    Architecture:
        - Inherits from BaseWindow (Phase 0) for common window functionality
        - Uses design system constants from styles module (Phase 0)
        - Integrates components from Phase 1-2 (dialogs, navigator, plotters)
        - Coordinates core processing engines (integrator, processors)

    State Management:
        - Core components: integrator, file_manager, config_manager
        - Processing state: processing_active, integration_active
        - File paths: current_nmr_file, current_peak_file
        - Selection state: selected_peak_info
        - Workflow mode: peak_list vs sn_threshold
    """

    # Signals for inter-component communication
    spectrum_loaded = Signal(str)  # Emitted when spectrum is loaded
    peaks_updated = Signal()  # Emitted when peak list changes
    processing_started = Signal()  # Emitted when processing begins
    processing_finished = Signal()  # Emitted when processing completes

    def __init__(self):
        """Initialize the main window and all subsystems."""
        super().__init__(
            title="lunaNMR v1.0 - NMR Peak Analysis",
            default_size=(1400, 900),
            min_size=(1000, 700),
            enable_status_bar=True
        )

        # Initialize state variables first (before GUI setup)
        self.init_variables()

        # Setup the user interface
        self.setup_ui()

        # Apply user preferences (window geometry, etc.)
        self.apply_user_preferences()

        # Sync all parameters to integrator
        self.sync_parameters_to_integrator()

        logger.info("lunaNMR v1.0 main window initialized successfully")
        self.update_status("Application ready")

    def init_variables(self):
        """Initialize all state variables and core components.

        This method sets up:
        1. Core processing components (integrator, file_manager, etc.)
        2. Configuration management (config_manager, user_prefs, proc_params)
        3. Processing state flags (processing_active, integration_active, etc.)
        4. File path tracking (current_nmr_file, current_peak_file, etc.)
        5. Selection state (selected_peak_info, selected_peak_number)
        6. Workflow configuration (workflow_mode, sn_threshold, etc.)

        Note: This is called BEFORE setup_ui() to ensure all state is ready
        when GUI components are created.
        """
        logger.debug("Initializing state variables and core components")

        # ===== Core Processing Components =====
        self.integrator = EnhancedVoigtIntegrator()
        self.file_manager = NMRFileManager()
        self.param_manager = NMRParameterManager()

        # Configuration management
        self.config_manager = ConfigurationManager("lunaNMR_v1")
        self.user_prefs = UserPreferences(self.config_manager)
        self.proc_params = ProcessingParameters(self.config_manager)

        # ===== Processing State Flags =====
        self.processing_active = False
        self.integration_active = False
        self.integration_paused = False
        self.integration_start_time = None

        # Processors (created on-demand)
        self.single_spectrum_processor = None
        self.multi_spectrum_processor = None

        # ===== File Path Tracking =====
        self.current_nmr_file = None
        self.current_peak_file = None
        self.current_nmr_folder = None
        self.current_peak_folder = None

        # Results tracking
        self.current_voigt_result = None
        self.last_fitting_results = None  # Store last fitting results for series integration
        self.batch_results = None
        self.series_output_folder = None

        # Project management
        self.current_project_path = None  # Path to currently loaded .lunaNMR project
        self.series_all_results = None  # For reopening multi-spectrum viewer
        self.project_manager = None  # Initialized after self is fully set up

        # DynamiXs state (for project save/load)
        self.dynamixs_state = None  # DynamiXs dialog parameters
        self.dynamixs_file_refs = None  # DynamiXs input file paths
        self.dynamixs_dialog = None  # Reference to open DynamiXs dialog

        # ===== Selection State =====
        self.selected_peak_info = None  # {type, index, data}
        self.selected_peak_number = 1  # For peak navigation

        # ===== Workflow Configuration =====
        # Workflow mode: "peak_list" (default) or "sn_threshold"
        self.workflow_mode = "peak_list"
        self.sn_threshold = 3.0
        self.expected_peak_count = 50

        # ===== Display Options =====
        self.show_detected = True
        self.show_assigned = True
        self.show_fitted_curves = True
        self.show_ellipses = False  # Debug tool for PS2D

        # ===== Contour Display Parameters =====
        self.contour_levels = 24  # Number of contour levels (5-50)
        self.contour_min = 0.05  # Minimum contour level (0.01-1.0)
        self.contour_increment = 1.1  # Contour increment factor (0.01-10.0)

        # ===== Zoom Parameters =====
        self.zoom_x_range = 4.0  # X axis zoom range in ppm (default: 4.0)
        self.zoom_y_range = 30.0  # Y axis zoom range in ppm (default: 30.0)
        self.saved_xlim = None  # Saved X axis limits for zoom preservation
        self.saved_ylim = None  # Saved Y axis limits for zoom preservation

        # ===== Processing Parameters =====
        # These will be populated from config_manager by load_parameters_from_config()
        # For now, just initialize placeholders
        self.use_voigt_fitting = True
        self.use_parallel_processing = True
        self.use_adaptive_optimization = True  # Grid search for optimal radF1/radF2
        self.use_ps2d_multi_peak = True
        self.fix_linewidths = False
        self.fix_positions = False

        # ===== Peak Detection Parameters (Expert Mode) =====
        # Detection parameters from parameter_manager defaults
        self.noise_threshold = 3.0
        self.search_window_x = 0.08  # 1H dimension (ppm)
        self.search_window_y = 0.8   # 15N dimension (ppm)
        self.detection_square_size = 3  # X-dimension/1H (pixels)
        self.detection_rectangle_y = 1  # Y-dimension/15N (pixels)
        self.detection_square_ppm_x = "0.000"  # Auto-calculated ppm conversion
        self.detection_square_ppm_y = "0.000"  # Auto-calculated ppm conversion

        # Peak Centroid Detection parameters
        self.centroid_window_x_ppm = 0.01  # X-window size in ppm (1H dimension)
        self.centroid_window_y_ppm = 0.1   # Y-window size in ppm (15N dimension)
        self.use_centroid_refinement = True  # Enable centroid refinement
        self.centroid_noise_multiplier = 2.0  # Noise threshold for centroid

        # Auto-add dummy peaks (for peak_list mode only)
        self.auto_add_dummy_peaks = True  # Auto-add nearby unmatched peaks

        # ===== Peak Detection Parameters (for series integration) =====
        self.height_threshold = 0.1  # Peak height threshold
        self.distance_factor = 2.0  # Peak distance factor
        self.prominence_threshold = 0.05  # Peak prominence threshold
        self.smoothing_sigma = 1.0  # Smoothing sigma
        self.max_peaks_fit = 50  # Maximum peaks to fit
        self.max_optimization_iterations = 50  # Optimization iterations

        # ===== Fitting Quality Parameters =====
        self.min_r_squared = 0.5  # Minimum R² quality threshold
        self.max_iterations = 1000  # Levenberg-Marquardt max iterations
        self.fitting_window_x = 0.2  # 1H fitting window (ppm)
        self.fitting_window_y = 2.0  # 15N fitting window (ppm)
        self.use_global_optimization = False  # Global optimization flag

        # ===== PS2D Algorithm Configuration (Expert Mode) =====
        self.nucleus_type = '15N'  # Default to 15N-HSQC
        # PS2D parameters initialized from ps2d_config
        ps2d_config = get_ps2d_config()
        self.ps2d_radF1 = ps2d_config.radF1
        self.ps2d_radF2 = ps2d_config.radF2
        self.ps2d_max_iterations = ps2d_config.max_iterations
        self.ps2d_max_cluster_size = ps2d_config.max_cluster_size
        self.ps2d_overlap_x = ps2d_config.overlap_threshold_x
        self.ps2d_overlap_y = ps2d_config.overlap_threshold_y

        # ===== Custom Linewidth Parameters (Expert Mode) =====
        self.use_custom_linewidths = False
        self.lw_lorentz_1h = 0.001
        self.lw_gauss_1h = 0.03    # 1H Gaussian (ppm)
        self.lw_lorentz_15n = 0.0001
        self.lw_gauss_15n = 0.3    # 15N Gaussian (ppm)

        # ===== Series Integration Parameters (Expert Mode) =====
        self.use_ps2d_linewidth_reuse = False  # PS2D Linewidth Reuse for ~40% speedup
        self.series_peak_source = "cascade"    # "detected" or "cascade"
        self.enable_cascade_drift_limit = True  # Enforce absolute drift limits in cascade mode (ON by default)

        # ===== Expert Mode Dialog State =====
        self.params_visible = False           # Peak Detection Parameters visibility
        self.ps2d_config_params_visible = False  # PS2D Config Parameters visibility
        self.ps2d_fit_params_visible = False  # PS2D Fit Parameters visibility
        self.custom_lw_visible = False        # Custom linewidth inputs visibility

        # ===== Shift Peak List Parameters =====
        self.shift_1h_value = 0.0  # 1H offset in ppm
        self.shift_15n_value = 0.0  # 15N/13C offset in ppm

        # ===== Peak Edition Parameters =====
        self.peak_edition_visible = False  # Collapsible section state
        self.peak_edit_mode = False  # Whether edit mode is active
        self.edit_reference_peaks = False
        self.edit_detected_peaks = False
        self.peak_deletion_mode = False
        self.peak_addition_mode = False
        self.click_tolerance = 0.15  # Click distance tolerance in ppm (v0o9 default)
        self.edit_connection_id = None  # Matplotlib event connection ID

        # ===== 3D Voigt Analysis Control Variables (v0o9 line 1840-1848) =====
        self.show_exp_3d = True  # Show experimental data
        self.show_fit_3d = True  # Show fitted data
        self.show_individual_3d = False  # Show individual peaks (hidden by default)
        self.show_peak_labels_3d = False  # Show peak labels (hidden by default)
        self.show_resid_3d = False  # Show residuals (hidden by default)
        self.limit_peak_display_3d = True  # Limit peak extent (ON by default)
        self.residual_mode_3d = 'overlay'  # Residual mode: 'overlay' or 'separate'
        self.color_scheme_3d = 'Clean'  # Color scheme: 'Classic', 'Clean', 'Dark', 'Warm'
        self.intensity_scale_3d = 100.0  # Intensity scale: 50-200%

        # ===== GUI Components (to be created in setup_ui) =====
        # Left panel components
        self.file_panel = None
        self.controls_panel = None
        self.contour_params_widget = None  # Collapsible contour settings widget

        # Workflow Selection section GUI components
        self.workflow_peak_list_radio = None
        self.workflow_sn_threshold_radio = None
        self.sn_params_widget = None  # Collapsible S/N parameters widget
        self.sn_threshold_spin = None
        self.expected_peaks_spin = None

        # Data Loading section GUI components
        self.nmr_file_label = None
        self.peak_file_label = None

        # Center panel components
        self.spectrum_plotter = None

        # Right panel components
        self.peak_navigator = None

        # Bottom panel components
        self.voigt_analysis_plotter = None
        self.voigt_2d_plotter = None  # 2D Voigt analysis plotter (2x2 grid)
        self.peak_params_text = None  # Text widget for peak parameters
        self.series_plotter = None

        # Splitters for layout
        self.main_splitter = None

        # Shift Peak List GUI components
        self.shift_1h_spin = None
        self.shift_15n_spin = None

        logger.debug("State variables initialized")

    def setup_ui(self):
        """Setup the complete user interface.

        This creates the main window layout with:
        - Menu bar at top
        - Three-panel layout (left controls, center spectrum, right navigator)
        - Status bar at bottom
        - Optional bottom panel for Voigt Analysis

        Layout hierarchy:
            QMainWindow
            ├── Menu Bar (File, View, Process, Help)
            ├── Central Widget
            │   └── Main Splitter (horizontal)
            │       ├── Left Panel (file browser + controls)
            │       ├── Center Panel (spectrum display)
            │       └── Right Panel (peak navigator)
            └── Status Bar

        Note: Detailed component creation will be done by specialized agents
        in later phases.
        """
        logger.debug("Setting up user interface")

        # Create menu bar
        self.setup_menu_bar()

        # Create central widget and main layout
        self.setup_central_widget()

        logger.debug("User interface setup complete")

    def setup_central_widget(self):
        """Create and configure the 3-panel layout.

        Layout structure:
            - Main horizontal splitter containing 3 panels
            - Left: Scrollable controls panel (min 350px)
            - Center: Tabbed plot display (stretches 3x)
            - Right: Peak Navigator (min 200px)

        Initial sizes follow 1:3:1 ratio for 1400px window:
            - Left: 280px (1 unit)
            - Center: 840px (3 units)
            - Right: 280px (1 unit)
        """
        # Create main horizontal splitter
        main_splitter = QSplitter(Qt.Horizontal)
        main_splitter.setHandleWidth(1)  # Thin handle for cleaner look
        main_splitter.setChildrenCollapsible(False)  # Prevent panels from collapsing

        # Create the three panels
        left_panel = self.create_left_panel()
        center_panel = self.create_center_panel()
        right_panel = self.create_right_panel()

        # Add panels to splitter
        main_splitter.addWidget(left_panel)
        main_splitter.addWidget(center_panel)
        main_splitter.addWidget(right_panel)

        # Calculate proportional sizes based on available screen width
        from PySide6.QtWidgets import QApplication
        screen = QApplication.primaryScreen()
        if screen:
            available_width = screen.availableGeometry().width()
        else:
            available_width = 1400  # Fallback

        # Use 1:3:1 ratio, scaled to screen width (with some margin for window chrome)
        usable_width = int(available_width * 0.85)  # 85% of screen width
        left_width = int(usable_width * 0.20) + 50   # 20% + 50px for left panel
        center_width = int(usable_width * 0.55) - 25 # 55% - 25px for center
        right_width = int(usable_width * 0.25) - 25  # 25% - 25px for right panel

        main_splitter.setSizes([left_width, center_width, right_width])

        # Set stretch factors (how much each panel grows/shrinks)
        main_splitter.setStretchFactor(0, 1)  # Left panel - standard growth
        main_splitter.setStretchFactor(1, 3)  # Center panel - 3x growth
        main_splitter.setStretchFactor(2, 1)  # Right panel - standard growth

        # Set minimum widths
        left_panel.setMinimumWidth(250)  # Increased by 50px
        right_panel.setMinimumWidth(180)
        center_panel.setMinimumWidth(250)

        # Set central widget
        self.setCentralWidget(main_splitter)

        # Store references for later use
        self.main_splitter = main_splitter
        self.left_panel = left_panel
        self.center_panel = center_panel
        self.right_panel = right_panel

    def create_left_panel(self) -> QWidget:
        """Create the left panel with scrollable controls area.

        Returns:
            QScrollArea containing control sections including:
            - Data Loading section
            - Spectrum Display Controls section (contour settings, zoom reset)
            - Control Center with action buttons
            - Peak Operations section (placeholder)
            - Parameters section (placeholder)
        """
        # Create scroll area for controls
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll_area.setFrameShape(QScrollArea.NoFrame)

        # Set background color
        scroll_area.setStyleSheet(f"background-color: {PANEL_BG_COLOR};")

        # Create content widget for scrolling
        content_widget = QWidget()
        content_layout = QVBoxLayout(content_widget)
        content_layout.setContentsMargins(SPACING_SM, SPACING_MD, SPACING_SM, SPACING_MD)  # Reduced horizontal margins
        content_layout.setSpacing(SPACING_MD)

        # Add Workflow Selection section (FIRST - at top)
        workflow_section = self.create_workflow_selection_section()
        content_layout.addWidget(workflow_section)

        # Add Data Loading section
        data_loading_section = self.create_data_loading_section()
        content_layout.addWidget(data_loading_section)

        # Add Spectrum Display Controls section
        spectrum_display_section = self.create_spectrum_display_section()
        content_layout.addWidget(spectrum_display_section)

        # Add Shift Peak List section
        shift_peak_list_section = self.create_shift_peak_list_section()
        content_layout.addWidget(shift_peak_list_section)

        # Add Control Center section
        control_center_section = self.create_control_center_section()
        content_layout.addWidget(control_center_section)

        # Add Processing Statistics section
        statistics_section = self.create_statistics_panel_section()
        content_layout.addWidget(statistics_section)

        # Add Advanced Options section
        advanced_options_section = self.create_advanced_options_section()
        content_layout.addWidget(advanced_options_section)

        # Add Voigt Fitting section
        voigt_fitting_section = self.create_voigt_fitting_section()
        content_layout.addWidget(voigt_fitting_section)

        # Add stretch to push content to top
        content_layout.addStretch()

        # Set content widget to scroll area
        scroll_area.setWidget(content_widget)

        return scroll_area

    def create_data_loading_section(self) -> QGroupBox:
        """Create the Data Loading section with Load Data button and status labels.

        Returns:
            QGroupBox containing:
            - Large "Load Data" button (primary action)
            - Status labels showing loaded NMR and peak files
        """
        group = QGroupBox("Data Loading")
        layout = QVBoxLayout(group)
        layout.setSpacing(SPACING_SM)
        layout.setContentsMargins(SPACING_SM, SPACING_MD, SPACING_SM, SPACING_SM)  # Reduced horizontal margins

        # Large Load Data button (PRIMARY ACTION)
        load_btn = QPushButton("Load Data")
        load_btn.setToolTip("Load NMR Spectrum and Peak List files")
        load_btn.setProperty("class", "primary")
        load_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)  # 40px - prominent
        load_btn.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        load_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {PRIMARY_BUTTON_BG};
                color: {PRIMARY_BUTTON_TEXT};
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                font-size: {FONT_SIZE_BODY}px;
                font-weight: bold;
                padding: 6px 8px;
            }}
            QPushButton:hover {{
                background-color: {PRIMARY_BUTTON_HOVER};
            }}
            QPushButton:pressed {{
                background-color: {PRIMARY_BUTTON_HOVER};
            }}
        """)
        load_btn.clicked.connect(self.on_load_data)
        layout.addWidget(load_btn)

        # Current Data label (section header)
        current_label = QLabel("Current Data:")
        current_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_SMALL}px;
                color: {PRIMARY_TEXT};
                font-weight: bold;
                padding-top: {SPACING_SM}px;
            }}
        """)
        layout.addWidget(current_label)

        # NMR file status label
        self.nmr_file_label = QLabel("NMR Spectrum: Not loaded")
        self.nmr_file_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_SMALL}px;
                color: {SECONDARY_TEXT};
                padding: 2px 0px;
            }}
        """)
        self.nmr_file_label.setWordWrap(True)
        layout.addWidget(self.nmr_file_label)

        # Peak file status label
        self.peak_file_label = QLabel("Peak List: Not loaded")
        self.peak_file_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_SMALL}px;
                color: {SECONDARY_TEXT};
                padding: 2px 0px;
            }}
        """)
        self.peak_file_label.setWordWrap(True)
        layout.addWidget(self.peak_file_label)

        return group

    def create_workflow_selection_section(self) -> QGroupBox:
        """Create the Peak Detection Workflow section with radio buttons.

        This section allows users to choose between two peak detection workflows:
        1. With Peak List (Standard): Uses a reference peak list to guide detection
        2. S/N Threshold (Native Detection): Detects peaks using signal-to-noise threshold

        When "S/N Threshold" is selected, additional parameters are shown:
        - S/N Threshold: Minimum signal-to-noise ratio (1.0-100.0, default 3.0)
        - Expected Peak Count: Approximate number of peaks expected (10-1000, default 50)

        Returns:
            QGroupBox containing workflow selection controls
        """
        group = QGroupBox("Peak Detection Workflow")
        group.setStyleSheet(f"""
            QGroupBox {{
                font-size: {FONT_SIZE_SECTION_LABEL}px;
                font-weight: bold;
                color: {PRIMARY_TEXT};
                border: 1px solid #D1D1D6;
                border-radius: 8px;
                margin-top: 8px;
                padding-top: 12px;
                background-color: {FRAME_BG_COLOR};
            }}
            QGroupBox::title {{
                subcontrol-origin: margin;
                subcontrol-position: top left;
                left: 10px;
                padding: 0 5px;
            }}
        """)

        layout = QVBoxLayout(group)
        layout.setContentsMargins(SPACING_SM, SPACING_MD, SPACING_SM, SPACING_SM)  # Reduced horizontal margins
        layout.setSpacing(SPACING_SM)

        # Create radio buttons
        self.workflow_peak_list_radio = QRadioButton("With Peak List (Standard)")
        self.workflow_sn_threshold_radio = QRadioButton("S/N Threshold (Native Detection)")

        # Style radio buttons
        radio_style = f"""
            QRadioButton {{
                font-size: {FONT_SIZE_BODY}px;
                color: {PRIMARY_TEXT};
                padding: 2px;
            }}
        """
        self.workflow_peak_list_radio.setStyleSheet(radio_style)
        self.workflow_sn_threshold_radio.setStyleSheet(radio_style)

        # Create button group for mutual exclusivity
        workflow_button_group = QButtonGroup(self)
        workflow_button_group.addButton(self.workflow_peak_list_radio)
        workflow_button_group.addButton(self.workflow_sn_threshold_radio)

        # Set default: Peak List
        self.workflow_peak_list_radio.setChecked(True)

        # Add radio buttons to layout
        layout.addWidget(self.workflow_peak_list_radio)
        layout.addWidget(self.workflow_sn_threshold_radio)

        # Create collapsible S/N parameters widget
        self.sn_params_widget = QWidget()
        params_layout = QFormLayout(self.sn_params_widget)
        params_layout.setContentsMargins(SPACING_MD, SPACING_SM, SPACING_MD, SPACING_SM)
        params_layout.setSpacing(SPACING_XS)
        params_layout.setLabelAlignment(Qt.AlignRight)

        # S/N Threshold spinbox
        self.sn_threshold_spin = QDoubleSpinBox()
        self.sn_threshold_spin.setRange(1.0, 100.0)
        self.sn_threshold_spin.setValue(3.0)
        self.sn_threshold_spin.setSingleStep(0.5)
        self.sn_threshold_spin.setDecimals(1)
        self.sn_threshold_spin.setStyleSheet(f"""
            QDoubleSpinBox {{
                background-color: {FRAME_BG_COLOR};
                color: {PRIMARY_TEXT};
                border: 1px solid #D1D1D6;
                border-radius: 4px;
                padding: {SPACING_XS}px;
                font-size: {FONT_SIZE_BODY}px;
            }}
        """)
        params_layout.addRow("S/N Threshold:", self.sn_threshold_spin)

        # Expected peak count spinbox
        self.expected_peaks_spin = QSpinBox()
        self.expected_peaks_spin.setRange(10, 1000)
        self.expected_peaks_spin.setValue(50)
        self.expected_peaks_spin.setSingleStep(10)
        self.expected_peaks_spin.setStyleSheet(f"""
            QSpinBox {{
                background-color: {FRAME_BG_COLOR};
                color: {PRIMARY_TEXT};
                border: 1px solid #D1D1D6;
                border-radius: 4px;
                padding: {SPACING_XS}px;
                font-size: {FONT_SIZE_BODY}px;
            }}
        """)
        params_layout.addRow("Expected Peaks:", self.expected_peaks_spin)

        # Initially hide S/N parameters
        self.sn_params_widget.setVisible(False)
        layout.addWidget(self.sn_params_widget)

        # Connect radio button toggle signals
        # Note: Only need to connect one since they're mutually exclusive
        # But connecting both makes the behavior more explicit
        self.workflow_peak_list_radio.toggled.connect(self.on_workflow_changed)
        self.workflow_sn_threshold_radio.toggled.connect(self.on_workflow_sn_changed)

        return group

    def on_workflow_changed(self, checked):
        """Handle workflow mode radio button changes.

        Args:
            checked: True if "With Peak List" radio is checked, False otherwise

        When "With Peak List" is selected:
        - Hides S/N parameter controls
        - Sets workflow_mode to "peak_list"

        When "S/N Threshold" is selected:
        - Shows S/N parameter controls
        - Sets workflow_mode to "sn_threshold"
        """
        if checked:
            # "With Peak List" selected
            self.sn_params_widget.setVisible(False)
            self.workflow_mode = "peak_list"
            logger.debug(f"Workflow mode: peak_list (checked={checked}, S/N params hidden, visible={self.sn_params_widget.isVisible()})")
        else:
            # "S/N Threshold" selected
            self.sn_params_widget.setVisible(True)
            self.workflow_mode = "sn_threshold"
            logger.debug(f"Workflow mode: sn_threshold (checked={checked}, S/N params visible, visible={self.sn_params_widget.isVisible()})")

        # Enable/disable auto-add dummy peaks checkbox based on workflow mode
        self._update_auto_add_dummy_enabled()

        self.update_status(f"Workflow mode: {self.workflow_mode}")

    def on_workflow_sn_changed(self, checked):
        """Handle S/N Threshold radio button changes.

        Args:
            checked: True if "S/N Threshold" radio is checked, False otherwise

        This is the complementary handler to on_workflow_changed.
        When S/N Threshold is selected, show the S/N parameters.
        When it's deselected (Peak List is selected), hide them.
        """
        if checked:
            # "S/N Threshold" selected
            self.sn_params_widget.setVisible(True)
            self.workflow_mode = "sn_threshold"
            logger.debug(f"Workflow mode: sn_threshold (checked={checked}, S/N params visible)")
        else:
            # "Peak List" selected (S/N deselected)
            self.sn_params_widget.setVisible(False)
            self.workflow_mode = "peak_list"
            logger.debug(f"Workflow mode: peak_list (checked={checked}, S/N params hidden)")

        # Enable/disable auto-add dummy peaks checkbox based on workflow mode
        self._update_auto_add_dummy_enabled()

        self.update_status(f"Workflow mode: {self.workflow_mode}")

    def _update_auto_add_dummy_enabled(self):
        """Enable/disable auto-add dummy peaks checkbox based on workflow mode.

        Only enabled in peak_list mode (grayed out in S/N threshold mode).
        """
        if hasattr(self, 'expert_auto_add_dummy_cb'):
            enabled = (self.workflow_mode == "peak_list")
            self.expert_auto_add_dummy_cb.setEnabled(enabled)

    # ===== Button Style Helpers =====

    def _get_green_button_style(self) -> str:
        """Get green (success) button style."""
        return f"""
            QPushButton {{
                background-color: {SUCCESS_GREEN};
                color: white;
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                font-size: {FONT_SIZE_BODY}px;
                font-weight: bold;
                padding: 6px 8px;
            }}
            QPushButton:hover {{
                background-color: #2AA64A;
            }}
            QPushButton:pressed {{
                background-color: #228B3C;
            }}
            QPushButton:disabled {{
                background-color: #C7C7CC;
                color: #8E8E93;
            }}
        """

    def _get_blue_button_style(self) -> str:
        """Get blue (primary) button style."""
        return f"""
            QPushButton {{
                background-color: {PRIMARY_BUTTON_BG};
                color: {PRIMARY_BUTTON_TEXT};
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                font-size: {FONT_SIZE_BODY}px;
                font-weight: bold;
                padding: 6px 8px;
            }}
            QPushButton:hover {{
                background-color: {PRIMARY_BUTTON_HOVER};
            }}
            QPushButton:pressed {{
                background-color: #3A7CC3;
            }}
            QPushButton:disabled {{
                background-color: #C7C7CC;
                color: #8E8E93;
            }}
        """

    def _update_button_colors_for_detection_complete(self):
        """Update button colors after detection completes.

        - Detect button → blue (action completed)
        - Fit Spectrum & Fit Series → green (ready for next action)
        """
        self.detect_peaks_button.setStyleSheet(self._get_blue_button_style())
        self.fit_all_peaks_button.setStyleSheet(self._get_green_button_style())
        self.series_button.setStyleSheet(self._get_green_button_style())

    def _reset_button_colors(self):
        """Reset button colors to initial state.

        - Detect button → green (primary action)
        - Fit Spectrum & Fit Series → blue (secondary action)
        """
        self.detect_peaks_button.setStyleSheet(self._get_green_button_style())
        self.fit_all_peaks_button.setStyleSheet(self._get_blue_button_style())
        self.series_button.setStyleSheet(self._get_blue_button_style())

    def create_control_center_section(self) -> QGroupBox:
        """Create the Control Center section with main action buttons.

        This section contains:
        - Three main action buttons (Detect, Fit Spectrum, Fit Series)
        - Collapsible Peak Edition subsection with:
          * Mode status label
          * Edit mode checkboxes (Reference/Detected peaks)
          * Delete/Add mode checkboxes
          * Selected peak info label

        Returns:
            QGroupBox containing all control center elements
        """
        # Create main group box
        group_box = QGroupBox("Control Center")
        group_box.setStyleSheet(f"""
            QGroupBox {{
                font-size: {FONT_SIZE_SECTION_LABEL}px;
                font-weight: bold;
                color: {PRIMARY_TEXT};
                border: 1px solid #D1D1D6;
                border-radius: 8px;
                margin-top: 8px;
                padding-top: 12px;
                background-color: {FRAME_BG_COLOR};
            }}
            QGroupBox::title {{
                subcontrol-origin: margin;
                subcontrol-position: top left;
                left: 10px;
                padding: 0 5px;
            }}
        """)

        layout = QVBoxLayout(group_box)
        layout.setContentsMargins(SPACING_SM, SPACING_MD, SPACING_SM, SPACING_SM)  # Reduced horizontal margins
        layout.setSpacing(SPACING_SM)

        # ===== Main Action Buttons =====
        # Create horizontal layout for three main buttons
        button_layout = QHBoxLayout()
        button_layout.setSpacing(SPACING_XS)

        # Detect button (SUCCESS action - green)
        self.detect_peaks_button = QPushButton("Detect")
        self.detect_peaks_button.setToolTip("Detect peaks in spectrum")
        self.detect_peaks_button.setMinimumHeight(36)
        self.detect_peaks_button.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        self.detect_peaks_button.clicked.connect(self.on_detect_peaks)
        self.detect_peaks_button.setStyleSheet(f"""
            QPushButton {{
                background-color: {SUCCESS_GREEN};
                color: white;
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                font-size: {FONT_SIZE_BODY}px;
                font-weight: bold;
                padding: 6px 8px;
            }}
            QPushButton:hover {{
                background-color: #2AA64A;
            }}
            QPushButton:pressed {{
                background-color: #228B3C;
            }}
            QPushButton:disabled {{
                background-color: #C7C7CC;
                color: #8E8E93;
            }}
        """)
        button_layout.addWidget(self.detect_peaks_button)

        # Fit Spectrum button (PRIMARY action - blue)
        self.fit_all_peaks_button = QPushButton("Fit Spectrum")
        self.fit_all_peaks_button.setToolTip("Fit all peaks in current spectrum")
        self.fit_all_peaks_button.setMinimumHeight(36)
        self.fit_all_peaks_button.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        self.fit_all_peaks_button.clicked.connect(self.on_fit_spectrum)
        self.fit_all_peaks_button.setStyleSheet(f"""
            QPushButton {{
                background-color: {PRIMARY_BUTTON_BG};
                color: {PRIMARY_BUTTON_TEXT};
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                font-size: {FONT_SIZE_BODY}px;
                font-weight: bold;
                padding: 6px 8px;
            }}
            QPushButton:hover {{
                background-color: {PRIMARY_BUTTON_HOVER};
            }}
            QPushButton:pressed {{
                background-color: #3A7CC3;
            }}
            QPushButton:disabled {{
                background-color: #C7C7CC;
                color: #8E8E93;
            }}
        """)
        button_layout.addWidget(self.fit_all_peaks_button)

        # Fit Series button (PRIMARY action - blue)
        self.series_button = QPushButton("Fit Series")
        self.series_button.setToolTip("Run series integration on multiple spectra")
        self.series_button.setMinimumHeight(36)
        self.series_button.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        self.series_button.clicked.connect(self.on_fit_series)
        self.series_button.setStyleSheet(f"""
            QPushButton {{
                background-color: {PRIMARY_BUTTON_BG};
                color: {PRIMARY_BUTTON_TEXT};
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                font-size: {FONT_SIZE_BODY}px;
                font-weight: bold;
                padding: 6px 8px;
            }}
            QPushButton:hover {{
                background-color: {PRIMARY_BUTTON_HOVER};
            }}
            QPushButton:pressed {{
                background-color: #3A7CC3;
            }}
            QPushButton:disabled {{
                background-color: #C7C7CC;
                color: #8E8E93;
            }}
        """)
        button_layout.addWidget(self.series_button)

        layout.addLayout(button_layout)

        # ===== Peak Edition Toggle Button =====
        self.peak_edit_toggle_button = QPushButton("Peak Edition ▼")
        self.peak_edit_toggle_button.setToolTip("Show/hide peak edition controls")
        self.peak_edit_toggle_button.setMinimumHeight(32)
        self.peak_edit_toggle_button.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        self.peak_edit_toggle_button.clicked.connect(self.toggle_peak_edition)
        self.peak_edit_toggle_button.setStyleSheet(f"""
            QPushButton {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {SECONDARY_BUTTON_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                font-size: {FONT_SIZE_BODY}px;
                padding: 4px 6px;
                text-align: left;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
            QPushButton:pressed {{
                background-color: #C0C0C5;
            }}
        """)
        layout.addWidget(self.peak_edit_toggle_button)

        # ===== Collapsible Peak Edition Controls =====
        self.peak_edition_widget = QWidget()
        peak_edition_layout = QVBoxLayout(self.peak_edition_widget)
        peak_edition_layout.setContentsMargins(SPACING_MD, SPACING_SM, SPACING_MD, SPACING_SM)
        peak_edition_layout.setSpacing(SPACING_SM)

        # Mode status label
        self.edit_mode_status_label = QLabel("Mode: View Only")
        self.edit_mode_status_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_SMALL}px;
                color: {SECONDARY_TEXT};
            }}
        """)
        peak_edition_layout.addWidget(self.edit_mode_status_label)

        # Peak list selection checkboxes
        peak_list_layout = QHBoxLayout()
        peak_list_layout.setSpacing(SPACING_SM)

        self.edit_reference_checkbox = QCheckBox("Edit Reference peaks")
        self.edit_reference_checkbox.setStyleSheet(f"""
            QCheckBox {{
                font-size: {FONT_SIZE_BODY}px;
                color: {PRIMARY_TEXT};
            }}
        """)
        self.edit_reference_checkbox.stateChanged.connect(self.on_edit_reference_changed)
        peak_list_layout.addWidget(self.edit_reference_checkbox)

        self.edit_detected_checkbox = QCheckBox("Edit Detected peaks")
        self.edit_detected_checkbox.setStyleSheet(f"""
            QCheckBox {{
                font-size: {FONT_SIZE_BODY}px;
                color: {PRIMARY_TEXT};
            }}
        """)
        self.edit_detected_checkbox.stateChanged.connect(self.on_edit_detected_changed)
        peak_list_layout.addWidget(self.edit_detected_checkbox)

        peak_edition_layout.addLayout(peak_list_layout)

        # Edit mode checkboxes
        mode_layout = QHBoxLayout()
        mode_layout.setSpacing(SPACING_SM)

        self.delete_mode_checkbox = QCheckBox("Delete Mode")
        self.delete_mode_checkbox.setStyleSheet(f"""
            QCheckBox {{
                font-size: {FONT_SIZE_BODY}px;
                color: {PRIMARY_TEXT};
            }}
        """)
        self.delete_mode_checkbox.stateChanged.connect(self.on_deletion_mode_changed)
        mode_layout.addWidget(self.delete_mode_checkbox)

        self.add_mode_checkbox = QCheckBox("Add Mode")
        self.add_mode_checkbox.setStyleSheet(f"""
            QCheckBox {{
                font-size: {FONT_SIZE_BODY}px;
                color: {PRIMARY_TEXT};
            }}
        """)
        self.add_mode_checkbox.stateChanged.connect(self.on_addition_mode_changed)
        mode_layout.addWidget(self.add_mode_checkbox)

        peak_edition_layout.addLayout(mode_layout)

        # Selected peak info label
        self.selected_peak_label = QLabel("No peak selected")
        self.selected_peak_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_SMALL}px;
                color: {SECONDARY_TEXT};
            }}
        """)
        peak_edition_layout.addWidget(self.selected_peak_label)

        # Initially hidden
        self.peak_edition_widget.setVisible(False)

        # Add Peak Edition widget to parent layout
        layout.addWidget(self.peak_edition_widget)

        return group_box

    def create_shift_peak_list_section(self) -> QGroupBox:
        """Create the Shift Peak List section with two-row layout.

        This section allows users to apply coordinate offsets to all peaks
        in the loaded peak list. Common use case: correcting systematic shifts
        between experimental and reference peak positions.

        Layout:
            - Row 1: 1H spinbox | 15N/13C spinbox
            - Row 2: Apply button | Reset button

        Returns:
            QGroupBox containing the shift controls
        """
        group = QGroupBox("Shift Peak List")
        group.setStyleSheet(f"""
            QGroupBox {{
                font-size: {FONT_SIZE_SECTION_LABEL}px;
                font-weight: bold;
                color: {PRIMARY_TEXT};
                background-color: {FRAME_BG_COLOR};
                border: 1px solid #ddd;
                border-radius: 6px;
                margin-top: 8px;
                padding-top: 10px;
            }}
            QGroupBox::title {{
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 5px;
            }}
        """)

        # Main vertical layout
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)
        main_layout.setSpacing(SPACING_XS)

        # Row 1: Spinboxes
        spinbox_layout = QHBoxLayout()
        spinbox_layout.setSpacing(SPACING_XS)

        # 1H offset control
        h_label = QLabel("1H:")
        h_label.setStyleSheet(f"font-size: {FONT_SIZE_BODY}px; color: {PRIMARY_TEXT};")
        spinbox_layout.addWidget(h_label)

        self.shift_1h_spin = QDoubleSpinBox()
        self.shift_1h_spin.setRange(-10.0, 10.0)
        self.shift_1h_spin.setValue(self.shift_1h_value)
        self.shift_1h_spin.setSingleStep(0.001)
        self.shift_1h_spin.setDecimals(3)
        self.shift_1h_spin.setSuffix(" ppm")
        self.shift_1h_spin.setStyleSheet(f"""
            QDoubleSpinBox {{
                font-size: {FONT_SIZE_BODY}px;
                padding: 2px;
                border: 1px solid #ccc;
                border-radius: 4px;
            }}
            QDoubleSpinBox:focus {{
                border: 1px solid #0078d4;
            }}
        """)
        self.shift_1h_spin.valueChanged.connect(self._on_shift_value_changed)
        spinbox_layout.addWidget(self.shift_1h_spin)

        # 15N/13C offset control
        n_label = QLabel("15N/13C:")
        n_label.setStyleSheet(f"font-size: {FONT_SIZE_BODY}px; color: {PRIMARY_TEXT};")
        spinbox_layout.addWidget(n_label)

        self.shift_15n_spin = QDoubleSpinBox()
        self.shift_15n_spin.setRange(-50.0, 50.0)
        self.shift_15n_spin.setValue(self.shift_15n_value)
        self.shift_15n_spin.setSingleStep(0.01)
        self.shift_15n_spin.setDecimals(2)
        self.shift_15n_spin.setSuffix(" ppm")
        self.shift_15n_spin.setStyleSheet(f"""
            QDoubleSpinBox {{
                font-size: {FONT_SIZE_BODY}px;
                padding: 2px;
                border: 1px solid #ccc;
                border-radius: 4px;
            }}
            QDoubleSpinBox:focus {{
                border: 1px solid #0078d4;
            }}
        """)
        self.shift_15n_spin.valueChanged.connect(self._on_shift_value_changed)
        spinbox_layout.addWidget(self.shift_15n_spin)

        spinbox_layout.addStretch()
        main_layout.addLayout(spinbox_layout)

        # Row 2: Buttons
        button_layout = QHBoxLayout()
        button_layout.setSpacing(SPACING_XS)

        # Apply button (primary style - blue)
        apply_btn = QPushButton("Apply")
        apply_btn.setToolTip("Apply reference shift to peak positions")
        apply_btn.setProperty("class", "primary")
        apply_btn.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        apply_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {PRIMARY_BUTTON_BG};
                color: {PRIMARY_BUTTON_TEXT};
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: 4px 6px;
                font-size: {FONT_SIZE_BODY}px;
                font-weight: bold;
            }}
            QPushButton:hover {{
                background-color: {PRIMARY_BUTTON_HOVER};
            }}
            QPushButton:pressed {{
                background-color: #005a9e;
            }}
            QPushButton:disabled {{
                background-color: #cccccc;
                color: #666666;
            }}
        """)
        apply_btn.clicked.connect(self.on_apply_shift)
        button_layout.addWidget(apply_btn)

        # Reset button (secondary style - grey)
        reset_btn = QPushButton("Reset")
        reset_btn.setToolTip("Reset shift values to zero")
        reset_btn.setProperty("class", "secondary")
        reset_btn.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        reset_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {SECONDARY_BUTTON_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: 4px 6px;
                font-size: {FONT_SIZE_BODY}px;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
            QPushButton:pressed {{
                background-color: #d0d0d0;
            }}
        """)
        reset_btn.clicked.connect(self.on_reset_shift)
        button_layout.addWidget(reset_btn)

        button_layout.addStretch()
        main_layout.addLayout(button_layout)

        group.setLayout(main_layout)
        return group

    def _on_shift_value_changed(self, value: float):
        """Handle changes to shift spinbox values.

        Args:
            value: New spinbox value (not used directly, read from spinboxes)
        """
        # Update internal state
        if self.shift_1h_spin:
            self.shift_1h_value = self.shift_1h_spin.value()
        if self.shift_15n_spin:
            self.shift_15n_value = self.shift_15n_spin.value()

    def _on_intensity_scale_change_3d(self, value: int):
        """Handle intensity scale slider change for 3D Voigt plot.

        Based on v0o9 line 815-827.

        Args:
            value: Slider value (50-200)
        """
        # Update label
        self.intensity_scale_3d_label.setText(f"{value}%")

        # Debounce: Only update plot after slider stops moving
        # Qt QTimer.singleShot provides a simple debounce mechanism
        if hasattr(self, '_intensity_scale_timer_3d'):
            # Cancel previous timer if exists
            try:
                self._intensity_scale_timer_3d.stop()
            except:
                pass

        # Schedule update after 100ms delay
        from PySide6.QtCore import QTimer
        self._intensity_scale_timer_3d = QTimer()
        self._intensity_scale_timer_3d.setSingleShot(True)
        self._intensity_scale_timer_3d.timeout.connect(
            lambda: self._apply_intensity_scale_3d(float(value))
        )
        self._intensity_scale_timer_3d.start(100)

        # Update internal state
        self.intensity_scale_3d = float(value)

    def _on_color_scheme_change_3d(self, scheme: str):
        """Handle color scheme dropdown change for 3D Voigt plot.

        Based on v0o9 line 829-831.

        Args:
            scheme: Color scheme name ('Classic', 'Clean', 'Dark', 'Warm')
        """
        self.color_scheme_3d = scheme
        if hasattr(self, 'voigt_3d_plotter') and self.voigt_3d_plotter:
            self.voigt_3d_plotter.set_color_scheme(scheme)

    def _apply_intensity_scale_3d(self, value: float):
        """Apply intensity scale to 3D Voigt plotter (called after debounce)."""
        if hasattr(self, 'voigt_3d_plotter') and self.voigt_3d_plotter:
            self.voigt_3d_plotter.set_intensity_scale(value)
        logger.debug(f"3D Voigt: Intensity Scale = {value}%")

    # =================== 3D VOIGT LAYER VISIBILITY HANDLERS ===================
    # Use toggled(bool) signal for checkboxes - matches v0.9 behavior exactly

    def _on_show_exp_3d_change(self, checked: bool):
        """Handle Experimental checkbox toggle for 3D Voigt plot (v0o9 line 1855)."""
        self.show_exp_3d = checked
        if hasattr(self, 'voigt_3d_plotter') and self.voigt_3d_plotter:
            self.voigt_3d_plotter.toggle_experimental(checked)
        logger.debug(f"3D Voigt: Show Experimental = {checked}")

    def _on_show_fit_3d_change(self, checked: bool):
        """Handle Fitted checkbox toggle for 3D Voigt plot (v0o9 line 1858)."""
        self.show_fit_3d = checked
        if hasattr(self, 'voigt_3d_plotter') and self.voigt_3d_plotter:
            self.voigt_3d_plotter.toggle_fitted(checked)
        logger.debug(f"3D Voigt: Show Fitted = {checked}")

    def _on_show_individual_3d_change(self, checked: bool):
        """Handle Individual Peaks checkbox toggle for 3D Voigt plot (v0o9 line 1861)."""
        self.show_individual_3d = checked
        if hasattr(self, 'voigt_3d_plotter') and self.voigt_3d_plotter:
            self.voigt_3d_plotter.toggle_individual_peaks(checked)
        logger.debug(f"3D Voigt: Show Individual Peaks = {checked}")

    def _on_show_resid_3d_change(self, checked: bool):
        """Handle Residuals checkbox toggle for 3D Voigt plot (v0o9 line 1867)."""
        self.show_resid_3d = checked
        if hasattr(self, 'voigt_3d_plotter') and self.voigt_3d_plotter:
            self.voigt_3d_plotter.toggle_residuals(checked)
        logger.debug(f"3D Voigt: Show Residuals = {checked}")

    def _on_limit_peak_display_3d_change(self, checked: bool):
        """Handle Limit Peak Extent checkbox toggle for 3D Voigt plot (v0o9 line 1870)."""
        self.limit_peak_display_3d = checked
        if hasattr(self, 'voigt_3d_plotter') and self.voigt_3d_plotter:
            self.voigt_3d_plotter.toggle_peak_clipping(checked)
        logger.debug(f"3D Voigt: Limit Peak Extent = {checked}")

    def on_apply_shift(self):
        """Apply coordinate offsets to all peaks in the current peak list.

        This method:
        1. Validates that a peak list is loaded
        2. Checks that at least one offset is non-zero
        3. Applies offsets to all peaks (both DataFrame and list formats)
        4. Updates the visualization
        5. Shows status message

        Based on original implementation in main_gui.py:3506-3576
        """
        logger.debug("Apply shift button clicked")

        # Check if integrator has a peak list
        if not hasattr(self.integrator, 'peak_list') or self.integrator.peak_list is None:
            QMessageBox.critical(
                self,
                "Error",
                "No peak list loaded. Please load a peak list first."
            )
            return

        # Get current offset values
        h_shift = self.shift_1h_spin.value()
        n_shift = self.shift_15n_spin.value()

        # Check if offsets are zero
        if h_shift == 0.0 and n_shift == 0.0:
            QMessageBox.warning(
                self,
                "Warning",
                "Both offsets are zero. No adjustment needed."
            )
            return

        try:
            # Handle different peak list formats
            peak_list = self.integrator.peak_list

            # Check if peak list is empty
            if hasattr(peak_list, 'empty') and peak_list.empty:
                QMessageBox.critical(
                    self,
                    "Error",
                    "Peak list is empty. Please load a peak list first."
                )
                return
            elif isinstance(peak_list, list) and len(peak_list) == 0:
                QMessageBox.critical(
                    self,
                    "Error",
                    "Peak list is empty. Please load a peak list first."
                )
                return

            # Apply offsets based on format
            if hasattr(peak_list, 'columns'):
                # DataFrame format (reference peaks)
                peak_list = peak_list.copy()
                original_count = len(peak_list)

                # Apply offsets to columns
                if 'Position_X' in peak_list.columns:
                    peak_list['Position_X'] += h_shift
                if 'Position_Y' in peak_list.columns:
                    peak_list['Position_Y'] += n_shift

                # Update the integrator's peak list
                self.integrator.peak_list = peak_list
            else:
                # List format (detected peaks)
                import copy
                peak_list = copy.deepcopy(peak_list)
                original_count = len(peak_list)

                # Apply offsets to list of dictionaries
                for peak in peak_list:
                    if isinstance(peak, dict):
                        if 'Position_X' in peak:
                            peak['Position_X'] += h_shift
                        if 'Position_Y' in peak:
                            peak['Position_Y'] += n_shift

                # Update the integrator's peak list
                self.integrator.peak_list = peak_list

            # Update status
            status_msg = f"Peak list shifted by {h_shift:+.3f} ppm (1H), {n_shift:+.2f} ppm (15N/13C) - {original_count} peaks"
            self.update_status(status_msg)

            # Update visualization after peak list shift (v0.9 equivalent: update_main_plot())
            self.update_spectrum_plot()
            self.update_statistics_panel()

            logger.info(f"Applied coordinate offsets: ΔX={h_shift:+.3f} ppm, ΔY={n_shift:+.2f} ppm to {original_count} peaks")

        except Exception as e:
            error_msg = f"Failed to apply coordinate offsets: {str(e)}"
            logger.error(error_msg)
            QMessageBox.critical(
                self,
                "Adjustment Error",
                error_msg
            )
            self.update_status(f"Error: {error_msg}")

    def on_reset_shift(self):
        """Reset offset values to zero.

        This simply resets the spinbox values without modifying the peak list.
        Users can reapply the original offsets if needed by reloading the peak list.

        Based on original implementation in main_gui.py:3577-3581
        """
        logger.debug("Reset shift button clicked")

        # Reset spinbox values
        self.shift_1h_spin.setValue(0.0)
        self.shift_15n_spin.setValue(0.0)

        # Update internal state
        self.shift_1h_value = 0.0
        self.shift_15n_value = 0.0

        # Update status
        self.update_status("Peak list shift reset to 0")
        logger.info("Reset coordinate offsets to zero")

        peak_list_layout.addWidget(self.edit_reference_checkbox)

    def create_spectrum_display_section(self) -> QGroupBox:
        """Create the Spectrum Display Controls section.

        This section contains:
        - Contour Settings button (collapsible)
        - Reset Zoom button
        - Collapsible contour parameters (Levels, Min Level, Increment)

        Returns:
            QGroupBox containing spectrum display controls
        """
        # Create group box
        group = QGroupBox("Spectrum Display Controls")
        group.setStyleSheet(f"""
            QGroupBox {{
                font-size: {FONT_SIZE_SECTION_LABEL}px;
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
        """)

        # Main layout
        layout = QVBoxLayout(group)
        layout.setContentsMargins(SPACING_SM, SPACING_MD, SPACING_SM, SPACING_SM)  # Reduced horizontal margins
        layout.setSpacing(SPACING_SM)

        # Button row: Contour Settings and Reset Zoom
        button_layout = QHBoxLayout()
        button_layout.setSpacing(SPACING_XS)

        # Contour Settings button (toggle)
        self.contour_toggle_button = QPushButton("Contour Settings ▼")
        self.contour_toggle_button.setToolTip("Show/hide contour display settings")
        self.contour_toggle_button.setProperty("class", "secondary")
        self.contour_toggle_button.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        self.contour_toggle_button.setStyleSheet(f"""
            QPushButton {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {SECONDARY_BUTTON_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: 4px 6px;
                font-size: {FONT_SIZE_BODY}px;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
            QPushButton:pressed {{
                background-color: {SECONDARY_BUTTON_BORDER};
            }}
        """)
        self.contour_toggle_button.clicked.connect(self.toggle_contour_settings)
        button_layout.addWidget(self.contour_toggle_button)

        # Reset Zoom button
        reset_zoom_button = QPushButton("Reset Zoom")
        reset_zoom_button.setToolTip("Reset view to full spectrum")
        reset_zoom_button.setProperty("class", "secondary")
        reset_zoom_button.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        reset_zoom_button.setStyleSheet(f"""
            QPushButton {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {SECONDARY_BUTTON_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: 4px 6px;
                font-size: {FONT_SIZE_BODY}px;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
            QPushButton:pressed {{
                background-color: {SECONDARY_BUTTON_BORDER};
            }}
        """)
        reset_zoom_button.clicked.connect(self.on_reset_zoom)
        button_layout.addWidget(reset_zoom_button)

        layout.addLayout(button_layout)

        # Collapsible contour parameters widget
        self.contour_params_widget = QWidget()
        params_layout = QFormLayout(self.contour_params_widget)
        params_layout.setContentsMargins(0, SPACING_SM, 0, 0)
        params_layout.setSpacing(SPACING_XS)
        params_layout.setLabelAlignment(Qt.AlignRight)

        # Levels spinbox (5-50, default 24)
        self.contour_levels_spin = QSpinBox()
        self.contour_levels_spin.setRange(5, 50)
        self.contour_levels_spin.setValue(self.contour_levels)
        self.contour_levels_spin.valueChanged.connect(self.on_contour_levels_changed)
        self.contour_levels_spin.setStyleSheet(f"""
            QSpinBox {{
                background-color: {FRAME_BG_COLOR};
                color: {PRIMARY_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS//2}px;
                padding: {SPACING_XS}px;
                font-size: {FONT_SIZE_BODY}px;
            }}
        """)
        params_layout.addRow("Levels:", self.contour_levels_spin)

        # Min Level spinbox (0.01-1.0, default 0.05)
        self.contour_min_spin = QDoubleSpinBox()
        self.contour_min_spin.setRange(0.01, 1.0)
        self.contour_min_spin.setValue(self.contour_min)
        self.contour_min_spin.setSingleStep(0.01)
        self.contour_min_spin.setDecimals(3)
        self.contour_min_spin.valueChanged.connect(self.on_contour_min_changed)
        self.contour_min_spin.setStyleSheet(f"""
            QDoubleSpinBox {{
                background-color: {FRAME_BG_COLOR};
                color: {PRIMARY_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS//2}px;
                padding: {SPACING_XS}px;
                font-size: {FONT_SIZE_BODY}px;
            }}
        """)
        params_layout.addRow("Min Level:", self.contour_min_spin)

        # Increment spinbox (0.01-10.0, default 1.1)
        self.contour_inc_spin = QDoubleSpinBox()
        self.contour_inc_spin.setRange(0.01, 10.0)
        self.contour_inc_spin.setValue(self.contour_increment)
        self.contour_inc_spin.setSingleStep(0.1)
        self.contour_inc_spin.setDecimals(2)
        self.contour_inc_spin.valueChanged.connect(self.on_contour_increment_changed)
        self.contour_inc_spin.setStyleSheet(f"""
            QDoubleSpinBox {{
                background-color: {FRAME_BG_COLOR};
                color: {PRIMARY_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS//2}px;
                padding: {SPACING_XS}px;
                font-size: {FONT_SIZE_BODY}px;
            }}
        """)
        params_layout.addRow("Increment:", self.contour_inc_spin)

        # Initially hidden
        self.contour_params_widget.setVisible(False)
        layout.addWidget(self.contour_params_widget)

        return group

    def toggle_contour_settings(self):
        """Toggle visibility of contour parameters."""
        visible = self.contour_params_widget.isVisible()

        if visible:
            self.contour_params_widget.hide()
            arrow = "▼"
            logger.debug("Contour settings hidden")
        else:
            self.contour_params_widget.show()
            arrow = "▶"
            logger.debug("Contour settings shown")

        self.contour_toggle_button.setText(f"Contour Settings {arrow}")

    def on_reset_zoom(self):
        """Handle Reset Zoom button click.

        Resets the spectrum view to show the full data range.
        Based on v0.9 reset_view() (main_gui.py:6243-6249).
        """
        logger.info("Reset Zoom requested")

        if hasattr(self, 'spectrum_plotter') and self.spectrum_plotter:
            self.spectrum_plotter.reset_zoom()
            self.update_status("Zoom reset to full view")
        else:
            self.update_status("No spectrum loaded")

    def on_contour_levels_changed(self, value: int):
        """Handle contour levels value change.

        Updates the spectrum plot with new contour levels setting.
        Based on v0.9 main_gui.py:1183-1186.
        """
        self.contour_levels = value
        logger.debug(f"Contour levels changed to {value}")
        self.update_spectrum_plot()

    def on_contour_min_changed(self, value: float):
        """Handle contour min level value change.

        Updates the spectrum plot with new minimum contour level.
        Based on v0.9 main_gui.py:1190-1193.
        """
        self.contour_min = value
        logger.debug(f"Contour min level changed to {value}")
        self.update_spectrum_plot()

    def on_contour_increment_changed(self, value: float):
        """Handle contour increment value change.

        Updates the spectrum plot with new contour increment.
        Based on v0.9 main_gui.py:1197-1200.
        """
        self.contour_increment = value
        logger.debug(f"Contour increment changed to {value}")
        self.update_spectrum_plot()

    def create_statistics_panel_section(self) -> QGroupBox:
        """Create the Processing Statistics section.

        This section displays statistics from peak detection and fitting operations,
        including peak counts, fitting success rates, and processing time.

        Returns:
            QGroupBox containing statistics display
        """
        group = QGroupBox("Processing Statistics")
        group.setStyleSheet(f"""
            QGroupBox {{
                font-size: {FONT_SIZE_SECTION_LABEL}px;
                font-weight: bold;
                color: {PRIMARY_TEXT};
                border: 1px solid #D1D1D6;
                border-radius: 8px;
                margin-top: 8px;
                padding-top: 12px;
                background-color: {FRAME_BG_COLOR};
            }}
            QGroupBox::title {{
                subcontrol-origin: margin;
                subcontrol-position: top left;
                left: 10px;
                padding: 0 5px;
            }}
        """)

        layout = QVBoxLayout(group)
        layout.setContentsMargins(SPACING_SM, SPACING_MD, SPACING_SM, SPACING_SM)  # Reduced horizontal margins
        layout.setSpacing(SPACING_SM)

        # Statistics label
        self.stats_label = QLabel("No processing performed yet")
        self.stats_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_BODY}px;
                color: {PRIMARY_TEXT};
                padding: {SPACING_SM}px;
            }}
        """)
        self.stats_label.setWordWrap(True)
        layout.addWidget(self.stats_label)

        return group

    def create_advanced_options_section(self) -> QGroupBox:
        """Create the Advanced Options section with Expert Mode button.

        This section provides access to advanced configuration options
        for expert users who need fine-grained control over fitting parameters.

        Returns:
            QGroupBox containing advanced options controls
        """
        group = QGroupBox("Advanced Options")
        group.setStyleSheet(f"""
            QGroupBox {{
                font-size: {FONT_SIZE_SECTION_LABEL}px;
                font-weight: bold;
                color: {PRIMARY_TEXT};
                border: 1px solid #D1D1D6;
                border-radius: 8px;
                margin-top: 8px;
                padding-top: 12px;
                background-color: {FRAME_BG_COLOR};
            }}
            QGroupBox::title {{
                subcontrol-origin: margin;
                subcontrol-position: top left;
                left: 10px;
                padding: 0 5px;
            }}
        """)

        layout = QVBoxLayout(group)
        layout.setContentsMargins(SPACING_SM, SPACING_MD, SPACING_SM, SPACING_SM)  # Reduced horizontal margins
        layout.setSpacing(SPACING_SM)

        # Expert Mode button
        expert_btn = QPushButton("Expert Mode")
        expert_btn.setToolTip("Open advanced parameter settings")
        expert_btn.setProperty("class", "secondary")
        expert_btn.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        expert_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {SECONDARY_BUTTON_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: 4px 6px;
                font-size: {FONT_SIZE_BODY}px;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
            QPushButton:pressed {{
                background-color: {SECONDARY_BUTTON_BORDER};
            }}
        """)
        expert_btn.clicked.connect(self.on_expert_mode)
        layout.addWidget(expert_btn)

        return group

    def create_voigt_fitting_section(self) -> QGroupBox:
        """Create the Voigt Fitting section with Reset Results button.

        This section contains operations related to Voigt profile fitting results,
        including clearing results and resetting fitting parameters.

        Returns:
            QGroupBox containing Voigt fitting controls
        """
        group = QGroupBox("Voigt Fitting")
        group.setStyleSheet(f"""
            QGroupBox {{
                font-size: {FONT_SIZE_SECTION_LABEL}px;
                font-weight: bold;
                color: {PRIMARY_TEXT};
                border: 1px solid #D1D1D6;
                border-radius: 8px;
                margin-top: 8px;
                padding-top: 12px;
                background-color: {FRAME_BG_COLOR};
            }}
            QGroupBox::title {{
                subcontrol-origin: margin;
                subcontrol-position: top left;
                left: 10px;
                padding: 0 5px;
            }}
        """)

        layout = QVBoxLayout(group)
        layout.setContentsMargins(SPACING_SM, SPACING_MD, SPACING_SM, SPACING_SM)  # Reduced horizontal margins
        layout.setSpacing(SPACING_SM)

        # Reset Results button (DESTRUCTIVE - red)
        reset_btn = QPushButton("Reset Results")
        reset_btn.setToolTip("Clear all fitting results")
        reset_btn.setProperty("class", "destructive")
        reset_btn.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        reset_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {DESTRUCTIVE_BUTTON_BG};
                color: {DESTRUCTIVE_BUTTON_TEXT};
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: 6px 8px;
                font-size: {FONT_SIZE_BODY}px;
                font-weight: bold;
            }}
            QPushButton:hover {{
                background-color: {DESTRUCTIVE_BUTTON_HOVER};
            }}
            QPushButton:pressed {{
                background-color: #C44539;
            }}
            QPushButton:disabled {{
                background-color: #C7C7CC;
                color: #8E8E93;
            }}
        """)
        reset_btn.clicked.connect(self.on_reset_results)
        layout.addWidget(reset_btn)

        return group

    def on_expert_mode(self):
        """Handle Expert Mode button click (v0o9 pattern line 2832).

        Opens a full Expert Mode dialog with all advanced parameters matching v0.9.
        All parameters are interactive and immediately affect processing.
        """
        logger.info("Expert Mode button clicked")

        # Reuse existing dialog if already created (prevents duplicate signal connections)
        if hasattr(self, '_expert_dialog') and self._expert_dialog is not None:
            # Update checkbox states to reflect current values before showing
            self._sync_expert_dialog_values()
            self._expert_dialog.show()
            self._expert_dialog.raise_()
            self._expert_dialog.activateWindow()
            return

        # Create parameter dialog (only once)
        dialog = QDialog(self)
        self._expert_dialog = dialog  # Cache the dialog
        dialog.setWindowTitle("Expert Mode - Advanced Parameters")
        dialog.setMinimumSize(700, 600)
        dialog.setModal(True)

        # Main layout
        main_layout = QVBoxLayout(dialog)
        main_layout.setSpacing(SPACING_SM)
        main_layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)

        # Scroll area for parameters
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        # Content widget
        content = QWidget()
        content_layout = QVBoxLayout(content)
        content_layout.setSpacing(SPACING_MD)

        # =================== PEAK DETECTION PARAMETERS SECTION ===================
        detection_group = QGroupBox("Peak Detection Parameters")
        detection_layout = QVBoxLayout(detection_group)

        # Parameters toggle button
        self.expert_params_toggle_btn = QPushButton("Parameters ▼")
        self.expert_params_toggle_btn.setToolTip("Show/hide detection parameter settings")
        self.expert_params_toggle_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {SECONDARY_BUTTON_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: 6px;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
        """)
        detection_layout.addWidget(self.expert_params_toggle_btn)

        # Collapsible detection parameters frame
        self.expert_params_frame = QWidget()
        params_frame_layout = QGridLayout(self.expert_params_frame)
        params_frame_layout.setSpacing(SPACING_SM)

        # Row 0: 1H/15N search window
        params_frame_layout.addWidget(QLabel("1H/15N (ppm):"), 0, 0)
        self.expert_search_x_spin = QDoubleSpinBox()
        self.expert_search_x_spin.setRange(0.01, 0.2)
        self.expert_search_x_spin.setSingleStep(0.01)
        self.expert_search_x_spin.setDecimals(3)
        self.expert_search_x_spin.setValue(self.search_window_x)
        self.expert_search_x_spin.valueChanged.connect(self._on_expert_param_change)
        params_frame_layout.addWidget(self.expert_search_x_spin, 0, 1)

        self.expert_search_y_spin = QDoubleSpinBox()
        self.expert_search_y_spin.setRange(0.05, 4.0)
        self.expert_search_y_spin.setSingleStep(0.05)
        self.expert_search_y_spin.setDecimals(2)
        self.expert_search_y_spin.setValue(self.search_window_y)
        self.expert_search_y_spin.valueChanged.connect(self._on_expert_param_change)
        params_frame_layout.addWidget(self.expert_search_y_spin, 0, 2)

        # Row 1: Noise Threshold
        params_frame_layout.addWidget(QLabel("Noise Threshold:"), 1, 0)
        self.expert_noise_spin = QDoubleSpinBox()
        self.expert_noise_spin.setRange(0.01, 10.0)
        self.expert_noise_spin.setSingleStep(0.1)
        self.expert_noise_spin.setDecimals(2)
        self.expert_noise_spin.setValue(self.noise_threshold)
        self.expert_noise_spin.valueChanged.connect(self._on_expert_param_change)
        params_frame_layout.addWidget(self.expert_noise_spin, 1, 1)

        # Row 2: Detection rectangle size (pixels)
        params_frame_layout.addWidget(QLabel("X×Y (pixels):"), 2, 0)
        self.expert_det_x_spin = QSpinBox()
        self.expert_det_x_spin.setRange(1, 9)
        self.expert_det_x_spin.setSingleStep(2)
        self.expert_det_x_spin.setValue(self.detection_square_size)
        self.expert_det_x_spin.valueChanged.connect(self._on_expert_detection_size_change)
        params_frame_layout.addWidget(self.expert_det_x_spin, 2, 1)

        self.expert_det_y_spin = QSpinBox()
        self.expert_det_y_spin.setRange(1, 5)
        self.expert_det_y_spin.setSingleStep(1)
        self.expert_det_y_spin.setValue(self.detection_rectangle_y)
        self.expert_det_y_spin.valueChanged.connect(self._on_expert_detection_size_change)
        params_frame_layout.addWidget(self.expert_det_y_spin, 2, 2)

        self.expert_det_ppm_label = QLabel(self.detection_square_ppm_x)
        self.expert_det_ppm_label.setStyleSheet("font-size: 10px; color: gray;")
        params_frame_layout.addWidget(self.expert_det_ppm_label, 2, 3)

        # Row 3: Peak Centroid Detection header
        centroid_label = QLabel("🎯 Peak Centroid Detection:")
        centroid_label.setStyleSheet("font-weight: bold;")
        params_frame_layout.addWidget(centroid_label, 3, 0, 1, 4)

        # Row 4: Centroid window ppm
        params_frame_layout.addWidget(QLabel("Window X ppm:"), 4, 0)
        self.expert_centroid_x_spin = QDoubleSpinBox()
        self.expert_centroid_x_spin.setRange(0.01, 0.2)
        self.expert_centroid_x_spin.setSingleStep(0.01)
        self.expert_centroid_x_spin.setDecimals(3)
        self.expert_centroid_x_spin.setValue(self.centroid_window_x_ppm)
        self.expert_centroid_x_spin.valueChanged.connect(self._on_expert_param_change)
        params_frame_layout.addWidget(self.expert_centroid_x_spin, 4, 1)

        self.expert_centroid_y_spin = QDoubleSpinBox()
        self.expert_centroid_y_spin.setRange(0.01, 0.5)
        self.expert_centroid_y_spin.setSingleStep(0.02)
        self.expert_centroid_y_spin.setDecimals(3)
        self.expert_centroid_y_spin.setValue(self.centroid_window_y_ppm)
        self.expert_centroid_y_spin.valueChanged.connect(self._on_expert_param_change)
        params_frame_layout.addWidget(self.expert_centroid_y_spin, 4, 2)

        # Row 5: Auto-add dummy peaks checkbox (peak_list mode only)
        self.expert_auto_add_dummy_cb = QCheckBox("Auto-add nearby unmatched peaks")
        self.expert_auto_add_dummy_cb.setChecked(self.auto_add_dummy_peaks)
        self.expert_auto_add_dummy_cb.setToolTip(
            "Automatically add peaks detected in spectrum but not in user peak list.\n"
            "Uses 3rd percentile intensity threshold and 3× overlap proximity.\n"
            "Only applies to Peak List mode (grayed out in S/N Threshold mode)."
        )
        self.expert_auto_add_dummy_cb.stateChanged.connect(self._on_auto_add_dummy_change)
        params_frame_layout.addWidget(self.expert_auto_add_dummy_cb, 5, 0, 1, 4)

        self.expert_params_frame.setVisible(False)
        detection_layout.addWidget(self.expert_params_frame)

        # Connect toggle button
        self.expert_params_toggle_btn.clicked.connect(self._toggle_expert_params)

        content_layout.addWidget(detection_group)

        # =================== PS2D ALGORITHM CONFIGURATION SECTION ===================
        ps2d_config_group = QGroupBox("PS2D Algorithm Configuration")
        ps2d_config_layout = QVBoxLayout(ps2d_config_group)

        # Nucleus Type - always visible
        nucleus_row = QHBoxLayout()
        nucleus_label = QLabel("Nucleus Type:")
        nucleus_label.setStyleSheet("font-weight: bold;")
        nucleus_row.addWidget(nucleus_label)

        self.expert_nucleus_combo = QComboBox()
        self.expert_nucleus_combo.addItems(['15N', '13C'])
        self.expert_nucleus_combo.setCurrentText(self.nucleus_type)
        self.expert_nucleus_combo.currentTextChanged.connect(self._on_nucleus_type_change)
        nucleus_row.addWidget(self.expert_nucleus_combo)

        nucleus_hint = QLabel("(for 2D overlap fitting)")
        nucleus_hint.setStyleSheet("font-size: 10px; color: gray;")
        nucleus_row.addWidget(nucleus_hint)
        nucleus_row.addStretch()
        ps2d_config_layout.addLayout(nucleus_row)

        # Button row: Apply Changes and Parameters toggle
        button_row = QHBoxLayout()
        apply_btn = QPushButton("Apply Changes")
        apply_btn.setToolTip("Apply PS2D configuration changes")
        apply_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {PRIMARY_BUTTON_BG};
                color: {PRIMARY_BUTTON_TEXT};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: 6px 12px;
            }}
            QPushButton:hover {{
                background-color: {PRIMARY_BUTTON_HOVER};
            }}
        """)
        apply_btn.clicked.connect(self._on_ps2d_params_apply)
        button_row.addWidget(apply_btn)

        self.expert_ps2d_config_toggle_btn = QPushButton("Parameters ▼")
        self.expert_ps2d_config_toggle_btn.setToolTip("Show/hide PS2D configuration parameters")
        self.expert_ps2d_config_toggle_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {SECONDARY_BUTTON_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: 6px;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
        """)
        self.expert_ps2d_config_toggle_btn.clicked.connect(self._toggle_ps2d_config_params)
        button_row.addWidget(self.expert_ps2d_config_toggle_btn)
        ps2d_config_layout.addLayout(button_row)

        # Collapsible PS2D config parameters frame
        self.expert_ps2d_config_frame = QWidget()
        ps2d_params_layout = QGridLayout(self.expert_ps2d_config_frame)
        ps2d_params_layout.setSpacing(SPACING_SM)

        # Row 0: Ellipse Radii
        radii_label = QLabel("Ellipse Radii:")
        radii_label.setStyleSheet("font-weight: bold;")
        ps2d_params_layout.addWidget(radii_label, 0, 0)

        ps2d_params_layout.addWidget(QLabel("F1:"), 0, 1)
        self.expert_radF1_spin = QDoubleSpinBox()
        self.expert_radF1_spin.setRange(0.01, 2.0)
        self.expert_radF1_spin.setSingleStep(0.01)
        self.expert_radF1_spin.setDecimals(3)
        self.expert_radF1_spin.setValue(self.ps2d_radF1)
        ps2d_params_layout.addWidget(self.expert_radF1_spin, 0, 2)

        ps2d_params_layout.addWidget(QLabel("F2:"), 0, 3)
        self.expert_radF2_spin = QDoubleSpinBox()
        self.expert_radF2_spin.setRange(0.001, 0.5)
        self.expert_radF2_spin.setSingleStep(0.005)
        self.expert_radF2_spin.setDecimals(4)
        self.expert_radF2_spin.setValue(self.ps2d_radF2)
        ps2d_params_layout.addWidget(self.expert_radF2_spin, 0, 4)

        ps2d_params_layout.addWidget(QLabel("(ppm)"), 0, 5)

        # Row 1: Show Ellipses checkbox
        self.expert_show_ellipses_check = QCheckBox("Show Ellipses")
        self.expert_show_ellipses_check.setChecked(self.show_ellipses)
        self.expert_show_ellipses_check.stateChanged.connect(self._on_show_ellipses_change)
        ps2d_params_layout.addWidget(self.expert_show_ellipses_check, 1, 0, 1, 3)

        # Row 2: Overlap Thresholds
        overlap_label = QLabel("Overlap Thresholds:")
        overlap_label.setStyleSheet("font-weight: bold;")
        ps2d_params_layout.addWidget(overlap_label, 2, 0)

        ps2d_params_layout.addWidget(QLabel("X:"), 2, 1)
        self.expert_overlap_x_spin = QDoubleSpinBox()
        self.expert_overlap_x_spin.setRange(0.01, 1.0)
        self.expert_overlap_x_spin.setSingleStep(0.01)
        self.expert_overlap_x_spin.setDecimals(3)
        self.expert_overlap_x_spin.setValue(self.ps2d_overlap_x)
        ps2d_params_layout.addWidget(self.expert_overlap_x_spin, 2, 2)

        ps2d_params_layout.addWidget(QLabel("Y:"), 2, 3)
        self.expert_overlap_y_spin = QDoubleSpinBox()
        self.expert_overlap_y_spin.setRange(0.01, 2.0)
        self.expert_overlap_y_spin.setSingleStep(0.05)
        self.expert_overlap_y_spin.setDecimals(3)
        self.expert_overlap_y_spin.setValue(self.ps2d_overlap_y)
        ps2d_params_layout.addWidget(self.expert_overlap_y_spin, 2, 4)

        ps2d_params_layout.addWidget(QLabel("(ppm)"), 2, 5)

        # Row 3: Max Iterations
        max_iter_label = QLabel("Max Iterations:")
        max_iter_label.setStyleSheet("font-weight: bold;")
        ps2d_params_layout.addWidget(max_iter_label, 3, 0)

        self.expert_max_iter_spin = QSpinBox()
        self.expert_max_iter_spin.setRange(50, 2000)
        self.expert_max_iter_spin.setSingleStep(50)
        self.expert_max_iter_spin.setValue(self.ps2d_max_iterations)
        ps2d_params_layout.addWidget(self.expert_max_iter_spin, 3, 2)

        # Row 4: Max Cluster Size
        cluster_size_label = QLabel("Max Cluster Size:")
        cluster_size_label.setStyleSheet("font-weight: bold;")
        ps2d_params_layout.addWidget(cluster_size_label, 4, 0)

        self.expert_max_cluster_spin = QSpinBox()
        self.expert_max_cluster_spin.setRange(2, 15)
        self.expert_max_cluster_spin.setSingleStep(1)
        self.expert_max_cluster_spin.setValue(self.ps2d_max_cluster_size)
        self.expert_max_cluster_spin.setToolTip("Maximum peaks in overlap cluster for simultaneous 2D fitting")
        ps2d_params_layout.addWidget(self.expert_max_cluster_spin, 4, 2)

        self.expert_ps2d_config_frame.setVisible(False)
        ps2d_config_layout.addWidget(self.expert_ps2d_config_frame)

        content_layout.addWidget(ps2d_config_group)

        # =================== PS2D FIT PARAMETERS SECTION ===================
        ps2d_fit_group = QGroupBox("PS2D Fit Parameters")
        ps2d_fit_layout = QVBoxLayout(ps2d_fit_group)

        # PS2D fit parameters toggle button
        self.expert_ps2d_fit_toggle_btn = QPushButton("PS2D fit parameters ▼")
        self.expert_ps2d_fit_toggle_btn.setToolTip("Show/hide PS2D fit parameters")
        self.expert_ps2d_fit_toggle_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {SECONDARY_BUTTON_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: 6px;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
        """)
        self.expert_ps2d_fit_toggle_btn.clicked.connect(self._toggle_ps2d_fit_params)
        ps2d_fit_layout.addWidget(self.expert_ps2d_fit_toggle_btn)

        # Collapsible PS2D fit parameters frame
        self.expert_ps2d_fit_frame = QWidget()
        fit_params_layout = QVBoxLayout(self.expert_ps2d_fit_frame)
        fit_params_layout.setSpacing(SPACING_SM)

        # Fix Linewidths checkbox
        fix_lw_row = QHBoxLayout()
        self.expert_fix_lw_check = QCheckBox("🔒 Fix Linewidths (fixLW flag)")
        self.expert_fix_lw_check.setChecked(self.fix_linewidths)
        self.expert_fix_lw_check.stateChanged.connect(self._on_fix_linewidths_change)
        fix_lw_row.addWidget(self.expert_fix_lw_check)
        fix_lw_hint = QLabel("(Hold linewidths constant during multi-peak fitting)")
        fix_lw_hint.setStyleSheet("font-size: 10px; color: gray;")
        fix_lw_row.addWidget(fix_lw_hint)
        fix_lw_row.addStretch()
        fit_params_layout.addLayout(fix_lw_row)

        # Fix Positions checkbox
        fix_pos_row = QHBoxLayout()
        self.expert_fix_pos_check = QCheckBox("🔒 Fix Positions (fixPos flag)")
        self.expert_fix_pos_check.setChecked(self.fix_positions)
        self.expert_fix_pos_check.stateChanged.connect(self._on_fix_positions_change)
        fix_pos_row.addWidget(self.expert_fix_pos_check)
        fix_pos_hint = QLabel("(Hold peak positions constant during multi-peak fitting)")
        fix_pos_hint.setStyleSheet("font-size: 10px; color: gray;")
        fix_pos_row.addWidget(fix_pos_hint)
        fix_pos_row.addStretch()
        fit_params_layout.addLayout(fix_pos_row)

        # Custom Initial Linewidths checkbox
        custom_lw_row = QHBoxLayout()
        self.expert_custom_lw_check = QCheckBox("Custom Initial Linewidths")
        self.expert_custom_lw_check.setChecked(self.use_custom_linewidths)
        self.expert_custom_lw_check.stateChanged.connect(self._on_custom_lw_toggle)
        custom_lw_row.addWidget(self.expert_custom_lw_check)
        custom_lw_hint = QLabel("(Override defaults: Lor/Gauss 1H=0.001/0.03, 15N=0.0001/0.3)")
        custom_lw_hint.setStyleSheet("font-size: 10px; color: gray;")
        custom_lw_row.addWidget(custom_lw_hint)
        custom_lw_row.addStretch()
        fit_params_layout.addLayout(custom_lw_row)

        # Custom linewidth input fields (collapsible)
        self.expert_custom_lw_frame = QWidget()
        custom_lw_layout = QGridLayout(self.expert_custom_lw_frame)
        custom_lw_layout.setSpacing(SPACING_SM)

        # 1H linewidths
        custom_lw_layout.addWidget(QLabel("1H:"), 0, 0)
        custom_lw_layout.addWidget(QLabel("Lor:"), 0, 1)
        self.expert_lw_lor_1h_spin = QDoubleSpinBox()
        self.expert_lw_lor_1h_spin.setRange(0.0001, 0.1)
        self.expert_lw_lor_1h_spin.setSingleStep(0.001)
        self.expert_lw_lor_1h_spin.setDecimals(4)
        self.expert_lw_lor_1h_spin.setValue(self.lw_lorentz_1h)
        self.expert_lw_lor_1h_spin.valueChanged.connect(self._on_custom_lw_change)
        custom_lw_layout.addWidget(self.expert_lw_lor_1h_spin, 0, 2)

        custom_lw_layout.addWidget(QLabel("Gauss:"), 0, 3)
        self.expert_lw_gauss_1h_spin = QDoubleSpinBox()
        self.expert_lw_gauss_1h_spin.setRange(0.001, 1.0)
        self.expert_lw_gauss_1h_spin.setSingleStep(0.01)
        self.expert_lw_gauss_1h_spin.setDecimals(3)
        self.expert_lw_gauss_1h_spin.setValue(self.lw_gauss_1h)
        self.expert_lw_gauss_1h_spin.valueChanged.connect(self._on_custom_lw_change)
        custom_lw_layout.addWidget(self.expert_lw_gauss_1h_spin, 0, 4)

        # 15N linewidths
        custom_lw_layout.addWidget(QLabel("15N:"), 1, 0)
        custom_lw_layout.addWidget(QLabel("Lor:"), 1, 1)
        self.expert_lw_lor_15n_spin = QDoubleSpinBox()
        self.expert_lw_lor_15n_spin.setRange(0.00001, 0.01)
        self.expert_lw_lor_15n_spin.setSingleStep(0.0001)
        self.expert_lw_lor_15n_spin.setDecimals(5)
        self.expert_lw_lor_15n_spin.setValue(self.lw_lorentz_15n)
        self.expert_lw_lor_15n_spin.valueChanged.connect(self._on_custom_lw_change)
        custom_lw_layout.addWidget(self.expert_lw_lor_15n_spin, 1, 2)

        custom_lw_layout.addWidget(QLabel("Gauss:"), 1, 3)
        self.expert_lw_gauss_15n_spin = QDoubleSpinBox()
        self.expert_lw_gauss_15n_spin.setRange(0.001, 0.5)
        self.expert_lw_gauss_15n_spin.setSingleStep(0.005)
        self.expert_lw_gauss_15n_spin.setDecimals(3)
        self.expert_lw_gauss_15n_spin.setValue(self.lw_gauss_15n)
        self.expert_lw_gauss_15n_spin.valueChanged.connect(self._on_custom_lw_change)
        custom_lw_layout.addWidget(self.expert_lw_gauss_15n_spin, 1, 4)

        self.expert_custom_lw_frame.setVisible(False)
        fit_params_layout.addWidget(self.expert_custom_lw_frame)

        # Parallel processing checkbox
        parallel_row = QHBoxLayout()
        self.expert_parallel_check = QCheckBox("🚀 Use Parallel Processing (75% cores)")
        self.expert_parallel_check.setChecked(self.use_parallel_processing)
        self.expert_parallel_check.stateChanged.connect(self._on_parallel_change)
        parallel_row.addWidget(self.expert_parallel_check)
        parallel_row.addStretch()
        fit_params_layout.addLayout(parallel_row)

        # Adaptive optimization checkbox (grid search for optimal radF1/radF2)
        adaptive_row = QHBoxLayout()
        self.expert_adaptive_check = QCheckBox("🎯 Adaptive Parameter Optimization")
        self.expert_adaptive_check.setChecked(self.use_adaptive_optimization)
        self.expert_adaptive_check.stateChanged.connect(self._on_adaptive_optimization_change)
        adaptive_row.addWidget(self.expert_adaptive_check)
        adaptive_hint = QLabel("(Grid search for optimal radF1/radF2 after PASS1)")
        adaptive_hint.setStyleSheet("font-size: 10px; color: gray;")
        adaptive_row.addWidget(adaptive_hint)
        adaptive_row.addStretch()
        fit_params_layout.addLayout(adaptive_row)

        self.expert_ps2d_fit_frame.setVisible(False)
        ps2d_fit_layout.addWidget(self.expert_ps2d_fit_frame)

        content_layout.addWidget(ps2d_fit_group)

        # =================== SERIES INTEGRATION SECTION ===================
        series_group = QGroupBox("🚀 Series Integration")
        series_layout = QVBoxLayout(series_group)

        # PS2D Linewidth Reuse checkbox
        lw_reuse_row = QHBoxLayout()
        self.expert_lw_reuse_check = QCheckBox("🔒 PS2D Linewidth Reuse (Fix LW from reference, ~40% speedup)")
        self.expert_lw_reuse_check.setChecked(self.use_ps2d_linewidth_reuse)
        self.expert_lw_reuse_check.stateChanged.connect(self._on_lw_reuse_toggle)
        lw_reuse_row.addWidget(self.expert_lw_reuse_check)
        lw_reuse_row.addStretch()
        series_layout.addLayout(lw_reuse_row)

        # Separator
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        series_layout.addWidget(separator)

        # Peak List Source label
        peak_source_label = QLabel("🎯 Peak List Source:")
        peak_source_label.setStyleSheet("font-weight: bold;")
        series_layout.addWidget(peak_source_label)

        # Peak source radio buttons
        self.expert_peak_source_group = QButtonGroup(dialog)
        detected_radio = QRadioButton("Use detected peaks")
        cascade_radio = QRadioButton("Cascade mode (propagate positions: n→n+1→n+2)")
        self.expert_peak_source_group.addButton(detected_radio, 0)
        self.expert_peak_source_group.addButton(cascade_radio, 1)

        if self.series_peak_source == "detected":
            detected_radio.setChecked(True)
        else:
            cascade_radio.setChecked(True)

        self.expert_peak_source_group.buttonClicked.connect(self._on_peak_source_change)
        series_layout.addWidget(detected_radio)
        series_layout.addWidget(cascade_radio)

        # Cascade drift limit checkbox (indented under cascade mode)
        drift_limit_row = QHBoxLayout()
        drift_limit_row.addSpacing(20)  # Indent to show it's related to cascade mode
        self.expert_cascade_drift_check = QCheckBox("🔒 Enforce max drift from reference")
        self.expert_cascade_drift_check.setChecked(self.enable_cascade_drift_limit)
        self.expert_cascade_drift_check.stateChanged.connect(self._on_cascade_drift_limit_change)
        self.expert_cascade_drift_check.setToolTip(
            "When enabled, peak positions are constrained to stay within max drift\n"
            "of their ORIGINAL reference positions (15N: ±0.10/±0.05, 13C: ±0.05/±0.05 ppm).\n"
            "This prevents unbounded drift accumulation across many spectra.\n"
            "Disable if you expect large genuine chemical shift changes."
        )
        drift_limit_row.addWidget(self.expert_cascade_drift_check)
        drift_limit_hint = QLabel("(15N: ±0.10/±0.05, 13C: ±0.05/±0.05 ppm)")
        drift_limit_hint.setStyleSheet("font-size: 10px; color: gray;")
        drift_limit_row.addWidget(drift_limit_hint)
        drift_limit_row.addStretch()
        series_layout.addLayout(drift_limit_row)

        content_layout.addWidget(series_group)

        # =================== ADDITIONAL EXPERT OPTIONS SECTION ===================
        future_group = QGroupBox("Additional Expert Options")
        future_layout = QVBoxLayout(future_group)
        future_label = QLabel("More expert options will be added here...")
        future_label.setStyleSheet("color: gray;")
        future_layout.addWidget(future_label)
        content_layout.addWidget(future_group)

        # Add stretch
        content_layout.addStretch()

        # Set content to scroll area
        scroll.setWidget(content)
        main_layout.addWidget(scroll)

        # Close button at bottom
        close_btn = QPushButton("Close")
        close_btn.setToolTip("Close expert mode dialog")
        close_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {SECONDARY_BUTTON_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: 8px 16px;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
        """)
        close_btn.clicked.connect(dialog.accept)
        main_layout.addWidget(close_btn, alignment=Qt.AlignCenter)

        # Show dialog
        dialog.exec()
        self.update_status("Expert Mode parameters updated")

    def _sync_expert_dialog_values(self):
        """Sync Expert Mode dialog widgets to current parameter values.

        Called when reopening a cached dialog to ensure UI matches internal state.
        """
        # Detection parameters
        if hasattr(self, 'expert_search_x_spin'):
            self.expert_search_x_spin.setValue(self.search_window_x)
        if hasattr(self, 'expert_search_y_spin'):
            self.expert_search_y_spin.setValue(self.search_window_y)
        if hasattr(self, 'expert_noise_spin'):
            self.expert_noise_spin.setValue(self.noise_threshold)

        # PS2D config parameters
        if hasattr(self, 'expert_rad_x_spin'):
            self.expert_rad_x_spin.setValue(self.ps2d_rad_x)
        if hasattr(self, 'expert_rad_y_spin'):
            self.expert_rad_y_spin.setValue(self.ps2d_rad_y)
        if hasattr(self, 'expert_overlap_x_spin'):
            self.expert_overlap_x_spin.setValue(self.ps2d_overlap_x)
        if hasattr(self, 'expert_overlap_y_spin'):
            self.expert_overlap_y_spin.setValue(self.ps2d_overlap_y)
        if hasattr(self, 'expert_max_iter_spin'):
            self.expert_max_iter_spin.setValue(self.ps2d_max_iterations)
        if hasattr(self, 'expert_max_cluster_spin'):
            self.expert_max_cluster_spin.setValue(self.ps2d_max_cluster_size)

        # PS2D fit parameters (checkboxes) - block signals during sync to avoid duplicate logs
        if hasattr(self, 'expert_fix_lw_check'):
            self.expert_fix_lw_check.blockSignals(True)
            self.expert_fix_lw_check.setChecked(self.fix_linewidths)
            self.expert_fix_lw_check.blockSignals(False)
        if hasattr(self, 'expert_fix_pos_check'):
            self.expert_fix_pos_check.blockSignals(True)
            self.expert_fix_pos_check.setChecked(self.fix_positions)
            self.expert_fix_pos_check.blockSignals(False)
        if hasattr(self, 'expert_parallel_check'):
            self.expert_parallel_check.blockSignals(True)
            self.expert_parallel_check.setChecked(self.use_parallel_processing)
            self.expert_parallel_check.blockSignals(False)
        if hasattr(self, 'expert_adaptive_check'):
            self.expert_adaptive_check.blockSignals(True)
            self.expert_adaptive_check.setChecked(self.use_adaptive_optimization)
            self.expert_adaptive_check.blockSignals(False)
        if hasattr(self, 'expert_custom_lw_check'):
            self.expert_custom_lw_check.blockSignals(True)
            self.expert_custom_lw_check.setChecked(self.use_custom_linewidths)
            self.expert_custom_lw_check.blockSignals(False)
        if hasattr(self, 'expert_lw_reuse_check'):
            self.expert_lw_reuse_check.blockSignals(True)
            self.expert_lw_reuse_check.setChecked(self.use_ps2d_linewidth_reuse)
            self.expert_lw_reuse_check.blockSignals(False)

        # Custom linewidth values
        if hasattr(self, 'expert_lw_lor_1h_spin'):
            self.expert_lw_lor_1h_spin.setValue(self.lw_lorentz_1h)
        if hasattr(self, 'expert_lw_gauss_1h_spin'):
            self.expert_lw_gauss_1h_spin.setValue(self.lw_gauss_1h)
        if hasattr(self, 'expert_lw_lor_15n_spin'):
            self.expert_lw_lor_15n_spin.setValue(self.lw_lorentz_15n)
        if hasattr(self, 'expert_lw_gauss_15n_spin'):
            self.expert_lw_gauss_15n_spin.setValue(self.lw_gauss_15n)

    # =================== EXPERT MODE HANDLER METHODS ===================

    def _toggle_expert_params(self):
        """Toggle visibility of detection parameters in Expert Mode dialog."""
        visible = self.expert_params_frame.isVisible()
        self.expert_params_frame.setVisible(not visible)
        self.expert_params_toggle_btn.setText("Parameters ▲" if not visible else "Parameters ▼")

    def _toggle_ps2d_config_params(self):
        """Toggle visibility of PS2D config parameters in Expert Mode dialog."""
        visible = self.expert_ps2d_config_frame.isVisible()
        self.expert_ps2d_config_frame.setVisible(not visible)
        self.expert_ps2d_config_toggle_btn.setText("Parameters ▲" if not visible else "Parameters ▼")

    def _toggle_ps2d_fit_params(self):
        """Toggle visibility of PS2D fit parameters in Expert Mode dialog."""
        visible = self.expert_ps2d_fit_frame.isVisible()
        self.expert_ps2d_fit_frame.setVisible(not visible)
        self.expert_ps2d_fit_toggle_btn.setText("PS2D fit parameters ▲" if not visible else "PS2D fit parameters ▼")

    def _on_expert_param_change(self):
        """Handle detection parameter changes from Expert Mode spinboxes."""
        self.search_window_x = self.expert_search_x_spin.value()
        self.search_window_y = self.expert_search_y_spin.value()
        self.noise_threshold = self.expert_noise_spin.value()
        self.centroid_window_x_ppm = self.expert_centroid_x_spin.value()
        self.centroid_window_y_ppm = self.expert_centroid_y_spin.value()

        # Update integrator if available
        if hasattr(self, 'integrator') and self.integrator:
            self.integrator.set_search_window(self.search_window_x, self.search_window_y)
            self.integrator.set_threshold_multiplier(self.noise_threshold)

        logger.info(f"Expert params updated: search_x={self.search_window_x}, search_y={self.search_window_y}, noise={self.noise_threshold}")

    def _on_auto_add_dummy_change(self, state):
        """Handle auto-add dummy peaks checkbox change."""
        self.auto_add_dummy_peaks = (state == Qt.Checked)
        logger.info(f"Auto-add dummy peaks: {self.auto_add_dummy_peaks}")

    def _on_expert_detection_size_change(self):
        """Handle detection rectangle size changes."""
        self.detection_square_size = self.expert_det_x_spin.value()
        self.detection_rectangle_y = self.expert_det_y_spin.value()
        self._update_detection_square_ppm()

    def _update_detection_square_ppm(self):
        """Update ppm conversion display for detection square size."""
        if hasattr(self, 'integrator') and self.integrator:
            if hasattr(self.integrator, 'ppm_x_axis') and hasattr(self.integrator, 'ppm_y_axis'):
                try:
                    x_axis = self.integrator.ppm_x_axis
                    y_axis = self.integrator.ppm_y_axis
                    x_ppm_per_pixel = abs(x_axis[-1] - x_axis[0]) / len(x_axis)
                    y_ppm_per_pixel = abs(y_axis[-1] - y_axis[0]) / len(y_axis)

                    standard_x_ppm = x_ppm_per_pixel * self.detection_square_size
                    standard_y_ppm = y_ppm_per_pixel * self.detection_rectangle_y

                    self.detection_square_ppm_x = f"(≈{standard_x_ppm:.3f}×{standard_y_ppm:.2f} ppm)"
                    if hasattr(self, 'expert_det_ppm_label'):
                        self.expert_det_ppm_label.setText(self.detection_square_ppm_x)
                    return
                except:
                    pass
        self.detection_square_ppm_x = "(load data to see ppm)"
        if hasattr(self, 'expert_det_ppm_label'):
            self.expert_det_ppm_label.setText(self.detection_square_ppm_x)

    def _on_nucleus_type_change(self, nucleus: str):
        """Handle nucleus type selection change for PS2D configuration."""
        self.nucleus_type = nucleus
        set_ps2d_config(nucleus)
        config = get_ps2d_config()

        # Update instance variables
        self.ps2d_radF1 = config.radF1
        self.ps2d_radF2 = config.radF2
        self.ps2d_max_iterations = config.max_iterations
        self.ps2d_max_cluster_size = config.max_cluster_size
        self.ps2d_overlap_x = config.overlap_threshold_x
        self.ps2d_overlap_y = config.overlap_threshold_y

        # Update spinboxes if they exist
        if hasattr(self, 'expert_radF1_spin'):
            self.expert_radF1_spin.setValue(config.radF1)
            self.expert_radF2_spin.setValue(config.radF2)
            self.expert_max_iter_spin.setValue(config.max_iterations)
            self.expert_overlap_x_spin.setValue(config.overlap_threshold_x)
            self.expert_overlap_y_spin.setValue(config.overlap_threshold_y)
        if hasattr(self, 'expert_max_cluster_spin'):
            self.expert_max_cluster_spin.setValue(config.max_cluster_size)

        logger.info(f"PS2D config switched to: {nucleus} (radF1={config.radF1:.3f}, radF2={config.radF2:.4f})")
        self.update_status(f"PS2D configuration: {nucleus}-HSQC (radF1={config.radF1:.2f}, radF2={config.radF2:.3f})")

        # Update plot if ellipses are shown (v0.9 main_gui.py:3228-3229)
        if self.show_ellipses:
            self.update_spectrum_plot()

    def _apply_auto_detected_nucleus(self):
        """Apply auto-detected nucleus type from integrator to GUI.

        This method checks if the integrator has auto-detected a nucleus type
        (15N or 13C) based on spectral dimensions and updates the Expert Mode
        PS2D configuration accordingly.

        Based on v0.9 main_gui.py lines 3346-3358.
        """
        if not hasattr(self.integrator, 'auto_detected_nucleus'):
            return

        detected_nucleus = self.integrator.auto_detected_nucleus
        if not detected_nucleus:
            logger.info("Could not auto-detect nucleus type from spectral dimensions")
            return

        # Check if nucleus type changed
        current_setting = self.nucleus_type
        if current_setting != detected_nucleus:
            # Update internal state
            self.nucleus_type = detected_nucleus

            # Update Expert Mode combo box if it exists
            if hasattr(self, 'expert_nucleus_combo'):
                self.expert_nucleus_combo.setCurrentText(detected_nucleus)

            # Trigger PS2D configuration update
            self._on_nucleus_type_change(detected_nucleus)

            logger.info(f"PS2D configuration automatically switched from {current_setting} to {detected_nucleus}")
            self.update_status(f"Auto-detected: {detected_nucleus}-HSQC spectrum")
        else:
            logger.info(f"PS2D configuration confirmed: {detected_nucleus} (already set)")

    def _on_ps2d_params_apply(self):
        """Apply PS2D parameter changes from Expert Mode dialog."""
        try:
            radF1 = self.expert_radF1_spin.value()
            radF2 = self.expert_radF2_spin.value()
            max_iterations = self.expert_max_iter_spin.value()
            max_cluster_size = self.expert_max_cluster_spin.value()
            overlap_x = self.expert_overlap_x_spin.value()
            overlap_y = self.expert_overlap_y_spin.value()

            # Validate
            if radF1 <= 0 or radF2 <= 0:
                self.update_status("PS2D radii must be positive")
                return
            if max_iterations < 50:
                self.update_status("PS2D max iterations must be >= 50")
                return
            if max_cluster_size < 2:
                self.update_status("PS2D max cluster size must be >= 2")
                return
            if overlap_x <= 0 or overlap_y <= 0:
                self.update_status("PS2D overlap thresholds must be positive")
                return

            # Update instance variables
            self.ps2d_radF1 = radF1
            self.ps2d_radF2 = radF2
            self.ps2d_max_iterations = max_iterations
            self.ps2d_max_cluster_size = max_cluster_size
            self.ps2d_overlap_x = overlap_x
            self.ps2d_overlap_y = overlap_y

            # Update global ps2d_config
            config = get_ps2d_config()
            config.radF1 = radF1
            config.radF2 = radF2
            config.max_iterations = max_iterations
            config.max_cluster_size = max_cluster_size
            config.overlap_threshold_x = overlap_x
            config.overlap_threshold_y = overlap_y
            config.radF1_selector = radF1 * 1.5
            config.radF2_selector = radF2 * 1.5

            logger.info(f"PS2D params applied: radF1={radF1:.4f}, radF2={radF2:.4f}, max_iter={max_iterations}, max_cluster={max_cluster_size}")
            self.update_status(f"PS2D params: radF1={radF1:.3f}, radF2={radF2:.4f}, cluster_size={max_cluster_size}")

            # Update plot with new PS2D parameters
            self.update_spectrum_plot()

        except Exception as e:
            logger.error(f"Error applying PS2D params: {e}")
            self.update_status("Error updating PS2D parameters")

    def _on_show_ellipses_change(self, state):
        """Handle show ellipses checkbox change.

        Updates the spectrum plot to show/hide PS2D elliptical windows.
        Based on v0.9 main_gui.py:3007.
        """
        self.show_ellipses = (state == 2)  # Qt.CheckState.Checked = 2
        self.update_spectrum_plot()

    def _on_adaptive_optimization_change(self, state):
        """Handle adaptive optimization checkbox change.

        When enabled (default), the system learns optimal radF1/radF2 and
        overlap thresholds from isolated peaks using grid search optimization.
        When disabled, fixed config values from PS2D config are used.

        Note: Future improvement could use Bayesian Optimization (e.g., scikit-optimize)
        for larger parameter spaces, but grid search is optimal for current 9 combinations.
        """
        self.use_adaptive_optimization = (state == 2)  # Qt.CheckState.Checked = 2
        if self.use_adaptive_optimization:
            logger.info("Adaptive optimization ENABLED - grid search for optimal radF1/radF2")
        else:
            logger.info("Adaptive optimization DISABLED - using default parameters")

    def _on_fix_linewidths_change(self, state):
        """Handle fix linewidths checkbox change."""
        self.fix_linewidths = (state == 2)  # Qt.CheckState.Checked = 2
        logger.info(f"Fix linewidths: {self.fix_linewidths}")

    def _on_fix_positions_change(self, state):
        """Handle fix positions checkbox change."""
        self.fix_positions = (state == 2)  # Qt.CheckState.Checked = 2
        logger.info(f"Fix positions: {self.fix_positions}")

    def _on_custom_lw_toggle(self, state):
        """Toggle visibility of custom linewidth inputs."""
        self.use_custom_linewidths = (state == 2)  # Qt.CheckState.Checked = 2
        if hasattr(self, 'expert_custom_lw_frame'):
            self.expert_custom_lw_frame.setVisible(self.use_custom_linewidths)
        logger.info(f"Custom linewidths: {self.use_custom_linewidths}")

    def _on_custom_lw_change(self):
        """Handle custom linewidth value changes."""
        self.lw_lorentz_1h = self.expert_lw_lor_1h_spin.value()
        self.lw_gauss_1h = self.expert_lw_gauss_1h_spin.value()
        self.lw_lorentz_15n = self.expert_lw_lor_15n_spin.value()
        self.lw_gauss_15n = self.expert_lw_gauss_15n_spin.value()
        logger.info(f"Custom LW: 1H Lor={self.lw_lorentz_1h:.4f}, Gauss={self.lw_gauss_1h:.3f}")

    def _on_parallel_change(self, state):
        """Handle parallel processing checkbox change."""
        self.use_parallel_processing = (state == 2)  # Qt.CheckState.Checked = 2
        logger.info(f"Parallel processing: {self.use_parallel_processing}")

    def _on_lw_reuse_toggle(self, state):
        """Handle PS2D linewidth reuse toggle."""
        self.use_ps2d_linewidth_reuse = (state == 2)  # Qt.CheckState.Checked = 2
        if self.use_ps2d_linewidth_reuse:
            logger.info("PS2D Linewidth Reuse ENABLED - ~40% speedup for series")
        else:
            logger.info("PS2D Linewidth Reuse DISABLED")

    def _on_peak_source_change(self, button):
        """Handle peak source radio button change."""
        if self.expert_peak_source_group.id(button) == 0:
            self.series_peak_source = "detected"
        else:
            self.series_peak_source = "cascade"
        logger.info(f"Series peak source: {self.series_peak_source}")

    def _on_cascade_drift_limit_change(self, state):
        """Handle cascade drift limit toggle."""
        self.enable_cascade_drift_limit = (state == 2)  # Qt.CheckState.Checked = 2
        if self.enable_cascade_drift_limit:
            logger.info("Cascade drift limit ENABLED - positions bounded to nucleus-specific max drift from reference")
        else:
            logger.info("Cascade drift limit DISABLED - positions can drift freely in cascade mode")

    def on_reset_results(self):
        """Handle Reset Results button click.

        Resets the application to initial startup state:
        - No peaks loaded
        - Synthetic test spectrum loaded
        - All fitting results cleared
        """
        logger.info("Reset Results button clicked")

        # Confirm destructive action
        reply = QMessageBox.question(
            self,
            "Reset to Initial State",
            "This will reset the application to startup state:\n\n"
            "• Clear all peaks (reference and detected)\n"
            "• Clear all fitting results\n"
            "• Load synthetic test spectrum\n\n"
            "This action cannot be undone.",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )

        if reply == QMessageBox.Yes:
            logger.info("Resetting to initial state")

            # Reset the integrator to fresh instance
            from lunaNMR.core.enhanced_voigt_integrator import EnhancedVoigtIntegrator
            self.integrator = EnhancedVoigtIntegrator()

            # Clear fitting results
            self.current_voigt_result = None

            # Clear peak navigator (load empty lists)
            if hasattr(self, 'peak_navigator'):
                self.peak_navigator.load_reference_peaks([])
                self.peak_navigator.load_detected_peaks([])

            # Reset voigt plotters to placeholder state
            if hasattr(self, 'voigt_2d_plotter') and self.voigt_2d_plotter:
                self.voigt_2d_plotter.show_placeholder("Select a peak to view Voigt analysis")

            if hasattr(self, 'voigt_3d_plotter') and self.voigt_3d_plotter:
                self.voigt_3d_plotter.show_placeholder("Select a peak to view 3D Voigt analysis")

            # Load synthetic test spectrum
            self.create_test_spectrum()

            # Reset button colors to initial state
            self._reset_button_colors()

            logger.info("Application reset to initial state")
            self.update_status("Reset complete - synthetic spectrum loaded")

    def create_center_panel(self) -> QWidget:
        """Create the center panel with tabbed plot display.

        Returns:
            QTabWidget with Main Spectrum tab populated with SpectrumPlotter
            and placeholder tabs for other plot types.
        """
        # Import SpectrumPlotter
        from lunaNMR.gui.components.spectrum_plotter import SpectrumPlotter

        # Create tab widget
        tab_widget = QTabWidget()
        tab_widget.setStyleSheet(f"""
            QTabWidget::pane {{
                background-color: {FRAME_BG_COLOR};
                border: none;
            }}
            QTabBar::tab {{
                background-color: {PANEL_BG_COLOR};
                color: {PRIMARY_TEXT};
                padding: {SPACING_SM}px {SPACING_MD}px;
                margin-right: 2px;
            }}
            QTabBar::tab:selected {{
                background-color: {FRAME_BG_COLOR};
                font-weight: bold;
            }}
        """)

        # Create Main Spectrum tab with SpectrumPlotter
        main_spectrum_tab = QWidget()
        main_spectrum_layout = QVBoxLayout(main_spectrum_tab)
        main_spectrum_layout.setContentsMargins(0, 0, 0, 0)  # Zero margins for full space
        main_spectrum_layout.setSpacing(0)

        # Create SpectrumPlotter (no toolbar - uses custom navigation)
        self.spectrum_plotter = SpectrumPlotter(parent=main_spectrum_tab)
        main_spectrum_layout.addWidget(self.spectrum_plotter)

        # Connect peak edit signal (Shift+click = add, Ctrl+click = delete)
        self.spectrum_plotter.peak_edit_requested.connect(self._on_peak_edit_with_modifier)

        # Add Main Spectrum tab first
        tab_widget.addTab(main_spectrum_tab, "Main Spectrum")

        # Create Voigt 2D Analysis tab with VoigtAnalysisPlotter
        voigt_2d_tab = self._create_voigt_2d_tab()
        tab_widget.addTab(voigt_2d_tab, "📈 Voigt Analysis")

        # Create 3D Voigt Analysis tab (placeholder for Phase 4)
        voigt_3d_tab = self._create_voigt_3d_tab()
        tab_widget.addTab(voigt_3d_tab, "3D Voigt Analysis")

        # Create Peak Parameters tab
        peak_params_tab = self._create_peak_parameters_tab()
        tab_widget.addTab(peak_params_tab, "Peak Parameters")

        # Store reference for context-aware peak selection
        self.viz_tabs = tab_widget

        return tab_widget

    def create_right_panel(self) -> QWidget:
        """Create the right panel with Peak Navigator.

        Returns:
            PeakNavigator widget from Phase 1 implementation.
        """
        # Create Peak Navigator and save reference
        self.peak_navigator = PeakNavigator()

        # Load sample data for testing (simple list format)
        # Format: [assignment, x_coord, y_coord, height, peak_id]
        sample_peaks = [
            ["A1", 8.456, 120.234, 1.2e6, 1],
            ["A2", 8.123, 119.876, 9.5e5, 2],
            ["A3", 7.890, 118.543, 1.5e6, 3],
            ["A4", 7.654, 117.210, 7.8e5, 4],
            ["A5", 7.321, 115.877, 1.1e6, 5],
        ]
        self.peak_navigator.load_detected_peaks(sample_peaks)

        # Connect signals
        self.peak_navigator.peak_selected.connect(self._on_peak_selected)
        self.peak_navigator.navigation_requested.connect(self._on_navigation_requested)
        self.peak_navigator.peak_analysis_requested.connect(self._on_peak_analysis_requested)

        return self.peak_navigator

    def _create_voigt_2d_tab(self) -> QWidget:
        """Create the Voigt 2D Analysis tab with VoigtAnalysisPlotter.

        Returns:
            QWidget containing the VoigtAnalysisPlotter (2x2 grid layout)

        Features:
            - 2x2 grid showing cross-sections (F1/F2) and 2D/residuals
            - Layer toggling (experimental/fitted/residuals/peaks)
            - Interactive toolbar for navigation and zoom
            - Ready for data updates from fitting results

        Note:
            VoigtAnalysisPlotter is initialized but will remain empty until
            fitting results are loaded. Data is populated from core_integrator
            in Phase 4.
        """
        from lunaNMR.gui.components.voigt_analysis_plotter import VoigtAnalysisPlotter

        # Create tab widget
        tab_widget = QWidget()
        layout = QVBoxLayout(tab_widget)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        # Create VoigtAnalysisPlotter with toolbar
        self.voigt_2d_plotter = VoigtAnalysisPlotter(parent=tab_widget, toolbar=True)
        layout.addWidget(self.voigt_2d_plotter)

        logger.debug("Voigt 2D Analysis tab created with VoigtAnalysisPlotter")
        return tab_widget

    def _create_voigt_3d_tab(self) -> QWidget:
        """Create the 3D Voigt Analysis tab with VoigtAnalysisPlotter.

        Returns:
            QWidget containing the VoigtAnalysisPlotter (3D surface visualization)

        Features:
            - 3D surface plot of experimental and fitted Voigt profiles
            - Interactive rotation and zoom
            - Layer toggling (experimental/fitted/residuals/peaks)
            - Toolbar for navigation
            - Ready for data updates from fitting results

        Note:
            VoigtAnalysisPlotter is initialized FIRST before controls to ensure
            signal connections work properly. Data is populated after fitting completes.

        Based on v0o9 implementation (lines 1829-1978).
        """
        from lunaNMR.gui.components.voigt_analysis_plotter import VoigtAnalysisPlotter
        from PySide6.QtWidgets import (QGroupBox, QCheckBox, QRadioButton, QComboBox,
                                       QSlider, QLabel, QButtonGroup)
        from PySide6.QtCore import Qt

        # Create tab widget
        tab_widget = QWidget()
        main_layout = QVBoxLayout(tab_widget)
        main_layout.setContentsMargins(5, 5, 5, 5)
        main_layout.setSpacing(3)

        # ===== CREATE PLOTTER FIRST (before controls that reference it) =====
        self.voigt_3d_plotter = VoigtAnalysisPlotter(parent=tab_widget, toolbar=True)

        # ===== Control Bar (2 rows) - matching v0o9 lines 1851-1933 =====
        control_frame = QWidget()
        control_layout = QVBoxLayout(control_frame)
        control_layout.setContentsMargins(0, 0, 0, 0)
        control_layout.setSpacing(2)

        # === Row 1: Layer visibility checkboxes ===
        row1 = QWidget()
        row1_layout = QHBoxLayout(row1)
        row1_layout.setContentsMargins(0, 0, 0, 0)
        row1_layout.setSpacing(5)

        # Layer visibility group (bold text, 25% larger font)
        layer_group = QGroupBox("Layer Visibility")
        layer_group.setStyleSheet("""
            QGroupBox { font-weight: bold; font-size: 15px; }
            QCheckBox { font-weight: bold; font-size: 15px; }
        """)
        layer_layout = QHBoxLayout(layer_group)
        layer_layout.setContentsMargins(8, 3, 8, 3)
        layer_layout.setSpacing(15)

        self.show_exp_3d_checkbox = QCheckBox("Exp")
        self.show_exp_3d_checkbox.setChecked(self.show_exp_3d)
        self.show_exp_3d_checkbox.toggled.connect(self._on_show_exp_3d_change)
        layer_layout.addWidget(self.show_exp_3d_checkbox)

        self.show_fit_3d_checkbox = QCheckBox("Fit")
        self.show_fit_3d_checkbox.setChecked(self.show_fit_3d)
        self.show_fit_3d_checkbox.toggled.connect(self._on_show_fit_3d_change)
        layer_layout.addWidget(self.show_fit_3d_checkbox)

        self.show_individual_3d_checkbox = QCheckBox("Individual")
        self.show_individual_3d_checkbox.setChecked(self.show_individual_3d)
        self.show_individual_3d_checkbox.toggled.connect(self._on_show_individual_3d_change)
        layer_layout.addWidget(self.show_individual_3d_checkbox)

        self.show_resid_3d_checkbox = QCheckBox("Residuals")
        self.show_resid_3d_checkbox.setChecked(self.show_resid_3d)
        self.show_resid_3d_checkbox.toggled.connect(self._on_show_resid_3d_change)
        layer_layout.addWidget(self.show_resid_3d_checkbox)

        self.limit_peak_display_3d_checkbox = QCheckBox("Limit Peak")
        self.limit_peak_display_3d_checkbox.setChecked(self.limit_peak_display_3d)
        self.limit_peak_display_3d_checkbox.toggled.connect(self._on_limit_peak_display_3d_change)
        layer_layout.addWidget(self.limit_peak_display_3d_checkbox)

        # Keep peak labels state but don't show control (default off)
        self.show_peak_labels_3d = False

        row1_layout.addWidget(layer_group)

        # Color scheme group (next to Layer Visibility)
        color_scheme_group = QGroupBox("Color Scheme")
        color_scheme_layout = QHBoxLayout(color_scheme_group)
        color_scheme_layout.setContentsMargins(3, 3, 3, 3)
        color_scheme_layout.setSpacing(2)

        self.color_scheme_3d_dropdown = QComboBox()
        self.color_scheme_3d_dropdown.addItems(['Classic', 'Clean', 'Dark', 'Warm'])
        self.color_scheme_3d_dropdown.setCurrentText(self.color_scheme_3d)
        self.color_scheme_3d_dropdown.currentTextChanged.connect(self._on_color_scheme_change_3d)
        # Make dropdown larger with larger text
        self.color_scheme_3d_dropdown.setMinimumWidth(120)
        self.color_scheme_3d_dropdown.setStyleSheet("QComboBox { font-size: 13px; padding: 4px 8px; }")
        color_scheme_layout.addWidget(self.color_scheme_3d_dropdown)

        row1_layout.addWidget(color_scheme_group)
        row1_layout.addStretch()
        control_layout.addWidget(row1)

        # Keep residual mode and intensity scale state variables functional (hidden controls)
        # Default residual mode to 'separate'
        self.residual_mode_3d = 'separate'
        # Intensity scale slider/label not shown but state maintained via self.intensity_scale_3d

        # Add control frame first (at top)
        main_layout.addWidget(control_frame)

        # Add plotter below controls
        main_layout.addWidget(self.voigt_3d_plotter)

        logger.debug("3D Voigt Analysis tab created with control bar and VoigtAnalysisPlotter")
        return tab_widget

    def _create_peak_parameters_tab(self) -> QWidget:
        """Create the Peak Parameters tab with QTextEdit for parameter display.

        Returns:
            QWidget containing a read-only QTextEdit showing peak fit results

        Features:
            - Monospace font (Courier, 9pt) for aligned parameter display
            - Read-only text widget for parameter viewing
            - Placeholder text guides users
            - Auto-populated when fitting completes (Phase 4)

        Note:
            The QTextEdit is accessible via self.peak_params_text for updates
            from the fitting engine.
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
            "Peak parameters will appear here after fitting completes.\n\n"
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

        Displays detailed parameters for the selected peak/group of peaks
        including position, intensity, linewidths, and quality metrics.

        Based on v0.9 main_gui.py lines 6400-6488.

        Args:
            voigt_result: Dictionary containing Voigt fitting results with keys:
                - assignment: Peak assignment
                - fitting_quality: Quality rating (Excellent/Good/Fair/Poor)
                - avg_r_squared or r_squared: R² value
                - method: Fitting method used
                - all_peaks: List of individual peak parameters
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

    def _on_peak_selected(self, peak_id: int):
        """Handle peak selection from navigator with context-aware behavior.

        Behavior depends on currently active tab (v0.9 set_selected_peak lines 7826-7863):
        - Main Spectrum tab: Zoom to peak
        - 3D Voigt Analysis tab: Show Voigt analysis for peak

        Args:
            peak_id: Row index of the selected peak in the navigator table
        """
        import logging
        logger = logging.getLogger(__name__)

        try:
            # Check current active tab for context-aware behavior (v0.9 line 7834)
            current_tab_name = ""
            if hasattr(self, 'viz_tabs'):
                current_tab_name = self.viz_tabs.tabText(self.viz_tabs.currentIndex())

            # Get current peak type from navigator (reference or detected)
            peak_type = "reference"
            if hasattr(self, 'peak_navigator'):
                peak_type = self.peak_navigator.selected_peak_type

            # Get peak coordinates based on peak type
            peak_x, peak_y, assignment = None, None, f'Peak {peak_id}'

            if peak_type == "detected":
                # Use detected peaks list from navigator
                if hasattr(self.peak_navigator, 'detected_peaks') and peak_id < len(self.peak_navigator.detected_peaks):
                    peak_data = self.peak_navigator.detected_peaks[peak_id]
                    assignment = str(peak_data[0])
                    peak_x = float(peak_data[1])
                    peak_y = float(peak_data[2])
            else:
                # Use reference peaks from integrator peak_list
                if hasattr(self.integrator, 'peak_list') and self.integrator.peak_list is not None:
                    if peak_id < len(self.integrator.peak_list):
                        peak_row = self.integrator.peak_list.iloc[peak_id]
                        peak_x = peak_row.get('Position_X') or peak_row.get('ppm_x') or peak_row.get('x_ppm')
                        peak_y = peak_row.get('Position_Y') or peak_row.get('ppm_y') or peak_row.get('y_ppm')
                        assignment = peak_row.get('Assignment', f'Peak {peak_id}')

            if peak_x is not None and peak_y is not None:
                # Always zoom to the selected peak (v0o9 lines 7806-7808)
                self.spectrum_plotter.set_zoom(
                    peak_x, peak_y,
                    0.3,   # Tight X zoom (1H dimension)
                    2.5    # Tight Y zoom (15N dimension)
                )

                # Context-aware: If Voigt tabs are active, show analysis (v0.9 line 7860-7863)
                if current_tab_name in ("3D Voigt Analysis", "📈 Voigt Analysis"):
                    if peak_type == "detected":
                        self.show_detected_peak_analysis_by_index(peak_id, auto_switch_tab=False)
                    else:
                        self.show_reference_peak_analysis_by_index(peak_id, auto_switch_tab=False)
                else:
                    # Not on Voigt tab - just show selection status
                    self.update_status(f"Selected: {assignment} at ({peak_x:.3f}, {peak_y:.2f}) ppm")
            else:
                logger.warning(f"Could not extract coordinates for peak {peak_id}")
                self.update_status(f"Peak {peak_id} selected")

        except Exception as e:
            logger.error(f"Error highlighting peak {peak_id}: {e}")
            self.update_status(f"Peak {peak_id} selected")

    def _on_navigation_requested(self, direction: str):
        """Handle navigation request from navigator.

        Args:
            direction: Navigation direction ('prev' or 'next')
        """
        self.update_status(f"Navigate {direction} requested")

    @Slot(str, int)
    def _on_peak_analysis_requested(self, peak_type: str, peak_index: int):
        """Handle peak analysis request from navigator (exact v0o9 logic).

        This is called when user clicks the analysis button (🔬) in peak navigator.

        Args:
            peak_type: Type of peak ('reference' or 'detected')
            peak_index: Index of peak in the respective list

        Based on v0.9 main_gui.py:6568-6589 (navigator_show_peak_analysis)
        """
        if peak_type == "reference":
            # For reference peaks - show Voigt analysis (v0.9 main_gui.py:6570-6582)
            self.show_reference_peak_analysis_by_index(peak_index)

        elif peak_type == "detected":
            # For detected peaks - show Voigt analysis directly (v0o9 line 6586)
            self.show_detected_peak_analysis_by_index(peak_index)

        else:
            QMessageBox.critical(
                self,
                "Error",
                f"Unknown peak type: {peak_type}"
            )

    def show_detected_peak_analysis_by_index(self, peak_index: int, auto_switch_tab: bool = True):
        """Show Voigt analysis for detected peak by index (exact v0o9 logic from line 6591).

        Args:
            peak_index: Index of peak in detected peaks list
            auto_switch_tab: If True, automatically switch to 3D Voigt tab.
                           If False, stay on current tab (used when clicking peak in navigator)
        """
        # Check if we have fitted results (same check as navigation button)
        if not hasattr(self.integrator, 'fitted_peaks') or not self.integrator.fitted_peaks:
            QMessageBox.information(
                self,
                "No Results",
                "No fitting results available for detected peaks.\n\n"
                "Run peak detection and fitting first."
            )
            return

        # Search by assignment instead of using direct index (v0o9 line 6624-6640)
        # Navigator and fitted_peaks may have different orders (especially with failed peaks)
        selected_result = None

        if hasattr(self, 'peak_navigator') and hasattr(self.peak_navigator, 'detected_peaks'):
            if peak_index < len(self.peak_navigator.detected_peaks):
                nav_assignment = str(self.peak_navigator.detected_peaks[peak_index][0])

                for result in self.integrator.fitted_peaks:
                    if str(result.get('assignment', '')) == nav_assignment:
                        selected_result = result
                        break

        # Display logic (v0o9 line 6652-6666)
        if selected_result:
            # Switch to 3D Voigt analysis tab (only if auto_switch_tab=True)
            if auto_switch_tab and hasattr(self, 'viz_tabs'):
                # Find and activate the "3D Voigt Analysis" tab
                for i in range(self.viz_tabs.count()):
                    if self.viz_tabs.tabText(i) == "3D Voigt Analysis":
                        self.viz_tabs.setCurrentIndex(i)
                        break

            # Show the analysis for the selected peak
            if hasattr(self, 'voigt_2d_plotter') and self.voigt_2d_plotter:
                self.voigt_2d_plotter.plot_voigt_analysis(selected_result)

            # Also update 3D view
            if hasattr(self, 'voigt_3d_plotter') and self.voigt_3d_plotter:
                self.voigt_3d_plotter.plot_voigt_analysis_3d(selected_result)

            # Update peak parameters tab (v0.9 line 6546)
            self.update_peak_parameters(selected_result)

            # Update status
            assignment = selected_result.get('assignment', f'Det_{peak_index+1}')
            quality = selected_result.get('fitting_quality', 'Unknown')
            self.update_status(f"Showing Voigt analysis for {assignment} (Quality: {quality})")

        else:
            # Peak not fitted yet - offer to fit it (same as navigation button, v0o9 line 6673-6685)
            if peak_index < len(self.peak_navigator.detected_peaks):
                assignment = str(self.peak_navigator.detected_peaks[peak_index][0])
            else:
                assignment = f'Det_{peak_index+1}'

            response = QMessageBox.question(
                self,
                "Peak Not Fitted",
                f"Peak {assignment} has not been fitted yet or fitting failed.\n\n"
                "Would you like to fit all peaks now?",
                QMessageBox.Yes | QMessageBox.No
            )

            if response == QMessageBox.Yes:
                self.on_fit_spectrum()

    def show_reference_peak_analysis_by_index(self, peak_index: int, auto_switch_tab: bool = True):
        """Show Voigt analysis for reference peak by index.

        This is the Qt port of v0.9's show_selected_peak_analysis (main_gui.py:6490-6567).

        Args:
            peak_index: Index of peak in reference peaks list (0-based)
            auto_switch_tab: If True, automatically switch to 3D Voigt tab.
        """
        # Check if we have a peak list (v0.9 line 6498-6502)
        if (not hasattr(self.integrator, 'peak_list') or
            self.integrator.peak_list is None or
            (hasattr(self.integrator.peak_list, 'empty') and self.integrator.peak_list.empty)):
            QMessageBox.information(
                self,
                "No Peak List",
                "No peak list loaded.\nPlease load peak data first."
            )
            return

        # Convert 0-based index to 1-based peak number (v0.9 line 6504-6508)
        peak_number = peak_index + 1
        max_peaks = len(self.integrator.peak_list)

        if peak_number < 1 or peak_number > max_peaks:
            QMessageBox.warning(
                self,
                "Invalid Peak",
                f"Peak number {peak_number} is invalid.\nValid range: 1-{max_peaks}"
            )
            return

        # Check if the selected peak has fitting results (v0.9 line 6511-6531)
        selected_result = None

        # Check fitted_peaks
        if (hasattr(self.integrator, 'fitted_peaks') and
            self.integrator.fitted_peaks is not None and
            len(self.integrator.fitted_peaks) > 0):
            for fitted_result in self.integrator.fitted_peaks:
                if fitted_result.get('peak_number') == peak_number:
                    selected_result = fitted_result
                    break

        # Fallback to voigt_fits (legacy location) for backward compatibility
        if selected_result is None and (hasattr(self.integrator, 'voigt_fits') and
            self.integrator.voigt_fits is not None and
            len(self.integrator.voigt_fits) > 0):
            for voigt_result in self.integrator.voigt_fits:
                if voigt_result.get('peak_number') == peak_number:
                    selected_result = voigt_result
                    break

        if selected_result:
            # Switch to 3D Voigt analysis tab (v0.9 line 6533-6535)
            if auto_switch_tab and hasattr(self, 'viz_tabs'):
                for i in range(self.viz_tabs.count()):
                    if self.viz_tabs.tabText(i) == "3D Voigt Analysis":
                        self.viz_tabs.setCurrentIndex(i)
                        break

            # Show the analysis for the selected peak (v0.9 line 6537-6546)
            if hasattr(self, 'voigt_2d_plotter') and self.voigt_2d_plotter:
                self.voigt_2d_plotter.plot_voigt_analysis(selected_result)

            # Also update 3D view
            if hasattr(self, 'voigt_3d_plotter') and self.voigt_3d_plotter:
                self.voigt_3d_plotter.plot_voigt_analysis_3d(selected_result)

            # Update peak parameters tab
            self.update_peak_parameters(selected_result)

            # Update status (v0.9 line 6548-6551)
            assignment = selected_result.get('assignment', f'Peak_{peak_number}')
            quality = selected_result.get('fitting_quality', 'Unknown')
            self.update_status(f"Showing Voigt analysis for {assignment} (Quality: {quality})")

        else:
            # Peak not fitted yet - offer to fit it (v0.9 line 6553-6566)
            peak_data = self.integrator.peak_list.iloc[peak_number - 1]
            assignment = peak_data.get('Assignment', f'Peak_{peak_number}')

            response = QMessageBox.question(
                self,
                "Peak Not Fitted",
                f"Peak {peak_number} ({assignment}) has not been fitted yet.\n\n"
                "Would you like to fit it now?",
                QMessageBox.Yes | QMessageBox.No
            )

            if response == QMessageBox.Yes:
                # Fit all peaks and then show analysis
                self.on_fit_spectrum()
                # After fitting completes, try to show analysis again
                # (use QTimer to delay the call until after processing)
                from PySide6.QtCore import QTimer
                QTimer.singleShot(100, lambda: self.show_reference_peak_analysis_by_index(peak_index, auto_switch_tab))

    def on_load_data(self):
        """Handle Load Data button click.

        Opens DataLoadingDialog for file selection, then loads the selected files
        into the integrator and updates the display.

        Process:
        1. Open DataLoadingDialog for file selection
        2. Load NMR spectrum file
        3. Load peak list file (if provided)
        4. Update file manager and integrator
        5. Refresh displays and navigator
        6. Show success/error messages
        """
        logger.info("Load Data button clicked")

        # Import dialog
        from lunaNMR.gui.dialogs import DataLoadingDialog

        # Create and show dialog with current workflow mode and folders
        dialog = DataLoadingDialog(
            parent=self,
            workflow_mode=self.workflow_mode,
            current_nmr_folder=self.current_nmr_folder,
            current_peak_folder=self.current_peak_folder
        )

        # Execute dialog (modal)
        if not dialog.exec():
            # User clicked Cancel
            logger.info("Data loading cancelled")
            self.update_status("Data loading cancelled")
            return

        # User clicked Load - get results
        result = dialog.get_result()
        nmr_file = result['nmr_file']
        peak_file = result['peak_file']

        # Update current folders for next time
        self.current_nmr_folder = result['nmr_folder']
        self.current_peak_folder = result['peak_folder']

        logger.info(f"Files selected - NMR: {nmr_file}, Peak: {peak_file}")

        # Load files into integrator
        try:
            self.load_nmr_file(nmr_file)

            if peak_file:
                self.load_peak_file(peak_file)
                logger.info(f"Successfully loaded NMR file: {nmr_file}")
                logger.info(f"Successfully loaded peak file: {peak_file}")
                self.update_status("Data loaded successfully")
            else:
                logger.info(f"Successfully loaded NMR file: {nmr_file} (no peak list)")
                self.update_status("NMR spectrum loaded (no peak list)")

        except Exception as e:
            from lunaNMR.gui.components import show_error
            show_error(self, "Load Error", f"Failed to load data files: {e}")
            logger.error(f"Failed to load data: {e}")

    def load_nmr_file(self, file_path: str):
        """Load NMR spectrum file into the integrator.

        This method:
        1. Loads NMR data using the integrator's load_nmr_file method
        2. Updates the current file path
        3. Updates status labels
        4. Updates the spectrum plot
        5. Shows success status in status bar

        Args:
            file_path: Path to NMR spectrum file (.ft, .ft2, .pipe, .fid)

        Raises:
            Exception: If file loading fails
        """
        logger.info(f"Loading NMR file: {file_path}")

        try:
            # Use integrator's load_nmr_file method
            success = self.integrator.load_nmr_file(file_path)

            if not success:
                raise RuntimeError("Integrator returned False - file loading failed")

            # Update current file path
            self.current_nmr_file = file_path

            # Update status labels
            self.update_file_status_labels(nmr_file=file_path, peak_file=self.current_peak_file)

            # Update spectrum plot
            self.update_spectrum_plot()

            # Update status bar
            filename = Path(file_path).name
            data_shape = self.integrator.nmr_data.shape if hasattr(self.integrator, 'nmr_data') else 'unknown'
            self.update_status(f"NMR spectrum loaded: {filename} (shape: {data_shape})")

            # Auto-detect and apply nucleus type based on spectral dimensions
            self._apply_auto_detected_nucleus()

            logger.info(f"NMR file loaded successfully: {filename}")

        except Exception as e:
            logger.error(f"Failed to load NMR file: {e}")
            raise RuntimeError(f"Failed to load NMR file: {e}")

    def load_peak_file(self, file_path: str):
        """Load peak list file into the integrator.

        This method:
        1. Loads peak list using the integrator's load_peak_list_file method
        2. Updates the current file path
        3. Updates status labels
        4. Updates peak navigator with loaded peaks
        5. Updates the spectrum plot to show peaks
        6. Shows success status in status bar

        Args:
            file_path: Path to peak list file (.txt, .csv, .peaks)

        Raises:
            Exception: If file loading fails
        """
        logger.info(f"Loading peak file: {file_path}")

        try:
            # Use integrator's load_peak_list_file method
            success = self.integrator.load_peak_list_file(file_path)

            if not success:
                raise RuntimeError("Integrator returned False - peak list loading failed")

            # Update current file path
            self.current_peak_file = file_path

            # Update status labels
            self.update_file_status_labels(nmr_file=self.current_nmr_file, peak_file=file_path)

            # Get peak data for navigator
            if hasattr(self.integrator, 'peak_list') and self.integrator.peak_list is not None:
                # Convert DataFrame to list format for navigator
                # Format: [assignment, x_coord, y_coord, height, peak_id]
                peak_data = []
                for idx, row in self.integrator.peak_list.iterrows():
                    assignment = row.get('Assignment', f'Peak_{idx+1}')
                    x_coord = row.get('Position_X', 0.0)
                    y_coord = row.get('Position_Y', 0.0)
                    height = row.get('Height', row.get('Intensity', 0.0))
                    peak_id = idx + 1

                    peak_data.append([assignment, x_coord, y_coord, height, peak_id])

                # Update peak navigator (if initialized)
                if self.peak_navigator is not None:
                    self.peak_navigator.load_reference_peaks(peak_data)
                else:
                    logger.warning("Peak navigator not initialized yet - peaks loaded but not displayed")

                peak_count = len(peak_data)
            else:
                peak_count = 0

            # Update plot with peaks
            self.update_spectrum_plot()

            # Update status bar
            filename = Path(file_path).name
            self.update_status(f"Peak list loaded: {filename} ({peak_count} peaks)")

            logger.info(f"Peak file loaded successfully: {filename} with {peak_count} peaks")

        except Exception as e:
            logger.error(f"Failed to load peak file: {e}")
            raise RuntimeError(f"Failed to load peak file: {e}")

    def save_current_zoom(self):
        """Save the current zoom limits for later restoration (v0o9 line 7775)."""
        if hasattr(self, 'spectrum_plotter') and self.spectrum_plotter:
            ax = self.spectrum_plotter.ax
            if ax:
                self.saved_xlim = ax.get_xlim()
                self.saved_ylim = ax.get_ylim()
                logger.debug(f"Saved zoom: xlim={self.saved_xlim}, ylim={self.saved_ylim}")

    def restore_zoom(self):
        """Restore the previously saved zoom limits (v0o9 line 7781)."""
        if self.saved_xlim and self.saved_ylim:
            if hasattr(self, 'spectrum_plotter') and self.spectrum_plotter:
                ax = self.spectrum_plotter.ax
                if ax:
                    ax.set_xlim(self.saved_xlim)
                    ax.set_ylim(self.saved_ylim)
                    logger.debug(f"Restored zoom: xlim={self.saved_xlim}, ylim={self.saved_ylim}")

    def update_spectrum_plot(self, preserve_zoom: bool = True):
        """Update the spectrum plot with current data.

        This method:
        1. Saves current zoom state (if preserve_zoom=True)
        2. Checks if NMR data is loaded
        3. Plots the spectrum using the spectrum plotter
        4. Overlays peaks if they are loaded
        5. Uses current contour parameters from GUI
        6. Restores zoom state (if preserve_zoom=True)
        7. Refreshes the display

        Args:
            preserve_zoom: If True, saves and restores zoom state during update.
                          Set to False when you want the view to reset.

        Note: This is called automatically after loading data,
        or can be called manually to refresh the plot.
        """
        # Check if NMR data is available
        if not hasattr(self.integrator, 'nmr_data') or self.integrator.nmr_data is None:
            logger.debug("No NMR data available for plotting")
            return

        try:
            logger.debug("Updating spectrum plot")

            # Save current zoom state BEFORE replotting (v0o9 line 5291)
            if preserve_zoom:
                self.save_current_zoom()

            # Clear and replot spectrum
            self.spectrum_plotter.plot_spectrum(
                self.integrator,
                contour_min_level=self.contour_min,
                contour_levels=self.contour_levels,
                contour_increment=self.contour_increment
            )

            # Add peaks if loaded
            has_detected = hasattr(self.integrator, 'fitted_peaks') and self.integrator.fitted_peaks
            has_reference = hasattr(self.integrator, 'peak_list') and self.integrator.peak_list is not None

            if has_detected or has_reference:
                # Plot both detected and reference peaks
                self.spectrum_plotter.plot_peaks(
                    self.integrator,
                    show_detected=has_detected,
                    show_assigned=has_reference
                )

                if has_detected:
                    logger.debug(f"Plotted {len(self.integrator.fitted_peaks)} detected peaks")
                if has_reference:
                    logger.debug(f"Plotted {len(self.integrator.peak_list)} reference peaks")

            # Restore zoom state AFTER replotting (v0o9 line 5313)
            if preserve_zoom:
                self.restore_zoom()

            # Refresh the canvas
            self.spectrum_plotter.refresh()

            logger.debug("Spectrum plot updated successfully")

        except Exception as e:
            logger.error(f"Failed to update spectrum plot: {e}")
            # Don't raise - allow the app to continue even if plotting fails
            self.update_status(f"Warning: Plot update failed - {e}")

    def update_file_status_labels(self, nmr_file: Optional[str] = None,
                                   peak_file: Optional[str] = None):
        """Update the file status labels in the Data Loading section.

        Args:
            nmr_file: Path to NMR spectrum file (None if not loaded)
            peak_file: Path to peak list file (None if not loaded)

        Note:
            - Long paths are abbreviated to show just filename
            - "Not loaded" message shown when file is None
            - Labels use secondary text color
        """
        # Update NMR file label
        if self.nmr_file_label:
            if nmr_file:
                # Abbreviate path if too long (show only filename)
                filename = Path(nmr_file).name if nmr_file else "Unknown"
                self.nmr_file_label.setText(f"NMR Spectrum: {filename}")
                self.nmr_file_label.setToolTip(str(nmr_file))  # Full path in tooltip
            else:
                self.nmr_file_label.setText("NMR Spectrum: Not loaded")
                self.nmr_file_label.setToolTip("")

        # Update peak file label
        if self.peak_file_label:
            if peak_file:
                # Abbreviate path if too long (show only filename)
                filename = Path(peak_file).name if peak_file else "Unknown"
                self.peak_file_label.setText(f"Peak List: {filename}")
                self.peak_file_label.setToolTip(str(peak_file))  # Full path in tooltip
            else:
                self.peak_file_label.setText("Peak List: Not loaded")
                self.peak_file_label.setToolTip("")

    def setup_menu_bar(self):
        """Create and configure the menu bar.

        Exact copy of v0o9 create_menu() (lines 975-1047) using Qt instead of Tkinter.

        Menu structure matches v0o9:
            File: Export Peak List, Export Results, Save Configuration, Exit
            Processing: Detect Peaks, Integrate Peaks, Fit Selected/All, Series Integration
            View: Reset Zoom, Center on Peak, Show toggles (Detected/Reference/Fitted)
            Results Analysis: Browsers, Export Matrix/Report/Batch
            Modules: DynamiXs
            Tools: Configuration, Statistics, Validate Files
            Help: User Guide, About
        """
        from PySide6.QtGui import QAction
        from PySide6.QtWidgets import QMenuBar

        logger.debug("Creating menu bar (v0o9 exact copy)")

        menu_bar = self.create_menu_bar()

        # === File Menu (v0o9 lines 980-987) ===
        file_menu = menu_bar.addMenu("&File")

        # Project management actions
        open_project_action = QAction("Open Project...", self)
        open_project_action.triggered.connect(self.open_project)
        file_menu.addAction(open_project_action)

        save_project_action = QAction("Save Project", self)
        save_project_action.triggered.connect(self.save_project)
        file_menu.addAction(save_project_action)

        save_project_as_action = QAction("Save Project As...", self)
        save_project_as_action.triggered.connect(self.save_project_as)
        file_menu.addAction(save_project_as_action)

        file_menu.addSeparator()

        export_peak_list_action = QAction("Export Peak List...", self)
        export_peak_list_action.triggered.connect(self.export_peak_list)
        file_menu.addAction(export_peak_list_action)

        export_results_action = QAction("Export Results...", self)
        export_results_action.triggered.connect(self.export_results)
        file_menu.addAction(export_results_action)

        file_menu.addSeparator()

        exit_action = QAction("Exit", self)
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        # === Results Analysis Menu (v0o9 lines 1020-1029) ===
        results_menu = menu_bar.addMenu("&Results Analysis")

        browse_series_action = QAction("Browse Series Results...", self)
        browse_series_action.triggered.connect(self.open_results_browser)
        results_menu.addAction(browse_series_action)

        # Browse Individual Spectra hidden but code kept for future use
        browse_spectra_action = QAction("Browse Individual Spectra...", self)
        browse_spectra_action.triggered.connect(self.open_spectrum_browser)
        browse_spectra_action.setVisible(False)  # Hidden but functional
        results_menu.addAction(browse_spectra_action)

        multi_viewer_action = QAction("Multi-Spectrum Overlay Viewer...", self)
        multi_viewer_action.triggered.connect(self.open_multi_spectrum_viewer)
        results_menu.addAction(multi_viewer_action)

        series_manager_action = QAction("Manage Saved Series...", self)
        series_manager_action.triggered.connect(self.open_series_manager)
        results_menu.addAction(series_manager_action)

        results_menu.addSeparator()

        export_matrix_action = QAction("Export Data Matrix...", self)
        export_matrix_action.triggered.connect(self.export_data_matrix)
        results_menu.addAction(export_matrix_action)

        export_report_action = QAction("Export Analysis Report...", self)
        export_report_action.triggered.connect(self.export_analysis_report)
        results_menu.addAction(export_report_action)

        batch_export_action = QAction("Batch Export Results...", self)
        batch_export_action.triggered.connect(self.batch_export_results)
        results_menu.addAction(batch_export_action)

        # === Modules Menu (v0o9 lines 1031-1034) ===
        modules_menu = menu_bar.addMenu("&Modules")

        dynamixs_action = QAction("DynamiXs Relaxation Analysis", self)
        dynamixs_action.triggered.connect(self.launch_dynamixs)
        modules_menu.addAction(dynamixs_action)

        # === Help Menu (v0o9 lines 1043-1047) ===
        help_menu = menu_bar.addMenu("&Help")

        user_guide_action = QAction("User Guide", self)
        user_guide_action.triggered.connect(self.show_help)
        help_menu.addAction(user_guide_action)

        about_action = QAction("About", self)
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)

        logger.debug("Menu bar created with all actions connected")

    # =================== MENU ACTION HANDLERS (v0o9 exact copies) ===================

    def _on_show_detected_toggled(self, checked: bool):
        """Handle Show Detected Peaks menu toggle (v0o9 line 1006)."""
        self.show_detected = checked
        self.update_spectrum_plot()
        logger.debug(f"Show Detected Peaks: {checked}")

    def _on_show_assigned_toggled(self, checked: bool):
        """Handle Show Reference Peaks menu toggle (v0o9 lines 1008-1014)."""
        self.show_assigned = checked
        self.update_spectrum_plot()
        logger.debug(f"Show Reference Peaks: {checked}")

    def _on_show_fitted_toggled(self, checked: bool):
        """Handle Show Fitted Curves menu toggle (v0o9 line 1016)."""
        self.show_fitted_curves = checked
        self.update_spectrum_plot()
        logger.debug(f"Show Fitted Curves: {checked}")

    # =================== FILE MENU HANDLERS ===================

    def export_peak_list(self):
        """Export current peak list (v0o9 lines 6022-6062).

        Exports fitted_peaks to CSV or TXT file with columns:
        Assignment, Position_X, Position_Y, Detected, SNR, Quality
        """
        from PySide6.QtWidgets import QFileDialog, QMessageBox
        import pandas as pd

        # Validate data exists
        if not hasattr(self.integrator, 'fitted_peaks') or not self.integrator.fitted_peaks:
            QMessageBox.critical(self, "Error", "No peaks to export")
            return

        # Get filename
        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Export Peak List",
            "",
            "CSV files (*.csv);;Text files (*.txt);;All files (*.*)"
        )

        if filename:
            try:
                peak_data = []
                for peak in self.integrator.fitted_peaks:
                    peak_data.append({
                        'Assignment': peak.get('assignment', ''),
                        'Position_X': peak.get('ppm_x', 0),
                        'Position_Y': peak.get('ppm_y', 0),
                        'Detected': peak.get('detected', False),
                        'SNR': peak.get('snr', 0),
                        'Quality': peak.get('detection_quality', 'Unknown')
                    })

                df = pd.DataFrame(peak_data)
                if filename.endswith('.csv'):
                    df.to_csv(filename, index=False, float_format='%.6f')
                else:
                    df.to_csv(filename, sep='\t', index=False, float_format='%.6f')

                self.update_status(f"Peak list exported: {len(peak_data)} peaks")
                QMessageBox.information(self, "Export Successful",
                    f"Peak list exported to:\n{filename}")

            except Exception as e:
                self.update_status(f"Export failed: {str(e)}")
                QMessageBox.critical(self, "Export Error", str(e))

    def export_results(self):
        """Export current processing results (v0o9 lines 6064-6094).

        Exports integration_results to CSV or Excel file.
        """
        from PySide6.QtWidgets import QFileDialog, QMessageBox
        import pandas as pd

        # Validate data exists
        if not hasattr(self.integrator, 'integration_results') or not self.integrator.integration_results:
            QMessageBox.critical(self, "Error", "No integration results to export")
            return

        # Get filename
        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Export Integration Results",
            "",
            "CSV files (*.csv);;Excel files (*.xlsx);;All files (*.*)"
        )

        if filename:
            try:
                df = pd.DataFrame(self.integrator.integration_results)

                if filename.endswith('.xlsx'):
                    df.to_excel(filename, index=False)
                else:
                    df.to_csv(filename, index=False, float_format='%.6f')

                self.update_status(f"Results exported: {len(df)} integrations")
                QMessageBox.information(self, "Export Successful",
                    f"Results exported to:\n{filename}")

            except Exception as e:
                self.update_status(f"Export failed: {str(e)}")
                QMessageBox.critical(self, "Export Error", str(e))

    def save_config(self):
        """Save configuration to file (v0o9 lines 6128-6140)."""
        from PySide6.QtWidgets import QFileDialog, QMessageBox

        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Save Configuration",
            "",
            "JSON files (*.json);;All files (*.*)"
        )

        if filename:
            if self.config_manager.export_config(filename):
                QMessageBox.information(self, "Success",
                    f"Configuration saved to:\n{filename}")
            else:
                QMessageBox.critical(self, "Error", "Failed to save configuration")

    # =================== PROCESSING MENU HANDLERS ===================

    def integrate_peaks(self):
        """Integrate detected peaks (v0o9 lines 3804-3839).

        Performs peak integration after detection using standard mode.
        """
        from PySide6.QtWidgets import QMessageBox

        # CRITICAL: Sync GUI parameters BEFORE integration
        self.sync_parameters_to_integrator()

        # Validate data exists
        if not hasattr(self.integrator, 'fitted_peaks') or not self.integrator.fitted_peaks:
            QMessageBox.critical(self, "Error", "Please detect peaks first")
            return

        try:
            self.update_status("Integrating peaks (standard mode)...")
            QApplication.setOverrideCursor(Qt.WaitCursor)

            results = self.integrator.integrate_peaks()
            if results is not None and len(results) > 0:
                detected_count = sum(1 for r in results if r.get('Integration_Method') != 'Reference')
                good_quality = sum(1 for r in results if r.get('Quality') in ['Excellent', 'Good'])
                self.update_status(f"Integration: {detected_count}/{len(results)} peaks, {good_quality} good quality")
                self.update_statistics_panel()
            else:
                self.update_status("Integration failed")

        except Exception as e:
            self.update_status(f"Integration error: {str(e)}")
            QMessageBox.critical(self, "Integration Error", str(e))

        finally:
            QApplication.restoreOverrideCursor()

    def fit_selected_peak(self):
        """Fit Voigt profile to selected peak (v0o9 lines 3841-3917).

        Performs enhanced Voigt fitting on the currently selected peak.
        """
        from PySide6.QtWidgets import QMessageBox

        # CRITICAL: Sync GUI parameters BEFORE fitting
        self.sync_parameters_to_integrator()

        # Validate data exists
        if not hasattr(self.integrator, 'peak_list') or self.integrator.peak_list is None:
            QMessageBox.critical(self, "Error", "No peaks available")
            return

        # Get selected peak from navigator
        if not self.selected_peak_info:
            QMessageBox.critical(self, "Error", "No peak selected. Click a peak first.")
            return

        peak_index = self.selected_peak_info.get('index', -1)
        if peak_index < 0 or peak_index >= len(self.integrator.peak_list):
            QMessageBox.critical(self, "Error", f"Invalid peak index: {peak_index}")
            return

        try:
            self.update_status(f"Fitting Voigt profile for peak {peak_index + 1}...")
            QApplication.setOverrideCursor(Qt.WaitCursor)

            # Get peak position
            peak_row = self.integrator.peak_list.iloc[peak_index]
            peak_x = float(peak_row.get('Position_X', peak_row.get('ppm_x', 0)))
            peak_y = float(peak_row.get('Position_Y', peak_row.get('ppm_y', 0)))
            assignment = peak_row.get('Assignment', f'Peak_{peak_index + 1}')

            # Perform enhanced Voigt fitting
            result = self.integrator.enhanced_peak_fitting(peak_x, peak_y, assignment)

            if result:
                result['peak_number'] = peak_index + 1

                # Store in fitted_peaks
                if not hasattr(self.integrator, 'fitted_peaks'):
                    self.integrator.fitted_peaks = []

                # Update or add to fitted_peaks
                found = False
                for i, existing_peak in enumerate(self.integrator.fitted_peaks):
                    if existing_peak.get('peak_number') == peak_index + 1:
                        self.integrator.fitted_peaks[i] = result
                        found = True
                        break

                if not found:
                    self.integrator.fitted_peaks.append(result)

                # Update 3D Voigt plot
                self.current_voigt_result = result
                if hasattr(self, 'voigt_3d_plotter') and self.voigt_3d_plotter:
                    self.voigt_3d_plotter.plot_voigt_analysis_3d(result)

                quality = result.get('fitting_quality', 'Unknown')
                self.update_status(f"Voigt fit: {assignment} - Quality: {quality}")
                self.update_statistics_panel()
            else:
                self.update_status(f"Voigt fitting failed for {assignment}")

        except Exception as e:
            self.update_status(f"Fitting error: {str(e)}")
            QMessageBox.critical(self, "Fitting Error", str(e))

        finally:
            QApplication.restoreOverrideCursor()

    def start_series_integration(self):
        """Start series integration (v0o9 lines 4452-4600).

        Processes multiple spectra using MultiSpectrumProcessor.
        Opens SeriesIntegrationDialog for configuration and progress tracking.
        """
        # CRITICAL: Sync GUI parameters BEFORE series integration
        self.sync_parameters_to_integrator()
        logger.info(f"Parameters synced for series: fix_positions={self.fix_positions}, fix_linewidths={self.fix_linewidths}")

        from lunaNMR.gui.dialogs import SeriesIntegrationDialog

        # Open series integration dialog
        dialog = SeriesIntegrationDialog(parent=self, main_window=self)

        # Connect signal to handle results
        dialog.processing_complete.connect(self._on_series_complete)

        dialog.exec()

        self.update_status("Series integration dialog closed")

    def _on_series_complete(self, batch_results):
        """Handle series integration completion.

        Args:
            batch_results: BatchResults object from series processing
        """
        self.batch_results = batch_results

        if batch_results and hasattr(batch_results, 'get_summary'):
            summary = batch_results.get_summary()
            self.update_status(
                f"Series complete: {summary.get('successful', 0)}/{summary.get('total_spectra', 0)} spectra"
            )
        else:
            self.update_status("Series integration complete")

    # =================== VIEW MENU HANDLERS ===================

    def reset_view(self):
        """Reset main spectrum view (v0o9 lines 6243-6248)."""
        if hasattr(self, 'spectrum_plotter') and self.spectrum_plotter:
            self.spectrum_plotter.reset_zoom()
        self.update_status("View reset")

    def center_on_selected_peak(self):
        """Center view on selected peak (v0o9 lines 5417-5434)."""
        if not hasattr(self.integrator, 'peak_list') or self.integrator.peak_list is None:
            return

        if not self.selected_peak_info:
            self.update_status("No peak selected")
            return

        peak_index = self.selected_peak_info.get('index', -1)
        if peak_index < 0 or peak_index >= len(self.integrator.peak_list):
            return

        # Get peak position
        peak_row = self.integrator.peak_list.iloc[peak_index]
        peak_x = float(peak_row.get('Position_X', peak_row.get('ppm_x', 0)))
        peak_y = float(peak_row.get('Position_Y', peak_row.get('ppm_y', 0)))

        # Highlight and center on peak
        if hasattr(self, 'spectrum_plotter') and self.spectrum_plotter:
            self.spectrum_plotter.highlight_peak(peak_x, peak_y)
            self.spectrum_plotter.set_zoom(peak_x, peak_y, self.zoom_x_range, self.zoom_y_range)

        assignment = peak_row.get('Assignment', f'Peak_{peak_index + 1}')
        self.update_status(f"Centered on peak {peak_index + 1}: {assignment}")

    # =================== RESULTS ANALYSIS MENU HANDLERS ===================

    def open_results_browser(self):
        """Open comprehensive results browser window (v0o9 lines 6785-7080).

        Displays series results in a multi-tab browser window.
        """
        from PySide6.QtWidgets import QMessageBox
        from lunaNMR.gui.dialogs import ResultsBrowserDialog

        # Validate data exists
        if not hasattr(self, 'batch_results') or not self.batch_results:
            QMessageBox.information(self, "No Data",
                "No series results available to browse.\n\n"
                "Please run series integration first.")
            return

        dialog = ResultsBrowserDialog(
            parent=self,
            batch_results=self.batch_results,
            main_window=self
        )
        dialog.exec()

    def open_spectrum_browser(self):
        """Open spectrum browser (v0o9 lines 7446-7520).

        Displays individual spectrum results with navigation.
        """
        from PySide6.QtWidgets import QMessageBox
        from lunaNMR.gui.dialogs import SpectrumBrowserDialog

        # Validate data exists
        if not hasattr(self, 'batch_results') or not self.batch_results:
            QMessageBox.warning(self, "No Data",
                "No series results available. Run series integration first.")
            return

        # Get series processor if available
        series_processor = getattr(self, 'series_processor', None)
        original_data_folder = getattr(self, 'current_nmr_folder', None)

        # Open spectrum browser dialog
        dialog = SpectrumBrowserDialog(
            parent=self,
            batch_results=self.batch_results,
            series_processor=series_processor,
            original_data_folder=original_data_folder
        )
        dialog.show()  # Non-modal

        self.update_status("Spectrum browser opened")

    def open_multi_spectrum_viewer(self):
        """Open multi-spectrum overlay viewer (v0o9 lines 7525-7575).

        Displays multiple spectra overlaid for comparison.
        If multiple series exist, prompts user to select which one to view.
        """
        import os
        from PySide6.QtWidgets import QMessageBox, QInputDialog
        from lunaNMR.gui.dialogs import MultiSpectrumViewerDialog

        # Check for multiple saved series
        saved_series = getattr(self, 'saved_series', {}) or {}

        # If multiple series exist, prompt user to select one
        if len(saved_series) > 1:
            # Build list of options with spectrum counts
            options = []
            for name, batch in saved_series.items():
                count = len(batch.results) if hasattr(batch, 'results') else 0
                options.append(f"{name} ({count} spectra)")

            selected, ok = QInputDialog.getItem(
                self,
                "Select Series",
                "Multiple series integrations found.\nSelect which one to view:",
                options,
                0,  # Default to first item
                False  # Not editable
            )

            if not ok or not selected:
                return

            # Extract series name from selection
            series_name = selected.rsplit(' (', 1)[0]
            batch_results = saved_series.get(series_name)

            if not batch_results:
                QMessageBox.warning(self, "Error", f"Could not load series '{series_name}'")
                return
        elif len(saved_series) == 1:
            # Only one series, use it directly
            series_name = list(saved_series.keys())[0]
            batch_results = saved_series[series_name]
        elif hasattr(self, 'batch_results') and self.batch_results:
            # Fallback to batch_results for backward compatibility
            batch_results = self.batch_results
        else:
            QMessageBox.warning(self, "No Data",
                "No series results available. Run series integration first.")
            return

        # Get data folder from batch_results metadata (preferred) or fallback to current_nmr_folder
        data_folder = None
        if hasattr(batch_results, 'metadata'):
            data_folder = batch_results.metadata.get('data_folder')
        elif isinstance(batch_results, dict):
            metadata = batch_results.get('metadata', {})
            data_folder = metadata.get('data_folder')

        # Fallback to current_nmr_folder if not in metadata
        if not data_folder:
            data_folder = getattr(self, 'current_nmr_folder', None)

        logger.debug(f"Data folder for spectrum paths: {data_folder}")

        # Get all results as list
        all_results = []
        if hasattr(batch_results, 'results'):
            for name, result in batch_results.results.items():
                result_copy = result.copy()
                result_copy['spectrum_name'] = name

                # CRITICAL: Construct full path for spectrum_file (v0.9 main_gui.py:7547-7553)
                # result may contain 'spectrum_file' but it's just a filename, not full path
                if 'spectrum_file' not in result_copy or not os.path.isabs(result_copy.get('spectrum_file', '')):
                    if data_folder:
                        result_copy['spectrum_file'] = os.path.join(data_folder, name)
                        logger.debug(f"Constructed full path for {name}: {result_copy['spectrum_file']}")
                    else:
                        logger.warning(f"No data folder - cannot construct path for {name}")

                all_results.append(result_copy)

        # Sort by spectrum name for consistent ordering
        all_results.sort(key=lambda x: x.get('spectrum_name', ''))

        # Open multi-spectrum viewer dialog
        dialog = MultiSpectrumViewerDialog(
            parent=self,
            all_results=all_results,
            file_manager=self.file_manager
        )
        dialog.show()  # Non-modal

        self.update_status(f"Multi-spectrum viewer opened with {len(all_results)} spectra")

    def open_series_manager(self):
        """Open Series Manager dialog for viewing/renaming/deleting saved series."""
        from lunaNMR.gui.dialogs import SeriesManagerDialog

        dialog = SeriesManagerDialog(parent=self, main_window=self)
        dialog.exec()

        self.update_status("Series manager closed")

    def export_data_matrix(self):
        """Export comprehensive data matrix (v0o9 lines 7100-7230).

        Exports peak data across all spectra in series to CSV or Excel.
        """
        import pandas as pd
        from PySide6.QtWidgets import QFileDialog, QMessageBox

        # Validate data exists
        if not hasattr(self, 'batch_results') or not self.batch_results:
            QMessageBox.warning(self, "No Data",
                "No series results to export. Run series integration first.")
            return

        # Get save filename
        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Export Comprehensive Data Matrix",
            "",
            "CSV files (*.csv);;Excel files (*.xlsx);;All files (*.*)"
        )

        if not filename:
            return

        try:
            # Build data matrix from batch results
            all_results = []

            if hasattr(self.batch_results, 'results'):
                for spectrum_name, spectrum_data in self.batch_results.results.items():
                    # Get fitted results from spectrum data
                    fitted_results = spectrum_data.get('fitted_results', [])
                    integration_results = spectrum_data.get('integration_results', [])

                    # Use integration_results if available, else fitted_results
                    results_to_export = integration_results if integration_results else fitted_results

                    for peak_result in results_to_export:
                        if isinstance(peak_result, dict):
                            result_row = peak_result.copy()

                            # Add spectrum metadata
                            result_row['spectrum'] = spectrum_name
                            result_row['spectrum_success'] = spectrum_data.get('status') == 'success'
                            result_row['detection_rate'] = spectrum_data.get('detection_rate', 0)
                            result_row['total_peaks'] = spectrum_data.get('total_peaks', 0)

                            all_results.append(result_row)

            if not all_results:
                QMessageBox.warning(self, "No Data",
                    "No peak results found in batch results.")
                return

            # Create DataFrame
            df = pd.DataFrame(all_results)

            # Reorder columns for better readability
            priority_cols = ['spectrum', 'assignment', 'Assignment',
                           'pos_f2', 'pos_f1', 'Position_X', 'Position_Y',
                           'intensity', 'volume', 'r_squared', 'snr']
            ordered_cols = [c for c in priority_cols if c in df.columns]
            other_cols = [c for c in df.columns if c not in ordered_cols]
            df = df[ordered_cols + other_cols]

            # Export based on file extension
            if filename.endswith('.xlsx'):
                df.to_excel(filename, index=False, sheet_name='Peak_Data')
            else:
                df.to_csv(filename, index=False)

            QMessageBox.information(self, "Export Complete",
                f"Data matrix exported successfully!\n\n"
                f"File: {filename}\n"
                f"Spectra: {len(self.batch_results.results)}\n"
                f"Peak records: {len(all_results)}")

            self.update_status(f"Data matrix exported: {len(all_results)} records")

        except Exception as e:
            logger.error(f"Data matrix export error: {e}")
            QMessageBox.critical(self, "Export Error",
                f"Failed to export data matrix:\n{str(e)}")

    def export_analysis_report(self):
        """Export comprehensive analysis report (v0o9 lines 7234-7255).

        Generates JSON analysis report from series results.
        """
        import json
        from datetime import datetime
        from PySide6.QtWidgets import QFileDialog, QMessageBox

        if not hasattr(self, 'batch_results') or not self.batch_results:
            QMessageBox.warning(self, "No Data",
                "No series results available for report generation.")
            return

        # Get save filename
        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Export Analysis Report",
            "",
            "JSON files (*.json);;Text files (*.txt);;All files (*.*)"
        )

        if not filename:
            return

        try:
            # Build comprehensive report
            report = {
                'report_metadata': {
                    'generated_at': datetime.now().isoformat(),
                    'software': 'lunaNMR v1.0',
                    'report_type': 'Series Integration Analysis'
                },
                'summary': {},
                'quality_statistics': {},
                'spectra_details': {}
            }

            # Get summary if available
            if hasattr(self.batch_results, 'get_summary'):
                report['summary'] = self.batch_results.get_summary()
            else:
                # Build summary manually
                total = len(self.batch_results.results) if hasattr(self.batch_results, 'results') else 0
                successful = sum(1 for r in self.batch_results.results.values()
                               if r.get('status') == 'success') if hasattr(self.batch_results, 'results') else 0
                report['summary'] = {
                    'total_spectra': total,
                    'successful': successful,
                    'failed': total - successful,
                    'success_rate': (successful / total * 100) if total > 0 else 0
                }

            # Quality statistics
            if hasattr(self.batch_results, 'results'):
                r_squared_values = []
                detection_rates = []

                for name, data in self.batch_results.results.items():
                    if data.get('avg_r_squared'):
                        r_squared_values.append(data['avg_r_squared'])
                    if data.get('detection_rate'):
                        detection_rates.append(data['detection_rate'])

                    # Add to spectra details
                    report['spectra_details'][name] = {
                        'status': data.get('status', 'unknown'),
                        'detected_peaks': data.get('detected_peaks', 0),
                        'total_peaks': data.get('total_peaks', 0),
                        'detection_rate': data.get('detection_rate', 0),
                        'avg_r_squared': data.get('avg_r_squared', 0),
                        'processing_time': data.get('processing_time', 0)
                    }

                if r_squared_values:
                    import numpy as np
                    report['quality_statistics']['r_squared'] = {
                        'mean': float(np.mean(r_squared_values)),
                        'std': float(np.std(r_squared_values)),
                        'min': float(np.min(r_squared_values)),
                        'max': float(np.max(r_squared_values))
                    }

                if detection_rates:
                    import numpy as np
                    report['quality_statistics']['detection_rate'] = {
                        'mean': float(np.mean(detection_rates)),
                        'std': float(np.std(detection_rates)),
                        'min': float(np.min(detection_rates)),
                        'max': float(np.max(detection_rates))
                    }

            # Write report
            with open(filename, 'w') as f:
                json.dump(report, f, indent=2)

            QMessageBox.information(self, "Export Complete",
                f"Analysis report exported successfully!\n\n"
                f"File: {filename}")

            self.update_status("Analysis report exported")

        except Exception as e:
            logger.error(f"Analysis report export error: {e}")
            QMessageBox.critical(self, "Export Error",
                f"Failed to export analysis report:\n{str(e)}")

    def batch_export_results(self):
        """Batch export results in multiple formats (v0o9 lines 7259-7290).

        Exports CSV, Excel, and JSON versions of results to selected directory.
        """
        import os
        import json
        import pandas as pd
        from datetime import datetime
        from PySide6.QtWidgets import (QMessageBox, QFileDialog, QDialog,
                                       QVBoxLayout, QHBoxLayout, QCheckBox,
                                       QLabel, QPushButton, QGroupBox)

        if not hasattr(self, 'batch_results') or not self.batch_results:
            QMessageBox.warning(self, "No Data",
                "No series results to export.")
            return

        # Create format selection dialog
        dialog = QDialog(self)
        dialog.setWindowTitle("Batch Export Options")
        dialog.setMinimumWidth(350)
        layout = QVBoxLayout(dialog)

        # Format options group
        format_group = QGroupBox("Export Formats")
        format_layout = QVBoxLayout(format_group)

        csv_check = QCheckBox("Individual CSV files (one per spectrum)")
        csv_check.setChecked(True)
        excel_check = QCheckBox("Combined Excel workbook (all spectra)")
        excel_check.setChecked(True)
        json_check = QCheckBox("JSON analysis report")
        json_check.setChecked(True)
        summary_check = QCheckBox("Summary CSV (one row per spectrum)")
        summary_check.setChecked(True)

        format_layout.addWidget(csv_check)
        format_layout.addWidget(excel_check)
        format_layout.addWidget(json_check)
        format_layout.addWidget(summary_check)
        layout.addWidget(format_group)

        # Button layout
        button_layout = QHBoxLayout()
        export_btn = QPushButton("Select Directory && Export")
        cancel_btn = QPushButton("Cancel")
        button_layout.addWidget(export_btn)
        button_layout.addWidget(cancel_btn)
        layout.addLayout(button_layout)

        export_btn.clicked.connect(dialog.accept)
        cancel_btn.clicked.connect(dialog.reject)

        if dialog.exec() != QDialog.Accepted:
            return

        # Check at least one format selected
        if not any([csv_check.isChecked(), excel_check.isChecked(),
                    json_check.isChecked(), summary_check.isChecked()]):
            QMessageBox.warning(self, "No Format Selected",
                "Please select at least one export format.")
            return

        # Select output directory
        output_dir = QFileDialog.getExistingDirectory(
            self, "Select Export Directory", "")

        if not output_dir:
            return

        try:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            exported_files = []

            # Build data from batch results
            all_results = []
            summary_data = []

            if hasattr(self.batch_results, 'results'):
                for spectrum_name, spectrum_data in self.batch_results.results.items():
                    fitted_results = spectrum_data.get('fitted_results', [])
                    integration_results = spectrum_data.get('integration_results', [])
                    results_to_export = integration_results if integration_results else fitted_results

                    # Summary row for this spectrum
                    summary_row = {
                        'spectrum': spectrum_name,
                        'status': spectrum_data.get('status', 'unknown'),
                        'total_peaks': spectrum_data.get('total_peaks', 0),
                        'detected_peaks': spectrum_data.get('detected_peaks', 0),
                        'detection_rate': spectrum_data.get('detection_rate', 0),
                        'avg_r_squared': 0,
                        'avg_snr': 0
                    }

                    # Calculate averages
                    r_squared_values = []
                    snr_values = []

                    for peak_result in results_to_export:
                        if isinstance(peak_result, dict):
                            result_row = peak_result.copy()
                            result_row['spectrum'] = spectrum_name
                            all_results.append(result_row)

                            if 'r_squared' in peak_result:
                                r_squared_values.append(peak_result['r_squared'])
                            if 'snr' in peak_result:
                                snr_values.append(peak_result['snr'])

                    if r_squared_values:
                        summary_row['avg_r_squared'] = sum(r_squared_values) / len(r_squared_values)
                    if snr_values:
                        summary_row['avg_snr'] = sum(snr_values) / len(snr_values)

                    summary_data.append(summary_row)

                    # Export individual CSV if requested
                    if csv_check.isChecked() and results_to_export:
                        safe_name = "".join(c if c.isalnum() or c in '-_' else '_' for c in spectrum_name)
                        csv_file = os.path.join(output_dir, f"{safe_name}_{timestamp}.csv")
                        spectrum_df = pd.DataFrame([r for r in results_to_export if isinstance(r, dict)])
                        spectrum_df.to_csv(csv_file, index=False)
                        exported_files.append(csv_file)

            # Export combined Excel workbook
            if excel_check.isChecked() and all_results:
                excel_file = os.path.join(output_dir, f"batch_results_{timestamp}.xlsx")
                with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
                    # All peaks sheet
                    all_df = pd.DataFrame(all_results)
                    all_df.to_excel(writer, sheet_name='All_Peaks', index=False)
                    # Summary sheet
                    if summary_data:
                        summary_df = pd.DataFrame(summary_data)
                        summary_df.to_excel(writer, sheet_name='Summary', index=False)
                exported_files.append(excel_file)

            # Export JSON report
            if json_check.isChecked():
                json_file = os.path.join(output_dir, f"analysis_report_{timestamp}.json")
                report = {
                    'metadata': {
                        'generated_at': datetime.now().isoformat(),
                        'software': 'lunaNMR v1.0',
                        'total_spectra': len(summary_data),
                        'total_peaks': len(all_results)
                    },
                    'summary': summary_data,
                    'peaks': all_results
                }
                with open(json_file, 'w') as f:
                    json.dump(report, f, indent=2, default=str)
                exported_files.append(json_file)

            # Export summary CSV
            if summary_check.isChecked() and summary_data:
                summary_file = os.path.join(output_dir, f"summary_{timestamp}.csv")
                summary_df = pd.DataFrame(summary_data)
                summary_df.to_csv(summary_file, index=False)
                exported_files.append(summary_file)

            QMessageBox.information(self, "Batch Export Complete",
                f"Successfully exported {len(exported_files)} files to:\n"
                f"{output_dir}\n\n"
                f"Files created:\n" +
                "\n".join([f"• {os.path.basename(f)}" for f in exported_files[:10]]) +
                (f"\n... and {len(exported_files) - 10} more" if len(exported_files) > 10 else ""))

            self.update_status(f"Batch export complete: {len(exported_files)} files")

        except Exception as e:
            logger.error(f"Batch export error: {e}")
            QMessageBox.critical(self, "Export Error",
                f"Failed to complete batch export:\n{str(e)}")

    # =================== MODULES MENU HANDLERS ===================

    def launch_dynamixs(self):
        """Launch DynamiXs Relaxation Analysis as embedded dialog.

        Opens DynamiXs as an embedded dialog within lunaNMR, allowing full
        state synchronization and project save/load integration.
        """
        from lunaNMR.gui.dialogs import DynamiXsDialog

        # Reuse existing dialog if already open
        if hasattr(self, 'dynamixs_dialog') and self.dynamixs_dialog and self.dynamixs_dialog.isVisible():
            self.dynamixs_dialog.raise_()
            self.dynamixs_dialog.activateWindow()
            return

        # Create and show new dialog
        self.dynamixs_dialog = DynamiXsDialog(parent=self, main_window=self)
        self.dynamixs_dialog.show()

    # =================== TOOLS MENU HANDLERS ===================

    def open_config_dialog(self):
        """Open configuration dialog (v0o9 lines 6207-6218).

        Opens ConfigurationDialog for viewing and managing settings.
        """
        from lunaNMR.gui.dialogs import ConfigurationDialog

        dialog = ConfigurationDialog(self, main_window=self)
        dialog.exec()

        # Update status after dialog closes
        self.update_status("Configuration dialog closed")

    def show_statistics(self):
        """Show statistics in separate window (v0o9 lines 6268-6302).

        Opens a statistics dialog showing processing metrics.
        """
        from PySide6.QtWidgets import QDialog, QVBoxLayout, QTextEdit, QMessageBox

        # Create statistics dialog
        dialog = QDialog(self)
        dialog.setWindowTitle("Processing Statistics")
        dialog.resize(500, 400)

        layout = QVBoxLayout(dialog)

        text_edit = QTextEdit()
        text_edit.setReadOnly(True)
        layout.addWidget(text_edit)

        # Generate statistics text
        stats_text = "=== Processing Statistics ===\n\n"

        if hasattr(self.integrator, 'fitted_peaks') and self.integrator.fitted_peaks:
            peaks = self.integrator.fitted_peaks
            detected = sum(1 for p in peaks if p.get('detected', False))
            stats_text += f"Fitted Peaks: {len(peaks)}\n"
            stats_text += f"Detected: {detected}\n"
            stats_text += f"Detection Rate: {detected / len(peaks) * 100:.1f}%\n\n"

            # Quality breakdown
            excellent = sum(1 for p in peaks if p.get('fitting_quality') == 'Excellent')
            good = sum(1 for p in peaks if p.get('fitting_quality') == 'Good')
            fair = sum(1 for p in peaks if p.get('fitting_quality') == 'Fair')
            poor = sum(1 for p in peaks if p.get('fitting_quality') == 'Poor')
            stats_text += f"Quality Breakdown:\n"
            stats_text += f"  Excellent: {excellent}\n"
            stats_text += f"  Good: {good}\n"
            stats_text += f"  Fair: {fair}\n"
            stats_text += f"  Poor: {poor}\n"
        else:
            stats_text += "No fitting results available.\n"
            stats_text += "Run 'Fit All Peaks' to generate statistics.\n"

        if hasattr(self, 'batch_results') and self.batch_results:
            stats_text += f"\n=== Series Results ===\n"
            stats_text += f"Series data available\n"

        text_edit.setPlainText(stats_text)
        dialog.exec()

    def validate_current_files(self):
        """Validate currently loaded files (v0o9 lines 6221-6240).

        Checks NMR and peak files for validity.
        """
        results = {'nmr_valid': False, 'peak_valid': False}
        messages = []

        # Validate NMR file
        if self.current_nmr_file:
            valid, message = self.file_manager.validate_nmr_file(self.current_nmr_file)
            results['nmr_valid'] = valid
            status = "OK" if valid else "FAILED"
            messages.append(f"NMR validation: {status} - {message}")

        # Validate peak file
        if self.current_peak_file:
            valid, message = self.file_manager.validate_peak_file(self.current_peak_file)
            results['peak_valid'] = valid
            status = "OK" if valid else "FAILED"
            messages.append(f"Peak list validation: {status} - {message}")

        # Show results
        from PySide6.QtWidgets import QMessageBox
        if messages:
            QMessageBox.information(self, "File Validation",
                "\n".join(messages))
        else:
            QMessageBox.information(self, "File Validation",
                "No files loaded to validate.")

        return results

    # =================== HELP MENU HANDLERS ===================

    def show_help(self):
        """Show help dialog (v0o9 lines 7579-7632).

        Displays user guide in a dialog window.
        """
        from PySide6.QtWidgets import QDialog, QVBoxLayout, QTextEdit

        help_text = """
NMR Peaks Series Analysis - User Guide

PROCESSING MODES:
• Full Detection: Complete spectrum analysis with validation
• In-Place Fitting: Direct fitting to reference positions

BASIC WORKFLOW:
1. Select NMR spectra folder and peak list folder
2. Choose processing mode using radio buttons
3. Load individual files or start series processing
4. Adjust detection and fitting parameters
5. View results in multiple visualization tabs

KEY FEATURES:
• Enhanced file management with preview
• Dual-mode processing (Full/In-Place)
• Advanced Voigt profile fitting
• Comprehensive series integration
• Statistics and quality assessment
• Configuration management

KEYBOARD SHORTCUTS:
• F1: Show this help
• F5: Refresh current view
• Ctrl+O: Open file
• Ctrl+S: Save results
• Ctrl+Q: Quit

For more information, visit the project documentation.
        """

        dialog = QDialog(self)
        dialog.setWindowTitle("User Guide")
        dialog.resize(600, 500)

        layout = QVBoxLayout(dialog)

        text_edit = QTextEdit()
        text_edit.setReadOnly(True)
        text_edit.setPlainText(help_text)
        layout.addWidget(text_edit)

        dialog.exec()

    def show_about(self):
        """Show about dialog (v0o9 lines 7634-7654).

        Displays application information.
        """
        from PySide6.QtWidgets import QMessageBox

        about_text = f"""
LunaNMR - NMR Peaks Series Analysis
Version 1.0.0 (Qt Port)

Enhanced Multi-Mode NMR Peak Detection and Analysis

Features:
• Dual processing modes (Full Detection / In-Place Fitting)
• Advanced Voigt profile fitting with PS2D
• Comprehensive series integration
• Enhanced visualization and statistics
• Configuration management

Qt port of original Tkinter version.
Developed using Python with PySide6, matplotlib, numpy, pandas, and scipy.
        """

        QMessageBox.about(self, "About LunaNMR", about_text)

    def apply_user_preferences(self):
        """Apply saved user preferences to the window.

        This restores:
        - Window geometry (size and position)
        - Splitter positions
        - Display options
        - Last used directories

        Uses Qt's QSettings for cross-platform persistence.
        """
        logger.debug("Applying user preferences")

        try:
            from PySide6.QtCore import QSettings

            settings = QSettings("LunaNMR", "LunaNMR_v1")

            # Restore window geometry
            geometry = settings.value("window/geometry")
            if geometry:
                self.restoreGeometry(geometry)
                logger.debug("Restored window geometry")

            # Restore window state (dock positions, toolbars, etc.)
            state = settings.value("window/state")
            if state:
                self.restoreState(state)
                logger.debug("Restored window state")

            # Restore splitter positions
            if hasattr(self, 'main_splitter'):
                splitter_state = settings.value("splitter/main")
                if splitter_state:
                    self.main_splitter.restoreState(splitter_state)
                    logger.debug("Restored main splitter state")

            # Restore last used directories
            last_nmr_folder = settings.value("paths/last_nmr_folder", "")
            if last_nmr_folder:
                self.current_nmr_folder = last_nmr_folder

            last_peak_folder = settings.value("paths/last_peak_folder", "")
            if last_peak_folder:
                self.current_peak_folder = last_peak_folder

            # Restore display options
            self.show_detected = settings.value("display/show_detected", True, type=bool)
            self.show_assigned = settings.value("display/show_assigned", True, type=bool)
            self.show_fitted_curves = settings.value("display/show_fitted_curves", True, type=bool)
            self.show_ellipses = settings.value("display/show_ellipses", False, type=bool)

            # Restore contour settings
            self.contour_levels = settings.value("contour/levels", 24, type=int)
            self.contour_min = settings.value("contour/min", 0.03, type=float)
            self.contour_increment = settings.value("contour/increment", 1.2, type=float)

            logger.debug("User preferences applied successfully")

        except Exception as e:
            logger.warning(f"Failed to apply user preferences: {e}")
            # Continue - don't fail startup for preference issues

    def sync_parameters_to_integrator(self):
        """Synchronize all GUI parameters to the integrator engine.

        This ensures that:
        1. Parameter manager has latest values from GUI
        2. Parameters are validated
        3. Integrator is configured with effective parameters
        4. Special parameters (S/N threshold, workflow mode) are synced

        This is called:
        - Once during initialization
        - Whenever parameters are changed by the user
        - Before starting any processing operation

        Based on v0.9 main_gui.py:695-775 (sync_parameters)
        """
        logger.debug("Synchronizing parameters to integrator")

        if not self.integrator:
            logger.warning("Integrator not initialized, skipping parameter sync")
            return

        try:
            # Update parameter manager from GUI attributes
            # v1.0 uses plain Python attributes instead of tkinter variables
            if hasattr(self, 'param_manager') and self.param_manager:
                # Update detection parameters
                self.param_manager.current_params['search_window_x'] = self.search_window_x
                self.param_manager.current_params['search_window_y'] = self.search_window_y
                self.param_manager.current_params['noise_threshold'] = self.noise_threshold

                # Update centroid parameters if available
                if hasattr(self, 'centroid_window_x_ppm'):
                    self.param_manager.current_params['centroid_window_x_ppm'] = self.centroid_window_x_ppm
                if hasattr(self, 'centroid_window_y_ppm'):
                    self.param_manager.current_params['centroid_window_y_ppm'] = self.centroid_window_y_ppm
                if hasattr(self, 'centroid_noise_multiplier'):
                    self.param_manager.current_params['centroid_noise_multiplier'] = self.centroid_noise_multiplier

                # Update PS2D parameters
                if hasattr(self, 'fix_linewidths'):
                    self.param_manager.current_params['fix_linewidths'] = self.fix_linewidths
                if hasattr(self, 'fix_positions'):
                    self.param_manager.current_params['fix_positions'] = self.fix_positions
                if hasattr(self, 'use_parallel_processing'):
                    self.param_manager.current_params['use_parallel_processing'] = self.use_parallel_processing
                if hasattr(self, 'use_adaptive_optimization'):
                    self.param_manager.current_params['use_adaptive_optimization'] = self.use_adaptive_optimization

                # Update custom linewidth parameters
                if hasattr(self, 'use_custom_linewidths') and self.use_custom_linewidths:
                    # Only pass custom linewidths if the feature is enabled
                    if hasattr(self, 'lw_lorentz_1h'):
                        self.param_manager.current_params['lw_lorentz_1h'] = self.lw_lorentz_1h
                    if hasattr(self, 'lw_gauss_1h'):
                        self.param_manager.current_params['lw_gauss_1h'] = self.lw_gauss_1h
                    if hasattr(self, 'lw_lorentz_15n'):
                        self.param_manager.current_params['lw_lorentz_15n'] = self.lw_lorentz_15n
                    if hasattr(self, 'lw_gauss_15n'):
                        self.param_manager.current_params['lw_gauss_15n'] = self.lw_gauss_15n
                else:
                    # Clear custom linewidths when feature is disabled (use defaults)
                    self.param_manager.current_params['lw_lorentz_1h'] = None
                    self.param_manager.current_params['lw_gauss_1h'] = None
                    self.param_manager.current_params['lw_lorentz_15n'] = None
                    self.param_manager.current_params['lw_gauss_15n'] = None

                # Handle simplified mode flag (v0.9 main_gui.py:705-716)
                if hasattr(self, 'use_simplified_parameters'):
                    self.param_manager.use_simplified_mode = self.use_simplified_parameters
                    if self.use_simplified_parameters:
                        self.param_manager.update_simplified_parameters(
                            window_scale=getattr(self, 'simplified_window_scale', 1.0),
                            quality_target=getattr(self, 'simplified_quality_target', 'standard'),
                            max_peaks_fit=getattr(self, 'simplified_max_peaks_fit', 50),
                            max_iterations=getattr(self, 'simplified_max_iterations', 1000),
                            noise_estimation_method=getattr(self, 'simplified_noise_method', 'auto'),
                            baseline_method=getattr(self, 'simplified_baseline_method', 'auto')
                        )

                # Validate parameters
                validation_errors = self.param_manager.validate_all_parameters()
                if validation_errors:
                    logger.warning(f"Parameter validation warnings: {', '.join(validation_errors[:3])}")

                # Get effective parameters based on mode
                integrator_params = self.param_manager.get_effective_parameters()

                # Apply parameters to integrator (v0.9 main_gui.py:727-735)
                detection_params = integrator_params['detection_params']
                self.integrator.set_search_window(
                    detection_params['search_window_x'],
                    detection_params['search_window_y']
                )
                self.integrator.set_threshold_multiplier(detection_params['noise_threshold'])

                # Update fitting parameters
                if hasattr(self.integrator, 'fitting_parameters'):
                    self.integrator.fitting_parameters.update(integrator_params['fitting_params'])

                # Update GUI params in integrator
                if hasattr(self.integrator, 'gui_params'):
                    self.integrator.gui_params = integrator_params['gui_params']

                logger.debug(f"Parameter manager updated with {len(integrator_params)} parameter groups")

            else:
                # Direct sync without parameter manager (fallback)
                self.integrator.set_search_window(self.search_window_x, self.search_window_y)
                self.integrator.set_threshold_multiplier(self.noise_threshold)

            # Sync S/N threshold and workflow mode (v0.9 main_gui.py:739-760)
            # Read from spinbox widgets directly (not from plain attributes)
            if hasattr(self, 'sn_threshold_spin') and self.sn_threshold_spin is not None:
                try:
                    sn_threshold_value = float(self.sn_threshold_spin.value())
                    self.integrator.sn_threshold = sn_threshold_value
                except (ValueError, AttributeError):
                    pass

            if hasattr(self, 'expected_peaks_spin') and self.expected_peaks_spin is not None:
                try:
                    expected_count_value = int(self.expected_peaks_spin.value())
                    self.integrator.expected_peak_count = expected_count_value
                except (ValueError, AttributeError):
                    pass

            # Sync PS2D expert parameters from spinboxes to global config
            ps2d_config = get_ps2d_config()
            if hasattr(self, 'expert_radF1_spin') and self.expert_radF1_spin is not None:
                try:
                    ps2d_config.radF1 = float(self.expert_radF1_spin.value())
                except (ValueError, AttributeError):
                    pass
            if hasattr(self, 'expert_radF2_spin') and self.expert_radF2_spin is not None:
                try:
                    ps2d_config.radF2 = float(self.expert_radF2_spin.value())
                except (ValueError, AttributeError):
                    pass
            if hasattr(self, 'expert_overlap_x_spin') and self.expert_overlap_x_spin is not None:
                try:
                    ps2d_config.overlap_threshold_x = float(self.expert_overlap_x_spin.value())
                except (ValueError, AttributeError):
                    pass
            if hasattr(self, 'expert_overlap_y_spin') and self.expert_overlap_y_spin is not None:
                try:
                    ps2d_config.overlap_threshold_y = float(self.expert_overlap_y_spin.value())
                except (ValueError, AttributeError):
                    pass
            if hasattr(self, 'expert_max_iter_spin') and self.expert_max_iter_spin is not None:
                try:
                    ps2d_config.max_iterations = int(self.expert_max_iter_spin.value())
                except (ValueError, AttributeError):
                    pass

            # Sync auto-add dummy peaks setting
            if hasattr(self, 'auto_add_dummy_peaks'):
                self.integrator.auto_add_dummy_peaks = self.auto_add_dummy_peaks

            # Sync workflow mode
            if hasattr(self, 'workflow_mode'):
                workflow_mode = self.workflow_mode
                if workflow_mode == "sn_threshold":
                    self.integrator.set_processing_mode('sn_native')
                elif workflow_mode == "peak_list":
                    self.integrator.set_processing_mode('in_place')
                else:
                    self.integrator.set_processing_mode('full_detection')

            logger.debug("Parameters synchronized to integrator successfully")

        except Exception as e:
            logger.error(f"Failed to sync parameters to integrator: {e}")
            # Continue - don't fail the whole operation for parameter sync issues

    # ===== Control Center Button Handlers =====

    def on_detect_peaks(self):
        """Handle Detect Peaks button click.

        This initiates peak detection on the loaded spectrum using
        the current detection parameters in a background thread.

        Process:
        1. Validate NMR data is loaded
        2. Create PeakDetectionWorker thread
        3. Show progress dialog
        4. Run detection in background
        5. Update GUI with results (peaks, plots, statistics)
        """
        logger.info("Detect Peaks button clicked")

        # CRITICAL: Sync GUI parameters BEFORE detection
        # This ensures Expert Mode settings (search windows, noise threshold, etc.) are passed
        self.sync_parameters_to_integrator()
        logger.info(f"Parameters synced for detection: search_window=({self.search_window_x}, {self.search_window_y}), noise_threshold={self.noise_threshold}")

        # Validate NMR data is loaded
        if not hasattr(self.integrator, 'nmr_data') or self.integrator.nmr_data is None:
            from lunaNMR.gui.components import show_error
            show_error(self, "No Data", "Please load an NMR spectrum before detecting peaks.")
            return

        # For peak_list mode, validate peak list is loaded
        if self.workflow_mode == 'peak_list':
            if not hasattr(self.integrator, 'peak_list') or self.integrator.peak_list is None:
                from lunaNMR.gui.components import show_error
                show_error(self, "No Peak List", "Peak list mode requires a reference peak list. Please load one or switch to S/N threshold mode.")
                return

        # Import progress dialog
        from lunaNMR.gui.components import ProgressDialog

        # Create progress dialog
        progress_dialog = ProgressDialog(
            parent=self,
            title="Peak Detection"
        )

        # Create and configure worker
        self.detection_worker = PeakDetectionWorker(self.integrator, self.workflow_mode)

        # Connect worker signals to dialog and handlers
        self.detection_worker.progress_updated.connect(progress_dialog.update_progress)
        self.detection_worker.detection_complete.connect(
            lambda peaks: self.on_detection_complete(peaks, progress_dialog)
        )
        self.detection_worker.detection_failed.connect(
            lambda error: self.on_detection_failed(error, progress_dialog)
        )
        # Do not auto-close - let user review log and close manually
        self.detection_worker.finished.connect(
            lambda: progress_dialog.complete("Detection complete - click Close to continue")
        )

        # Start worker and show progress
        logger.info(f"Starting peak detection (mode={self.workflow_mode})")
        self.update_status("Detecting peaks...")
        self.detection_worker.start()
        progress_dialog.exec()

    @Slot(list)
    def on_detection_complete(self, detected_peaks, progress_dialog):
        """Handle successful peak detection completion.

        Args:
            detected_peaks: List of detected peak dictionaries
            progress_dialog: Progress dialog to close
        """
        try:
            logger.info(f"Detection complete: {len(detected_peaks)} peaks detected")

            # Convert fitted_peaks to DataFrame for peak_list
            # This ensures consistency between detection and fitting steps
            if self.integrator.fitted_peaks:
                detected_df = self._convert_fitted_to_dataframe(self.integrator.fitted_peaks)
                if detected_df is not None:
                    # In peak_list mode, backup original reference list
                    if self.workflow_mode == 'peak_list':
                        if not hasattr(self.integrator, 'peak_list_original'):
                            self.integrator.peak_list_original = self.integrator.peak_list.copy()

                    # Set peak_list to detected positions
                    self.integrator.peak_list = detected_df

            # Update peak navigator with detected peaks (if initialized)
            if self.peak_navigator is not None:
                self.peak_navigator.load_detected_peaks(self.integrator.fitted_peaks)
            else:
                logger.warning("Peak navigator not initialized - detected peaks not displayed")

            # Update spectrum plot with detected peaks
            self.update_spectrum_plot()

            # Update statistics
            stats = self.integrator.get_detection_statistics()
            status_msg = f"Detection complete: {stats['detected_peaks']}/{stats.get('total_peaks', len(detected_peaks))} peaks"
            if 'detection_rate' in stats:
                status_msg += f" ({stats['detection_rate']:.1f}%)"
            self.update_status(status_msg)

            # Update statistics panel
            self.update_statistics_panel()

            # Update button colors to indicate detection is complete
            self._update_button_colors_for_detection_complete()

            logger.info("Detection results updated in GUI")

        except Exception as e:
            logger.error(f"Error handling detection results: {e}", exc_info=True)
            from lunaNMR.gui.components import show_error
            show_error(self, "Update Error", f"Detection succeeded but GUI update failed: {e}")

    @Slot(str)
    def on_detection_failed(self, error_message, progress_dialog):
        """Handle peak detection failure.

        Args:
            error_message: Error description
            progress_dialog: Progress dialog to close
        """
        logger.error(f"Peak detection failed: {error_message}")
        self.update_status(f"Detection failed: {error_message}")

        from lunaNMR.gui.components import show_error
        show_error(self, "Detection Failed", error_message)

    # ===== Peak Detection Helper Methods =====

    def _convert_fitted_to_dataframe(self, fitted_peaks):
        """Convert fitted_peaks list to peak_list DataFrame.

        This ensures fitting uses detected positions, not reference positions.
        Called after detection to replace reference peak_list with detected peaks.

        Args:
            fitted_peaks: List of peak dictionaries with keys:
                - assignment: Peak label (e.g., 'G10')
                - ppm_x: X-axis position in ppm
                - ppm_y: Y-axis position in ppm
                - intensity: Peak intensity/volume
                - snr: Signal-to-noise ratio
                - detected: Boolean detection flag
                - detection_quality: Quality descriptor

        Returns:
            pandas.DataFrame with columns: Assignment, Position_X, Position_Y,
            Height, Intensity, SNR, Detected, Detection_Quality
            Returns None if fitted_peaks is empty.
        """
        import pandas as pd

        if not fitted_peaks:
            return None

        # Extract data from fitted_peaks
        data = {
            'Assignment': [],
            'Position_X': [],
            'Position_Y': [],
            'Height': [],
            'Intensity': [],
            'SNR': [],
            'Detected': [],
            'Detection_Quality': []
        }

        for peak in fitted_peaks:
            data['Assignment'].append(peak.get('assignment', 'Unknown'))
            # Check multiple key variants: fitting uses peak_x/peak_y, detection uses ppm_x/ppm_y
            data['Position_X'].append(peak.get('ppm_x', peak.get('peak_x', 0.0)))
            data['Position_Y'].append(peak.get('ppm_y', peak.get('peak_y', 0.0)))
            data['Height'].append(peak.get('intensity', 0.0))
            data['Intensity'].append(peak.get('intensity', 0.0))
            data['SNR'].append(peak.get('snr', 0.0))
            data['Detected'].append(peak.get('detected', False))
            data['Detection_Quality'].append(peak.get('detection_quality', 'Unknown'))

        df = pd.DataFrame(data)
        logger.debug(f"Converted {len(fitted_peaks)} fitted_peaks to DataFrame")
        return df

    def update_statistics_panel(self):
        """Update the Statistics panel with detection/fitting results.

        Retrieves statistics from the integrator and updates the statistics
        display labels in the left panel.
        """
        try:
            stats = self.integrator.get_detection_statistics()

            # Format statistics text
            stats_text = []
            if 'detected_peaks' in stats:
                stats_text.append(f"Detected: {stats['detected_peaks']}")
            if 'total_peaks' in stats:
                stats_text.append(f"Total: {stats['total_peaks']}")
            if 'detection_rate' in stats:
                stats_text.append(f"Rate: {stats['detection_rate']:.1f}%")
            if 'failed_peaks' in stats:
                stats_text.append(f"Failed: {stats['failed_peaks']}")

            # Update statistics display (if it exists)
            if hasattr(self, 'statistics_display'):
                self.statistics_display.setText('\n'.join(stats_text))

            logger.debug(f"Statistics panel updated: {stats}")

        except Exception as e:
            logger.warning(f"Could not update statistics panel: {e}")

    # ===== Fitting Button Handlers =====

    def on_fit_spectrum(self):
        """Handle Fit Spectrum button click.

        This initiates Voigt fitting on all detected/reference peaks using
        SingleSpectrumProcessor with PS2D multi-peak fitting support.

        This follows the exact logic from lunaNMR_v0o9.

        Process:
        1. Validate peak_list is loaded
        2. Create VoigtFittingWorker thread with SingleSpectrumProcessor
        3. Show progress dialog
        4. Run fitting in background
        5. Update GUI with results (peaks, plots, statistics)
        """
        logger.info("Fit Spectrum button clicked")

        # CRITICAL: Sync GUI parameters to param_manager BEFORE fitting
        # This ensures Expert Mode settings (fix_positions, fix_linewidths, etc.) are passed
        self.sync_parameters_to_integrator()
        logger.info(f"Parameters synced: fix_positions={self.fix_positions}, fix_linewidths={self.fix_linewidths}")

        # Validate peak_list is loaded (required for fitting)
        if not hasattr(self.integrator, 'peak_list') or self.integrator.peak_list is None:
            from lunaNMR.gui.components import show_error
            show_error(self, "No Peak List",
                      "Please load a peak list or run peak detection first.")
            return

        if len(self.integrator.peak_list) == 0:
            from lunaNMR.gui.components import show_error
            show_error(self, "Empty Peak List",
                      "Peak list is empty. Please load peaks or run detection first.")
            return

        # Import progress dialog
        from lunaNMR.gui.components import ProgressDialog

        # Create progress dialog
        progress_dialog = ProgressDialog(
            parent=self,
            title="Voigt Fitting"
        )

        # Build processing options from GUI parameters (v0o9 pattern)
        processing_options = {
            'use_parallel': self.use_parallel_processing,
            'use_global_optimization': getattr(self, 'use_global_optimization', False),
            'use_voigt_fitting': self.use_voigt_fitting
        }

        logger.info(f"Processing options: {processing_options}")

        # Create and configure worker with SingleSpectrumProcessor
        self.fitting_worker = VoigtFittingWorker(
            self.integrator,
            self.param_manager,
            processing_options
        )

        # Connect worker signals to dialog and handlers
        # Track fitting phase for statistics counting
        # PASS 1 = calibration (learn linewidths) - don't count
        # OPTIMIZE = grid search - don't count
        # PASS 1-bis = re-fit isolated with optimal params - don't count (intermediate)
        # PASS 2 = final multi-peak cluster fitting - COUNT
        fitting_phase = [None]  # Track current phase: None, 'pass1bis', 'pass2'

        def on_fitting_progress(p, t, log):
            """Handle progress update and add to detailed log."""
            progress_dialog.update_progress(p, f"{t}\n{log}" if log else t)

            # Track phase transitions (don't reset timer - keep cumulative)
            if t:
                if "PASS 2" in t and fitting_phase[0] != 'pass2':
                    # Entering PASS 2 - reset counters, start counting
                    fitting_phase[0] = 'pass2'
                    progress_dialog.reset_statistics(reset_timer=False)
                elif "PASS 1-bis" in t and fitting_phase[0] is None:
                    # Entering PASS 1-bis - track but don't count yet
                    fitting_phase[0] = 'pass1bis'

            if log:
                # Only count in statistics during PASS 2
                count_this = (fitting_phase[0] == 'pass2')
                progress_dialog.add_log_message(log, is_error=False, count_in_stats=count_this)

        self.fitting_worker.progress_updated.connect(on_fitting_progress)
        self.fitting_worker.fitting_complete.connect(
            lambda summary: self.on_fitting_complete(summary, progress_dialog)
        )
        self.fitting_worker.fitting_failed.connect(
            lambda error: self.on_fitting_failed(error, progress_dialog)
        )
        # Do not auto-close - let user review log and close manually
        self.fitting_worker.finished.connect(
            lambda: progress_dialog.complete("Fitting complete - click Close to continue")
        )

        # Start worker and show progress
        logger.info(f"Starting Voigt fitting for {len(self.integrator.peak_list)} peaks")
        self.update_status(f"Fitting {len(self.integrator.peak_list)} peaks...")
        self.fitting_worker.start()
        progress_dialog.exec()

    @Slot(dict)
    def on_fitting_complete(self, summary, progress_dialog):
        """Handle successful fitting completion (exact v0o9 logic).

        Args:
            summary: Dictionary with keys:
                - total_peaks: int
                - successful_peaks: int
                - failed_peaks: int
                - average_quality: float
                - fitted_results: list of fit results
            progress_dialog: Progress dialog to close
        """
        try:
            logger.info(f"Fitting complete: {summary['successful_peaks']}/{summary['total_peaks']} successful")

            # Store results (v0o9 pattern line 4204)
            if 'fitted_results' in summary and summary['fitted_results']:
                self.last_fitting_results = summary['fitted_results']

                # Standardize result keys with bidirectional mapping (v0o9 lines 4207-4233)
                standardized_results = []
                for result in summary['fitted_results']:
                    # Create standardized structure matching v0o9 detection workflow
                    # Check multiple key variants: fitting uses peak_x/peak_y, detection uses ppm_x/ppm_y
                    standardized_result = {
                        'assignment': result.get('Assignment', result.get('assignment', '')),
                        'ppm_x': float(result.get('Position_X', result.get('ppm_x', result.get('peak_x', 0)))),
                        'ppm_y': float(result.get('Position_Y', result.get('ppm_y', result.get('peak_y', 0)))),
                        'detected': True,
                        'fitted': True
                    }
                    # Preserve additional fitting data
                    for key, value in result.items():
                        if key not in ['Assignment', 'Position_X', 'Position_Y']:
                            standardized_result[key] = value

                    # Bidirectional quality field mapping for visualization compatibility (v0o9 lines 4222-4228)
                    if 'fitting_quality' in standardized_result:
                        standardized_result['Quality'] = standardized_result['fitting_quality']
                        standardized_result['quality'] = standardized_result['fitting_quality']
                    elif 'Quality' in standardized_result:
                        standardized_result['fitting_quality'] = standardized_result['Quality']
                        standardized_result['quality'] = standardized_result['Quality']

                    standardized_results.append(standardized_result)

                # Update integrator.fitted_peaks (v0o9 line 4232)
                self.integrator.fitted_peaks = standardized_results
                logger.info(f"Standardized {len(standardized_results)} fitting results for series integration compatibility")

                # Update peak navigator with heights (v0o9 lines 4236-4286)
                if hasattr(self, 'peak_navigator') and self.peak_navigator:
                    # Standardize result keys for Peak Navigator compatibility (bidirectional)
                    standardized_for_navigator = []
                    for result in summary['fitted_results']:
                        nav_result = result.copy()

                        # Bidirectional assignment key mapping
                        if 'Assignment' in nav_result and 'assignment' not in nav_result:
                            nav_result['assignment'] = nav_result['Assignment']
                        elif 'assignment' in nav_result and 'Assignment' not in nav_result:
                            nav_result['Assignment'] = nav_result['assignment']

                        # Bidirectional height/amplitude mapping
                        if 'amplitude' in nav_result and 'height' not in nav_result:
                            nav_result['height'] = nav_result['amplitude']
                        elif 'height' in nav_result and 'amplitude' not in nav_result:
                            nav_result['amplitude'] = nav_result['height']

                        # Bidirectional r_squared mapping
                        if 'avg_r_squared' in nav_result and 'r_squared' not in nav_result:
                            nav_result['r_squared'] = nav_result['avg_r_squared']
                        elif 'r_squared' in nav_result and 'avg_r_squared' not in nav_result:
                            nav_result['avg_r_squared'] = nav_result['r_squared']

                        # Bidirectional quality mapping (for color-coded visualization)
                        if 'fitting_quality' in nav_result and 'Quality' not in nav_result:
                            nav_result['Quality'] = nav_result['fitting_quality']
                        elif 'Quality' in nav_result and 'fitting_quality' not in nav_result:
                            nav_result['fitting_quality'] = nav_result['Quality']
                        # Also add lowercase 'quality' for maximum compatibility
                        if 'fitting_quality' in nav_result and 'quality' not in nav_result:
                            nav_result['quality'] = nav_result['fitting_quality']
                        elif 'Quality' in nav_result and 'quality' not in nav_result:
                            nav_result['quality'] = nav_result['Quality']

                        # Ensure fitted flag is set
                        if 'fitted' not in nav_result:
                            nav_result['fitted'] = (nav_result.get('success', False) or
                                                     nav_result.get('height', 0) > 0 or
                                                     nav_result.get('amplitude', 0) > 0)

                        standardized_for_navigator.append(nav_result)

                    logger.info(f"Standardized {len(standardized_for_navigator)} results for Peak Navigator")
                    if standardized_for_navigator:
                        logger.debug(f"Sample result keys: {list(standardized_for_navigator[0].keys())}")
                        logger.debug(f"Sample assignment: {standardized_for_navigator[0].get('assignment', 'N/A')}")
                        logger.debug(f"Sample height: {standardized_for_navigator[0].get('height', 'N/A')}")
                        logger.debug(f"Sample r_squared: {standardized_for_navigator[0].get('r_squared', 'N/A')}")
                        logger.debug(f"Sample fitted: {standardized_for_navigator[0].get('fitted', 'N/A')}")
                    # Use update_heights_from_results (v0o9 pattern line 4286)
                    self.peak_navigator.update_heights_from_results(standardized_for_navigator)

            # Print detailed summary (v0o9 lines 4289-4295)
            logger.info(f"\nSingle Spectrum Processing Summary:")
            logger.info(f"   Total peaks: {summary['total_peaks']}")
            logger.info(f"   Successful: {summary['successful_peaks']} ({summary.get('success_rate', 0):.1f}%)")
            if 'average_r_squared' in summary:
                logger.info(f"   Average R²: {summary['average_r_squared']:.3f}")

            if 'quality_distribution' in summary:
                quality_dist = summary['quality_distribution']
                logger.info(f"   Quality distribution: {quality_dist['excellent']} excellent, {quality_dist['good']} good, {quality_dist['poor']} poor")

            # Update main spectrum plot to show color-coded peaks (v0o9 line 4298)
            self.update_spectrum_plot()

            # Update statistics panel
            self.update_statistics_panel()

            # Update Voigt Analysis plots with first successful fit (v0o9 lines 3886-3892)
            if self.integrator.fitted_peaks and len(self.integrator.fitted_peaks) > 0:
                for peak_result in self.integrator.fitted_peaks:
                    quality = peak_result.get('quality', 0)
                    if isinstance(quality, str):
                        try:
                            quality = float(quality)
                        except (ValueError, TypeError):
                            quality = 0
                    if quality > 0:
                        self.current_voigt_result = peak_result
                        # Update 2D Voigt Analysis plot
                        if hasattr(self, 'voigt_2d_plotter') and self.voigt_2d_plotter:
                            self.voigt_2d_plotter.plot_voigt_analysis(peak_result)
                            logger.debug(f"Updated 2D Voigt plot with peak: {peak_result.get('peak_id', 'unknown')}")
                        # Update 3D Voigt Analysis plot
                        if hasattr(self, 'voigt_3d_plotter') and self.voigt_3d_plotter:
                            self.voigt_3d_plotter.plot_voigt_analysis_3d(peak_result)
                            logger.debug(f"Updated 3D Voigt plot with peak: {peak_result.get('peak_id', 'unknown')}")
                        break

            # Update status bar with summary
            status_msg = f"Fitting complete: {summary['successful_peaks']}/{summary['total_peaks']} peaks"
            if 'average_quality' in summary and summary['average_quality'] > 0:
                status_msg += f" (avg R²: {summary['average_quality']:.3f})"
            self.update_status(status_msg)

            # Save fitting results to CSV file (matching series integration output format)
            if 'fitted_results' in summary and summary['fitted_results']:
                output_file = self.save_fitting_results_to_csv(summary['fitted_results'])
                if output_file:
                    status_msg += f" - Results saved"
                    self.update_status(status_msg)

            logger.info("Fitting results updated in GUI")

        except Exception as e:
            logger.error(f"Error handling fitting results: {e}", exc_info=True)
            from lunaNMR.gui.components import show_error
            show_error(self, "Update Error", f"Fitting succeeded but GUI update failed: {e}")

    @Slot(str)
    def on_fitting_failed(self, error_message, progress_dialog):
        """Handle fitting failure.

        Args:
            error_message: Error description
            progress_dialog: Progress dialog to close
        """
        logger.error(f"Voigt fitting failed: {error_message}")
        self.update_status(f"Fitting failed: {error_message}")

        from lunaNMR.gui.components import show_error
        show_error(self, "Fitting Failed", error_message)

    def save_fitting_results_to_csv(self, fitted_results: list) -> str:
        """Save single spectrum fitting results to CSV file.

        Creates a CSV file with columns matching series integration output:
        Assignment, Position_X, Position_Y, Height, Volume, LW_X, LW_Y, R_Squared, Quality

        Formatting matches series integration:
        - Position_X, Position_Y: 3 decimal places
        - LW_X, LW_Y: 3 decimal places
        - R_Squared: 2 decimal places
        - Height, Volume: scientific notation

        Args:
            fitted_results: List of fitting result dictionaries

        Returns:
            Path to saved CSV file, or empty string if save failed
        """
        import os
        import pandas as pd
        import numpy as np

        if not fitted_results:
            logger.warning("No fitting results to save")
            return ""

        # Derive output filename from current NMR file
        if not self.current_nmr_file:
            logger.warning("No current NMR file - cannot determine output path")
            return ""

        # Create output path: same directory as NMR file, with _fitting_results.csv suffix
        nmr_dir = os.path.dirname(self.current_nmr_file)
        nmr_basename = os.path.basename(self.current_nmr_file)
        # Remove extension and add suffix
        base_name = os.path.splitext(nmr_basename)[0]
        output_file = os.path.join(nmr_dir, f"{base_name}_fitting_results.csv")

        try:
            # Build DataFrame from results
            rows = []
            for i, result in enumerate(fitted_results):
                # Extract values with fallbacks for different key naming conventions
                assignment = result.get('assignment', result.get('Assignment', f'Peak_{i+1}'))

                # Position: try multiple key names
                pos_x = (result.get('pos_f2') or result.get('ppm_x') or
                        result.get('center_x') or result.get('Position_X') or
                        result.get('peak_x') or 0.0)
                pos_y = (result.get('pos_f1') or result.get('ppm_y') or
                        result.get('center_y') or result.get('Position_Y') or
                        result.get('peak_y') or 0.0)

                # Height and Volume
                height = result.get('height', result.get('Height', result.get('amplitude', 0.0)))
                volume = result.get('volume', result.get('Volume', result.get('intensity', 0.0)))

                # Linewidths - sum of Gaussian and Lorentzian FWHM components
                # NOTE: sigma/gamma in x_fit/y_fit are ALREADY FWHM values (not true sigma/gamma)
                # So we just add them directly (matching Peak Parameters display)
                lw_x = 0.0
                lw_y = 0.0

                # Try to get from x_fit/y_fit dicts first
                x_fit = result.get('x_fit', {})
                y_fit = result.get('y_fit', {})

                if x_fit:
                    x_sigma = x_fit.get('sigma', 0)  # Actually Gaussian FWHM
                    x_gamma = x_fit.get('gamma', 0)  # Actually Lorentzian FWHM
                    if x_sigma or x_gamma:
                        lw_x = x_sigma + x_gamma  # Simple sum of FWHM components

                if y_fit:
                    y_sigma = y_fit.get('sigma', 0)  # Actually Gaussian FWHM
                    y_gamma = y_fit.get('gamma', 0)  # Actually Lorentzian FWHM
                    if y_sigma or y_gamma:
                        lw_y = y_sigma + y_gamma  # Simple sum of FWHM components

                # Fallback to direct lw_x/lw_y if x_fit/y_fit not available
                if lw_x == 0.0:
                    lw_x = result.get('lw_x', 0.0)
                if lw_y == 0.0:
                    lw_y = result.get('lw_y', 0.0)

                # Final fallback: compute from sigma_x/gamma_x (also FWHM values)
                if lw_x == 0.0:
                    sigma_x = result.get('sigma_x', 0)  # Gaussian FWHM
                    gamma_x = result.get('gamma_x', 0)  # Lorentzian FWHM
                    if sigma_x or gamma_x:
                        lw_x = sigma_x + gamma_x

                if lw_y == 0.0:
                    sigma_y = result.get('sigma_y', 0)  # Gaussian FWHM
                    gamma_y = result.get('gamma_y', 0)  # Lorentzian FWHM
                    if sigma_y or gamma_y:
                        lw_y = sigma_y + gamma_y

                # R-squared
                r_squared = (result.get('r_squared') or result.get('R_Squared') or
                            result.get('avg_r_squared') or 0.0)

                # Quality
                quality = (result.get('quality') or result.get('Quality') or
                          result.get('fitting_quality') or 'Unknown')

                rows.append({
                    'Assignment': assignment,
                    'Position_X': float(pos_x) if pos_x else 0.0,
                    'Position_Y': float(pos_y) if pos_y else 0.0,
                    'Height': float(height) if height else 0.0,
                    'Volume': float(volume) if volume else 0.0,
                    'LW_X': float(lw_x) if lw_x else 0.0,
                    'LW_Y': float(lw_y) if lw_y else 0.0,
                    'R_Squared': float(r_squared) if r_squared else 0.0,
                    'Quality': quality
                })

            df = pd.DataFrame(rows)

            # Format columns with specific decimal places (matching series integration)
            df_formatted = df.copy()
            df_formatted['Position_X'] = df['Position_X'].apply(lambda x: f'{x:.3f}')
            df_formatted['Position_Y'] = df['Position_Y'].apply(lambda x: f'{x:.3f}')
            df_formatted['Height'] = df['Height'].apply(lambda x: f'{x:.2e}')
            df_formatted['Volume'] = df['Volume'].apply(lambda x: f'{x:.2e}')
            df_formatted['LW_X'] = df['LW_X'].apply(lambda x: f'{x:.3f}')
            df_formatted['LW_Y'] = df['LW_Y'].apply(lambda x: f'{x:.3f}')
            df_formatted['R_Squared'] = df['R_Squared'].apply(lambda x: f'{x:.2f}')

            df_formatted.to_csv(output_file, index=False)

            logger.info(f"Saved fitting results to: {output_file}")
            print(f"💾 Saved fitting results: {len(rows)} peaks → {output_file}")

            return output_file

        except Exception as e:
            logger.error(f"Failed to save fitting results: {e}")
            return ""

    def on_fit_series(self):
        """Handle Fit Series button click.

        This opens the series integration dialog for batch processing
        multiple spectra. Exact copy of v0.9 start_series_integration flow.

        Based on v0.9 main_gui.py lines 4452-4510.
        """
        logger.info("Fit Series button clicked")

        # CRITICAL: Sync GUI parameters BEFORE series integration
        # This ensures Expert Mode settings are up-to-date for the dialog to read
        self.sync_parameters_to_integrator()
        logger.info(f"Parameters synced: fix_positions={self.fix_positions}, fix_linewidths={self.fix_linewidths}")

        # Import the dialog
        from lunaNMR.gui.dialogs import SeriesIntegrationDialog

        # Validate NMR folder is set
        if not hasattr(self, 'current_nmr_folder') or not self.current_nmr_folder:
            QMessageBox.warning(
                self,
                "No NMR Folder",
                "No NMR data folder selected.\n\n"
                "Please load NMR data first using 'Load Data'."
            )
            return

        # Create and show the dialog
        dialog = SeriesIntegrationDialog(parent=self, main_window=self)

        # Connect processing complete signal
        dialog.processing_complete.connect(self._on_series_processing_complete)

        # Show dialog (non-blocking)
        dialog.exec()

    def _on_series_processing_complete(self, batch_results):
        """Handle series processing completion.

        Updates GUI components with batch results including:
        - Statistics display
        - Results browser availability
        - Voigt visualization

        Based on v0.9 _complete_new_multi_spectrum_processing() (main_gui.py:4717-4823)

        Args:
            batch_results: BatchResults object or dict containing processing results
        """
        logger.info("Series processing complete, updating GUI")

        if not batch_results:
            logger.warning("No batch results returned")
            return

        # Store results
        self.batch_results = batch_results

        # Get summary - handle both dict format and BatchResults object
        if hasattr(batch_results, 'get_summary'):
            summary = batch_results.get_summary()
        elif isinstance(batch_results, dict):
            summary = batch_results.get('summary', {})
        else:
            summary = {}

        successful = summary.get('successful', 0)
        total = summary.get('total_spectra', 0)
        success_rate = summary.get('success_rate', 0)

        # Update status
        status_msg = f"Series processing complete: {successful}/{total} successful ({success_rate:.1f}%)"
        self.update_status(status_msg)

        # Switch to 3D Voigt Analysis tab (v0.9 line 4791-4792)
        if hasattr(self, 'viz_tabs'):
            for i in range(self.viz_tabs.count()):
                if self.viz_tabs.tabText(i) == "3D Voigt Analysis":
                    self.viz_tabs.setCurrentIndex(i)
                    break

        # Show completion message
        QMessageBox.information(
            self,
            "Series Processing Complete",
            f"Series integration complete!\n\n"
            f"Processed: {total} spectra\n"
            f"Successful: {successful}\n"
            f"Failed: {total - successful}\n"
            f"Success rate: {success_rate:.1f}%\n\n"
            "Use View > Results Browser to analyze results."
        )

        logger.info(f"Series processing: {successful}/{total} successful")

    def toggle_peak_edition(self):
        """Toggle peak editing mode on/off (v0o9 line 5558 logic)."""
        # Check current state by visibility of controls
        is_visible = self.peak_edition_widget.isVisible()

        if not is_visible:
            # Enable edit mode
            self.peak_edit_mode = True
            self.edit_connection_id = self.spectrum_plotter.canvas.mpl_connect('button_press_event', self.on_peak_edit_click)
            self.spectrum_plotter.canvas.setCursor(Qt.CrossCursor)

            # Show editing controls
            self.peak_edition_widget.setVisible(True)

            # Update button text
            self.peak_edit_toggle_button.setText("Peak Edition ▲")

            logger.info("Peak editing mode enabled")
        else:
            # Disable edit mode
            self.peak_edit_mode = False
            if self.edit_connection_id is not None:
                self.spectrum_plotter.canvas.mpl_disconnect(self.edit_connection_id)
                self.edit_connection_id = None
            self.spectrum_plotter.canvas.setCursor(Qt.ArrowCursor)
            self.selected_peak_info = None
            if hasattr(self, 'selected_peak_label'):
                self.selected_peak_label.setText("No peak selected")

            # Hide editing controls
            self.peak_edition_widget.setVisible(False)

            # Update button text
            self.peak_edit_toggle_button.setText("Peak Edition ▼")

            logger.info("Peak editing mode disabled")

        # Update status display
        self.update_edit_mode_status()

    def on_edit_reference_changed(self, state):
        """Handle Edit Reference peaks checkbox state change.

        Args:
            state: int - Qt.CheckState value (0=Unchecked, 2=Checked)
        """
        self.edit_reference_peaks = (state == 2)  # Qt.Checked = 2
        logger.debug(f"Edit Reference peaks: {self.edit_reference_peaks}")
        self.update_edit_mode_status()

    def on_edit_detected_changed(self, state):
        """Handle Edit Detected peaks checkbox state change.

        Args:
            state: int - Qt.CheckState value (0=Unchecked, 2=Checked)
        """
        self.edit_detected_peaks = (state == 2)  # Qt.Checked = 2
        logger.debug(f"Edit Detected peaks: {self.edit_detected_peaks}")
        self.update_edit_mode_status()

    def on_deletion_mode_changed(self, state):
        """Handle Delete Mode checkbox state change.

        Args:
            state: int - Qt.CheckState value (0=Unchecked, 2=Checked)
        """
        self.peak_deletion_mode = (state == 2)  # Qt.Checked = 2
        logger.debug(f"Delete Mode: {self.peak_deletion_mode}")

        # If deletion mode is enabled, disable addition mode
        if self.peak_deletion_mode:
            self.add_mode_checkbox.setChecked(False)
            self.peak_addition_mode = False

        self.update_edit_mode_status()

    def on_addition_mode_changed(self, state):
        """Handle Add Mode checkbox state change.

        Args:
            state: int - Qt.CheckState value (0=Unchecked, 2=Checked)
        """
        self.peak_addition_mode = (state == 2)  # Qt.Checked = 2
        logger.debug(f"Add Mode: {self.peak_addition_mode}")

        # If addition mode is enabled, disable deletion mode
        if self.peak_addition_mode:
            self.delete_mode_checkbox.setChecked(False)
            self.peak_deletion_mode = False

        self.update_edit_mode_status()

    def update_edit_mode_status(self):
        """Update the edit mode status label based on current settings.

        This reflects the active editing modes (view only, delete, add)
        and which peak lists are being edited.
        """
        modes = []

        # Check which lists are being edited
        if self.edit_reference_peaks:
            modes.append("Reference")
        if self.edit_detected_peaks:
            modes.append("Detected")

        # Check the active mode
        if self.peak_deletion_mode:
            action = "Delete"
        elif self.peak_addition_mode:
            action = "Add"
        else:
            action = "View Only"

        # Build status message
        if modes:
            peaks_str = " + ".join(modes)
            status = f"Mode: {action} ({peaks_str} peaks)"
        else:
            status = f"Mode: {action}"

        self.edit_mode_status_label.setText(status)
        logger.debug(f"Edit mode status updated: {status}")

    # ===== Peak Editing Helper Methods (v0o9 lines 5632-6001) =====

    def _on_peak_edit_with_modifier(self, x_ppm: float, y_ppm: float, modifiers):
        """Handle peak editing via modifier keys (Shift/Ctrl + click).

        This is the new modifier-based peak editing that works with the
        NMRNavigationHandler. It allows adding and deleting peaks without
        entering a separate edit mode.

        Args:
            x_ppm: X coordinate in ppm where clicked
            y_ppm: Y coordinate in ppm where clicked
            modifiers: Qt.KeyboardModifiers indicating which keys are held
        """
        if modifiers & Qt.ShiftModifier:
            # Shift+Click: Add new peak at position
            self.add_new_peak(x_ppm, y_ppm)
            logger.info(f"Shift+Click: Added peak at ({x_ppm:.4f}, {y_ppm:.4f})")
        elif modifiers & Qt.ControlModifier:
            # Ctrl+Click: Delete nearest peak
            peak_info = self.find_nearest_peak(x_ppm, y_ppm)
            if peak_info:
                self.delete_selected_peak(peak_info)
                logger.info(f"Ctrl+Click: Deleted peak near ({x_ppm:.4f}, {y_ppm:.4f})")
            else:
                logger.warning(f"Ctrl+Click: No peak found near ({x_ppm:.4f}, {y_ppm:.4f})")

    def on_peak_edit_click(self, event):
        """Handle mouse clicks in peak edit mode (v0o9 line 5632)."""
        if not self.peak_edit_mode or event.inaxes != self.spectrum_plotter.ax:
            return

        click_x, click_y = event.xdata, event.ydata

        # Check if addition mode is enabled (highest priority)
        if self.peak_addition_mode:
            # Addition mode: add new peak at click position
            self.add_new_peak(click_x, click_y)
            return

        # Check if deletion mode is enabled
        if self.peak_deletion_mode:
            # Deletion mode: find and delete peak immediately
            peak_info = self.find_nearest_peak(click_x, click_y)
            if peak_info:
                self.delete_selected_peak(peak_info)
            return

        # Regular editing mode (move peaks)
        if self.selected_peak_info is None:
            # First click: select peak
            peak_info = self.find_nearest_peak(click_x, click_y)
            if peak_info:
                self.selected_peak_info = peak_info
                self.update_selected_peak_display()
        else:
            # Second click: move selected peak
            self.move_selected_peak(click_x, click_y)
            self.selected_peak_info = None
            if hasattr(self, 'selected_peak_label'):
                self.selected_peak_label.setText("Peak moved! No peak selected")

    def find_nearest_peak(self, click_x, click_y):
        """Find the nearest peak to click position (v0o9 line 5666)."""
        import math
        candidates = []

        # Check reference peaks (only if enabled for editing)
        if (self.edit_reference_peaks and
            hasattr(self.integrator, 'peak_list') and self.integrator.peak_list is not None):
            for idx, row in self.integrator.peak_list.iterrows():
                peak_x = float(row['Position_X'])
                peak_y = float(row['Position_Y'])
                distance = math.sqrt((click_x - peak_x)**2 + (click_y - peak_y)**2)
                if distance < self.click_tolerance:
                    candidates.append({
                        'type': 'reference',
                        'index': idx,
                        'distance': distance,
                        'x': peak_x,
                        'y': peak_y,
                        'assignment': row.get('Assignment', f'Peak_{idx+1}')
                    })

        # Check detected peaks (only if enabled for editing)
        if (self.edit_detected_peaks and
            hasattr(self.integrator, 'fitted_peaks') and self.integrator.fitted_peaks):
            for idx, peak in enumerate(self.integrator.fitted_peaks):
                peak_x = float(peak.get('ppm_x', 0))
                peak_y = float(peak.get('ppm_y', 0))
                distance = math.sqrt((click_x - peak_x)**2 + (click_y - peak_y)**2)
                if distance < self.click_tolerance:
                    candidates.append({
                        'type': 'detected',
                        'index': idx,
                        'distance': distance,
                        'x': peak_x,
                        'y': peak_y,
                        'assignment': peak.get('assignment', f'Det_{idx+1}')
                    })

        # Return closest peak
        if candidates:
            closest = min(candidates, key=lambda x: x['distance'])
            logger.info(f"Selected {closest['type']} peak {closest['index']}: {closest['assignment']} at ({closest['x']:.3f}, {closest['y']:.1f})")
            return closest

        return None

    def update_selected_peak_display(self):
        """Update display to show selected peak info (v0o9 line 5713)."""
        if self.selected_peak_info and hasattr(self, 'selected_peak_label'):
            peak_info = self.selected_peak_info
            if self.peak_deletion_mode:
                instruction_text = "In DELETE mode - click will delete this peak"
            else:
                instruction_text = "Click new position to move"

            self.selected_peak_label.setText(
                f"Selected: {peak_info['type']} peak '{peak_info['assignment']}' at ({peak_info['x']:.3f}, {peak_info['y']:.1f}) - {instruction_text}"
            )

    def move_selected_peak(self, new_x, new_y):
        """Move selected peak to new position (v0o9 line 5726)."""
        if not self.selected_peak_info:
            return

        peak_info = self.selected_peak_info
        old_x, old_y = peak_info['x'], peak_info['y']

        logger.info(f"Moving {peak_info['type']} peak '{peak_info['assignment']}' from ({old_x:.3f}, {old_y:.1f}) to ({new_x:.3f}, {new_y:.1f})")

        if peak_info['type'] == 'reference':
            # Update DataFrame
            peak_idx = peak_info['index']
            self.integrator.peak_list.loc[peak_idx, 'Position_X'] = new_x
            self.integrator.peak_list.loc[peak_idx, 'Position_Y'] = new_y

            # Update peak navigator if showing reference peaks
            if hasattr(self, 'peak_navigator') and hasattr(self.peak_navigator, 'selected_peak_type') and self.peak_navigator.selected_peak_type == 'reference':
                self.peak_navigator.load_reference_peaks(self.integrator.peak_list)

        elif peak_info['type'] == 'detected':
            # Update list of dictionaries
            peak_idx = peak_info['index']
            self.integrator.fitted_peaks[peak_idx]['ppm_x'] = new_x
            self.integrator.fitted_peaks[peak_idx]['ppm_y'] = new_y

            # Update peak navigator if showing detected peaks
            if hasattr(self, 'peak_navigator') and hasattr(self.peak_navigator, 'selected_peak_type') and self.peak_navigator.selected_peak_type == 'detected':
                self.peak_navigator.load_detected_peaks(self.integrator.fitted_peaks)

        # Refresh the main plot
        self.update_spectrum_plot()

        logger.info(f"✅ Peak position updated successfully")

    def delete_selected_peak(self, peak_info):
        """Delete the selected peak from the appropriate peak list (v0o9 line 5761)."""
        # Confirm deletion
        assignment = peak_info['assignment']
        peak_type = peak_info['type']

        result = QMessageBox.question(
            self,
            "Confirm Peak Deletion",
            f"Are you sure you want to delete {peak_type} peak '{assignment}'?\n\n"
            f"Position: ({peak_info['x']:.3f}, {peak_info['y']:.1f})\n"
            f"This action cannot be undone.",
            QMessageBox.Yes | QMessageBox.No
        )

        if result != QMessageBox.Yes:
            return

        try:
            if peak_info['type'] == 'reference':
                # Delete from DataFrame
                peak_idx = peak_info['index']

                # Remove the row from the DataFrame
                self.integrator.peak_list = self.integrator.peak_list.drop(index=peak_idx).reset_index(drop=True)

                logger.info(f"✅ Deleted reference peak '{assignment}' at index {peak_idx}")

                # Update peak navigator if showing reference peaks
                if hasattr(self, 'peak_navigator') and hasattr(self.peak_navigator, 'selected_peak_type') and self.peak_navigator.selected_peak_type == 'reference':
                    self.peak_navigator.load_reference_peaks(self.integrator.peak_list)

            elif peak_info['type'] == 'detected':
                # Delete from list of dictionaries
                peak_idx = peak_info['index']

                # Remove the peak from the list
                if 0 <= peak_idx < len(self.integrator.fitted_peaks):
                    del self.integrator.fitted_peaks[peak_idx]

                    logger.info(f"✅ Deleted detected peak '{assignment}' at index {peak_idx}")

                    # Update peak navigator if showing detected peaks
                    if hasattr(self, 'peak_navigator') and hasattr(self.peak_navigator, 'selected_peak_type') and self.peak_navigator.selected_peak_type == 'detected':
                        self.peak_navigator.load_detected_peaks(self.integrator.fitted_peaks)
                else:
                    logger.error(f"❌ Invalid peak index: {peak_idx}")
                    return

            # Clear selected peak
            self.selected_peak_info = None
            if hasattr(self, 'selected_peak_label'):
                self.selected_peak_label.setText(f"Deleted {peak_type} peak '{assignment}'")

            # Refresh the main plot
            self.update_spectrum_plot()

            # Update statistics
            self.update_statistics_panel()

            logger.info(f"✅ Peak deletion completed successfully")

        except Exception as e:
            QMessageBox.critical(self, "Deletion Error", f"Failed to delete peak: {str(e)}")
            logger.error(f"❌ Peak deletion failed: {e}")

    def _find_max_assignment_number(self, peaks_list, assignment_key='Assignment'):
        """Helper: Find highest assignment number from peak list (v0o9 line 5824)."""
        import re
        max_assignment = 0

        # Handle DataFrame (reference peaks)
        if hasattr(peaks_list, 'iterrows'):
            for idx, row in peaks_list.iterrows():
                assignment = row.get(assignment_key, '')
                max_assignment = max(max_assignment, self._extract_number_from_assignment(assignment))

        # Handle list of dicts (detected peaks)
        elif isinstance(peaks_list, list):
            for peak in peaks_list:
                assignment = peak.get(assignment_key, '') if isinstance(peak, dict) else ''
                max_assignment = max(max_assignment, self._extract_number_from_assignment(assignment))

        return max_assignment

    def _extract_number_from_assignment(self, assignment):
        """Helper: Extract numeric value from assignment string (v0o9 line 5855)."""
        import re
        try:
            if isinstance(assignment, (int, float)):
                return int(assignment)
            elif isinstance(assignment, str):
                # Extract digits from string (e.g., "Peak_123" -> 123, "45" -> 45)
                numbers = re.findall(r'\d+', str(assignment))
                if numbers:
                    return int(numbers[-1])  # Take the last number found
        except (ValueError, TypeError):
            pass
        return 0

    def _extract_intensity_at_position(self, x_ppm, y_ppm, window=3):
        """Extract intensity at specified position using local maximum (v0o9 line 5870)."""
        try:
            if not hasattr(self.integrator, 'nmr_data') or self.integrator.nmr_data is None:
                return 0.0

            # Find indices closest to clicked position
            x_idx = np.argmin(np.abs(self.integrator.ppm_x_axis - x_ppm))
            y_idx = np.argmin(np.abs(self.integrator.ppm_y_axis - y_ppm))

            # Extract local region (window around click point)
            x_slice = slice(max(0, x_idx - window), min(len(self.integrator.ppm_x_axis), x_idx + window + 1))
            y_slice = slice(max(0, y_idx - window), min(len(self.integrator.ppm_y_axis), y_idx + window + 1))

            local_data = self.integrator.nmr_data[y_slice, x_slice]

            # Return maximum absolute intensity in local region
            return float(np.max(np.abs(local_data)))

        except Exception as e:
            logger.warning(f"⚠️ Could not extract intensity: {e}")
            return 0.0

    def add_new_peak(self, click_x, click_y):
        """Add a new peak at the clicked position (v0o9 line 5906)."""
        import pandas as pd

        try:
            # Extract intensity from spectrum at clicked position
            intensity = self._extract_intensity_at_position(click_x, click_y)
            logger.info(f"📍 Adding peak at ({click_x:.3f}, {click_y:.1f})")
            logger.info(f"   Extracted intensity: {intensity:.2e}")

            # Find the highest assignment number from BOTH lists
            max_assignment_ref = 0
            max_assignment_det = 0

            if hasattr(self.integrator, 'peak_list') and self.integrator.peak_list is not None and len(self.integrator.peak_list) > 0:
                max_assignment_ref = self._find_max_assignment_number(self.integrator.peak_list, 'Assignment')

            if hasattr(self.integrator, 'fitted_peaks') and self.integrator.fitted_peaks:
                max_assignment_det = self._find_max_assignment_number(self.integrator.fitted_peaks, 'assignment')

            # Use the highest assignment number from either list
            max_assignment = max(max_assignment_ref, max_assignment_det)
            new_assignment = str(max_assignment + 1)

            # ============================================================
            # ADD TO REFERENCE PEAK LIST (DataFrame)
            # ============================================================
            new_peak_ref = pd.DataFrame([{
                'Assignment': new_assignment,
                'Position_X': click_x,
                'Position_Y': click_y,
                'Height': intensity
            }])

            if hasattr(self.integrator, 'peak_list') and self.integrator.peak_list is not None:
                self.integrator.peak_list = pd.concat([self.integrator.peak_list, new_peak_ref], ignore_index=True)
            else:
                self.integrator.peak_list = new_peak_ref

            logger.info(f"✅ Added to reference peak list: '{new_assignment}'")

            # ============================================================
            # ADD TO DETECTED PEAK LIST (list of dicts)
            # ============================================================
            new_peak_det = {
                'assignment': new_assignment,
                'ppm_x': click_x,
                'ppm_y': click_y,
                'height': intensity,
                'intensity': intensity,
                'r_squared': 0.0,
                'status': 'manual_add',
                'manual_add': True,
                'detected': True
            }

            if hasattr(self.integrator, 'fitted_peaks') and self.integrator.fitted_peaks:
                self.integrator.fitted_peaks.append(new_peak_det)
            else:
                self.integrator.fitted_peaks = [new_peak_det]

            logger.info(f"✅ Added to detected peak list: '{new_assignment}'")

            # Update peak navigator (refresh both lists)
            if hasattr(self, 'peak_navigator'):
                if hasattr(self.peak_navigator, 'selected_peak_type'):
                    if self.peak_navigator.selected_peak_type == 'reference':
                        self.peak_navigator.load_reference_peaks(self.integrator.peak_list)
                    elif self.peak_navigator.selected_peak_type == 'detected':
                        self.peak_navigator.load_detected_peaks(self.integrator.fitted_peaks)

            # Update status label
            if hasattr(self, 'selected_peak_label'):
                self.selected_peak_label.setText(
                    f"Added peak '{new_assignment}' at ({click_x:.3f}, {click_y:.1f}) | Intensity: {intensity:.2e}"
                )

            # Refresh the main plot
            self.update_spectrum_plot()

            # Update statistics
            self.update_statistics_panel()

            logger.info(f"✅ Peak addition completed: added to BOTH lists with intensity extraction")

        except Exception as e:
            QMessageBox.critical(self, "Addition Error", f"Failed to add peak: {str(e)}")
            logger.error(f"❌ Peak addition failed: {e}")
            import traceback
            traceback.print_exc()

    # ===== Project Management =====

    def _get_project_manager(self) -> ProjectManager:
        """Get or create ProjectManager instance."""
        if self.project_manager is None:
            self.project_manager = ProjectManager(self)
        return self.project_manager

    def open_project(self):
        """Open a saved project bundle."""
        from PySide6.QtWidgets import QFileDialog
        from lunaNMR.gui.dialogs.missing_files_dialog import MissingFilesDialog

        # Check for unsaved changes
        if not self._check_unsaved_changes():
            return

        # Get project path
        project_path = QFileDialog.getExistingDirectory(
            self,
            "Open Project",
            str(Path.home()),
        )

        if not project_path:
            return

        project_path = Path(project_path)

        # Validate it's a .lunaNMR bundle
        if not project_path.suffix == '.lunaNMR' and not project_path.name.endswith('.lunaNMR'):
            QMessageBox.warning(
                self,
                "Invalid Project",
                "Please select a .lunaNMR project folder."
            )
            return

        # Check for missing files before loading
        pm = self._get_project_manager()
        missing_files = pm.get_missing_files_structured(project_path)

        remapped_paths = {}
        skipped_files = set()

        # Show dialog if there are missing files
        if missing_files:
            dialog = MissingFilesDialog(self, missing_files=missing_files)
            if not dialog.exec():
                # User cancelled
                return

            remapped_paths = dialog.get_remapped_paths()
            skipped_files = dialog.get_skipped_files()

        # Load project
        success, error_messages, summary = pm.load_project(project_path)

        if success:
            self.current_project_path = project_path

            # Apply path remapping
            if remapped_paths:
                pm.apply_path_remapping(remapped_paths)

            # Log skipped files
            if skipped_files:
                logger.info(f"Skipped files: {skipped_files}")

            # Refresh display
            self._refresh_after_project_load()

            self.update_status(f"Project loaded: {project_path.name}")
            logger.info(f"Project loaded from {project_path}")

            # Show summary popup
            self._show_load_summary(project_path, summary)
        else:
            QMessageBox.critical(
                self,
                "Load Failed",
                f"Failed to load project:\n{', '.join(error_messages)}"
            )

    def save_project(self):
        """Save project to current path, or prompt for path if none."""
        if self.current_project_path is None:
            self.save_project_as()
        else:
            self._do_save_project(self.current_project_path)

    def save_project_as(self):
        """Save project to a new location."""
        from PySide6.QtWidgets import QFileDialog

        # Suggest a name based on NMR file
        suggested_name = "untitled.lunaNMR"
        if self.current_nmr_file:
            suggested_name = Path(self.current_nmr_file).stem + ".lunaNMR"

        project_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Project As",
            str(Path.home() / suggested_name),
            "LunaNMR Project (*.lunaNMR)"
        )

        if not project_path:
            return

        project_path = Path(project_path)

        # Ensure .lunaNMR extension
        if not project_path.suffix == '.lunaNMR':
            project_path = project_path.with_suffix('.lunaNMR')

        self._do_save_project(project_path)

    def _do_save_project(self, project_path: Path):
        """Actually save the project to the given path."""
        pm = self._get_project_manager()
        success, summary = pm.save_project(project_path)

        if success:
            self.current_project_path = project_path
            self.update_status(f"Project saved: {project_path.name}")
            logger.info(f"Project saved to {project_path}")

            # Show summary popup
            self._show_save_summary(project_path, summary)
        else:
            QMessageBox.critical(
                self,
                "Save Failed",
                "Failed to save project. Check logs for details."
            )

    def _show_save_summary(self, project_path: Path, summary: dict):
        """Show a popup summarizing what was saved."""
        lines = [f"Project saved: {project_path.name}", ""]

        # General items
        if summary.get('saved_items'):
            lines.append("Saved items:")
            for item in summary['saved_items']:
                lines.append(f"  • {item}")

        # DynamiXs summary
        dynamixs = summary.get('dynamixs', {})
        if dynamixs:
            lines.append("")
            lines.append("DynamiXs Analysis:")

            if 't1t2_fitting' in dynamixs:
                t1t2 = dynamixs['t1t2_fitting']
                experiments = t1t2.get('fitted_experiments', [])
                if experiments:
                    exp_str = ', '.join(e.replace('_', ' ').title() for e in experiments)
                    lines.append(f"  • T1/T2 Fitting: {exp_str}")

            if 'spectral_density' in dynamixs:
                lines.append("  • Spectral Density: Analysis complete")

            if 'model_free' in dynamixs:
                lines.append("  • Model Free: Analysis complete")

        if not dynamixs and len(summary.get('saved_items', [])) <= 2:
            lines.append("")
            lines.append("No analysis results to save.")
            lines.append("Complete analyses in DynamiXs to save results.")

        logger.info(f"Showing save summary popup")
        msg = QMessageBox(self)
        msg.setWindowTitle("Project Saved")
        msg.setText("\n".join(lines))
        msg.setIcon(QMessageBox.Information)
        msg.setWindowModality(Qt.ApplicationModal)
        msg.raise_()
        msg.activateWindow()
        msg.exec()

    def _show_load_summary(self, project_path: Path, summary: dict):
        """Show a popup summarizing what was loaded."""
        lines = [f"Project loaded: {project_path.name}", ""]

        # General items
        if summary.get('loaded_items'):
            lines.append("Loaded items:")
            for item in summary['loaded_items']:
                lines.append(f"  • {item}")

        # DynamiXs summary
        dynamixs = summary.get('dynamixs', {})
        if dynamixs:
            lines.append("")
            lines.append("DynamiXs Analysis:")

            if 't1t2_fitting' in dynamixs:
                t1t2 = dynamixs['t1t2_fitting']
                experiments = t1t2.get('fitted_experiments', [])
                if experiments:
                    exp_str = ', '.join(e.replace('_', ' ').title() for e in experiments)
                    lines.append(f"  • T1/T2 Fitting: {exp_str}")

            if 'spectral_density' in dynamixs:
                lines.append("  • Spectral Density: Analysis complete")

            if 'model_free' in dynamixs:
                lines.append("  • Model Free: Analysis complete")

            lines.append("")
            lines.append("Open DynamiXs to view results.")
        else:
            lines.append("")
            lines.append("No DynamiXs analysis results in project.")

        logger.info(f"Showing load summary popup")
        msg = QMessageBox(self)
        msg.setWindowTitle("Project Loaded")
        msg.setText("\n".join(lines))
        msg.setIcon(QMessageBox.Information)
        msg.setWindowModality(Qt.ApplicationModal)
        msg.raise_()
        msg.activateWindow()
        msg.exec()

    def _check_unsaved_changes(self) -> bool:
        """Check for unsaved changes and prompt user.

        Returns:
            True if OK to proceed (no changes or user chose to discard)
            False if user cancelled
        """
        # For now, always return True - we can add change tracking later
        # TODO: Implement change tracking to detect unsaved modifications
        return True

    def _refresh_after_project_load(self):
        """Refresh UI after loading a project."""
        # Update file labels
        self.update_file_status_labels(
            nmr_file=self.current_nmr_file,
            peak_file=self.current_peak_file
        )

        # Sync fit results to integrator.fitted_peaks BEFORE loading spectrum
        # (load_nmr_file uses integrator.fitted_peaks for peak display)
        if self.last_fitting_results:
            self.integrator.fitted_peaks = self.last_fitting_results
            logger.info(f"Synced {len(self.last_fitting_results)} fit results to integrator.fitted_peaks")

        # Refresh spectrum display if NMR file is loaded
        if self.current_nmr_file and Path(self.current_nmr_file).exists():
            try:
                self.load_nmr_file(self.current_nmr_file)
            except Exception as e:
                logger.warning(f"Could not reload spectrum: {e}")

        # Update peak navigator with reference peaks and fit results
        if hasattr(self, 'peak_navigator'):
            # Load reference peaks from peak list
            if self.integrator.peak_list is not None:
                self.peak_navigator.load_reference_peaks(self.integrator.peak_list)

            # Load detected/fitted peaks from fit results
            if self.last_fitting_results:
                self.peak_navigator.load_detected_peaks(self.last_fitting_results)
            else:
                self.peak_navigator.refresh_peak_list()

        # Update Voigt analysis plotter if fit results exist
        if hasattr(self, 'voigt_plotter') and self.last_fitting_results:
            try:
                # Select first peak to display
                if len(self.last_fitting_results) > 0:
                    self.voigt_plotter.display_results(self.last_fitting_results[0])
                    logger.info(f"Loaded {len(self.last_fitting_results)} fit results into Voigt plotter")
            except Exception as e:
                logger.warning(f"Could not update Voigt plotter: {e}")

    # ===== Test Spectrum Creation =====

    def create_test_spectrum(self):
        """Create and display a synthetic 2D NMR spectrum for testing.

        This generates:
        - Synthetic 2D data with Gaussian peaks
        - Coordinate axes (1H: 6-10 ppm, 15N: 100-130 ppm)
        - Multiple test peaks at known positions

        This is used for:
        - Verification of SpectrumPlotter integration
        - Testing visualization features (zoom, pan, toolbar)
        - Development without requiring real NMR data
        """
        logger.info("Creating synthetic test spectrum")

        import numpy as np

        # Create coordinate axes
        # 1H dimension: 6-10 ppm (high to low)
        ppm_x = np.linspace(10.0, 6.0, 256)
        # 15N dimension: 100-130 ppm (high to low)
        ppm_y = np.linspace(130.0, 100.0, 256)

        # Create 2D data grid
        X, Y = np.meshgrid(ppm_x, ppm_y)

        # Initialize data with noise
        data = np.random.normal(0, 0.05, (256, 256))

        # Add synthetic peaks (Gaussian profiles) - Smiley face
        # Peak parameters: (x_pos, y_pos, amplitude, width_x, width_y)
        # Note: NMR Y-axis is inverted (high ppm at top), so lower Y = higher on plot
        test_peaks = [
            # === Left eye (lower ppm = appears at top) ===
            (7.3, 111.0, 1.5, 0.07, 0.7),
            # === Right eye ===
            (8.7, 111.0, 1.5, 0.07, 0.7),
            # === Smile (arc curving up in ppm = appears curving down on plot) ===
            (7.0, 118.0, 1.0, 0.05, 0.5),    # Left corner
            (7.25, 119.5, 1.0, 0.05, 0.5),
            (7.5, 120.5, 1.0, 0.05, 0.5),
            (7.75, 121.0, 1.0, 0.05, 0.5),
            (8.0, 121.5, 1.0, 0.05, 0.5),    # Bottom center (highest ppm)
            (8.25, 121.0, 1.0, 0.05, 0.5),
            (8.5, 120.5, 1.0, 0.05, 0.5),
            (8.75, 119.5, 1.0, 0.05, 0.5),
            (9.0, 118.0, 1.0, 0.05, 0.5),    # Right corner
        ]

        # Add Gaussian peaks
        for x_pos, y_pos, amp, width_x, width_y in test_peaks:
            peak = amp * np.exp(
                -((X - x_pos)**2 / (2 * width_x**2) +
                  (Y - y_pos)**2 / (2 * width_y**2))
            )
            data += peak

        # Create mock integrator object with NMR data
        class MockIntegrator:
            def __init__(self):
                self.nmr_data = data
                self.ppm_x_axis = ppm_x
                self.ppm_y_axis = ppm_y
                self.threshold = 0.1

        mock_integrator = MockIntegrator()

        # Plot the test spectrum
        try:
            self.spectrum_plotter.plot_spectrum(
                mock_integrator,
                contour_min_level=0.1,
                contour_levels=15,
                contour_increment=1.2,
                colormap='viridis',
                show_colorbar=False
            )
            logger.info("Test spectrum plotted successfully")
            self.update_status("Test spectrum loaded - smiley face (11 peaks)")
        except Exception as e:
            logger.error(f"Failed to plot test spectrum: {e}")
            import traceback
            traceback.print_exc()
            self.update_status(f"Error plotting test spectrum: {e}")

# ===== Module-level functions for testing =====

def main():
    """Launch the main window for testing."""
    import sys
    from PySide6.QtWidgets import QApplication

    app = QApplication(sys.argv)

    # Create and show main window
    window = LunaNMRMainWindow()
    window.show()

    # Create test spectrum for visualization testing
    window.create_test_spectrum()

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
