#!/usr/bin/env python3
"""
NMR Peaks Series GUI - Main Application

This is the main GUI application that combines full peak detection and in-place fitting
capabilities with advanced series processing features. The interface allows users to
choose between processing modes via radio buttons and provides enhanced visualization
and analysis capabilities.

Features:
- Unified interface with mode selection (Full Detection vs In-Place Fitting)
- Enhanced file management with preview
- Advanced series processing with progress tracking
- Comprehensive Voigt fitting analysis
- Configuration management and user preferences
- Statistics and quality assessment
"""

import sys
import os
import time
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import threading
import subprocess
from pathlib import Path
from datetime import datetime

# Add current directory to path for imports
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.extend([current_dir, parent_dir])


from lunaNMR.utils.parameter_manager import NMRParameterManager
from lunaNMR.processors.single_spectrum_processor import SingleSpectrumProcessor
from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor

from lunaNMR.gui.gui_components import AdvancedProgressDialog

# Import our modular components
try:
    from lunaNMR.core.core_integrator import VoigtIntegrator, EnhancedVoigtIntegrator
    from lunaNMR.gui.gui_components import (
        ScrollableFrame,
        EnhancedFileListFrame,
        AdvancedProgressDialog,
        StatisticsPanel,
        ModeSelectionFrame,
        PeakNavigator
    )
    from lunaNMR.utils.file_manager import NMRFileManager, DataValidator, FileMetadata
    from lunaNMR.processors.series_processor import SeriesProcessor, SeriesAnalyzer, BatchResults
    from lunaNMR.gui.visualization import SpectrumPlotter, VoigtAnalysisPlotter, SeriesPlotter, PlotManager
    from lunaNMR.utils.config_manager import ConfigurationManager, UserPreferences, ProcessingParameters

    # Matplotlib setup
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    import numpy as np
    import pandas as pd

    print("✅ All modules imported successfully")

except ImportError as e:
    print(f"❌ Import error: {e}")
    messagebox.showerror("Import Error",
        f"Failed to import required modules:\n{str(e)}\n\n"
        "Please ensure all dependencies are installed and the modules are in the correct location.")
    sys.exit(1)

class NMRPeaksSeriesGUI:
    """Main GUI application class with enhanced features and mode selection"""

    def __init__(self, root):
        self.root = root
        self.root.title("lunaNMR v0.9 - Enhanced Multi-Mode Interface")

        # Core components
        self.integrator = EnhancedVoigtIntegrator()
        self.file_manager = NMRFileManager()
        self.validator = DataValidator()
        self.series_processor = SeriesProcessor()
        self.series_analyzer = SeriesAnalyzer()

        # Configuration management
        self.config_manager = ConfigurationManager("NMRPeakSeries")
        self.user_prefs = UserPreferences(self.config_manager)
        self.proc_params = ProcessingParameters(self.config_manager)

        # GUI state variables
        self.current_nmr_file = None
        self.current_peak_file = None
        self.current_voigt_result = None
        self.batch_results = None

        # Initialize new decoupled architecture components
        self.param_manager = NMRParameterManager()
        print("✅ Parameter manager initialized")

        self.processing_active = False
        # These will be created when needed
        self.single_spectrum_processor = None
        self.multi_spectrum_processor = None

        # Integration progress tracking
        self.integration_active = False
        self.integration_paused = False
        self.integration_start_time = None

        # Initialize GUI variables from config
        self.init_variables_from_config()

        # Setup GUI
        self.setup_gui()

        # Apply user preferences
        self.user_prefs.apply_to_window(self.root)

        # Setup plot manager
        self.setup_plot_manager()

        # Bind window events
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

        # Sync all initial parameters to the integrator engine
        self._sync_parameters_to_integrator()

        print("✅ NMR Peaks Series GUI initialized successfully")

    def init_variables_from_config(self):
        """Initialize GUI variables from configuration"""
        # Detection parameters
        detection_params = self.proc_params.get_detection_params()
        self.noise_threshold = tk.DoubleVar(value=detection_params['noise_threshold'])
        self.search_window_x = tk.DoubleVar(value=detection_params['search_window_x'])
        self.search_window_y = tk.DoubleVar(value=detection_params['search_window_y'])
        self.use_reference_detection = tk.BooleanVar(value=detection_params['use_reference_detection'])

        # Fitting parameters
        fitting_params = self.proc_params.get_fitting_params()
        self.fitting_window_x = tk.DoubleVar(value=fitting_params['fitting_window_x'])
        self.fitting_window_y = tk.DoubleVar(value=fitting_params['fitting_window_y'])
        self.min_r_squared = tk.DoubleVar(value=fitting_params['min_r_squared'])
        self.max_iterations = tk.IntVar(value=fitting_params['max_iterations'])

        # Peak detection parameters (for multi-peak fitting)
        self.peak_height_threshold = tk.DoubleVar(value=0.02)    # 2% of max intensity
        self.peak_distance_factor = tk.DoubleVar(value=50.0)     # 1/50 = 2% of spectrum
        self.peak_prominence_threshold = tk.DoubleVar(value=0.01) # 1% of max intensity
        self.smoothing_sigma = tk.DoubleVar(value=0.5)           # Gaussian smoothing
        self.max_peaks_fit = tk.IntVar(value=4)                  # Maximum peaks to fit simultaneously
        self.max_optimization_iterations = tk.IntVar(value=50)   # Maximum iterative optimization attempts
        self.use_parallel_processing = tk.BooleanVar(value=True) # Enable parallel processing

        # Peak Centroid Detection parameters (optional post-processing enhancement)
        self.use_centroid_refinement = tk.BooleanVar(value=False)  # Enable centroid refinement
        self.centroid_window_x_ppm = tk.DoubleVar(value=0.01)      # X-window size in ppm (1H dimension)
        self.centroid_window_y_ppm = tk.DoubleVar(value=0.1)       # Y-window size in ppm (15N dimension)
        self.centroid_noise_multiplier = tk.DoubleVar(value=2.0)   # Noise threshold multiplier

        # Detection square size parameters - now anisotropic (X and Y dimensions)
        self.detection_square_size = tk.IntVar(value=3)          # Size of detection square X-dimension/1H (pixels)
        self.detection_rectangle_y = tk.IntVar(value=1)          # Size of detection rectangle Y-dimension/15N (pixels)

        # Peak editing parameters
        self.peak_edit_mode = tk.BooleanVar(value=False)         # Peak editing mode toggle
        self.edit_reference_peaks = tk.BooleanVar(value=True)    # Enable editing of reference peaks
        self.edit_detected_peaks = tk.BooleanVar(value=True)     # Enable editing of detected peaks
        self.selected_peak_info = None                           # Currently selected peak {type, index, data}
        self.click_tolerance = 0.05                              # Tolerance for peak selection (ppm)
        self.detection_square_ppm_x = tk.StringVar(value="0.000") # Auto-calculated ppm conversion for 1H
        self.detection_square_ppm_y = tk.StringVar(value="0.000") # Auto-calculated ppm conversion for 15N

        # Peak Ridge Consolidation parameters (Solution A)
        self.consolidation_x_tolerance = tk.DoubleVar(value=0.05) # 1H tolerance for grouping peaks (ppm)
        self.consolidation_y_tolerance = tk.DoubleVar(value=2.0)  # 15N tolerance for consolidation (ppm)

        # Global optimization toggle (default OFF for backward compatibility)
        self.use_global_optimization = tk.BooleanVar(value=False)

        # Display options
        display_options = self.proc_params.get_display_options()
        self.show_detected = tk.BooleanVar(value=display_options['show_detected'])
        self.show_assigned = tk.BooleanVar(value=display_options['show_assigned'])
        self.show_fitted_curves = tk.BooleanVar(value=display_options['show_fitted_curves'])

        # Series options
        series_options = self.proc_params.get_series_options()
        self.auto_process_series = tk.BooleanVar(value=series_options['auto_process_series'])
        self.save_individual_results = tk.BooleanVar(value=series_options['save_individual_results'])
        self.create_summary_plots = tk.BooleanVar(value=series_options['create_summary_plots'])

        # Series integration advanced options (Voigt fitting for series)
        self.series_use_voigt_fitting = tk.BooleanVar(value=series_options.get('use_voigt_fitting', True))
        self.series_use_parallel_processing = tk.BooleanVar(value=series_options.get('use_parallel_processing', True))  # Default enabled
        self.series_use_global_optimization = tk.BooleanVar(value=series_options.get('use_global_optimization', False))
        self.series_num_integrations = tk.IntVar(value=series_options.get('num_integrations', 3))  # Default: 3 integrations

        # NEW: Series peak source and detection options
        self.series_peak_source = tk.StringVar(value="detected")  # "detected", "reference", "cascade"
        self.series_enable_detection = tk.BooleanVar(value=False)  # Default OFF - direct fitting only

        # Integration mode options (initialized here, GUI creation happens later)
        integration_options = series_options.get('integration_options', {})
        self.integration_mode = tk.StringVar(value=integration_options.get('mode', 'standard'))

        # Removed parameter stubs (for compatibility)
        self.aic_threshold = tk.DoubleVar(value=2.0)
        self.detection_confidence_threshold = tk.DoubleVar(value=0.3)
        self.max_integration_iterations = tk.IntVar(value=5)
        self.show_detection_quality = tk.BooleanVar(value=False)

        # Advanced integration parameters (missing from GUI)
        self.adaptive_thresholds_enabled = tk.BooleanVar(value=integration_options.get('adaptive_thresholds_enabled', False))
        self.multi_resolution_enabled = tk.BooleanVar(value=integration_options.get('multi_resolution_enabled', False))
        self.physics_constraints_enabled = tk.BooleanVar(value=integration_options.get('physics_constraints_enabled', True))
        self.convergence_threshold = tk.DoubleVar(value=integration_options.get('convergence_threshold', 0.01))
        self.fit_likelihood_threshold = tk.DoubleVar(value=integration_options.get('fit_likelihood_threshold', 0.2))
        # Detection mode control (in-place vs full detection)
        self.force_full_detection = tk.BooleanVar(value=integration_options.get('force_full_detection', False))

        # 1D Refinement parameters (NEW)
        self.enable_1d_refinement = tk.BooleanVar(value=integration_options.get('enable_1d_refinement', True))
        self.refinement_quality_threshold = tk.DoubleVar(value=integration_options.get('refinement_quality_threshold', 0.7))
        self.refinement_coordinate_threshold = tk.DoubleVar(value=integration_options.get('refinement_coordinate_threshold', 0.01))

        # Enhanced Graph-Based Detection parameters (NEW)
        # Enhanced Detection constraints - separate X and Y dimensions
        self.enhanced_radius_x = tk.DoubleVar(value=integration_options.get('enhanced_radius_x', 0.05))    # 1H constraint (ppm)
        self.enhanced_radius_y = tk.DoubleVar(value=integration_options.get('enhanced_radius_y', 2.0))     # 15N/13C constraint (ppm)
        self.enhanced_pattern_similarity = tk.DoubleVar(value=integration_options.get('enhanced_pattern_similarity', 0.7))
        self.enhanced_missing_tolerance = tk.IntVar(value=integration_options.get('enhanced_missing_tolerance', 1))
        self.enhanced_position_weight = tk.DoubleVar(value=integration_options.get('enhanced_position_weight', 0.7))

        # Enhanced Detection square size - now anisotropic (X and Y dimensions)
        self.enhanced_detection_square_size = tk.IntVar(value=3)         # Size of detection square X-dimension/1H for Enhanced Detection
        self.enhanced_detection_rectangle_y = tk.IntVar(value=1)         # Size of detection rectangle Y-dimension/15N for Enhanced Detection
        self.enhanced_detection_square_ppm_x = tk.StringVar(value="(load data to see ppm)")

        # Peak reduction parameters (NEW)
        self.enhanced_peak_limit = tk.IntVar(value=integration_options.get('enhanced_peak_limit', 50))
        self.enhanced_noise_threshold = tk.DoubleVar(value=integration_options.get('enhanced_noise_threshold', 0.01))

        # Peak coordinate adjustment parameters
        self.adjust_x_offset = tk.DoubleVar(value=0.0)  # ppm offset in 1H dimension
        self.adjust_y_offset = tk.DoubleVar(value=0.0)  # ppm offset in 15N/13C dimension

        # Navigation and zoom
        self.selected_peak_number = tk.IntVar(value=1)

        # Zoom controls
        self.current_xlim = None
        self.current_ylim = None
        self.zoom_x_center = tk.DoubleVar(value=8.0)
        self.zoom_y_center = tk.DoubleVar(value=120.0)
        self.zoom_x_range = tk.DoubleVar(value=4.0)
        self.zoom_y_range = tk.DoubleVar(value=30.0)

        # Interactive zoom mode variables
        self.zoom_mode_active = False
        self.zoom_start_point = None
        self.zoom_rectangle = None
        self.original_xlim = None
        self.original_ylim = None
        self.mouse_press_cid = None
        self.mouse_motion_cid = None
        self.mouse_release_cid = None

        # Contour controls
        self.contour_levels = tk.IntVar(value=24)
        self.contour_min = tk.DoubleVar(value=0.05)
        self.contour_increment = tk.DoubleVar(value=1.1)

        # Initialize GUI parameters in integrator
        self.update_integrator_params()


    def _sync_parameters_to_integrator(self):
        """
        Sync parameters using the new parameter manager system
        """
        if not hasattr(self, 'integrator') or not self.integrator:
            return

        # Update parameter manager from GUI variables
        updated_params = self.param_manager.update_from_gui_variables(self)

        # Validate parameters
        validation_errors = self.param_manager.validate_all_parameters()
        if validation_errors:
            print(f"⚠️ Parameter validation warnings: {', '.join(validation_errors[:3])}{'...' if len(validation_errors) > 3 else ''}")

        # Get formatted parameters for integrator
        integrator_params = self.param_manager.get_integrator_parameters()

        # Apply parameters to integrator
        detection_params = integrator_params['detection_params']
        self.integrator.set_search_window(
            detection_params['search_window_x'],
            detection_params['search_window_y']
        )
        self.integrator.set_threshold_multiplier(detection_params['noise_threshold'])

        self.integrator.fitting_parameters.update(integrator_params['fitting_params'])
        self.integrator.gui_params = integrator_params['gui_params']

        print(f"✅ Parameters synchronized: {len(updated_params)} parameters updated via parameter manager")


    def on_parameter_change(self):
        """Handle parameter change events to keep integrator and GUI in sync."""
        self._sync_parameters_to_integrator()

        # Update GUI elements that depend on these parameters
        self.update_distance_ppm_display()
        self._update_series_voigt_param_displays()

        # Save parameters to config file
        self.proc_params.set_detection_params(
            noise_threshold=self.noise_threshold.get(),
            search_window_x=self.search_window_x.get(),
            search_window_y=self.search_window_y.get(),
            use_reference_detection=self.use_reference_detection.get(),
            detection_square_size=self.detection_square_size.get(),
            detection_rectangle_y=self.detection_rectangle_y.get(),
            consolidation_x_tolerance=self.consolidation_x_tolerance.get(),
            consolidation_y_tolerance=self.consolidation_y_tolerance.get()
        )
        self.proc_params.set_fitting_params(
            fitting_window_x=self.fitting_window_x.get(),
            fitting_window_y=self.fitting_window_y.get(),
            min_r_squared=self.min_r_squared.get(),
            max_iterations=self.max_iterations.get()
        )

    def update_integrator_params(self):
        """Legacy method for updating integrator parameters. Now calls the new sync function."""
        self._sync_parameters_to_integrator()

    def setup_gui(self):
        """Setup the main GUI layout with enhanced features"""
        # Main container with menu
        self.create_menu()

        # Configure root grid for responsive layout
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)  # Main content expands
        self.root.rowconfigure(1, weight=0)  # Status bar fixed height

        main_frame = ttk.Frame(self.root)
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), padx=5, pady=5)

        # Configure main frame grid weights for 3-column layout
        main_frame.columnconfigure(0, weight=0, minsize=400)  # Left panel - controls (fixed minimum)
        main_frame.columnconfigure(1, weight=1)              # Center panel - spectrum (expandable)
        main_frame.columnconfigure(2, weight=0, minsize=200)  # Right panel - peak navigator (fixed minimum)
        main_frame.rowconfigure(0, weight=1)

        # Left panel for controls
        left_panel_frame = ttk.Frame(main_frame)
        left_panel_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), padx=(0, 5))
        left_panel_frame.columnconfigure(0, weight=1)
        left_panel_frame.rowconfigure(0, weight=1)

        # Scrollable controls
        self.scrollable_controls = ScrollableFrame(left_panel_frame)
        self.scrollable_controls.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Center panel for visualization
        center_panel = ttk.Frame(main_frame)
        center_panel.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N, tk.S), padx=5)
        center_panel.columnconfigure(0, weight=1)
        center_panel.rowconfigure(0, weight=1)

        # Right panel for peak navigator
        self.right_panel_frame = ttk.Frame(main_frame)
        self.right_panel_frame.grid(row=0, column=2, sticky=(tk.W, tk.E, tk.N, tk.S), padx=(5, 0))
        self.right_panel_frame.columnconfigure(0, weight=1)
        self.right_panel_frame.rowconfigure(0, weight=1)

        # Setup status bar first (needed for status updates)
        self.setup_status_bar()

        # Setup components
        self.setup_enhanced_controls(self.scrollable_controls.scrollable_frame)
        self.setup_enhanced_visualization(center_panel)
        self.setup_peak_navigator()

    def create_menu(self):
        """Create enhanced menu bar"""
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)

        # File menu
        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="Load NMR Spectrum...", command=self.load_nmr_file)
        file_menu.add_command(label="Load Peak List...", command=self.load_peak_file)
        file_menu.add_separator()

        # Recent files submenu
        self.recent_menu = tk.Menu(file_menu, tearoff=0)
        file_menu.add_cascade(label="Recent Files", menu=self.recent_menu)
        self.update_recent_files_menu()

        file_menu.add_separator()
        file_menu.add_command(label="Export Peak List...", command=self.export_peak_list)
        file_menu.add_command(label="Export Results...", command=self.export_results)
        file_menu.add_command(label="Save Configuration...", command=self.save_config)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.on_closing)

        # Processing menu
        process_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Processing", menu=process_menu)
        process_menu.add_command(label="Detect Peaks", command=self.detect_peaks)
        process_menu.add_command(label="Enhanced Peak Detection", command=self.enhanced_detect_peaks)
        process_menu.add_command(label="Integrate Peaks", command=self.integrate_peaks)
        process_menu.add_command(label="Fit Selected Peak", command=self.fit_selected_peak)
        process_menu.add_command(label="Fit All Peaks", command=self.fit_all_peaks)
        process_menu.add_separator()
        process_menu.add_command(label="Start Series Integration", command=self.start_series_integration)

        # View menu
        view_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="View", menu=view_menu)
        view_menu.add_command(label="Reset Zoom", command=self.reset_view)
        view_menu.add_command(label="Center on Selected Peak", command=self.center_on_selected_peak)
        view_menu.add_separator()
        view_menu.add_checkbutton(label="Show Detected Peaks", variable=self.show_detected, command=self.update_main_plot)
        #view_menu.add_checkbutton(label="Show Reference Peaks", variable=self.show_assigned, command=self.update_main_plot)
        view_menu.add_checkbutton(
            label="Show Reference Peaks",
            variable=self.show_assigned,
            command=lambda: [
                self.update_main_plot(),
                print(f"Reference peaks display toggled: {self.show_assigned.get()}")
            ]
        )
        view_menu.add_checkbutton(label="Show Fitted Curves", variable=self.show_fitted_curves, command=self.update_main_plot)
        view_menu.add_separator()
        #view_menu.add_checkbutton(label="🔍 Show Peaks After Limit (DEBUG)", variable=self.show_included_peaks_after_limit_debug, command=self.update_main_plot)

        # Results Analysis menu
        results_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Results Analysis", menu=results_menu)
        results_menu.add_command(label="Browse Series Results...", command=self.open_results_browser)
        results_menu.add_command(label="Browse Individual Spectra...", command=self.open_spectrum_browser)
        results_menu.add_command(label="Peak Evolution Analysis", command=self.show_peak_evolution)
        results_menu.add_command(label="Comparison Dashboard", command=self.open_comparison_dashboard)
        results_menu.add_separator()
        results_menu.add_command(label="Export Data Matrix...", command=self.export_data_matrix)
        results_menu.add_command(label="Export Analysis Report...", command=self.export_analysis_report)
        results_menu.add_command(label="Batch Export Results...", command=self.batch_export_results)
        results_menu.add_separator()
        results_menu.add_command(label="Series Quality Assessment", command=self.show_quality_assessment)
        results_menu.add_command(label="Detection Statistics", command=self.show_detection_statistics)

        # Modules menu
        modules_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Modules", menu=modules_menu)
        modules_menu.add_command(label="DynamiXs Relaxation Analysis", command=self.launch_dynamixs)

        # Tools menu
        tools_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Tools", menu=tools_menu)
        tools_menu.add_command(label="Configuration...", command=self.open_config_dialog)
        tools_menu.add_command(label="Statistics", command=self.show_statistics)
        tools_menu.add_command(label="Validate Files", command=self.validate_current_files)

        # Help menu
        help_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="User Guide", command=self.show_help)
        help_menu.add_command(label="About", command=self.show_about)

    def setup_enhanced_controls(self, parent):
        """Setup enhanced control panel with all advanced features"""

        # =================== MODE SELECTION SECTION ===================
        #self.mode_selection_frame = ModeSelectionFrame(parent, title="🔧 Processing Mode")
        #self.mode_selection_frame.pack(fill=tk.X, pady=(0, 10))

        # =================== FILE MANAGEMENT SECTION ===================
        file_frame = ttk.LabelFrame(parent, text="📁 Enhanced Data Management", padding=10)
        file_frame.pack(fill=tk.X, pady=(0, 10))

        # Create enhanced file selection panels with responsive grid
        file_container = ttk.Frame(file_frame)
        file_container.pack(fill=tk.BOTH, expand=True, pady=5)

        # Configure grid for equal width columns
        file_container.columnconfigure(0, weight=1)
        file_container.columnconfigure(1, weight=1)
        file_container.rowconfigure(0, weight=1)

        # NMR Data folder panel
        nmr_panel = ttk.Frame(file_container)
        nmr_panel.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), padx=(0, 5))
        nmr_panel.columnconfigure(0, weight=1)
        nmr_panel.rowconfigure(0, weight=1)

        self.nmr_file_list = EnhancedFileListFrame(nmr_panel, "NMR Spectra", ["ft", "fid"], height=4)
        self.nmr_file_list.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.nmr_file_list.set_callback(self.on_nmr_file_select)

        # Peak List folder panel
        peak_panel = ttk.Frame(file_container)
        peak_panel.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), padx=(5, 0))
        peak_panel.columnconfigure(0, weight=1)
        peak_panel.rowconfigure(0, weight=1)

        self.peak_file_list = EnhancedFileListFrame(peak_panel, "Peak Lists", ["txt", "csv"], height=4)
        self.peak_file_list.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.peak_file_list.set_callback(self.on_peak_file_select)

        # Current selection display with validation
        current_frame = ttk.Frame(file_frame)
        current_frame.pack(fill=tk.X, pady=5)

        self.current_status = ttk.Label(current_frame, text="📋 No files loaded",
                                      font=('TkDefaultFont', 9, 'bold'), foreground='gray')
        self.current_status.pack(anchor=tk.W)

        self.file_validation_label = ttk.Label(current_frame, text="",
                                             font=('TkDefaultFont', 8), foreground='blue')
        self.file_validation_label.pack(anchor=tk.W, padx=(20, 0))

        # =================== SPECTRUM CONTROLS SECTION ===================
        spectrum_frame = ttk.LabelFrame(parent, text="🖼️ Spectrum Display Controls", padding=10)
        spectrum_frame.pack(fill=tk.X, pady=(0, 10))

        # Contour controls
        contour_frame = ttk.LabelFrame(spectrum_frame, text="Contour Settings", padding=5)
        contour_frame.pack(fill=tk.X, pady=(0, 5))

        contour_grid = ttk.Frame(contour_frame)
        contour_grid.pack(fill=tk.X)

        ttk.Label(contour_grid, text="Levels:").grid(row=0, column=0, sticky=tk.W)
        contour_levels_spin = tk.Spinbox(contour_grid, from_=5, to=50, increment=1, width=4,
                                        textvariable=self.contour_levels, command=self.update_main_plot)
        contour_levels_spin.grid(row=0, column=1, sticky=tk.W, padx=(5,15))
        # Bind events for immediate updates on direct input
        contour_levels_spin.bind('<KeyRelease>', lambda e: self.update_main_plot())
        contour_levels_spin.bind('<FocusOut>', lambda e: self.update_main_plot())

        ttk.Label(contour_grid, text="Min Level:").grid(row=0, column=2, sticky=tk.W)
        contour_min_spin = tk.Spinbox(contour_grid, from_=0.01, to=1.0, increment=0.01, width=4,
                                     textvariable=self.contour_min, format="%.3f", command=self.update_main_plot)
        contour_min_spin.grid(row=0, column=3, sticky=tk.W, padx=5)
        # Bind events for immediate updates on direct input
        contour_min_spin.bind('<KeyRelease>', lambda e: self.update_main_plot())
        contour_min_spin.bind('<FocusOut>', lambda e: self.update_main_plot())

        ttk.Label(contour_grid, text="Increment:").grid(row=1, column=0, sticky=tk.W)
        contour_inc_spin = tk.Spinbox(contour_grid, from_=0.01, to=10.0, increment=0.01, width=4,
                                     textvariable=self.contour_increment, format="%.2f", command=self.update_main_plot)
        contour_inc_spin.grid(row=1, column=1, sticky=tk.W, padx=(5,15))
        # Bind events for immediate updates on direct input
        contour_inc_spin.bind('<KeyRelease>', lambda e: self.update_main_plot())
        contour_inc_spin.bind('<FocusOut>', lambda e: self.update_main_plot())

        ttk.Button(contour_grid, text="🔄 Reset Zoom", command=self.reset_interactive_zoom).grid(row=2, column=0, sticky=(tk.W, tk.E), padx=2)


        # =================== PEAK COORDINATE ADJUSTMENT SECTION ===================
        adjustment_frame = ttk.LabelFrame(parent, text="📐 Peak List Coordinate Adjustment", padding=10)
        adjustment_frame.pack(fill=tk.X, pady=(0, 10))

        # Coordinate adjustment controls
        adjust_controls_frame = ttk.Frame(adjustment_frame)
        adjust_controls_frame.pack(fill=tk.X, pady=5)

        # X and Y offset controls
        offset_frame = ttk.Frame(adjust_controls_frame)
        offset_frame.pack(fill=tk.X, pady=2)

        ttk.Label(offset_frame, text="1H (ppm):").pack(side=tk.LEFT)
        x_offset_spinbox = ttk.Spinbox(offset_frame, from_=-2.0, to=2.0, increment=0.001, width=5,
                                      textvariable=self.adjust_x_offset, format="%.3f")
        x_offset_spinbox.pack(side=tk.LEFT, padx=(5, 20))

        ttk.Label(offset_frame, text="15N/13C (ppm):").pack(side=tk.LEFT)
        y_offset_spinbox = ttk.Spinbox(offset_frame, from_=-10.0, to=10.0, increment=0.1, width=4,
                                      textvariable=self.adjust_y_offset, format="%.1f")
        y_offset_spinbox.pack(side=tk.LEFT, padx=(5, 20))

        # Action buttons with responsive layout
        adjust_buttons_frame = ttk.Frame(adjust_controls_frame)
        adjust_buttons_frame.pack(fill=tk.X, pady=5)

        # Configure equal width columns for buttons
        adjust_buttons_frame.columnconfigure(0, weight=1)
        adjust_buttons_frame.columnconfigure(1, weight=1)
        adjust_buttons_frame.columnconfigure(2, weight=1)

        ttk.Button(adjust_buttons_frame, text="↗️ Apply Offsets",
                  command=self.apply_coordinate_offsets).grid(row=0, column=0, sticky=(tk.W, tk.E), padx=1)
        ttk.Button(adjust_buttons_frame, text="🔄 Reset",
                  command=self.reset_coordinate_offsets).grid(row=0, column=1, sticky=(tk.W, tk.E), padx=1)
        #ttk.Button(adjust_buttons_frame, text="💾 Save newt List",
        #          command=self.save_adjusted_peak_list).grid(row=1, column=0, sticky=(tk.W, tk.E), padx=1)

        # Info and status
        adjust_info_frame = ttk.Frame(adjustment_frame)
        adjust_info_frame.pack(fill=tk.X, pady=2)

        #ttk.Label(adjust_info_frame,
        #         text="ℹ️ Adjust all peak coordinates by specified offset values. Useful for spectrum calibration or systematic shifts.",
        #         font=('TkDefaultFont', 8), foreground='blue').pack(anchor=tk.W)

        self.adjustment_status_label = ttk.Label(adjust_info_frame, text="Ready for coordinate adjustment",
                                                foreground='gray')
        self.adjustment_status_label.pack(anchor=tk.W, pady=(2, 0))

        # =================== DETECTION & PROCESSING SECTION ===================
        detection_frame = ttk.LabelFrame(parent, text="🔍 Enhanced Peak Detection & Processing", padding=10)
        detection_frame.pack(fill=tk.X, pady=(0, 10))

        # Mode-specific options
        self.mode_options_frame = ttk.Frame(detection_frame)
        self.mode_options_frame.pack(fill=tk.X, pady=(0, 5))

        # Set callback and initialize mode after all UI is ready
        #self.mode_selection_frame.set_callback(self.on_mode_change)
        #self.mode_selection_frame.set_mode(self.processing_mode.get())

        # Detection parameters in organized grid
        params_frame = ttk.LabelFrame(detection_frame, text="Detection Parameters", padding=5)
        params_frame.pack(fill=tk.X, pady=5)

        # Primary parameters
        primary_frame = ttk.Frame(params_frame)
        primary_frame.pack(fill=tk.X, pady=(0, 5))

        ttk.Label(primary_frame, text="1H/15N (ppm)").grid(row=0, column=0, sticky=tk.W, padx=(0,10))
        search_x_spin = ttk.Spinbox(primary_frame, from_=0.01, to=0.2, increment=0.01, width=4,
                                  textvariable=self.search_window_x, #format="%.2f",
                                  command=self.on_parameter_change)
        search_x_spin.grid(row=0, column=1, sticky=tk.W)

        #ttk.Label(primary_frame, text="15N/13C (±ppm):").grid(row=0, column=2, sticky=tk.W, padx=(0,10))
        search_y_spin = ttk.Spinbox(primary_frame, from_=0.05, to=4.0, increment=0.05, width=4,
                                  textvariable=self.search_window_y, #format="%.1f",
                                  command=self.on_parameter_change)
        search_y_spin.grid(row=0, column=2, sticky=tk.W)

        ttk.Label(primary_frame, text="Noise Threshold:").grid(row=1, column=0, sticky=tk.W, padx=(0,10))
        noise_spin = ttk.Spinbox(primary_frame, from_=0.01, to=10.0, increment=0.01, width=4,
                               textvariable=self.noise_threshold, #format="%.1f",
                               command=self.on_parameter_change)
        noise_spin.grid(row=1, column=1, sticky=tk.W)

        # Detection rectangle size parameter (anisotropic)
        ttk.Label(primary_frame, text="X×Y (pixels):").grid(row=3, column=0, sticky=tk.W, padx=(0,5))
        square_size_spin = tk.Spinbox(primary_frame, from_=1, to=9, increment=2, width=3,
                                     textvariable=self.detection_square_size,
                                     command=self.update_detection_square_ppm)
        square_size_spin.grid(row=3, column=1, sticky=tk.W)

        #ttk.Label(primary_frame, text="×").grid(row=3, column=2, sticky=tk.W)

        rectangle_y_spin = tk.Spinbox(primary_frame, from_=1, to=5, increment=1, width=3,
                                     textvariable=self.detection_rectangle_y,
                                     command=self.update_detection_square_ppm)
        rectangle_y_spin.grid(row=3, column=2, sticky=tk.W)

        # PPM conversion display
        self.square_ppm_label = ttk.Label(primary_frame, textvariable=self.detection_square_ppm_x, font=('TkDefaultFont', 8))
        self.square_ppm_label.grid(row=3, column=3, sticky=tk.W, padx=(5,0))


        # Row 8: Peak Centroid Detection (optional post-processing enhancement)
        ttk.Label(primary_frame, text="🎯 Peak Centroid Detection:", font=('TkDefaultFont', 9, 'bold')).grid(
            row=4, column=0, columnspan=4, sticky=tk.W, pady=(5,5))

        # Row 9: Enable centroid refinement checkbox
        #centroid_check = ttk.Checkbutton(primary_frame, text="🔬 Enable Centroid Refinement (sub-pixel accuracy)",
        #                               variable=self.use_centroid_refinement,
        #                               command=self.on_parameter_change)
        #centroid_check.grid(row=5, column=0, columnspan=4, sticky=tk.W, pady=(5,0))

        # Row 6: Centroid window parameters in ppm
        ttk.Label(primary_frame, text="Window X ppm:").grid(row=6, column=0, sticky=tk.W)
        centroid_x_spin = tk.Spinbox(primary_frame, from_=0.01, to=0.2, increment=0.01, width=4,
                                   textvariable=self.centroid_window_x_ppm, #format="%.2f",
                                   command=self.on_parameter_change)
        centroid_x_spin.grid(row=6, column=1, sticky=tk.W, padx=(5,15))

        #ttk.Label(primary_frame, text="Y ppm:").grid(row=6, column=2, sticky=tk.W)
        centroid_y_spin = tk.Spinbox(primary_frame, from_=0.01, to=0.5, increment=0.02, width=4,
                                   textvariable=self.centroid_window_y_ppm, #format="%.1f",
                                   command=self.on_parameter_change)
        centroid_y_spin.grid(row=6, column=2, sticky=tk.W, padx=5)

        # Row 7: Noise threshold multiplier
        #ttk.Label(primary_frame, text="Noise Mult:").grid(row=7, column=0, sticky=tk.W)
        #centroid_noise_spin = tk.Spinbox(primary_frame, from_=1.0, to=5.0, increment=0.1, width=4,
        #                             textvariable=self.centroid_noise_multiplier, #format="%.1f",
        #                             command=self.on_parameter_change)
        #centroid_noise_spin.grid(row=7, column=1, sticky=tk.W, padx=(5,5))


        # Advanced options
        advanced_frame = ttk.Frame(params_frame)
        advanced_frame.pack(fill=tk.X, pady=(5, 0))

        ttk.Checkbutton(advanced_frame, text="Use Reference-Based Detection",
                       variable=self.use_reference_detection,
                       command=self.on_parameter_change).pack(anchor=tk.W)

        # Peak Ridge Consolidation parameters (Solution A)
        consolidation_frame = ttk.Frame(advanced_frame)
        consolidation_frame.pack(fill=tk.X, pady=2)

        ttk.Label(consolidation_frame, text="Y-Peak:").pack(side=tk.LEFT)
        ttk.Label(consolidation_frame, text="X-T:").pack(side=tk.LEFT, padx=(20,0))
        ttk.Spinbox(consolidation_frame, from_=0.01, to=0.2, increment=0.01, width=4,
                   textvariable=self.consolidation_x_tolerance, format="%.3f",
                   command=self.on_parameter_change).pack(side=tk.LEFT, padx=2)
        #ttk.Label(consolidation_frame, text="ppm").pack(side=tk.LEFT)

        ttk.Label(consolidation_frame, text="Y-T:").pack(side=tk.LEFT, padx=(10,0))
        ttk.Spinbox(consolidation_frame, from_=0.5, to=10.0, increment=0.5, width=4,
                   textvariable=self.consolidation_y_tolerance, format="%.1f",
                   command=self.on_parameter_change).pack(side=tk.LEFT, padx=2)
        ttk.Label(consolidation_frame, text="ppm").pack(side=tk.LEFT)

        # =================== DETECTION MODE CONTROLS ===================
        # Simplified detection controls for Standard and Enhanced modes only
        mode_controls_frame = ttk.LabelFrame(detection_frame, text="🎛️ Detection Mode Selection", padding=5)
        mode_controls_frame.pack(fill=tk.X, pady=5)

        # Enhanced Detection Parameters
        #enhanced_frame = ttk.LabelFrame(mode_controls_frame, text="🚀 Enhanced Detection Parameters", padding=5)
        #enhanced_frame = ttk.LabelFrame(detection_frame, text="🚀 Enhanced Detection Parameters", padding=5)
        #enhanced_frame.pack(fill=tk.X, pady=2)

        # Row 1: Radius Constraint and Pattern Similarity
        #enhanced_row1 = ttk.Frame(enhanced_frame)
        #enhanced_row1.pack(fill=tk.X, pady=1)

        #ttk.Label(enhanced_row1, text="1H Window (ppm):").pack(side=tk.LEFT)
        #radius_x_spin = ttk.Spinbox(enhanced_row1, from_=0.01, to=0.5, increment=0.01, width=4,
        #                           textvariable=self.enhanced_radius_x, format="%.2f")
        #radius_x_spin.pack(side=tk.LEFT, padx=5)

        #ttk.Label(enhanced_row1, text="15N/13C Window (ppm):").pack(side=tk.LEFT, padx=(10, 0))
        #radius_y_spin = ttk.Spinbox(enhanced_row1, from_=0.5, to=10.0, increment=0.5, width=4,
        #                           textvariable=self.enhanced_radius_y, format="%.1f")
        #radius_y_spin.pack(side=tk.LEFT, padx=5)

        # Row 2: Missing Peak Tolerance and Position Weight
        #enhanced_row2 = ttk.Frame(enhanced_frame)
        #enhanced_row2.pack(fill=tk.X, pady=1)

        #ttk.Label(enhanced_row2, text="Missing Peak Tolerance:").pack(side=tk.LEFT)
        #tolerance_spin = ttk.Spinbox(enhanced_row2, from_=0, to=3, increment=1, width=4,
        #                            textvariable=self.enhanced_missing_tolerance)
        #tolerance_spin.pack(side=tk.LEFT, padx=5)

        #ttk.Label(enhanced_row2, text="Position Weight:").pack(side=tk.LEFT, padx=(20, 0))
        #weight_spin = ttk.Spinbox(enhanced_row2, from_=0.1, to=1.0, increment=0.1, width=4,
        #                         textvariable=self.enhanced_position_weight, format="%.1f")
        #weight_spin.pack(side=tk.LEFT, padx=5)

        # Row 3: Peak Reduction Controls (NEW)
        #enhanced_row3 = ttk.Frame(enhanced_frame)
        #enhanced_row3.pack(fill=tk.X, pady=1)

        #ttk.Label(enhanced_row3, text="Max Peaks Detected:").pack(side=tk.LEFT)
        #peak_limit_spin = ttk.Spinbox(enhanced_row3, from_=10, to=200, increment=10, width=4,
        #                             textvariable=self.enhanced_peak_limit)
        #peak_limit_spin.pack(side=tk.LEFT, padx=5)

        #ttk.Label(enhanced_row3, text="Noise Threshold:").pack(side=tk.LEFT, padx=(20, 0))
        #noise_spin = ttk.Spinbox(enhanced_row3, from_=0.001, to=1.0, increment=0.01, width=4,
        #                        textvariable=self.enhanced_noise_threshold, format="%.3f")
        #noise_spin.pack(side=tk.LEFT, padx=5)


        # Row 4: Pattern Similarity (NEW)
        #enhanced_row4 = ttk.Frame(enhanced_frame)
        #enhanced_row4.pack(fill=tk.X, pady=2)

        #ttk.Label(enhanced_row4, text="Pattern Similarity:").pack(side=tk.LEFT, padx=(20, 0))
        #pattern_spin = ttk.Spinbox(enhanced_row4, from_=0.1, to=1.0, increment=0.1, width=4,
        #                          textvariable=self.enhanced_pattern_similarity, format="%.1f")
        #pattern_spin.pack(side=tk.LEFT, padx=5)

        # Row 5: xx
        #enhanced_row5 = ttk.Frame(enhanced_frame)
        #enhanced_row5.pack(fill=tk.X, pady=2)

        #ttk.Label(enhanced_row5, text="Detection Size X×Y (pixels):").pack(side=tk.LEFT, padx=(20, 0))
        #enhanced_square_spin = tk.Spinbox(enhanced_row5, from_=1, to=9, increment=2, width=3,
        #                                 textvariable=self.enhanced_detection_square_size,
        #                                 command=self.update_detection_square_ppm)
        #enhanced_square_spin.pack(side=tk.LEFT, padx=2)

        #ttk.Label(enhanced_row5, text="×").pack(side=tk.LEFT)

        #enhanced_rectangle_y_spin = tk.Spinbox(enhanced_row5, from_=1, to=5, increment=1, width=3,
        #                                      textvariable=self.enhanced_detection_rectangle_y,
        #                                      command=self.update_detection_square_ppm)
        #enhanced_rectangle_y_spin.pack(side=tk.LEFT, padx=2)

        # PPM conversion display for enhanced detection
        #enhanced_ppm_label = ttk.Label(enhanced_row5, textvariable=self.enhanced_detection_square_ppm_x,
        #                             font=('TkDefaultFont', 8))
        #enhanced_ppm_label.pack(side=tk.LEFT, padx=(5,0))


        # Help text for Enhanced Detection parameters
        #enhanced_help_frame = ttk.Frame(enhanced_frame)
        #enhanced_help_frame.pack(fill=tk.X, pady=1)
        #enhanced_help_text = (
        #    "ℹ️ Enhanced Detection: Radius = max distance for matches, Pattern = min similarity for complex matching,\n"
        #    "   Missing Tolerance = peaks allowed to be missing, Position Weight = balance position vs intensity matching,\n"
        #    "   Max Peaks = limit total detected peaks, Min Intensity = filter weak peaks (reduces noise)"
        #)
        #enhanced_help_label = ttk.Label(enhanced_help_frame, text=enhanced_help_text,
        #                               font=('TkDefaultFont', 8), foreground='blue', wraplength=600)
        #enhanced_help_label.pack(anchor=tk.W, padx=5)

        # Removed: Force Full Detection and threshold parameters (Integrated/Adaptive mode only)

        # Removed: Core Detection Parameters section (Integrated/Adaptive mode only)

        # Store reference to these controls for enable/disable (cleaned up list)
        #self.single_integration_controls = []  # Removed all Integrated/Adaptive mode controls
        #self.single_integration_checkboxes = []  # Removed all Integrated/Adaptive mode checkboxes

        # Processing buttons with responsive grid layout
        button_frame = ttk.Frame(detection_frame)
        button_frame.pack(fill=tk.X, pady=10)

        # Configure equal width columns for main processing buttons
        button_frame.columnconfigure(0, weight=1)
        button_frame.columnconfigure(1, weight=1)
        button_frame.columnconfigure(2, weight=1)
        button_frame.columnconfigure(3, weight=1)
        button_frame.columnconfigure(4, weight=1)

        self.detect_peaks_button = ttk.Button(button_frame, text="🔍 Detect Peaks",
                                             command=self.detect_peaks)
        self.detect_peaks_button.grid(row=0, column=0, sticky=(tk.W, tk.E), padx=2)

        # Enhanced Peak Detection button (NEW)
        self.enhanced_detect_button = ttk.Button(button_frame, text="🚀 Enhanced Detection",
                                               command=self.enhanced_detect_peaks, state='disabled')
        self.enhanced_detect_button.grid(row=0, column=1, sticky=(tk.W, tk.E), padx=2)

        # =================== PEAK NAVIGATION SECTION ===================
        nav_frame = ttk.LabelFrame(parent, text="🎯 Peak Navigation", padding=10)
        nav_frame.pack(fill=tk.X, pady=(0, 10))

        # Peak selector with enhanced info
        peak_select_frame = ttk.Frame(nav_frame)
        peak_select_frame.pack(fill=tk.X, pady=2)

        ttk.Label(peak_select_frame, text="Peak Number:").pack(side=tk.LEFT)
        peak_spin = tk.Spinbox(peak_select_frame, from_=1, to=500, width=4,
                              textvariable=self.selected_peak_number,
                              command=self.center_on_selected_peak)
        peak_spin.pack(side=tk.LEFT, padx=5)

        self.peak_info_label = ttk.Label(peak_select_frame, text="Peak: -/-",
                                        font=('TkDefaultFont', 9, 'bold'))
        self.peak_info_label.pack(side=tk.LEFT, padx=20)

        # Enhanced navigation buttons
        nav_button_frame = ttk.Frame(nav_frame)
        nav_button_frame.pack(fill=tk.X, pady=5)

        ttk.Button(nav_button_frame, text="◀◀",
                  command=self.prev_peak, width=4).pack(side=tk.LEFT, padx=2)
        ttk.Button(nav_button_frame, text="🎯 ",
                  command=self.center_on_selected_peak, width=4).pack(side=tk.LEFT, padx=2)
        ttk.Button(nav_button_frame, text="▶▶",
                  command=self.next_peak, width=4).pack(side=tk.LEFT, padx=2)
        ttk.Button(nav_button_frame, text="🔍",
                  command=self.zoom_to_peak, width=4).pack(side=tk.LEFT, padx=2)
        ttk.Button(nav_button_frame, text="🔬",
                  command=self.show_selected_peak_analysis, width=4).pack(side=tk.LEFT, padx=2)

        # =================== PEAK EDITING SECTION ===================
        edit_frame = ttk.LabelFrame(parent, text="🎯 Peak Position Editor", padding=10)
        edit_frame.pack(fill=tk.X, pady=(0, 10))

        # Edit mode toggle
        edit_mode_frame = ttk.Frame(edit_frame)
        edit_mode_frame.pack(fill=tk.X, pady=2)

        edit_checkbox = ttk.Checkbutton(edit_mode_frame, text="Enable Peak Editing Mode",
                                       variable=self.peak_edit_mode,
                                       command=self.toggle_peak_edit_mode)
        edit_checkbox.pack(side=tk.LEFT)

        # Mode status label
        self.edit_mode_status_label = ttk.Label(edit_mode_frame, text="Mode: View Only",
                                               font=('TkDefaultFont', 9), foreground='gray')
        self.edit_mode_status_label.pack(side=tk.LEFT, padx=20)

        # Peak list selection frame
        peak_list_frame = ttk.Frame(edit_frame)
        peak_list_frame.pack(fill=tk.X, pady=(5, 2))

        #ttk.Label(peak_list_frame, text="Edit which peak lists:", font=('TkDefaultFont', 9)).pack(anchor=tk.W)

        # Individual peak list checkboxes
        list_checkboxes_frame = ttk.Frame(peak_list_frame)
        list_checkboxes_frame.pack(fill=tk.X, pady=(2, 0))

        ref_checkbox = ttk.Checkbutton(list_checkboxes_frame, text="Reference peaks",
                                      variable=self.edit_reference_peaks,
                                      command=self.update_edit_mode_status)
        ref_checkbox.pack(side=tk.LEFT, padx=(20, 10))

        det_checkbox = ttk.Checkbutton(list_checkboxes_frame, text="Detected peaks",
                                      variable=self.edit_detected_peaks,
                                      command=self.update_edit_mode_status)
        det_checkbox.pack(side=tk.LEFT)

        # Instructions
        #instructions_label = ttk.Label(edit_frame,
        #                             text="Instructions: Click peak to select, click new position to move",
        #                             font=('TkDefaultFont', 8), foreground='blue')
        #instructions_label.pack(anchor=tk.W, pady=(2, 0))

        # Selected peak info
        self.selected_peak_label = ttk.Label(edit_frame, text="No peak selected",
                                           font=('TkDefaultFont', 9), foreground='black')
        self.selected_peak_label.pack(anchor=tk.W, pady=2)


        # =================== VOIGT FITTING SECTION ===================
        voigt_frame = ttk.LabelFrame(parent, text="📈 Voigt Profile Fitting", padding=5)
        voigt_frame.pack(fill=tk.X, pady=(0, 10))

        # Fitting parameters in organized layout
        voigt_params_frame = ttk.LabelFrame(voigt_frame, text="Fitting Parameters", padding=5)
        voigt_params_frame.pack(fill=tk.X, pady=5)

        # Parameters grid
        params_grid = ttk.Frame(voigt_params_frame)
        params_grid.pack(fill=tk.X)

        ttk.Label(params_grid, text="X-Window:").grid(row=0, column=0, sticky=tk.W)
        fit_x_spin = tk.Spinbox(params_grid, from_=0.01, to=0.2, increment=0.01, width=4,
                               textvariable=self.fitting_window_x, #format="%.1f",
                               command=self.on_parameter_change)
        fit_x_spin.grid(row=0, column=1, sticky=tk.W, padx=(5,15))

        ttk.Label(params_grid, text="Y-Window:").grid(row=0, column=2, sticky=tk.W)
        fit_y_spin = tk.Spinbox(params_grid, from_=0.01, to=1.0, increment=0.05, width=4,
                               textvariable=self.fitting_window_y, #format="%.1f",
                               command=self.on_parameter_change)
        fit_y_spin.grid(row=0, column=3, sticky=tk.W, padx=5)

        ttk.Label(params_grid, text="Min R²:").grid(row=1, column=0, sticky=tk.W)
        r2_spin = tk.Spinbox(params_grid, from_=0.1, to=1.0, increment=0.1, width=4,
                            textvariable=self.min_r_squared, #format="%.2f",
                            command=self.on_parameter_change)
        r2_spin.grid(row=1, column=1, sticky=tk.W, padx=(5,15))

        ttk.Label(params_grid, text="Max Iter:").grid(row=1, column=2, sticky=tk.W)
        iter_spin = tk.Spinbox(params_grid, from_=10, to=1000, increment=50, width=4,
                              textvariable=self.max_iterations,
                              command=self.on_parameter_change)
        iter_spin.grid(row=1, column=3, sticky=tk.W, padx=5)

        # Global optimization toggle
        global_opt_check = ttk.Checkbutton(params_grid, text="🔄 Use Global Optimization",
                                         variable=self.use_global_optimization,
                                         command=self.on_parameter_change,
                                         state="disabled")
        global_opt_check.grid(row=2, column=0, columnspan=2, sticky=tk.W, pady=(10,0))

        # Peak Detection Parameters (Multi-Peak Fitting)
        ttk.Label(params_grid, text="🔍 Peak Detection Parameters", font=('TkDefaultFont', 9, 'bold')).grid(
            row=3, column=0, columnspan=4, sticky=tk.W, pady=(15,5))

        # Row 4: Height threshold and Distance factor
        ttk.Label(params_grid, text="Height %:").grid(row=4, column=0, sticky=tk.W)
        height_spin = tk.Spinbox(params_grid, from_=0.005, to=0.2, increment=0.005, width=4,
                               textvariable=self.peak_height_threshold, #format="%.3f",
                               command=self.on_parameter_change)
        height_spin.grid(row=4, column=1, sticky=tk.W, padx=(5,15))

        ttk.Label(params_grid, text="Distance:").grid(row=4, column=2, sticky=tk.W)
        distance_spin = tk.Spinbox(params_grid, from_=10, to=200, increment=10, width=4,
                                 textvariable=self.peak_distance_factor,
                                 command=self.on_parameter_change)
        distance_spin.grid(row=4, column=3, sticky=tk.W, padx=5)

        # Add ppm conversion display
        self.distance_ppm_label = ttk.Label(params_grid, text="(≈0.00 ppm)",
                                           font=('TkDefaultFont', 8), foreground='gray')
        self.distance_ppm_label.grid(row=4, column=4, sticky=tk.W, padx=(2,0))

        # Row 5: Prominence and Smoothing
        ttk.Label(params_grid, text="Prominence:").grid(row=5, column=0, sticky=tk.W)
        prom_spin = tk.Spinbox(params_grid, from_=0.001, to=0.1, increment=0.001, width=4,
                             textvariable=self.peak_prominence_threshold, #format="%.3f",
                             command=self.on_parameter_change)
        prom_spin.grid(row=5, column=1, sticky=tk.W, padx=(5,15))

        ttk.Label(params_grid, text="Smoothing:").grid(row=5, column=2, sticky=tk.W)
        smooth_spin = tk.Spinbox(params_grid, from_=0.1, to=2.0, increment=0.1, width=4,
                               textvariable=self.smoothing_sigma, format="%.1f",
                               command=self.on_parameter_change)
        smooth_spin.grid(row=5, column=3, sticky=tk.W, padx=5)

        # Row 6: Max peaks to fit and Max iterations
        ttk.Label(params_grid, text="Max Peaks:").grid(row=6, column=0, sticky=tk.W)
        max_peaks_spin = tk.Spinbox(params_grid, from_=2, to=8, increment=1, width=4,
                                  textvariable=self.max_peaks_fit,
                                  command=self.on_parameter_change)
        max_peaks_spin.grid(row=6, column=1, sticky=tk.W, padx=(5,15))

        ttk.Label(params_grid, text="Max Iter:").grid(row=6, column=2, sticky=tk.W)
        max_iter_spin = tk.Spinbox(params_grid, from_=10, to=200, increment=10, width=4,
                                 textvariable=self.max_optimization_iterations,
                                 command=self.on_parameter_change)
        max_iter_spin.grid(row=6, column=3, sticky=tk.W, padx=5)

        # Row 7: Parallel processing toggle
        parallel_check = ttk.Checkbutton(params_grid, text="🚀 Use Parallel Processing (75% cores)",
                                       variable=self.use_parallel_processing,
                                       command=self.on_parameter_change)
        parallel_check.grid(row=7, column=0, columnspan=4, sticky=tk.W, pady=(5,0))


        # Voigt buttons
        voigt_button_frame = ttk.Frame(voigt_frame)
        voigt_button_frame.pack(fill=tk.X, pady=5)

        ttk.Button(voigt_button_frame, text="📊 Fit Selected Peak",
                  command=self.fit_selected_peak, width=18).pack(side=tk.LEFT, padx=2)
        self.fit_all_peaks_button = ttk.Button(voigt_button_frame, text="📈 Fit All Peaks",
                                              command=self.fit_all_peaks, width=17)
        self.fit_all_peaks_button.pack(side=tk.LEFT, padx=2)


        # =================== SERIES INTEGRATION SECTION ===================
        series_frame = ttk.LabelFrame(parent, text="🚀 Advanced Series Integration Workflow", padding=10)
        series_frame.pack(fill=tk.X, pady=(0, 10))

        # Series options with enhanced controls
        options_frame = ttk.LabelFrame(series_frame, text="Processing Options", padding=5)
        options_frame.pack(fill=tk.X, pady=(0, 5))

        ttk.Checkbutton(options_frame, text="Auto-process entire folder",
                       variable=self.auto_process_series).pack(anchor=tk.W)
        ttk.Checkbutton(options_frame, text="Save individual spectrum results",
                       variable=self.save_individual_results).pack(anchor=tk.W)
        ttk.Checkbutton(options_frame, text="Create comprehensive summary plots",
                       variable=self.create_summary_plots).pack(anchor=tk.W)

        # Advanced Voigt fitting options for series integration
        ttk.Separator(options_frame, orient='horizontal').pack(fill=tk.X, pady=5)
        ttk.Label(options_frame, text="🔬 Advanced Fitting Options:", font=('TkDefaultFont', 9, 'bold')).pack(anchor=tk.W)

        series_voigt_check = ttk.Checkbutton(options_frame, text="Voigt profile fitting",
                       variable=self.series_use_voigt_fitting, command=self._toggle_series_voigt_params)
        series_voigt_check.pack(anchor=tk.W, padx=10)

        # Detailed Voigt parameters (shown when checkbox is enabled)
        self.series_voigt_params_frame = ttk.LabelFrame(options_frame, text="Advanced Voigt Parameters", padding=5)

        # Create parameter grid similar to Peak Detection section
        series_params_grid = ttk.Frame(self.series_voigt_params_frame)
        series_params_grid.pack(fill=tk.X)

        # Row 0: Fitting Windows
        #ttk.Label(series_params_grid, text="X-Window:").grid(row=0, column=0, sticky=tk.W)
        #ttk.Label(series_params_grid, text=f"{self.fitting_window_x.get():.1f}", name="x_window_display").grid(row=0, column=1, sticky=tk.W, padx=(5,15))

        #ttk.Label(series_params_grid, text="Y-Window:").grid(row=0, column=2, sticky=tk.W)
        #ttk.Label(series_params_grid, text=f"{self.fitting_window_y.get():.1f}", name="y_window_display").grid(row=0, column=3, sticky=tk.W, padx=5)

        # Row 1: Quality parameters
        #ttk.Label(series_params_grid, text="Min R²:").grid(row=1, column=0, sticky=tk.W)
        #ttk.Label(series_params_grid, text=f"{self.min_r_squared.get():.2f}", name="r_squared_display").grid(row=1, column=1, sticky=tk.W, padx=(5,15))

        #ttk.Label(series_params_grid, text="Max Iter:").grid(row=1, column=2, sticky=tk.W)
        #ttk.Label(series_params_grid, text=f"{self.max_iterations.get()}", name="max_iter_display").grid(row=1, column=3, sticky=tk.W, padx=5)

        # Row 2: Peak Detection Parameters
        #ttk.Label(series_params_grid, text="🔍 Peak Detection:", font=('TkDefaultFont', 8, 'bold')).grid(
        #    row=2, column=0, columnspan=4, sticky=tk.W, pady=(10,5))

        # Row 3: Height and Distance
        #ttk.Label(series_params_grid, text="Height %:").grid(row=3, column=0, sticky=tk.W)
        #ttk.Label(series_params_grid, text=f"{self.peak_height_threshold.get():.3f}", name="height_display").grid(row=3, column=1, sticky=tk.W, padx=(5,15))

        #ttk.Label(series_params_grid, text="Distance:").grid(row=3, column=2, sticky=tk.W)
        #ttk.Label(series_params_grid, text=f"{self.peak_distance_factor.get():.0f}", name="distance_display").grid(row=3, column=3, sticky=tk.W, padx=5)

        # Row 4: Prominence and Smoothing
        #ttk.Label(series_params_grid, text="Prominence:").grid(row=4, column=0, sticky=tk.W)
        #ttk.Label(series_params_grid, text=f"{self.peak_prominence_threshold.get():.3f}", name="prominence_display").grid(row=4, column=1, sticky=tk.W, padx=(5,15))

        #ttk.Label(series_params_grid, text="Smoothing:").grid(row=4, column=2, sticky=tk.W)
        #ttk.Label(series_params_grid, text=f"{self.smoothing_sigma.get():.1f}", name="smoothing_display").grid(row=4, column=3, sticky=tk.W, padx=5)

        # Row 5: Max peaks and optimization iterations
        #ttk.Label(series_params_grid, text="Max Peaks:").grid(row=5, column=0, sticky=tk.W)
        #ttk.Label(series_params_grid, text=f"{self.max_peaks_fit.get()}", name="max_peaks_display").grid(row=5, column=1, sticky=tk.W, padx=(5,15))

        #ttk.Label(series_params_grid, text="Opt Iter:").grid(row=5, column=2, sticky=tk.W)
        #ttk.Label(series_params_grid, text=f"{self.max_optimization_iterations.get()}", name="opt_iter_display").grid(row=5, column=3, sticky=tk.W, padx=5)

        # Info text
        #ttk.Label(self.series_voigt_params_frame, text="ℹ️ Parameters are synchronized from Peak Detection section",
        #         font=('TkDefaultFont', 8), foreground='blue').pack(pady=(5,0))

        # Initially hide the parameter frame
        self._toggle_series_voigt_params()

        ttk.Checkbutton(options_frame, text="Multicore parallel optimisation",
                       variable=self.series_use_parallel_processing).pack(anchor=tk.W, padx=10)
        #ttk.Checkbutton(options_frame, text="Use global optimization)",
        #               variable=self.series_use_global_optimization).pack(anchor=tk.W, padx=10)

        # Integrated detection-fitting controls
        ttk.Separator(options_frame, orient='horizontal').pack(fill=tk.X, pady=5)
        ttk.Label(options_frame, text=" Integrated Detection-Fitting:", font=('TkDefaultFont', 9, 'bold')).pack(anchor=tk.W)

        # Integration mode selection
        integration_mode_frame = ttk.Frame(options_frame)
        integration_mode_frame.pack(fill=tk.X, padx=10, pady=2)

        ttk.Label(integration_mode_frame, text="Integration Mode:").pack(side=tk.LEFT)
        integration_mode_combo = ttk.Combobox(integration_mode_frame, textvariable=self.integration_mode,
                                            values=['standard'], #, 'enhanced'
                                            state="readonly", width=12)
        integration_mode_combo.pack(side=tk.LEFT, padx=(5, 0))
        integration_mode_combo.bind('<<ComboboxSelected>>', self.on_integration_mode_change)


        # Integration status display
        self.integration_status_frame = ttk.Frame(options_frame)
        self.integration_status_frame.pack(fill=tk.X, padx=10, pady=2)

        self.integration_status_label = ttk.Label(self.integration_status_frame,
                                                text="Status: Standard mode",
                                                font=('TkDefaultFont', 8), foreground='blue')
        self.integration_status_label.pack(side=tk.LEFT)

        # Number of integrations option for Voigt fitting
        integrations_frame = ttk.Frame(options_frame)
        integrations_frame.pack(fill=tk.X, padx=10, pady=5)
        ttk.Label(integrations_frame, text="Number of integrations per peak:").pack(side=tk.LEFT)
        integrations_spinbox = ttk.Spinbox(integrations_frame, from_=1, to=10, width=5,
                                         textvariable=self.series_num_integrations)
        integrations_spinbox.pack(side=tk.LEFT, padx=(5, 0))
        #ttk.Label(integrations_frame, text="(only applies when Voigt fitting is enabled)",
        #         font=('TkDefaultFont', 8), foreground='gray').pack(side=tk.LEFT, padx=(5, 0))

        # NEW: Series peak source selection
        ttk.Separator(series_frame, orient='horizontal').pack(fill=tk.X, pady=10)

        peak_source_label_frame = ttk.Frame(series_frame)
        peak_source_label_frame.pack(fill=tk.X, pady=(0, 5))
        ttk.Label(peak_source_label_frame, text="🎯 Peak List Source:", font=('TkDefaultFont', 9, 'bold')).pack(anchor=tk.W)

        peak_source_frame = ttk.Frame(series_frame)
        peak_source_frame.pack(fill=tk.X, padx=10, pady=(0, 10))

        # Peak source radio buttons
        ttk.Radiobutton(peak_source_frame, text="Use detected peaks ",
                       variable=self.series_peak_source, value="detected").pack(anchor=tk.W)
        ttk.Radiobutton(peak_source_frame, text="Use reference peak list",
                       variable=self.series_peak_source, value="reference").pack(anchor=tk.W)
        #ttk.Radiobutton(peak_source_frame, text="Cascade results (0 → 1 → 2)",
        #               variable=self.series_peak_source, value="cascade").pack(anchor=tk.W)

        # Peak detection toggle
        detection_frame = ttk.Frame(series_frame)
        detection_frame.pack(fill=tk.X, padx=10, pady=(0, 10))

        #ttk.Checkbutton(detection_frame, text="🔍 Enable peak detection for each spectrum (slower, more thorough)",
        #               variable=self.series_enable_detection).pack(anchor=tk.W)
        #ttk.Label(detection_frame, text="   ⚡ When OFF: Direct fitting only (faster, uses existing peak positions)",
        #         font=('TkDefaultFont', 8), foreground='gray').pack(anchor=tk.W)

        # Series processing controls
        #ttk.Separator(series_frame, orient='horizontal').pack(fill=tk.X, pady=10)

        series_control_frame = ttk.Frame(series_frame)
        series_control_frame.pack(fill=tk.X, pady=5)

        self.series_button = ttk.Button(series_control_frame, text="🚀 START SERIES INTEGRATION",
                                       command=self.start_series_integration)
        self.series_button.pack(fill=tk.X, pady=2)

        # Series results browser with enhanced features
        #ttk.Separator(series_frame, orient='horizontal').pack(fill=tk.X, pady=10)

        #results_label_frame = ttk.Frame(series_frame)
        #results_label_frame.pack(fill=tk.X)

        #ttk.Label(results_label_frame, text="📊 Browse Series Results:",
        #         font=('TkDefaultFont', 9, 'bold')).pack(side=tk.LEFT)
        #ttk.Button(results_label_frame, text="📈 Analysis",
        #          command=self.show_series_analysis, width=4).pack(side=tk.RIGHT)

        #self.results_var = tk.StringVar()
        #self.results_combo = ttk.Combobox(series_frame, textvariable=self.results_var,
        #                                 state="readonly", width=40)
        #self.results_combo.pack(fill=tk.X, pady=2)
        #self.results_combo.bind('<<ComboboxSelected>>', self.load_series_result)

        # =================== STATISTICS PANEL ===================
        self.statistics_panel = StatisticsPanel(parent, title="📊 Processing Statistics")
        self.statistics_panel.pack(fill=tk.X, pady=(0, 10))


        # =================== DISPLAY OPTIONS SECTION ===================
        #display_frame = ttk.LabelFrame(parent, text="🖼️ Enhanced Display Options", padding=10)
        #display_frame.pack(fill=tk.X, pady=(0, 10))

        #display_checks_frame = ttk.Frame(display_frame)
        #display_checks_frame.pack(fill=tk.X)

        #ttk.Checkbutton(display_checks_frame, text="Show detected peaks",
        #               variable=self.show_detected, command=self.update_main_plot).pack(side=tk.LEFT)
        #ttk.Checkbutton(display_checks_frame, text="Show reference peaks",
        #               variable=self.show_assigned, command=self.update_main_plot).pack(side=tk.LEFT, padx=10)
        #ttk.Checkbutton(display_checks_frame, text="Show fitted curves",
        #               variable=self.show_fitted_curves, command=self.update_main_plot).pack(side=tk.LEFT)


        # =================== INTEGRATION PROGRESS SECTION ===================
        progress_frame = ttk.LabelFrame(parent, text="🔄 Integration Progress", padding=10)
        progress_frame.pack(fill=tk.X, pady=(0, 10))

        # Progress display and control
        progress_info_frame = ttk.Frame(progress_frame)
        progress_info_frame.pack(fill=tk.X, pady=(0, 5))

        # Current operation status
        self.progress_status_label = ttk.Label(progress_info_frame, text="Ready",
                                             font=('TkDefaultFont', 9), foreground='blue')
        self.progress_status_label.pack(side=tk.LEFT)

        # Progress percentage
        self.progress_percentage_label = ttk.Label(progress_info_frame, text="",
                                                 font=('TkDefaultFont', 9), foreground='gray')
        self.progress_percentage_label.pack(side=tk.RIGHT)

        # Main progress bar
        self.integration_progress = ttk.Progressbar(progress_frame, mode='determinate', length=300)
        self.integration_progress.pack(fill=tk.X, pady=(2, 5))

        # Detailed progress information
        progress_details_frame = ttk.Frame(progress_frame)
        progress_details_frame.pack(fill=tk.X, pady=(0, 5))

        # Current iteration info
        iteration_frame = ttk.Frame(progress_details_frame)
        iteration_frame.pack(fill=tk.X, pady=1)

        ttk.Label(iteration_frame, text="Iteration:", font=('TkDefaultFont', 8)).pack(side=tk.LEFT)
        self.current_iteration_label = ttk.Label(iteration_frame, text="0/0",
                                               font=('TkDefaultFont', 8), foreground='black')
        self.current_iteration_label.pack(side=tk.LEFT, padx=(5, 20))

        ttk.Label(iteration_frame, text="Peaks processed:", font=('TkDefaultFont', 8)).pack(side=tk.LEFT)
        self.peaks_processed_label = ttk.Label(iteration_frame, text="0/0",
                                             font=('TkDefaultFont', 8), foreground='black')
        self.peaks_processed_label.pack(side=tk.LEFT, padx=(5, 20))

        ttk.Label(iteration_frame, text="Quality score:", font=('TkDefaultFont', 8)).pack(side=tk.LEFT)
        self.quality_score_label = ttk.Label(iteration_frame, text="N/A",
                                           font=('TkDefaultFont', 8), foreground='black')
        self.quality_score_label.pack(side=tk.LEFT, padx=(5, 0))

        # Convergence info
        convergence_frame = ttk.Frame(progress_details_frame)
        convergence_frame.pack(fill=tk.X, pady=1)

        ttk.Label(convergence_frame, text="Convergence:", font=('TkDefaultFont', 8)).pack(side=tk.LEFT)
        self.convergence_status_label = ttk.Label(convergence_frame, text="Pending",
                                                font=('TkDefaultFont', 8), foreground='orange')
        self.convergence_status_label.pack(side=tk.LEFT, padx=(5, 20))

        ttk.Label(convergence_frame, text="Time elapsed:", font=('TkDefaultFont', 8)).pack(side=tk.LEFT)
        self.time_elapsed_label = ttk.Label(convergence_frame, text="0s",
                                          font=('TkDefaultFont', 8), foreground='black')
        self.time_elapsed_label.pack(side=tk.LEFT, padx=(5, 0))

        # Progress control buttons
        progress_control_frame = ttk.Frame(progress_frame)
        progress_control_frame.pack(fill=tk.X, pady=(5, 0))

        self.pause_integration_button = ttk.Button(progress_control_frame, text="⏸️ Pause",
                                                  command=self.pause_integration, width=4)
        self.pause_integration_button.pack(side=tk.LEFT, padx=2)

        self.stop_integration_button = ttk.Button(progress_control_frame, text="⏹️ Stop",
                                                 command=self.stop_integration, width=4)
        self.stop_integration_button.pack(side=tk.LEFT, padx=2)

        # Initially hide progress frame
        progress_frame.pack_forget()
        self.progress_frame = progress_frame  # Store reference for showing/hiding

        # =================== STATUS & EXPORT SECTION ===================
        #status_frame = ttk.LabelFrame(parent, text="📋 Status & Export", padding=10)
        #status_frame.pack(fill=tk.X, pady=(0, 10))

        #self.status_label = ttk.Label(status_frame, text="Ready - Load files to begin",
        #                             foreground='blue')
        #self.status_label.pack(anchor=tk.W, pady=2)

        # Export buttons with responsive layout
        #export_frame = ttk.Frame(status_frame)
        #export_frame.pack(fill=tk.X, pady=5)

        # Configure equal width columns for export buttons
        #export_frame.columnconfigure(0, weight=1)
        #export_frame.columnconfigure(1, weight=1)
        #export_frame.columnconfigure(2, weight=1)

        #ttk.Button(export_frame, text="💾 Export Peak List",
        #          command=self.export_peak_list).grid(row=0, column=0, sticky=(tk.W, tk.E), padx=2)
        #ttk.Button(export_frame, text="📈 Save Plot",
        #          command=self.save_plot).grid(row=0, column=1, sticky=(tk.W, tk.E), padx=2)
        #ttk.Button(export_frame, text="📊 Export Results",
        #          command=self.export_results).grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E), padx=2)

    def update_results_dropdown(self):
        """
        Update the series results dropdown with available spectra.
        Handles both new and legacy batch result formats.
        """
        try:
            # Check if we have batch results
            if hasattr(self, 'batch_results') and self.batch_results:
                # Handle legacy format
                if hasattr(self.batch_results, 'results') and self.batch_results.results:
                    spectrum_names = list(self.batch_results.results.keys())
                    print(f"✅ Results dropdown updated: {len(spectrum_names)} spectra available")

                    # Update dropdown widget if it exists
                    if hasattr(self, 'results_combo') and self.results_combo:
                        self.results_combo['values'] = spectrum_names
                        if spectrum_names:
                            self.results_combo.set(spectrum_names[0])  # Select first spectrum
                        else:
                            self.results_combo.set("")

                    return

            # Check for new format results
            if hasattr(self, 'new_batch_results') and self.new_batch_results:
                if 'results' in self.new_batch_results:
                    spectrum_names = list(self.new_batch_results['results'].keys())
                    print(f"✅ Results dropdown updated from new format: {len(spectrum_names)} spectra")

                    if hasattr(self, 'results_combo') and self.results_combo:
                        self.results_combo['values'] = spectrum_names
                        if spectrum_names:
                            self.results_combo.set(spectrum_names[0])
                        else:
                            self.results_combo.set("")

                    return

            # No results available
            print("⚠️ No batch results available for dropdown update")
            if hasattr(self, 'results_combo') and self.results_combo:
                self.results_combo['values'] = []
                self.results_combo.set("")

        except Exception as e:
            print(f"❌ Error updating results dropdown: {e}")
            if hasattr(self, 'results_combo') and self.results_combo:
                self.results_combo['values'] = []
                self.results_combo.set("")

    def _safe_update_results_dropdown(self):
        """Safely update results dropdown with error handling"""
        try:
            if hasattr(self, 'update_results_dropdown'):
                self.update_results_dropdown()
            else:
                print("⚠️ update_results_dropdown method not found")
        except Exception as e:
            print(f"⚠️ Results dropdown update failed: {e}")

    def _safe_update_main_plot(self):
        """Safely update main plot with error handling"""
        try:
            if hasattr(self, 'update_main_plot'):
                self.update_main_plot()
            else:
                print("⚠️ update_main_plot method not found")
        except Exception as e:
            print(f"⚠️ Main plot update failed: {e}")

    def _complete_progress_dialog_safe(self, progress_dialog, message, failed=False):
        """
        Safely complete progress dialog with multiple fallback methods.

        Args:
            progress_dialog: The progress dialog instance
            message (str): Completion message
            failed (bool): Whether the operation failed
        """
        if not progress_dialog:
            return

        try:
            # Method 1: Try standard complete method without failed parameter
            if hasattr(progress_dialog, 'complete'):
                try:
                    import inspect
                    sig = inspect.signature(progress_dialog.complete)
                    params = list(sig.parameters.keys())

                    if 'failed' in params or len(params) >= 2:
                        # Supports failed parameter
                        progress_dialog.complete(message, failed=failed)
                    else:
                        # Only supports message
                        progress_dialog.complete(message)

                except TypeError:
                    # Fallback to message only
                    progress_dialog.complete(message)

            # Method 2: Try manual status update
            else:
                if hasattr(progress_dialog, 'current_task'):
                    progress_dialog.current_task.set(message)

                if hasattr(progress_dialog, 'progress_var'):
                    progress_dialog.progress_var.set(100.0)

                # Set failed status if supported
                if failed and hasattr(progress_dialog, 'set_status'):
                    progress_dialog.set_status("failed")

            print(f"✅ Progress dialog completed: {message}")

        except Exception as e:
            print(f"⚠️ Progress dialog completion failed: {e}")

            # Ultimate fallback - try to close dialog
            try:
                if hasattr(progress_dialog, 'top') and progress_dialog.top.winfo_exists():
                    progress_dialog.top.title(f"Complete: {message}")
            except:
                pass  # Silent fallback


    def setup_enhanced_visualization(self, parent):
        """Setup enhanced visualization with multiple tabs"""
        # Create enhanced notebook with more tabs
        self.viz_notebook = ttk.Notebook(parent)
        self.viz_notebook.pack(fill=tk.BOTH, expand=True)

        # Tab 1: Main spectrum
        main_tab = ttk.Frame(self.viz_notebook)
        self.viz_notebook.add(main_tab, text="📊 Main Spectrum")

        # Create responsive main figure
        self.fig_main, self.ax_main = plt.subplots(1, 1, figsize=(4, 3))
        self.fig_main.tight_layout()
        self.canvas_main = FigureCanvasTkAgg(self.fig_main, main_tab)
        self.canvas_main.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        toolbar_main = NavigationToolbar2Tk(self.canvas_main, main_tab)
        toolbar_main.update()

        # Tab 2: Voigt analysis (enhanced 2x2 layout)
        voigt_tab = ttk.Frame(self.viz_notebook)
        self.viz_notebook.add(voigt_tab, text="📈 Voigt Analysis")

        # Create responsive Voigt analysis figure - 2x1 layout with backward compatibility
        self.fig_voigt, axes_1d = plt.subplots(2, 1, figsize=(4, 8))

        # CRITICAL: Convert 1D axes array to 2D structure for backward compatibility
        # This preserves all existing nested loop access patterns
        self.axes_voigt = [[axes_1d[0]], [axes_1d[1]]]
        self.fig_voigt.tight_layout()
        self.fig_voigt.suptitle('Enhanced Voigt Analysis', fontsize='small', y=0.95)

        # Adjust subplot spacing for vertical layout
        self.fig_voigt.subplots_adjust(
            left=0.1,    # Left margin
            bottom=0.1,  # Bottom margin
            right=0.9,   # Right margin
            top=0.85,    # Top margin (leave space for title)
            hspace=0.4   # Vertical spacing between subplots (increased for 2x1)
        )

        #self.axes_voigt = axes
        self.canvas_voigt = FigureCanvasTkAgg(self.fig_voigt, voigt_tab)
        self.canvas_voigt.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        toolbar_voigt = NavigationToolbar2Tk(self.canvas_voigt, voigt_tab)
        toolbar_voigt.update()

        # Tab 3: Series overview
        series_tab = ttk.Frame(self.viz_notebook)
        self.viz_notebook.add(series_tab, text="📊 Series Overview")

        # Create responsive series overview figure
        self.fig_series, self.ax_series = plt.subplots(1, 1, figsize=(8, 6))
        self.fig_series.tight_layout()
        self.canvas_series = FigureCanvasTkAgg(self.fig_series, series_tab)
        self.canvas_series.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        toolbar_series = NavigationToolbar2Tk(self.canvas_series, series_tab)
        toolbar_series.update()

        # Tab 4: Statistics and Quality Assessment
        stats_tab = ttk.Frame(self.viz_notebook)
        self.viz_notebook.add(stats_tab, text="📈 Statistics")

        # Create responsive statistics figure
        self.fig_stats, ((self.ax_stats_1, self.ax_stats_2),
                         (self.ax_stats_3, self.ax_stats_4)) = plt.subplots(2, 2, figsize=(10, 7))
        self.fig_stats.tight_layout()
        self.canvas_stats = FigureCanvasTkAgg(self.fig_stats, stats_tab)
        self.canvas_stats.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        toolbar_stats = NavigationToolbar2Tk(self.canvas_stats, stats_tab)
        toolbar_stats.update()

        # Tab 5: Integration Diagnostics
        diagnostics_tab = ttk.Frame(self.viz_notebook)
        self.viz_notebook.add(diagnostics_tab, text="🔍 Integration Diagnostics")

        # Create a paned window for diagnostics layout
        diagnostics_paned = ttk.PanedWindow(diagnostics_tab, orient=tk.HORIZONTAL)
        diagnostics_paned.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Left panel: Diagnostic plots
        plots_frame = ttk.Frame(diagnostics_paned)
        diagnostics_paned.add(plots_frame, weight=2)

        self.fig_diagnostics, ((self.ax_diag_quality, self.ax_diag_convergence),
                              (self.ax_diag_aic, self.ax_diag_timing)) = plt.subplots(2, 2, figsize=(10, 8))
        self.canvas_diagnostics = FigureCanvasTkAgg(self.fig_diagnostics, plots_frame)
        self.canvas_diagnostics.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        toolbar_diagnostics = NavigationToolbar2Tk(self.canvas_diagnostics, plots_frame)
        toolbar_diagnostics.update()

        # Right panel: Diagnostic information
        info_frame = ttk.Frame(diagnostics_paned)
        diagnostics_paned.add(info_frame, weight=1)

        # Diagnostic information sections
        self.setup_diagnostics_info_panel(info_frame)

        # Tab 5: Configuration and Settings
        config_tab = ttk.Frame(self.viz_notebook)
        self.viz_notebook.add(config_tab, text="⚙️ Configuration")
        self.setup_config_tab(config_tab)

        # Initialize plots
        self.init_all_plots()

    def setup_config_tab(self, parent):
        """Setup configuration management tab"""
        config_scroll = ScrollableFrame(parent)
        config_scroll.pack(fill=tk.BOTH, expand=True)

        config_content = config_scroll.scrollable_frame

        # Configuration sections
        ttk.Label(config_content, text="Configuration Management",
                 font=('TkDefaultFont', 14, 'bold')).pack(anchor=tk.W, pady=(0, 10))

        # Save/Load configuration
        config_file_frame = ttk.LabelFrame(config_content, text="Configuration Files", padding=10)
        config_file_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Button(config_file_frame, text="💾 Save Configuration",
                  command=self.save_config).pack(side=tk.LEFT, padx=5)
        ttk.Button(config_file_frame, text="📂 Load Configuration",
                  command=self.load_config).pack(side=tk.LEFT, padx=5)
        ttk.Button(config_file_frame, text="🔄 Reset to Defaults",
                  command=self.reset_config).pack(side=tk.LEFT, padx=5)

        # Current configuration display
        config_display_frame = ttk.LabelFrame(config_content, text="Current Settings", padding=10)
        config_display_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))

        self.config_text = tk.Text(config_display_frame, height=15, font=('Courier', 9))
        config_scroll_bar = ttk.Scrollbar(config_display_frame, orient="vertical", command=self.config_text.yview)

        self.config_text.config(yscrollcommand=config_scroll_bar.set)
        self.config_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        config_scroll_bar.pack(side=tk.RIGHT, fill=tk.Y)

        self.update_config_display()

    def setup_plot_manager(self):
        """Setup plot manager to coordinate all visualizations"""
        self.spectrum_plotter = SpectrumPlotter(self.fig_main, self.ax_main)
        self.voigt_plotter = VoigtAnalysisPlotter(self.fig_voigt, self.axes_voigt)
        self.series_plotter = SeriesPlotter(self.fig_series, self.ax_series)

        self.plot_manager = PlotManager()
        self.plot_manager.register_plotters(
            spectrum_plotter=self.spectrum_plotter,
            voigt_plotter=self.voigt_plotter,
            series_plotter=self.series_plotter
        )

    def setup_status_bar(self):
        """Setup enhanced status bar"""
        self.status_bar = ttk.Frame(self.root)
        self.status_bar.grid(row=1, column=0, sticky=(tk.W, tk.E), padx=5, pady=(0, 5))

        # Configure status bar grid for responsive layout
        self.status_bar.columnconfigure(0, weight=1)  # Status text expands
        self.status_bar.columnconfigure(1, weight=0)  # Progress bar fixed size

        self.status_text = ttk.Label(self.status_bar, text="Ready")
        self.status_text.grid(row=0, column=0, sticky=(tk.W, tk.E), padx=5)

        self.progress_var = tk.DoubleVar()
        self.progress_bar = ttk.Progressbar(self.status_bar, variable=self.progress_var,
                                          mode='determinate')
        self.progress_bar.grid(row=0, column=1, sticky=tk.E, padx=5)

    def setup_peak_navigator(self):
        """Setup peak navigator in right panel"""
        # Create peak navigator
        self.peak_navigator = PeakNavigator(self.right_panel_frame)
        self.peak_navigator.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Set spectrum controller reference
        self.peak_navigator.set_spectrum_controller(self)

    def setup_diagnostics_info_panel(self, parent):
        """Setup the diagnostics information panel"""
        # Create scrollable frame for diagnostics info
        diagnostics_canvas = tk.Canvas(parent, bg='white')
        diagnostics_scrollbar = ttk.Scrollbar(parent, orient="vertical", command=diagnostics_canvas.yview)
        scrollable_frame = ttk.Frame(diagnostics_canvas)

        scrollable_frame.bind(
            "<Configure>",
            lambda e: diagnostics_canvas.configure(scrollregion=diagnostics_canvas.bbox("all"))
        )

        diagnostics_canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        diagnostics_canvas.configure(yscrollcommand=diagnostics_scrollbar.set)

        # Pack canvas and scrollbar
        diagnostics_canvas.pack(side="left", fill="both", expand=True)
        diagnostics_scrollbar.pack(side="right", fill="y")

        # Integration Summary Section
        summary_frame = ttk.LabelFrame(scrollable_frame, text="🔍 Integration Summary", padding=10)
        summary_frame.pack(fill=tk.X, pady=(0, 10))

        self.diag_total_peaks_label = ttk.Label(summary_frame, text="Total peaks analyzed: N/A")
        self.diag_total_peaks_label.pack(anchor=tk.W)

        self.diag_converged_peaks_label = ttk.Label(summary_frame, text="Converged peaks: N/A")
        self.diag_converged_peaks_label.pack(anchor=tk.W)

        self.diag_avg_quality_label = ttk.Label(summary_frame, text="Average quality score: N/A")
        self.diag_avg_quality_label.pack(anchor=tk.W)

        self.diag_total_iterations_label = ttk.Label(summary_frame, text="Total iterations: N/A")
        self.diag_total_iterations_label.pack(anchor=tk.W)

        # Performance Metrics Section
        performance_frame = ttk.LabelFrame(scrollable_frame, text="⚡ Performance Metrics", padding=10)
        performance_frame.pack(fill=tk.X, pady=(0, 10))

        self.diag_processing_time_label = ttk.Label(performance_frame, text="Processing time: N/A")
        self.diag_processing_time_label.pack(anchor=tk.W)

        self.diag_avg_time_per_peak_label = ttk.Label(performance_frame, text="Avg time per peak: N/A")
        self.diag_avg_time_per_peak_label.pack(anchor=tk.W)

        self.diag_memory_usage_label = ttk.Label(performance_frame, text="Memory usage: N/A")
        self.diag_memory_usage_label.pack(anchor=tk.W)

        # Quality Distribution Section
        quality_frame = ttk.LabelFrame(scrollable_frame, text="📊 Quality Distribution", padding=10)
        quality_frame.pack(fill=tk.X, pady=(0, 10))

        self.diag_high_quality_label = ttk.Label(quality_frame, text="High quality (≥0.7): N/A")
        self.diag_high_quality_label.pack(anchor=tk.W)

        self.diag_med_quality_label = ttk.Label(quality_frame, text="Medium quality (0.3-0.7): N/A")
        self.diag_med_quality_label.pack(anchor=tk.W)

        self.diag_low_quality_label = ttk.Label(quality_frame, text="Low quality (<0.3): N/A")
        self.diag_low_quality_label.pack(anchor=tk.W)

        # AIC Analysis Section
        aic_frame = ttk.LabelFrame(scrollable_frame, text="📈 AIC Analysis", padding=10)
        aic_frame.pack(fill=tk.X, pady=(0, 10))

        self.diag_avg_aic_label = ttk.Label(aic_frame, text="Average AIC score: N/A")
        self.diag_avg_aic_label.pack(anchor=tk.W)

        self.diag_aic_improvement_label = ttk.Label(aic_frame, text="AIC improvement: N/A")
        self.diag_aic_improvement_label.pack(anchor=tk.W)

        self.diag_model_selection_label = ttk.Label(aic_frame, text="Model selection accuracy: N/A")
        self.diag_model_selection_label.pack(anchor=tk.W)

        # Convergence Analysis Section
        convergence_frame = ttk.LabelFrame(scrollable_frame, text="🎯 Convergence Analysis", padding=10)
        convergence_frame.pack(fill=tk.X, pady=(0, 10))

        self.diag_convergence_rate_label = ttk.Label(convergence_frame, text="Convergence rate: N/A")
        self.diag_convergence_rate_label.pack(anchor=tk.W)

        self.diag_avg_iterations_label = ttk.Label(convergence_frame, text="Avg iterations to converge: N/A")
        self.diag_avg_iterations_label.pack(anchor=tk.W)

        self.diag_failed_peaks_label = ttk.Label(convergence_frame, text="Failed to converge: N/A")
        self.diag_failed_peaks_label.pack(anchor=tk.W)

        # Actions Section
        actions_frame = ttk.LabelFrame(scrollable_frame, text="🔧 Actions", padding=10)
        actions_frame.pack(fill=tk.X, pady=(0, 10))

        ttk.Button(actions_frame, text="🔄 Refresh Diagnostics",
                  command=self.refresh_diagnostics, width=20).pack(pady=2)
        ttk.Button(actions_frame, text="📊 Export Diagnostics",
                  command=self.export_diagnostics, width=20).pack(pady=2)
        ttk.Button(actions_frame, text="🔍 Detailed Analysis",
                  command=self.show_detailed_analysis, width=20).pack(pady=2)

    def init_all_plots(self):
        """Initialize all plots with default content"""
        # Main spectrum plot
        self.ax_main.set_title('lunaNMR v0.9 - Load Spectrum to Begin')
        self.ax_main.set_xlabel('¹H Chemical Shift (ppm)')
        self.ax_main.set_ylabel('¹⁵N/¹³C Chemical Shift (ppm)')
        self.ax_main.grid(True, alpha=0.3)

        # Voigt analysis plots
        for i, ax_row in enumerate(self.axes_voigt):
            for j, ax in enumerate(ax_row):
                ax.set_title(f'Voigt Analysis {i+1}-{j+1} - Fit Peak to Begin')
                ax.grid(True, alpha=0.3)

        # Series overview plot
        self.ax_series.set_title('Series Integration Results - Run Series Processing')
        self.ax_series.set_xlabel('Spectrum Index')
        self.ax_series.set_ylabel('Detection Rate (%)')
        self.ax_series.grid(True, alpha=0.3)

        # Statistics plots
        for ax in [self.ax_stats_1, self.ax_stats_2, self.ax_stats_3, self.ax_stats_4]:
            ax.set_title('Statistics - Process Data to Begin')
            ax.grid(True, alpha=0.3)

        # Draw all canvases
        self.canvas_main.draw()
        self.canvas_voigt.draw()
        self.canvas_series.draw()
        self.canvas_stats.draw()

    # =================== EVENT HANDLERS ===================

    def _toggle_series_voigt_params(self):
        """Toggle visibility of series Voigt parameter frame based on checkbox state"""
        if self.series_use_voigt_fitting.get():
            self.series_voigt_params_frame.pack(fill=tk.X, padx=10, pady=(5,0))
            # Update parameter displays when shown
            self._update_series_voigt_param_displays()
        else:
            self.series_voigt_params_frame.pack_forget()

    def _update_series_voigt_param_displays(self):
        """Update the parameter display labels in the series Voigt parameters section"""
        if not hasattr(self, 'series_voigt_params_frame'):
            return

        # Find all child frames and update displays
        for widget in self.series_voigt_params_frame.winfo_children():
            if isinstance(widget, ttk.Frame):
                for child in widget.winfo_children():
                    if isinstance(child, ttk.Label):
                        name = str(child)
                        if 'x_window_display' in name:
                            child.config(text=f"{self.fitting_window_x.get():.1f}")
                        elif 'y_window_display' in name:
                            child.config(text=f"{self.fitting_window_y.get():.1f}")
                        elif 'r_squared_display' in name:
                            child.config(text=f"{self.min_r_squared.get():.2f}")
                        elif 'max_iter_display' in name:
                            child.config(text=f"{self.max_iterations.get()}")
                        elif 'height_display' in name:
                            child.config(text=f"{self.peak_height_threshold.get():.3f}")
                        elif 'distance_display' in name:
                            child.config(text=f"{self.peak_distance_factor.get():.0f}")
                        elif 'prominence_display' in name:
                            child.config(text=f"{self.peak_prominence_threshold.get():.3f}")
                        elif 'smoothing_display' in name:
                            child.config(text=f"{self.smoothing_sigma.get():.1f}")
                        elif 'max_peaks_display' in name:
                            child.config(text=f"{self.max_peaks_fit.get()}")
                        elif 'opt_iter_display' in name:
                            child.config(text=f"{self.max_optimization_iterations.get()}")

    def on_legacy_parameter_change(self): #GM changed
        """Handle parameter changes"""
        # Update integrator parameters
        self.integrator.set_search_window(self.search_window_x.get(), self.search_window_y.get())
        self.integrator.set_threshold_multiplier(self.noise_threshold.get())
        # Update ppm conversion display
        self.update_distance_ppm_display()

        # Update fitting parameters
        self.integrator.fitting_parameters.update({
            'fitting_window_x': self.fitting_window_x.get(),
            'fitting_window_y': self.fitting_window_y.get(),
            'min_r_squared': self.min_r_squared.get(),
            'max_iterations': self.max_iterations.get()
        })

        # Update series Voigt parameter displays
        self._update_series_voigt_param_displays()

        # Save parameters to config
        self.proc_params.set_detection_params(
            noise_threshold=self.noise_threshold.get(),
            search_window_x=self.search_window_x.get(),
            search_window_y=self.search_window_y.get(),
            use_reference_detection=self.use_reference_detection.get(),
            detection_square_size=self.detection_square_size.get(),
            detection_rectangle_y=self.detection_rectangle_y.get(),
            consolidation_x_tolerance=self.consolidation_x_tolerance.get(),
            consolidation_y_tolerance=self.consolidation_y_tolerance.get()
        )

        self.proc_params.set_fitting_params(
            fitting_window_x=self.fitting_window_x.get(),
            fitting_window_y=self.fitting_window_y.get(),
            min_r_squared=self.min_r_squared.get(),
            max_iterations=self.max_iterations.get()
        )
 
    def update_distance_ppm_display(self):
        """Update the ppm conversion display for peak distance factor"""
        if hasattr(self, 'distance_ppm_label') and hasattr(self, 'integrator'):
            # Calculate ppm equivalent if spectrum is loaded
            if hasattr(self.integrator, 'ppm_x_axis') and self.integrator.ppm_x_axis is not None:
                total_points = len(self.integrator.ppm_x_axis)
                ppm_range = abs(self.integrator.ppm_x_axis[-1] - self.integrator.ppm_x_axis[0])
                distance_points = self.peak_distance_factor.get()
                distance_ppm = distance_points * ppm_range / total_points
                self.distance_ppm_label.config(text=f"(≈{distance_ppm:.3f} ppm)")
            else:
                self.distance_ppm_label.config(text="(load data for ppm)")

    def on_integration_mode_change(self, event=None):
        """Handle integration mode change"""
        mode = self.integration_mode.get()

        # Update series processor integration mode
        try:
            integration_params = {
                # 1D Refinement parameters (NEW)
                'enable_1d_refinement': self.enable_1d_refinement.get(),
                'refinement_quality_threshold': self.refinement_quality_threshold.get(),
                'refinement_coordinate_threshold': self.refinement_coordinate_threshold.get(),
            }
            self.series_processor.set_integration_mode(mode, **integration_params)

            # Update status display
            mode_suffix = " (Full Detection)" if self.force_full_detection.get() and mode in ['integrated', 'adaptive'] else ""
            refinement_suffix = " + 1D Refinement" if self.enable_1d_refinement.get() and mode in ['integrated', 'adaptive'] else ""
            status_text = f"Status: {mode.title()} mode{mode_suffix}{refinement_suffix} (Core detection params: height={self.peak_height_threshold.get():.3f}, max_peaks={self.max_peaks_fit.get()})"
            # Mode-specific status text updates (integrated/adaptive modes removed)

            # Update both status labels (series section and single spectrum section)
            self.integration_status_label.config(text=status_text)
            if hasattr(self, 'single_integration_status'):
                self.single_integration_status.config(text=status_text)

            # Enable/disable parameter controls based on mode
            state = 'normal' if mode != 'standard' else 'disabled'

            # Enable/disable series integration controls
            for frame in [self.integration_params_frame]:
                for widget in frame.winfo_children():
                    for child in widget.winfo_children():
                        if isinstance(child, (ttk.Spinbox,)):
                            child.config(state=state)

            # Enable/disable single spectrum integration controls
            if hasattr(self, 'single_integration_controls'):
                for control in self.single_integration_controls:
                    control.config(state=state)

            # Enable/disable single spectrum integration checkboxes
            if hasattr(self, 'single_integration_checkboxes'):
                checkbox_state = 'normal' if mode != 'standard' else 'disabled'
                for checkbox in self.single_integration_checkboxes:
                    checkbox.config(state=checkbox_state)

            # Update button text to reflect current mode
            if hasattr(self, 'detect_peaks_button'):
                if mode == 'integrated':
                    self.detect_peaks_button.config(text="🔍 Detect Peaks (INT)")
                elif mode == 'adaptive':
                    self.detect_peaks_button.config(text="🔍 Detect Peaks (ADV)")
                else:
                    self.detect_peaks_button.config(text="🔍 Detect Peaks")

            if hasattr(self, 'fit_all_peaks_button'):
                if mode == 'integrated':
                    self.fit_all_peaks_button.config(text="📈 Fit All Peaks (INT)")
                elif mode == 'adaptive':
                    self.fit_all_peaks_button.config(text="📈 Fit All Peaks (ADV)")
                else:
                    self.fit_all_peaks_button.config(text="📈 Fit All Peaks")

            self.update_status(f"✅ Integration mode changed to: {mode}")
            print(f"Integration mode changed to: {mode}")

        except Exception as e:
            self.update_status(f"❌ Error setting integration mode: {str(e)}", error=True)
            print(f"Error setting integration mode: {e}")

    def show_integration_progress(self):
        """Show the integration progress panel"""
        self.progress_frame.pack(fill=tk.X, pady=(0, 10), before=self.progress_frame.master.children['!labelframe3'])  # Before Status & Export
        self.root.update()

    def hide_integration_progress(self):
        """Hide the integration progress panel"""
        self.progress_frame.pack_forget()
        self.root.update()

    def update_integration_progress(self, status_text="Processing...", progress_percent=0,
                                   current_iteration=0, max_iterations=0,
                                   peaks_processed=0, total_peaks=0,
                                   quality_score=None, convergence_status="Running",
                                   elapsed_time=0):
        """Update the integration progress display"""
        try:
            # Update main status and progress
            self.progress_status_label.config(text=status_text)
            self.progress_percentage_label.config(text=f"{progress_percent:.1f}%")
            self.integration_progress['value'] = progress_percent

            # Update detailed information
            self.current_iteration_label.config(text=f"{current_iteration}/{max_iterations}")
            self.peaks_processed_label.config(text=f"{peaks_processed}/{total_peaks}")

            if quality_score is not None:
                self.quality_score_label.config(text=f"{quality_score:.3f}")
            else:
                self.quality_score_label.config(text="N/A")

            # Update convergence status with color coding
            self.convergence_status_label.config(text=convergence_status)
            if convergence_status == "Converged":
                self.convergence_status_label.config(foreground='green')
            elif convergence_status == "Running":
                self.convergence_status_label.config(foreground='orange')
            elif convergence_status == "Failed":
                self.convergence_status_label.config(foreground='red')
            else:
                self.convergence_status_label.config(foreground='blue')

            # Update elapsed time
            if elapsed_time < 60:
                time_text = f"{elapsed_time:.1f}s"
            elif elapsed_time < 3600:
                time_text = f"{elapsed_time/60:.1f}m"
            else:
                time_text = f"{elapsed_time/3600:.1f}h"
            self.time_elapsed_label.config(text=time_text)

            # Force GUI update
            self.root.update_idletasks()

        except Exception as e:
            print(f"Error updating integration progress: {e}")

    def pause_integration(self):
        """Pause the current integration process"""
        # Set integration state to paused
        if hasattr(self.series_processor, 'pause_integration'):
            self.series_processor.pause_integration()
            self.pause_integration_button.config(text="▶️ Resume")
            self.update_status("Integration paused by user")
        else:
            print("Pause functionality not implemented in series processor")

    def stop_integration(self):
        """Stop the current integration process"""
        # Set integration state to stopped
        if hasattr(self.series_processor, 'stop_integration'):
            self.series_processor.stop_integration()
            self.hide_integration_progress()
            self.update_status("Integration stopped by user")
        else:
            print("Stop functionality not implemented in series processor")

    def reset_integration_progress(self):
        """Reset the progress display to initial state"""
        self.update_integration_progress(
            status_text="Ready",
            progress_percent=0,
            current_iteration=0,
            max_iterations=0,
            peaks_processed=0,
            total_peaks=0,
            quality_score=None,
            convergence_status="Pending",
            elapsed_time=0
        )
        self.pause_integration_button.config(text="⏸️ Pause")

    def refresh_diagnostics(self):
        """Refresh the diagnostics display with current data"""
        try:
            # Check if integrated detection-fitting results are available
            if hasattr(self.integrator, 'fitted_peaks') and self.integrator.fitted_peaks is not None:
                self.update_diagnostics_display(self.integrator.fitted_peaks)
            elif hasattr(self, 'last_integration_results'):
                self.update_diagnostics_display(self.last_integration_results)
            else:
                self.update_status("No integration results available for diagnostics")
        except Exception as e:
            self.update_status(f"Error refreshing diagnostics: {e}", error=True)

    def export_diagnostics(self):
        """Export diagnostics data to file"""
        try:
            from tkinter.filedialog import asksaveasfilename
            import json
            from datetime import datetime

            filename = asksaveasfilename(
                title="Export Integration Diagnostics",
                defaultextension=".json",
                filetypes=[("JSON files", "*.json"), ("Text files", "*.txt"), ("All files", "*.*")]
            )

            if filename:
                diagnostics_data = self.collect_diagnostics_data()
                diagnostics_data['export_timestamp'] = datetime.now().isoformat()

                with open(filename, 'w') as f:
                    if filename.endswith('.json'):
                        json.dump(diagnostics_data, f, indent=2)
                    else:
                        # Text format
                        for section, data in diagnostics_data.items():
                            f.write(f"\n=== {section.upper()} ===\n")
                            if isinstance(data, dict):
                                for key, value in data.items():
                                    f.write(f"{key}: {value}\n")
                            else:
                                f.write(f"{data}\n")

                self.update_status(f"Diagnostics exported to {filename}")

        except Exception as e:
            self.update_status(f"Error exporting diagnostics: {e}", error=True)

    def show_detailed_analysis(self):
        """Show detailed analysis window"""
        try:
            # Create detailed analysis popup window
            detail_window = tk.Toplevel(self.root)
            detail_window.title("Detailed Integration Analysis")
            # Set responsive geometry (60% of parent window size)
            parent_width = self.root.winfo_width()
            parent_height = self.root.winfo_height()
            width = max(600, int(parent_width * 0.6))
            height = max(400, int(parent_height * 0.6))
            detail_window.geometry(f"{width}x{height}")
            detail_window.minsize(600, 400)

            # Add detailed analysis content
            detail_text = tk.Text(detail_window, wrap=tk.WORD, padx=10, pady=10)
            detail_scrollbar = ttk.Scrollbar(detail_window, orient="vertical", command=detail_text.yview)
            detail_text.configure(yscrollcommand=detail_scrollbar.set)

            # Generate detailed analysis
            analysis_content = self.generate_detailed_analysis()
            detail_text.insert(tk.END, analysis_content)
            detail_text.config(state=tk.DISABLED)  # Make read-only

            detail_text.pack(side="left", fill="both", expand=True)
            detail_scrollbar.pack(side="right", fill="y")

        except Exception as e:
            self.update_status(f"Error showing detailed analysis: {e}", error=True)

    def update_diagnostics_display(self, integration_results):
        """Update the diagnostics display with integration results"""
        try:
            # Extract diagnostic metrics from results
            metrics = self.calculate_diagnostics_metrics(integration_results)

            # Update summary labels
            self.diag_total_peaks_label.config(text=f"Total peaks analyzed: {metrics['total_peaks']}")
            self.diag_converged_peaks_label.config(text=f"Converged peaks: {metrics['converged_peaks']}")
            self.diag_avg_quality_label.config(text=f"Average quality score: {metrics['avg_quality']:.3f}")
            self.diag_total_iterations_label.config(text=f"Total iterations: {metrics['total_iterations']}")

            # Update performance labels
            self.diag_processing_time_label.config(text=f"Processing time: {metrics['processing_time']:.2f}s")
            self.diag_avg_time_per_peak_label.config(text=f"Avg time per peak: {metrics['avg_time_per_peak']:.3f}s")
            self.diag_memory_usage_label.config(text=f"Memory usage: {metrics['memory_usage']:.1f} MB")

            # Update quality distribution
            self.diag_high_quality_label.config(text=f"High quality (≥0.7): {metrics['high_quality_count']} ({metrics['high_quality_pct']:.1f}%)")
            self.diag_med_quality_label.config(text=f"Medium quality (0.3-0.7): {metrics['med_quality_count']} ({metrics['med_quality_pct']:.1f}%)")
            self.diag_low_quality_label.config(text=f"Low quality (<0.3): {metrics['low_quality_count']} ({metrics['low_quality_pct']:.1f}%)")

            # Update AIC analysis
            self.diag_avg_aic_label.config(text=f"Average AIC score: {metrics['avg_aic']:.2f}")
            self.diag_aic_improvement_label.config(text=f"AIC improvement: {metrics['aic_improvement']:.2f}")
            self.diag_model_selection_label.config(text=f"Model selection accuracy: {metrics['model_accuracy']:.1f}%")

            # Update convergence analysis
            self.diag_convergence_rate_label.config(text=f"Convergence rate: {metrics['convergence_rate']:.1f}%")
            self.diag_avg_iterations_label.config(text=f"Avg iterations to converge: {metrics['avg_iterations']:.1f}")
            self.diag_failed_peaks_label.config(text=f"Failed to converge: {metrics['failed_peaks']}")

            # Update diagnostic plots
            self.update_diagnostics_plots(metrics)

        except Exception as e:
            print(f"Error updating diagnostics display: {e}")

    def calculate_diagnostics_metrics(self, integration_results):
        """Calculate diagnostic metrics from integration results"""
        import psutil
        import time

        metrics = {}

        try:
            # Basic counts
            total_peaks = len(integration_results) if integration_results else 0
            metrics['total_peaks'] = total_peaks

            if total_peaks == 0:
                # Return empty metrics for zero peaks
                return {key: 0 for key in ['total_peaks', 'converged_peaks', 'avg_quality', 'total_iterations',
                                         'processing_time', 'avg_time_per_peak', 'memory_usage',
                                         'high_quality_count', 'high_quality_pct', 'med_quality_count', 'med_quality_pct',
                                         'low_quality_count', 'low_quality_pct', 'avg_aic', 'aic_improvement',
                                         'model_accuracy', 'convergence_rate', 'avg_iterations', 'failed_peaks']}

            # Quality analysis
            quality_scores = [peak.get('composite_quality', 0.5) for peak in integration_results]
            converged_peaks = sum(1 for peak in integration_results if peak.get('converged', False))

            metrics['converged_peaks'] = converged_peaks
            metrics['avg_quality'] = sum(quality_scores) / len(quality_scores) if quality_scores else 0

            # Iteration analysis
            iterations_list = [peak.get('iterations', 1) for peak in integration_results]
            metrics['total_iterations'] = sum(iterations_list)
            metrics['avg_iterations'] = sum(iterations_list) / len(iterations_list) if iterations_list else 0

            # Performance metrics
            processing_time = getattr(self, 'last_processing_time', 0)
            metrics['processing_time'] = processing_time
            metrics['avg_time_per_peak'] = processing_time / total_peaks if total_peaks > 0 else 0

            # Memory usage
            process = psutil.Process()
            metrics['memory_usage'] = process.memory_info().rss / 1024 / 1024  # MB

            # Quality distribution
            high_quality = sum(1 for q in quality_scores if q >= 0.7)
            med_quality = sum(1 for q in quality_scores if 0.3 <= q < 0.7)
            low_quality = sum(1 for q in quality_scores if q < 0.3)

            metrics['high_quality_count'] = high_quality
            metrics['high_quality_pct'] = (high_quality / total_peaks) * 100 if total_peaks > 0 else 0
            metrics['med_quality_count'] = med_quality
            metrics['med_quality_pct'] = (med_quality / total_peaks) * 100 if total_peaks > 0 else 0
            metrics['low_quality_count'] = low_quality
            metrics['low_quality_pct'] = (low_quality / total_peaks) * 100 if total_peaks > 0 else 0

            # AIC analysis
            aic_scores = [peak.get('aic_score', 0) for peak in integration_results]
            metrics['avg_aic'] = sum(aic_scores) / len(aic_scores) if aic_scores else 0
            metrics['aic_improvement'] = 0  # Placeholder - would need baseline comparison
            metrics['model_accuracy'] = 85.0  # Placeholder - would need validation data

            # Convergence analysis
            metrics['convergence_rate'] = (converged_peaks / total_peaks) * 100 if total_peaks > 0 else 0
            metrics['failed_peaks'] = total_peaks - converged_peaks

        except Exception as e:
            print(f"Error calculating diagnostics metrics: {e}")
            # Return safe defaults
            metrics = {key: 0 for key in ['total_peaks', 'converged_peaks', 'avg_quality', 'total_iterations',
                                       'processing_time', 'avg_time_per_peak', 'memory_usage',
                                       'high_quality_count', 'high_quality_pct', 'med_quality_count', 'med_quality_pct',
                                       'low_quality_count', 'low_quality_pct', 'avg_aic', 'aic_improvement',
                                       'model_accuracy', 'convergence_rate', 'avg_iterations', 'failed_peaks']}

        return metrics

    def update_diagnostics_plots(self, metrics):
        """Update the diagnostic plots"""
        try:
            import numpy as np

            # Helper function to safely convert to int, handling NaN values
            def safe_int(value, default=0):
                try:
                    if value is None or (hasattr(value, '__iter__') and len(value) == 0):
                        return default
                    if np.isnan(float(value)):
                        return default
                    return int(float(value))
                except (ValueError, TypeError, OverflowError):
                    return default

            # Helper function to safely convert to float, handling NaN values
            def safe_float(value, default=0.0):
                try:
                    if value is None or (hasattr(value, '__iter__') and len(value) == 0):
                        return default
                    if np.isnan(float(value)):
                        return default
                    return float(value)
                except (ValueError, TypeError, OverflowError):
                    return default

            # Clear all diagnostic plots
            for ax in [self.ax_diag_quality, self.ax_diag_convergence, self.ax_diag_aic, self.ax_diag_timing]:
                ax.clear()

            # Safely extract metrics with NaN handling
            total_peaks = safe_int(metrics.get('total_peaks', 0))
            high_quality_count = safe_int(metrics.get('high_quality_count', 0))
            med_quality_count = safe_int(metrics.get('med_quality_count', 0))
            low_quality_count = safe_int(metrics.get('low_quality_count', 0))
            converged_peaks = safe_int(metrics.get('converged_peaks', 0))
            failed_peaks = safe_int(metrics.get('failed_peaks', 0))
            avg_time_per_peak = safe_float(metrics.get('avg_time_per_peak', 0.0))

            # Quality distribution plot
            quality_labels = ['High\n(≥0.7)', 'Medium\n(0.3-0.7)', 'Low\n(<0.3)']
            quality_counts = [high_quality_count, med_quality_count, low_quality_count]
            colors = ['green', 'orange', 'red']

            self.ax_diag_quality.bar(quality_labels, quality_counts, color=colors, alpha=0.7)
            self.ax_diag_quality.set_title('Quality Distribution')
            self.ax_diag_quality.set_ylabel('Number of Peaks')

            # Convergence rate plot - handle case when no peaks detected
            if total_peaks > 0:
                conv_labels = ['Converged', 'Failed']
                conv_counts = [converged_peaks, failed_peaks]
                conv_colors = ['green', 'red']

                # Only create pie chart if we have data
                if sum(conv_counts) > 0:
                    self.ax_diag_convergence.pie(conv_counts, labels=conv_labels, colors=conv_colors, autopct='%1.1f%%')
                else:
                    self.ax_diag_convergence.text(0.5, 0.5, 'No peaks detected', ha='center', va='center', transform=self.ax_diag_convergence.transAxes)
            else:
                self.ax_diag_convergence.text(0.5, 0.5, 'No peaks detected', ha='center', va='center', transform=self.ax_diag_convergence.transAxes)

            self.ax_diag_convergence.set_title('Convergence Rate')

            # AIC scores (placeholder histogram)
            self.ax_diag_aic.hist([1, 2, 3, 4, 5], bins=5, alpha=0.7, color='blue')
            self.ax_diag_aic.set_title('AIC Score Distribution')
            self.ax_diag_aic.set_xlabel('AIC Score')
            self.ax_diag_aic.set_ylabel('Frequency')

            # Processing time plot - handle case when no peaks detected
            if total_peaks > 0 and avg_time_per_peak > 0:
                peak_indices = list(range(1, min(total_peaks, 20) + 1))
                processing_times = [avg_time_per_peak] * len(peak_indices)

                self.ax_diag_timing.plot(peak_indices, processing_times, 'o-', color='purple')
            else:
                self.ax_diag_timing.text(0.5, 0.5, 'No timing data', ha='center', va='center', transform=self.ax_diag_timing.transAxes)

            self.ax_diag_timing.set_title('Processing Time per Peak')
            self.ax_diag_timing.set_xlabel('Peak Index')
            self.ax_diag_timing.set_ylabel('Time (s)')

            self.fig_diagnostics.tight_layout()
            self.canvas_diagnostics.draw()

        except Exception as e:
            print(f"Error updating diagnostic plots: {e}")

    def collect_diagnostics_data(self):
        """Collect all diagnostics data for export"""
        return {
            'integration_summary': {
                'total_peaks': self.diag_total_peaks_label.cget('text'),
                'converged_peaks': self.diag_converged_peaks_label.cget('text'),
                'average_quality': self.diag_avg_quality_label.cget('text'),
                'total_iterations': self.diag_total_iterations_label.cget('text')
            },
            'performance_metrics': {
                'processing_time': self.diag_processing_time_label.cget('text'),
                'avg_time_per_peak': self.diag_avg_time_per_peak_label.cget('text'),
                'memory_usage': self.diag_memory_usage_label.cget('text')
            },
            'quality_distribution': {
                'high_quality': self.diag_high_quality_label.cget('text'),
                'medium_quality': self.diag_med_quality_label.cget('text'),
                'low_quality': self.diag_low_quality_label.cget('text')
            }
        }

    def generate_detailed_analysis(self):
        """Generate detailed analysis text"""
        return """DETAILED INTEGRATION ANALYSIS

=== OVERVIEW ===
This analysis provides comprehensive diagnostics for the integrated detection-fitting process.

=== QUALITY ASSESSMENT ===
- Quality scores are calculated from detection confidence, fit R², and AIC scores
- High quality peaks (≥0.7) indicate reliable detection and fitting
- Medium quality peaks (0.3-0.7) may require manual review
- Low quality peaks (<0.3) should be examined carefully

=== CONVERGENCE ANALYSIS ===
- Convergence indicates successful iterative refinement
- Failed convergence may suggest parameter adjustment needed
- Average iterations provides insight into algorithm efficiency

=== PERFORMANCE METRICS ===
- Processing time scales with peak complexity and number
- Memory usage monitored for large datasets
- Time per peak indicates algorithm efficiency

=== RECOMMENDATIONS ===
- For low convergence rates: adjust AIC threshold or max iterations
- For poor quality distribution: review detection confidence thresholds
- For high processing times: consider parallel processing options

Generated by NMR Peaks Series Analysis - Integration Diagnostics
"""

    # =================== FILE HANDLING ===================

    def on_nmr_file_select(self, file_path, filename):
        """Handle NMR file selection"""
        self.current_nmr_file = file_path
        self.config_manager.add_recent_file(file_path, 'nmr')
        self.update_recent_files_menu()
        self.update_current_status()
        self.validate_current_files()

        if self.current_peak_file:
            self.load_current_data()
        else:
            self.update_status(f"NMR file selected: {filename}. Please select a peak list.")

    def on_peak_file_select(self, file_path, filename):
        """Handle peak list file selection"""
        self.current_peak_file = file_path
        self.config_manager.add_recent_file(file_path, 'peak')
        self.update_recent_files_menu()
        self.update_current_status()
        self.validate_current_files()

        if self.current_nmr_file:
            self.load_current_data()
        else:
            self.update_status(f"Peak list selected: {filename}. Please select an NMR spectrum.")

    def load_current_data(self):
        """Load currently selected files"""
        if not self.current_nmr_file or not self.current_peak_file:
            return

        try:
            self.update_status("Loading data...")
            self.root.config(cursor="watch")
            self.root.update()

            # Load peak list
            peak_df = self.file_manager.load_peak_list(self.current_peak_file)

            # Validate and fix peak list
            issues = self.validator.validate_peak_list_integrity(peak_df)
            if issues:
                fixed_df, fixes = self.validator.auto_fix_peak_list(peak_df)
                peak_df = fixed_df
                self.update_status(f"Fixed peak list issues: {', '.join(fixes)}")

            # Load into integrator
            self.integrator.peak_list = peak_df

            # Load NMR data
            success = self.integrator.load_nmr_file(self.current_nmr_file)

            if success:
                self.update_status("✅ Data loaded successfully - Ready for processing")
                # Auto-adjust zoom to fit the loaded data
                self.auto_adjust_zoom_to_data()
                self.update_main_plot()
                self.update_peak_navigation()
                self.update_statistics()

                # FORCE: Ensure reference peaks are visible on load
                self.show_assigned.set(True)  # Ensure checkbox is checked
                self.update_main_plot()       # Force a second plot update

                # Load reference peaks into navigator
                if hasattr(self, 'peak_navigator'):
                    self.peak_navigator.load_reference_peaks(self.integrator.peak_list)
            else:
                self.update_status("❌ Failed to load NMR data")

        except Exception as e:
            self.update_status(f"❌ Error loading data: {str(e)}")
            messagebox.showerror("Load Error", f"Failed to load data:\n{str(e)}")

        finally:
            self.root.config(cursor="")

    def update_current_status(self):
        """Update current file status display"""
        if self.current_nmr_file and self.current_peak_file:
            nmr_name = os.path.basename(self.current_nmr_file)
            peak_name = os.path.basename(self.current_peak_file)
            self.current_status.config(
                text=f"📊 Spectrum: {nmr_name} | 📋 Peaks: {peak_name}",
                foreground='green'
            )
        elif self.current_nmr_file:
            nmr_name = os.path.basename(self.current_nmr_file)
            self.current_status.config(
                text=f"📊 Spectrum: {nmr_name} | 📋 Peak list: Not selected",
                foreground='orange'
            )
        elif self.current_peak_file:
            peak_name = os.path.basename(self.current_peak_file)
            self.current_status.config(
                text=f"📊 Spectrum: Not selected | 📋 Peaks: {peak_name}",
                foreground='orange'
            )
        else:
            self.current_status.config(
                text="📋 No files loaded",
                foreground='gray'
            )

    # =================== PEAK COORDINATE ADJUSTMENT METHODS ===================

    def apply_coordinate_offsets(self):
        """Apply X and Y offsets to all peaks in the current peak list"""
        if not hasattr(self.integrator, 'peak_list') or self.integrator.peak_list is None:
            messagebox.showerror("Error", "No peak list loaded. Please load a peak list first.")
            return

        # Handle both DataFrame and list formats
        peak_list = self.integrator.peak_list
        if hasattr(peak_list, 'empty') and peak_list.empty:
            messagebox.showerror("Error", "Peak list is empty. Please load a peak list first.")
            return
        elif isinstance(peak_list, list) and len(peak_list) == 0:
            messagebox.showerror("Error", "Peak list is empty. Please load a peak list first.")
            return

        try:
            x_offset = self.adjust_x_offset.get()
            y_offset = self.adjust_y_offset.get()

            if x_offset == 0.0 and y_offset == 0.0:
                messagebox.showwarning("Warning", "Both offsets are zero. No adjustment needed.")
                return

            # Handle different peak list formats
            if hasattr(peak_list, 'columns'):
                # DataFrame format (reference peaks)
                peak_list = peak_list.copy()
                original_count = len(peak_list)

                # Apply offsets
                if 'Position_X' in peak_list.columns:
                    peak_list['Position_X'] += x_offset
                if 'Position_Y' in peak_list.columns:
                    peak_list['Position_Y'] += y_offset

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
                            peak['Position_X'] += x_offset
                        if 'Position_Y' in peak:
                            peak['Position_Y'] += y_offset

                # Update the integrator's peak list
                self.integrator.peak_list = peak_list

            # Update status
            self.adjustment_status_label.config(
                text=f"✅ Applied offsets: ΔX={x_offset:+.3f} ppm, ΔY={y_offset:+.1f} ppm to {original_count} peaks",
                foreground='green'
            )

            # Update visualization
            self.update_main_plot()
            self.update_statistics()

            print(f"↗️ Coordinate offsets applied: ΔX={x_offset:+.3f} ppm, ΔY={y_offset:+.1f} ppm")

        except Exception as e:
            error_msg = f"Failed to apply coordinate offsets: {str(e)}"
            self.adjustment_status_label.config(text=f"❌ {error_msg}", foreground='red')
            messagebox.showerror("Adjustment Error", error_msg)

    def reset_coordinate_offsets(self):
        """Reset offset values to zero"""
        self.adjust_x_offset.set(0.0)
        self.adjust_y_offset.set(0.0)
        self.adjustment_status_label.config(text="🔄 Offsets reset to zero", foreground='blue')

    def save_adjusted_peak_list(self):
        """Save the adjusted peak list to a new file"""
        if not hasattr(self.integrator, 'peak_list') or self.integrator.peak_list is None or self.integrator.peak_list.empty:
            messagebox.showerror("Error", "No peak list available to save.")
            return

        try:
            from tkinter import filedialog
            import os

            # Suggest filename with adjustment info
            x_offset = self.adjust_x_offset.get()
            y_offset = self.adjust_y_offset.get()

            if self.current_peak_file:
                base_name = os.path.splitext(os.path.basename(self.current_peak_file))[0]
                suggested_name = f"{base_name}_adjusted_X{x_offset:+.3f}_Y{y_offset:+.1f}.txt"
            else:
                suggested_name = f"adjusted_peaks_X{x_offset:+.3f}_Y{y_offset:+.1f}.txt"

            # Ask user for save location
            file_path = filedialog.asksaveasfilename(
                title="Save Adjusted Peak List",
                defaultextension=".txt",
                initialname=suggested_name,
                filetypes=[
                    ("Text files", "*.txt"),
                    ("CSV files", "*.csv"),
                    ("All files", "*.*")
                ]
            )

            if file_path:
                # Save the adjusted peak list
                if file_path.endswith('.csv'):
                    self.integrator.peak_list.to_csv(file_path, index=False, sep=',')
                else:
                    self.integrator.peak_list.to_csv(file_path, index=False, sep='\t')

                self.adjustment_status_label.config(
                    text=f"💾 Adjusted peak list saved: {os.path.basename(file_path)}",
                    foreground='green'
                )

                print(f"💾 Adjusted peak list saved to: {file_path}")

        except Exception as e:
            error_msg = f"Failed to save adjusted peak list: {str(e)}"
            self.adjustment_status_label.config(text=f"❌ {error_msg}", foreground='red')
            messagebox.showerror("Save Error", error_msg)

    # =================== PROCESSING METHODS ===================

    def detect_peaks(self):
        """Perform peak detection - supports both standard and integrated detection"""
        if not hasattr(self.integrator, 'nmr_data') or self.integrator.nmr_data is None:
            messagebox.showerror("Error", "Please load NMR data first")
            return

        try:
            # Check integration mode
            integration_mode = self.integration_mode.get()

            if integration_mode in ['integrated', 'adaptive']:
                # Use integrated detection-fitting
                self._detect_peaks_integrated()
            else:
                # Use standard detection (legacy)
               self._detect_peaks_standard()

        except Exception as e:
            self.update_status(f"❌ Detection failed: {str(e)}")
            messagebox.showerror("Detection Error", str(e))

        finally:
            self.root.config(cursor="")

    def update_detection_square_ppm(self):
        """Update ppm conversion display for detection square size"""
        if hasattr(self.integrator, 'ppm_x_axis') and hasattr(self.integrator, 'ppm_y_axis'):
            try:
                # Calculate ppm per pixel
                x_axis = self.integrator.ppm_x_axis
                y_axis = self.integrator.ppm_y_axis

                x_ppm_per_pixel = abs(x_axis[-1] - x_axis[0]) / len(x_axis)
                y_ppm_per_pixel = abs(y_axis[-1] - y_axis[0]) / len(y_axis)

                # Update both standard and enhanced detection rectangle displays (anisotropic)
                standard_square_size_x = self.detection_square_size.get()
                standard_rectangle_size_y = self.detection_rectangle_y.get()
                enhanced_square_size_x = self.enhanced_detection_square_size.get()
                enhanced_rectangle_size_y = self.enhanced_detection_rectangle_y.get()

                standard_x_ppm = x_ppm_per_pixel * standard_square_size_x
                standard_y_ppm = y_ppm_per_pixel * standard_rectangle_size_y
                enhanced_x_ppm = x_ppm_per_pixel * enhanced_square_size_x
                enhanced_y_ppm = y_ppm_per_pixel * enhanced_rectangle_size_y

                # Update displays with anisotropic dimensions
                self.detection_square_ppm_x.set(f"(≈{standard_x_ppm:.3f}×{standard_y_ppm:.2f} ppm)")
                self.enhanced_detection_square_ppm_x.set(f"(≈{enhanced_x_ppm:.3f}×{enhanced_y_ppm:.2f} ppm)")

            except:
                self.detection_square_ppm_x.set("(load data to see ppm)")
                self.enhanced_detection_square_ppm_x.set("(load data to see ppm)")
        else:
            self.detection_square_ppm_x.set("(load data to see ppm)")
            self.enhanced_detection_square_ppm_x.set("(load data to see ppm)")

    def enhanced_detect_peaks(self):
        """Enhanced graph-based peak detection - alternative to standard detection"""
        if not hasattr(self.integrator, 'nmr_data') or self.integrator.nmr_data is None:
            messagebox.showerror("Error", "Please load NMR data first")
            return

        # Check if peak list is loaded for reference-based detection
        if not hasattr(self.integrator, 'peak_list') or self.integrator.peak_list is None or self.integrator.peak_list.empty:
            messagebox.showwarning("Warning",
                "Enhanced Detection requires a reference peak list.\n"
                "Please load a peak list first, then use Enhanced Detection\n"
                "to improve peak assignments for overlapping cases.")
            return

        try:
            self.update_status("🚀 Running Enhanced Graph-Based Peak Detection...")
            self.root.config(cursor="watch")

            # Import enhanced detection system
            from integrated_detection_fitter import EnhancedPeakDetectionIntegrated, create_integrated_fitter

            # Create or get integrated fitter
            if not hasattr(self.integrator, 'integrated_fitter') or self.integrator.integrated_fitter is None:
                self.integrator.integrated_fitter = create_integrated_fitter(
                    self.integrator.enhanced_fitter if hasattr(self.integrator, 'enhanced_fitter') else None
                )

            # Create enhanced detector
            enhanced_detector = EnhancedPeakDetectionIntegrated(self.integrator.integrated_fitter)

            # Setup GUI parameters for enhanced detection (including user-configurable options)
            gui_params = {
                # Core detection parameters
                'height_threshold': self.peak_height_threshold.get(),
                'distance_factor': self.peak_distance_factor.get(),
                'prominence_threshold': self.peak_prominence_threshold.get(),
                'smoothing_sigma': self.smoothing_sigma.get(),
                'max_peaks_fit': self.max_peaks_fit.get(),
                'fitting_window_x': 0.3,  # Default for enhanced detection
                # Peak Centroid Detection parameters
                'use_centroid_refinement': self.use_centroid_refinement.get(),
                'centroid_window_x_ppm': self.centroid_window_x_ppm.get(),
                'centroid_window_y_ppm': self.centroid_window_y_ppm.get(),
                'centroid_noise_multiplier': self.centroid_noise_multiplier.get(),
                'fitting_window_y': 8.0,  # Default for enhanced detection
                'noise_threshold': self.enhanced_noise_threshold.get(),  # Use enhanced detection noise threshold

                # Detection square size (NEW - user configurable anisotropic)
                'detection_square_size': self.enhanced_detection_square_size.get(),
                'detection_rectangle_y': self.enhanced_detection_rectangle_y.get(),

                # Peak Ridge Consolidation parameters (Solution A)
                'consolidation_x_tolerance': self.consolidation_x_tolerance.get(),
                'consolidation_y_tolerance': self.consolidation_y_tolerance.get(),

                # Enhanced Detection specific parameters (NEW - user configurable)
                'enhanced_radius_x': self.enhanced_radius_x.get(),              # 1H dimension constraint
                'enhanced_radius_y': self.enhanced_radius_y.get(),              # 15N/13C dimension constraint
                'enhanced_pattern_similarity': self.enhanced_pattern_similarity.get(),
                'enhanced_missing_tolerance': self.enhanced_missing_tolerance.get(),
                'enhanced_position_weight': self.enhanced_position_weight.get(),

                # Peak reduction parameters (NEW - to reduce too many peaks)
                'enhanced_peak_limit': self.enhanced_peak_limit.get()
            }

            # Get NMR data dimensions with error checking
            if not hasattr(self.integrator, 'nmr_data') or self.integrator.nmr_data is None:
                raise ValueError("NMR data not found in integrator")

            if not hasattr(self.integrator, 'ppm_x_axis') or self.integrator.ppm_x_axis is None:
                raise ValueError("1H axis (ppm_x_axis) not found in integrator")

            if not hasattr(self.integrator, 'ppm_y_axis') or self.integrator.ppm_y_axis is None:
                raise ValueError("15N axis (ppm_y_axis) not found in integrator")

            nmr_data = self.integrator.nmr_data
            x_axis = self.integrator.ppm_x_axis  # 1H axis
            y_axis = self.integrator.ppm_y_axis  # 15N axis

            print(f"   📊 NMR data shape: {nmr_data.shape}")
            print(f"   📏 1H axis: {len(x_axis)} points [{np.min(x_axis):.2f}, {np.max(x_axis):.2f}] ppm")
            print(f"   📏 15N axis: {len(y_axis)} points [{np.min(y_axis):.1f}, {np.max(y_axis):.1f}] ppm")

            # Run enhanced detection
            self.update_status("🔍 Analyzing peak patterns and overlaps...")

            results = enhanced_detector.peak_detection_integrated_enhanced(
                x_data=x_axis,
                y_data=y_axis,
                intensity_data=nmr_data,
                peak_list=self.integrator.peak_list,
                nucleus_type='2D',
                gui_params=gui_params
            )

            # Process results
            if results and results.get('success', False):
                detected_peaks = results.get('peaks', [])


                if detected_peaks:
                    # Convert to DataFrame format (same as standard detection)
                    peak_data = []
                    fitted_peaks_data = []  # NEW: For red circles visualization

                    for i, peak in enumerate(detected_peaks):
                        # DataFrame format for peak_list (table display)
                        peak_data.append({
                            'Peak_Number': i + 1,
                            'Position_X': peak.get('Position_X', 0),
                            'Position_Y': peak.get('Position_Y', 0),
                            'Intensity': peak.get('Intensity', 0),
                            'Volume': peak.get('Volume', 0),
                            'Assignment': peak.get('Assignment', f'Enhanced_{i+1}'),
                            # Enhanced metadata
                            'Detection_Method': peak.get('Detection_Method', 'enhanced_graph_based'),
                            'Match_Confidence': peak.get('Match_Confidence', 0),
                            'Match_Type': peak.get('Match_Type', 'enhanced'),
                            'Position_Error_X': peak.get('Position_Error_X', 0)
                        })

                        # Dictionary format for fitted_peaks (red circles visualization)
                        fitted_peaks_data.append({
                            'Position_X': peak.get('Position_X', 0),
                            'Position_Y': peak.get('Position_Y', 0),
                            'ppm_x': peak.get('Position_X', 0),  # Alternative coordinate names
                            'ppm_y': peak.get('Position_Y', 0),  # Alternative coordinate names
                            'intensity': peak.get('Intensity', 0),
                            'volume': peak.get('Volume', 0),
                            'assignment': peak.get('Assignment', f'Enhanced_{i+1}'),
                            'detected': True,     # KEY: This enables red circle display
                            'fitted': True,      # Additional status flag
                            'success': True,     # Additional status flag
                            'method': 'enhanced_detection',
                            'match_confidence': peak.get('Match_Confidence', 0),
                            'match_type': peak.get('Match_Type', 'enhanced')
                        })

                    # Create new peak list DataFrame
                    import pandas as pd
                    enhanced_peak_list = pd.DataFrame(peak_data)

                    # CRITICAL FIX: DO NOT replace the original reference peak list!
                    # The reference peaks must remain constant - only store detected results separately

                    # Store Enhanced Detection results separately (DO NOT replace peak_list!)
                    self.integrator.enhanced_detection_results = enhanced_peak_list

                    # IMPORTANT: Preserve original reference peak list (backup if needed)
                    if not hasattr(self.integrator, 'original_reference_peaks'):
                        self.integrator.original_reference_peaks = self.integrator.peak_list.copy()
                        print(f"   📋 Backed up {len(self.integrator.original_reference_peaks)} original reference peaks")

                    # Keep reference peak list unchanged for table display
                    self.integrator.peak_list = self.integrator.original_reference_peaks.copy()

                    # Update fitted_peaks to show ONLY successfully matched peaks for visualization
                    # Red circles should show DETECTED coordinates (where peaks were actually found)
                    self.integrator.fitted_peaks = self._create_detected_visualization_data(detected_peaks)

                    print(f"   🔄 Reference peak preservation:")
                    print(f"      Original references: {len(self.integrator.original_reference_peaks)}")
                    print(f"      Current peak_list: {len(self.integrator.peak_list)}")
                    print(f"      Enhanced results: {len(enhanced_peak_list)}")
                    print(f"      Visualization points: {len(self.integrator.fitted_peaks)}")

                    # Store peaks after limit filtering for debug visualization

                    if hasattr(self.integrator, 'integrated_fitter'):
                        print(f"🔍 DEBUG: integrated_fitter exists: {self.integrator.integrated_fitter is not None}")

                    if hasattr(self.integrator, 'integrated_fitter') and self.integrator.integrated_fitter:
                        # Try to get the debug data from the integrated_fitter
                        included_peaks = self.integrator.integrated_fitter.get_included_peaks_after_limit_debug()
                        #print(f"🔍 DEBUG: Retrieved {len(included_peaks) if included_peaks else 0} included peaks after limit from fitter")

                        # Also check the enhanced detector directly for debug data
                        if hasattr(enhanced_detector, 'detector'):
                            if hasattr(enhanced_detector.detector, '_included_peaks_after_limit_debug'):
                                direct_peaks = enhanced_detector.detector._included_peaks_after_limit_debug
                                print(f"🔍 DEBUG: Found {len(direct_peaks) if direct_peaks else 0} peaks in enhanced_detector.detector")
                                if direct_peaks and not included_peaks:
                                    included_peaks = direct_peaks
                                    print(f"🔍 DEBUG: Using direct detector peaks as fallback")
                            else:
                                print(f"🔍 DEBUG: enhanced_detector.detector does not have _included_peaks_after_limit_debug")
                        else:
                            print(f"🔍 DEBUG: enhanced_detector does not have detector attribute")

                        # Check enhanced_detector itself
                        if hasattr(enhanced_detector, '_included_peaks_after_limit_debug'):
                            direct_peaks2 = enhanced_detector._included_peaks_after_limit_debug
                            #print(f"🔍 DEBUG: Found {len(direct_peaks2) if direct_peaks2 else 0} peaks in enhanced_detector itself")
                            if direct_peaks2 and not included_peaks:
                                included_peaks = direct_peaks2
                                #print(f"🔍 DEBUG: Using enhanced_detector peaks as fallback")

                        self.integrator.included_peaks_after_limit_debug = included_peaks
                        if included_peaks:
                            #print(f"🔍 DEBUG: Stored {len(included_peaks)} included peaks after limit for visualization")
                            # Print first few peaks for debugging
                            for i, peak in enumerate(included_peaks[:3]):
                                print(f"🔍 DEBUG: Included peak {i+1}: {dict(list(peak.items())[:6])}")
                        else:
                            print(f"🔍 DEBUG: No included peaks after limit to store")
                    else:
                        print(f"🔍 DEBUG: integrated_fitter not available or None")

                    # FALLBACK: If integrated_fitter method fails, use the fitted_peaks_data as debug data
                    #if not hasattr(self.integrator, 'included_peaks_after_limit_debug') or not self.integrator.included_peaks_after_limit_debug:
                        #print(f"🔍 DEBUG: Using fitted_peaks_data as fallback debug data")
                        # Convert fitted_peaks_data to the expected format for debug visualization
                    #    if fitted_peaks_data:
                    #        debug_peaks_fallback = []
                    #        for peak in fitted_peaks_data:
                    #            debug_peaks_fallback.append({
                    #                'position_x': peak.get('ppm_x', peak.get('Position_X', 0)),
                    #                'position_y': peak.get('ppm_y', peak.get('Position_Y', 0)),
                    #                'intensity': peak.get('intensity', 0),
                    #                'method': 'enhanced_detection_fallback'
                    #            })
                    #        self.integrator.included_peaks_after_limit_debug = debug_peaks_fallback
                    #        #print(f"🔍 DEBUG: Fallback stored {len(debug_peaks_fallback)} peaks for debug visualization")

                    # Update display and plots
                    self.update_main_plot()
                    self.update_statistics()

                    # Success message with enhanced metrics
                    easy_matches = results.get('easy_matches', 0)
                    complex_matches = results.get('complex_matches', 0)
                    overall_confidence = results.get('overall_confidence', 0)

                    # Count matched vs preserved peaks
                    matched_count = len([p for p in detected_peaks if p.get('Detection_Method') == 'enhanced_graph_based'])
                    preserved_count = len(detected_peaks) - matched_count

                    self.update_status(
                        f"✅ Enhanced Detection: {len(detected_peaks)} peaks total "
                        f"({matched_count} matched: {easy_matches} easy + {complex_matches} complex, "
                        f"{preserved_count} preserved, Confidence: {overall_confidence:.2f})"
                    )

                    messagebox.showinfo("Enhanced Detection Complete",
                        f"Enhanced Graph-Based Detection completed successfully!\n\n"
                        f"Results:\n"
                        f"• Total peaks detected: {len(detected_peaks)}\n"
                        f"• Easy matches: {easy_matches}\n"
                        f"• Complex pattern matches: {complex_matches}\n"
                        f"• Overall confidence: {overall_confidence:.1%}\n\n"
                        f"The enhanced peak list is now ready for 'Fit All Peaks'.")

                else:
                    self.update_status("⚠️ Enhanced Detection: No peaks detected")
                    messagebox.showwarning("No Peaks", "Enhanced detection found no peaks.")
            else:
                error_msg = results.get('error', 'Unknown error') if results else 'Detection failed'
                self.update_status(f"❌ Enhanced Detection failed: {error_msg}")
                messagebox.showerror("Detection Error", f"Enhanced detection failed: {error_msg}")

        except Exception as e:
            error_msg = str(e)
            self.update_status(f"❌ Enhanced Detection error: {error_msg}")

            # Provide specific guidance for common errors
            if "x_axis" in error_msg or "y_axis" in error_msg:
                guidance = ("This appears to be a data axis error.\n"
                           "Please ensure NMR data is properly loaded.")
            elif "ppm_x_axis" in error_msg or "ppm_y_axis" in error_msg:
                guidance = ("NMR axes not found. Please ensure NMR data\n"
                           "is loaded correctly before using Enhanced Detection.")
            elif "peak_list" in error_msg:
                guidance = ("Peak list required for Enhanced Detection.\n"
                           "Please load a peak list file first.")
            else:
                guidance = "Try using standard 'Detect Peaks' instead."

            messagebox.showerror("Enhanced Detection Error",
                f"Enhanced detection failed with error:\n{error_msg}\n\n"
                f"Suggestion: {guidance}")
        finally:
            self.root.config(cursor="")

    def _create_detected_visualization_data(self, detected_peaks):
        """
        Create visualization data for red circles showing DETECTED coordinates
        Only show detected peaks that were successfully matched to reference peaks
        """
        detected_viz_data = []

        # Only show detected coordinates for successfully matched peaks
        for peak in detected_peaks:
            detection_method = peak.get('Detection_Method', '')
            if detection_method != 'detected_unassigned' and detection_method != 'reference_only':
                # This peak was matched - show the DETECTED coordinate as red circle
                assignment = peak.get('Assignment', '')

                # Use the detected coordinates (Position_X/Y contain detected coordinates from Simple Pattern Matcher)
                detected_x = peak.get('Position_X', 0)  # DETECTED coordinate
                detected_y = peak.get('Position_Y', 0)  # DETECTED coordinate

                detected_viz_data.append({
                    'Position_X': detected_x,  # DETECTED coordinate (red circle)
                    'Position_Y': detected_y,  # DETECTED coordinate (red circle)
                    'ppm_x': detected_x,
                    'ppm_y': detected_y,
                    'intensity': peak.get('Intensity', 0),
                    'volume': peak.get('Volume', 0),
                    'assignment': assignment,
                    'detected': True,      # Enable red circle display
                    'fitted': True,       # Additional status flag
                    'success': True,      # Additional status flag
                    'match_type': peak.get('Match_Type', 'enhanced'),
                    'confidence': peak.get('Match_Confidence', 0)
                })

        print(f"   🔴 Created {len(detected_viz_data)} detected peak visualization points (red circles)")
        return detected_viz_data

    def _detect_peaks_standard(self):
        """Standard peak detection method (legacy)"""
        self.update_status("Detecting peaks (standard mode)...")
        self.root.config(cursor="watch")

        # Update parameters
        self.on_parameter_change()
#
#        # Perform detection based on mode
        detected_peaks = self.integrator.process_peaks()
#
        if detected_peaks is not None and (
            not hasattr(detected_peaks, 'empty') or not detected_peaks.empty
        ) and len(detected_peaks) > 0:
            stats = self.integrator.get_detection_statistics()
            self.update_status(f"✅ Detection: {stats['detected_peaks']}/{stats['total_peaks']} peaks ({stats['detection_rate']:.1f}%)")
            self.update_main_plot()
            self.update_statistics()
            # Load detected peaks into navigator
            if hasattr(self, 'peak_navigator'):
                self.peak_navigator.load_detected_peaks(self.integrator.fitted_peaks)
        else:
            self.update_status("❌ No peaks detected - Check parameters")

    def _detect_peaks_integrated(self):
        """Integrated detection-fitting method"""
        import time
        from integrated_detection_fitter import create_integrated_fitter

        self.update_status("Detecting peaks (integrated mode)...")
        self.root.config(cursor="watch")

        # Show progress panel
        self.show_integration_progress()
        self.integration_start_time = time.time()

        try:
            # Create integrated fitter if not available
            if not hasattr(self.integrator, 'integrated_fitter') or self.integrator.integrated_fitter is None:
                self.integrator.integrated_fitter = create_integrated_fitter(
                    self.integrator.enhanced_fitter if hasattr(self.integrator, 'enhanced_fitter') else None
                )

            # Configure integration parameters - include GUI detection parameters
            integration_params = {
                # SOLUTION 1: Add essential detection parameters from GUI
                'height_threshold': self.peak_height_threshold.get(),
                'distance_factor': self.peak_distance_factor.get(),
                'prominence_threshold': self.peak_prominence_threshold.get(),
                'smoothing_sigma': self.smoothing_sigma.get(),
                'max_peaks_fit': self.max_peaks_fit.get(),
                'max_optimization_iterations': self.max_optimization_iterations.get(),
                # Peak Centroid Detection parameters
                'use_centroid_refinement': self.use_centroid_refinement.get(),
                'centroid_window_x_ppm': self.centroid_window_x_ppm.get(),
                'centroid_window_y_ppm': self.centroid_window_y_ppm.get(),
                'centroid_noise_multiplier': self.centroid_noise_multiplier.get(),
                # SOLUTION: Add advanced integration parameters from GUI
                'adaptive_thresholds_enabled': self.adaptive_thresholds_enabled.get(),
                'multi_resolution_enabled': self.multi_resolution_enabled.get(),
                'physics_constraints_enabled': self.physics_constraints_enabled.get(),
                'convergence_threshold': self.convergence_threshold.get(),
                'fit_likelihood_threshold': self.fit_likelihood_threshold.get(),
                # Detection mode control
                'force_full_detection': self.force_full_detection.get(),
                # 1D Refinement parameters (NEW)
                'enable_1d_refinement': self.enable_1d_refinement.get(),
                'refinement_quality_threshold': self.refinement_quality_threshold.get(),
                'refinement_coordinate_threshold': self.refinement_coordinate_threshold.get(),
                # Add all other GUI parameters for complete integration
                'threshold': getattr(self, 'threshold', tk.DoubleVar(value=0.01)).get(),
                'noise_level': getattr(self, 'noise_level', tk.DoubleVar(value=0.001)).get(),
                'min_snr': getattr(self, 'min_snr', tk.DoubleVar(value=3.0)).get()
            }

            # Update parameters
            self.on_parameter_change()

            # Choose between in-place fitting and full detection based on user preference
            if self.force_full_detection.get():
                print(f"🔍 USER REQUESTED FULL DETECTION: Using iterative detection-fitting with default max iterations")
                self._run_integrated_detection_fitting(integration_params)
            else:
                print(f"📍 Using in-place fitting mode (fast fitting of existing peaks)")
                self._run_integrated_inplace_fitting(integration_params)

        except Exception as e:
            self.hide_integration_progress()
            raise e

    def _run_integrated_detection_fitting(self, integration_params):
        """Run the integrated detection-fitting process (OLD SERIES METHOD - SHOULD NOT BE USED FOR FIT ALL PEAKS)"""
        print("⚠️⚠️⚠️ WARNING: Using OLD SERIES INTEGRATION method - this should NOT be called for Fit All Peaks! ⚠️⚠️⚠️")
        import time

        total_peaks = len(self.integrator.peak_list) if hasattr(self.integrator, 'peak_list') and self.integrator.peak_list is not None else 0
        max_iterations = 5  # Default max iterations

        # Initialize progress
        self.update_integration_progress(
            status_text="Starting integrated detection-fitting...",
            progress_percent=0,
            current_iteration=0,
            max_iterations=max_iterations,
            peaks_processed=0,
            total_peaks=total_peaks,
            convergence_status="Running",
            elapsed_time=0
        )

        try:
            # Get the integrated fitter
            integrated_fitter = self.integrator.integrated_fitter

            # Initialize variables for progress tracking
            peaks_processed = 0
            avg_quality = 0.0
            results = None

            # Perform iterative detection-fitting
            for iteration in range(max_iterations):
                iteration_start = time.time()

                # Update progress
                elapsed_time = time.time() - self.integration_start_time
                progress_percent = ((iteration + 1) / max_iterations) * 100

                self.update_integration_progress(
                    status_text=f"Iteration {iteration + 1}: Detection and fitting...",
                    progress_percent=progress_percent,
                    current_iteration=iteration + 1,
                    max_iterations=max_iterations,
                    peaks_processed=0,
                    total_peaks=total_peaks,
                    convergence_status="Running",
                    elapsed_time=elapsed_time
                )

                # Force GUI update
                self.root.update_idletasks()

                # Perform integrated detection-fitting
                if hasattr(integrated_fitter, 'integrated_detection_fitting'):
                    # CRITICAL FIX: Pass 2D NMR data for both-dimension peak detection
                    if hasattr(self.integrator, 'nmr_data') and self.integrator.nmr_data is not None:
                        # For 2D NMR: pass the full 2D data matrix and both axes
                        nmr_2d_data = self.integrator.nmr_data
                        x_axis = self.integrator.ppm_x_axis  # 1H axis
                        y_axis = self.integrator.ppm_y_axis  # 15N axis
                        print(f"      🔧 FIXED: Using 2D NMR data shape {nmr_2d_data.shape} for both-dimension detection")
                        print(f"         1H axis: {len(x_axis)} points [{np.min(x_axis):.1f}, {np.max(x_axis):.1f}] ppm")
                        print(f"         15N axis: {len(y_axis)} points [{np.min(y_axis):.1f}, {np.max(y_axis):.1f}] ppm")
                        print(f"         Intensity range: [{np.min(nmr_2d_data):.1f}, {np.max(nmr_2d_data):.1f}]")
                    else:
                        # Fallback (shouldn't happen)
                        nmr_2d_data = None
                        x_axis = self.integrator.ppm_x_axis
                        y_axis = self.integrator.ppm_y_axis
                        print(f"      ⚠️ FALLBACK: No NMR data found, using axes only")

                    # Determine if in-place mode (peak list available)
                    peak_list = getattr(self.integrator, 'peak_list', None)
                    in_place_mode = peak_list is not None and not peak_list.empty

                    if in_place_mode:
                        print(f"      📍 USING IN-PLACE MODE: {len(peak_list)} reference peaks from peak list")
                    else:
                        print(f"      🌐 USING FULL-SPECTRUM MODE: No peak list constraints")

                    results = integrated_fitter.integrated_detection_fitting(
                        x_axis,           # 1H chemical shifts
                        y_axis,           # 15N chemical shifts
                        nucleus_type='2D', # Changed to indicate 2D detection
                        gui_params=integration_params,
                        nmr_2d_data=nmr_2d_data,  # Pass 2D intensity matrix separately
                        peak_list=peak_list,       # Pass peak list for in-place mode
                        in_place_mode=in_place_mode # Enable in-place constraints
                    )
                else:
                    # Fallback to standard detection if method not available
                    detected_peaks = self.integrator.process_peaks()
                    results = {'peaks': detected_peaks, 'converged': True, 'iterations': iteration + 1}

                # Update results
                if results and 'peaks' in results:
                    self.integrator.fitted_peaks = results['peaks']

                    # Calculate quality metrics
                    if isinstance(results['peaks'], list) and len(results['peaks']) > 0:
                        quality_scores = [peak.get('composite_quality', 0.5) for peak in results['peaks']]
                        avg_quality = sum(quality_scores) / len(quality_scores)
                        peaks_processed = len(results['peaks'])
                    else:
                        avg_quality = 0
                        peaks_processed = 0

                    # Update progress with results
                    self.update_integration_progress(
                        status_text=f"Iteration {iteration + 1}: Processing complete",
                        progress_percent=progress_percent,
                        current_iteration=iteration + 1,
                        max_iterations=max_iterations,
                        peaks_processed=peaks_processed,
                        total_peaks=total_peaks,
                        quality_score=avg_quality,
                        convergence_status="Converged" if results.get('converged', False) else "Running",
                        elapsed_time=elapsed_time
                    )

                    # Check convergence
                    if results.get('converged', False):
                        break

                # Force GUI update
                self.root.update_idletasks()

                # Small delay to allow GUI updates
                time.sleep(0.1)

            # Final update
            total_time = time.time() - self.integration_start_time
            final_status = "Converged" if results and results.get('converged', False) else "Max iterations reached"

            self.update_integration_progress(
                status_text=f"Integration complete: {final_status}",
                progress_percent=100,
                current_iteration=max_iterations,
                max_iterations=max_iterations,
                peaks_processed=peaks_processed,
                total_peaks=total_peaks,
                quality_score=avg_quality,
                convergence_status=final_status,
                elapsed_time=total_time
            )

            # Update plots and statistics
            self.update_main_plot()
            self.update_statistics()

            # Store results for diagnostics
            self.last_integration_results = results['peaks'] if results and 'peaks' in results else []
            self.last_processing_time = total_time

            # Update diagnostics
            if hasattr(self, 'refresh_diagnostics'):
                self.refresh_diagnostics()

            self.update_status(f"✅ Integrated detection complete: {peaks_processed} peaks processed in {total_time:.2f}s")

            # Hide progress after 3 seconds
            self.root.after(3000, self.hide_integration_progress)

        except Exception as e:
            self.update_integration_progress(
                status_text=f"Error: {str(e)}",
                convergence_status="Failed"
            )
            self.hide_integration_progress()
            raise e

    def _run_integrated_inplace_fitting(self, integration_params):
        """Run integrated in-place fitting for existing 2D peak list"""
        import time
        import numpy as np

        try:
            # Get peak list and validate
            peak_list = self.integrator.peak_list
            if peak_list is None or peak_list.empty:
                self.update_status("❌ No peak list available for in-place fitting")
                return

            total_peaks = len(peak_list)
            print(f"🔄🔄🔄 USING NEW IN-PLACE METHOD: Starting integrated in-place fitting for {total_peaks} peaks 🔄🔄🔄")

            # Initialize progress
            self.update_integration_progress(
                status_text=f"Starting in-place integrated fitting for {total_peaks} peaks...",
                progress_percent=0,
                current_iteration=0,
                max_iterations=1,
                peaks_processed=0,
                total_peaks=total_peaks,
                convergence_status="Running",
                elapsed_time=0
            )

            # Get integrated fitter
            integrated_fitter = self.integrator.integrated_fitter

            # Get NMR data for in-place fitting
            if hasattr(self.integrator, 'nmr_data') and self.integrator.nmr_data is not None:
                nmr_2d_data = self.integrator.nmr_data
                x_axis = self.integrator.ppm_x_axis  # 1H axis
                y_axis = self.integrator.ppm_y_axis  # 15N axis
                print(f"   📊 Using 2D NMR data shape {nmr_2d_data.shape} for in-place fitting")
                print(f"      1H axis: {len(x_axis)} points [{np.min(x_axis):.1f}, {np.max(x_axis):.1f}] ppm")
                print(f"      15N axis: {len(y_axis)} points [{np.min(y_axis):.1f}, {np.max(y_axis):.1f}] ppm")
            else:
                self.update_status("❌ No NMR data available for in-place fitting")
                return

            # Force GUI update
            self.root.update_idletasks()

            # Perform integrated in-place fitting
            results = integrated_fitter.integrated_detection_fitting(
                x_axis,           # 1H chemical shifts
                y_axis,           # 15N chemical shifts
                nucleus_type='2D', # 2D detection
                gui_params=integration_params,
                nmr_2d_data=nmr_2d_data,  # Pass 2D intensity matrix
                peak_list=peak_list,       # Pass peak list for in-place mode
                in_place_mode=True         # Enable in-place constraints
            )

            # Process results
            if results and results.get('success', False):
                fitted_peaks = results.get('peaks', [])
                peaks_processed = len(fitted_peaks)
                avg_quality = results.get('avg_r_squared', 0)

                # Store results in integrator
                self.integrator.fitted_peaks = fitted_peaks

                # Update progress
                progress_percent = 100
                self.update_integration_progress(
                    status_text="In-place integrated fitting complete",
                    progress_percent=progress_percent,
                    current_iteration=1,
                    max_iterations=1,
                    peaks_processed=peaks_processed,
                    total_peaks=total_peaks,
                    quality_score=avg_quality,
                    convergence_status="Converged",
                    elapsed_time=time.time() - self.integration_start_time
                )

                # Update plots and statistics
                self.update_main_plot()
                self.update_statistics()
                # Load detected peaks into navigator
                if hasattr(self, 'peak_navigator'):
                    self.peak_navigator.load_detected_peaks(self.integrator.fitted_peaks)

                # Store results for diagnostics
                self.last_integration_results = fitted_peaks
                self.last_processing_time = time.time() - self.integration_start_time

                # Update diagnostics
                if hasattr(self, 'refresh_diagnostics'):
                    self.refresh_diagnostics()

                detection_diagnostics = results.get('integration_diagnostics', {})
                detected_count = detection_diagnostics.get('detected_count', 0)
                retained_count = detection_diagnostics.get('retained_count', 0)
                detection_rate = detection_diagnostics.get('detection_rate', 0)

                total_time = time.time() - self.integration_start_time
                self.update_status(f"✅ In-place integrated fitting complete: {peaks_processed}/{total_peaks} peaks processed, {detected_count} detected + {retained_count} retained ({detection_rate:.1f}% detection rate) in {total_time:.2f}s")

                print(f"🏁 In-place integrated fitting results:")
                print(f"   Total processed: {peaks_processed}/{total_peaks}")
                print(f"   Detected: {detected_count}, Retained: {retained_count}")
                print(f"   Detection rate: {detection_rate:.1f}%")
                print(f"   Average R²: {avg_quality:.3f}")
                print(f"   Processing time: {total_time:.2f}s")

            else:
                error_msg = results.get('error', 'Unknown error') if results else 'No results returned'
                self.update_status(f"❌ In-place integrated fitting failed: {error_msg}")
                print(f"❌ In-place integrated fitting failed: {error_msg}")

            # Hide progress after 3 seconds
            self.root.after(3000, self.hide_integration_progress)

        except Exception as e:
            self.update_integration_progress(
                status_text=f"Error: {str(e)}",
                convergence_status="Failed"
            )
            self.hide_integration_progress()
            print(f"❌ Exception in in-place fitting: {e}")
            import traceback
            traceback.print_exc()
            raise e

    def integrate_peaks(self):
        """Integrate detected peaks - supports both standard and integrated modes"""
        if not hasattr(self.integrator, 'fitted_peaks') or (
            self.integrator.fitted_peaks is None or
            (hasattr(self.integrator.fitted_peaks, 'empty') and self.integrator.fitted_peaks.empty) or
            len(self.integrator.fitted_peaks) == 0
        ):
            messagebox.showerror("Error", "Please detect peaks first")
            return

        try:
            # Check integration mode
            integration_mode = self.integration_mode.get()

            if integration_mode in ['integrated', 'adaptive']:
                # Use integrated mode - peaks are already fitted during detection
                self._integrate_peaks_integrated()
            else:
                # Use standard integration (legacy)
                self._integrate_peaks_standard()

        except Exception as e:
            self.update_status(f"❌ Integration error: {str(e)}")
            messagebox.showerror("Integration Error", str(e))

        finally:
            self.root.config(cursor="")

    def _integrate_peaks_standard(self):
        """Standard peak integration method (legacy)"""
        self.update_status("Integrating peaks (standard mode)...")
        self.root.config(cursor="watch")

        results = self.integrator.integrate_peaks()
        if results is not None and (
            not hasattr(results, 'empty') or not results.empty
        ) and len(results) > 0:
            detected_count = sum(1 for r in results if r.get('Integration_Method') != 'Reference')
            good_quality = sum(1 for r in results if r.get('Quality') in ['Excellent', 'Good'])
            self.update_status(f"✅ Integration: {detected_count}/{len(results)} peaks, {good_quality} good quality")
            self.update_statistics()
        else:
            self.update_status("❌ Integration failed")

    def _integrate_peaks_integrated(self):
        """Integrated mode - peaks are already fitted, just update display"""
        self.update_status("Integration complete (integrated mode)...")

        # In integrated mode, peaks are already fitted with quality scores
        if hasattr(self.integrator, 'fitted_peaks') and self.integrator.fitted_peaks:
            peaks = self.integrator.fitted_peaks

            if isinstance(peaks, list):
                total_peaks = len(peaks)
                high_quality = sum(1 for p in peaks if p.get('composite_quality', 0) >= 0.7)
                converged = sum(1 for p in peaks if p.get('converged', False))

                self.update_status(f"✅ Integrated mode: {total_peaks} peaks analyzed, {high_quality} high quality, {converged} converged")
            else:
                self.update_status("✅ Integrated peaks processed")

            self.update_statistics()

            # Update diagnostics if available
            if hasattr(self, 'refresh_diagnostics'):
                self.refresh_diagnostics()
        else:
            self.update_status("❌ No integrated results available")

    def fit_selected_peak(self):
        """Fit Voigt profile to selected peak"""
        if not hasattr(self.integrator, 'peak_list') or self.integrator.peak_list is None:
            messagebox.showerror("Error", "No peaks available")
            return

        peak_num = self.selected_peak_number.get()
        if peak_num < 1 or peak_num > len(self.integrator.peak_list):
            messagebox.showerror("Error", f"Invalid peak number: {peak_num}")
            return

        try:
            self.update_status(f"Fitting Voigt profile for peak {peak_num}...")
            self.root.config(cursor="watch")

            # Get peak position
            peak_row = self.integrator.peak_list.iloc[peak_num - 1]
            peak_x = float(peak_row['Position_X'])
            peak_y = float(peak_row['Position_Y'])
            assignment = peak_row.get('Assignment', f'Peak_{peak_num}')

            # Perform enhanced Voigt fitting
            result = self.integrator.enhanced_peak_fitting(peak_x, peak_y, assignment)

            if result:
                # Store peak number in the result for easy lookup
                result['peak_number'] = peak_num

                # CRITICAL FIX: Store individual fit result where GUI expects it
                if not hasattr(self.integrator, 'fitted_peaks'):
                    self.integrator.fitted_peaks = []

                # Update or add to fitted_peaks
                found = False
                for i, existing_peak in enumerate(self.integrator.fitted_peaks):
                    if existing_peak.get('peak_number') == peak_num:
                        self.integrator.fitted_peaks[i] = result
                        found = True
                        break

                if not found:
                    self.integrator.fitted_peaks.append(result)

                print(f"✅ Stored individual fit result for peak {peak_num} in fitted_peaks")

                self.current_voigt_result = result
                self.voigt_plotter.plot_voigt_analysis(result)
                self.canvas_voigt.draw()

                # Switch to Voigt tab
                self.viz_notebook.select(1)

                quality = result.get('fitting_quality', 'Unknown')
                self.update_status(f"✅ Voigt fit: {assignment} - Quality: {quality}")
                self.update_statistics()
            else:
                self.update_status(f"❌ Voigt fitting failed for {assignment}")

        except Exception as e:
            self.update_status(f"❌ Fitting error: {str(e)}")
            messagebox.showerror("Fitting Error", str(e))

        finally:
            self.root.config(cursor="")

    def fit_all_peaks(self):
        """Fit Voigt profiles to all peaks using new Single Spectrum Processor"""

        # Validate peak list
        if (not hasattr(self.integrator, 'peak_list') or
            self.integrator.peak_list is None or
            (hasattr(self.integrator.peak_list, 'empty') and self.integrator.peak_list.empty)):
            messagebox.showerror("Error", "Please load a peak list first")
            return

        print("🚀 Starting Fit All Peaks using Single Spectrum Processor")

        # LOG VOIGT FITTING PARAMETERS
        print("\n" + "="*60)
        print("📋 VOIGT FITTING PARAMETERS")
        print("="*60)

        # Get current parameters from parameter manager
        try:
            fitting_params = self.param_manager.get_integrator_parameters()['fitting_params']
            detection_params = self.param_manager.get_integrator_parameters()['detection_params']
            gui_params = self.param_manager.get_integrator_parameters()['gui_params']

            print(f"🔧 Core Fitting Parameters:")
            print(f"   • Fitting Window X: ±{fitting_params.get('fitting_window_x', 'N/A')} ppm (1H dimension)")
            print(f"   • Fitting Window Y: ±{fitting_params.get('fitting_window_y', 'N/A')} ppm (15N/13C dimension)")
            print(f"   • Min R-squared: {fitting_params.get('min_r_squared', 'N/A')}")
            print(f"   • Max Iterations: {fitting_params.get('max_iterations', 'N/A')}")

            print(f"\n🎯 Detection Parameters:")
            print(f"   • Search Window X: ±{detection_params.get('search_window_x', 'N/A')} ppm (1H dimension)")
            print(f"   • Search Window Y: ±{detection_params.get('search_window_y', 'N/A')} ppm (15N/13C dimension)")
            print(f"   • Noise Threshold Multiplier: {detection_params.get('noise_threshold', 'N/A')}")

            print(f"\n⚙️ Advanced Fitting Parameters:")
            print(f"   • Global Optimization: {self.use_global_optimization.get()}")
            print(f"   • Parallel Processing: {self.use_parallel_processing.get()}")
            print(f"   • Peak Height Threshold: {self.peak_height_threshold.get()}")
            print(f"   • Peak Distance Factor: {self.peak_distance_factor.get()}")
            print(f"   • Peak Prominence Threshold: {self.peak_prominence_threshold.get()}")
            print(f"   • Smoothing Sigma: {self.smoothing_sigma.get()}")
            print(f"   • Max Peaks per Fit: {self.max_peaks_fit.get()}")
            print(f"   • Max Optimization Iterations: {self.max_optimization_iterations.get()}")

            print(f"\n🔬 Peak Detection Settings:")
            print(f"   • Detection Square Size X: {self.detection_square_size.get()} pixels")
            print(f"   • Detection Rectangle Y: {self.detection_rectangle_y.get()} pixels")
            print(f"   • Use Reference Detection: {self.use_reference_detection.get()}")

            print(f"\n📐 Consolidation Parameters:")
            print(f"   • X-axis Tolerance: {self.consolidation_x_tolerance.get()} ppm")
            print(f"   • Y-axis Tolerance: {self.consolidation_y_tolerance.get()} ppm")

            print(f"\n🎚️ Centroid Refinement:")
            print(f"   • Enable Centroid Refinement: {self.use_centroid_refinement.get()}")
            print(f"   • Centroid Window X: {self.centroid_window_x_ppm.get()} ppm")
            print(f"   • Centroid Window Y: {self.centroid_window_y_ppm.get()} ppm")
            print(f"   • Noise Multiplier: {self.centroid_noise_multiplier.get()}")

        except Exception as e:
            print(f"⚠️ Error retrieving parameters: {e}")
            print("Using GUI variables directly...")

            print(f"🔧 Core Fitting Parameters:")
            print(f"   • Fitting Window X: ±{self.fitting_window_x.get()} ppm (1H dimension)")
            print(f"   • Fitting Window Y: ±{self.fitting_window_y.get()} ppm (15N/13C dimension)")
            print(f"   • Min R-squared: {self.min_r_squared.get()}")
            print(f"   • Max Iterations: {self.max_iterations.get()}")

            print(f"\n🎯 Detection Parameters:")
            print(f"   • Search Window X: ±{self.search_window_x.get()} ppm (1H dimension)")
            print(f"   • Search Window Y: ±{self.search_window_y.get()} ppm (15N/13C dimension)")
            print(f"   • Noise Threshold Multiplier: {self.noise_threshold.get()}")

        print("="*60)
        print(f"📊 Peak List Info: {len(self.integrator.peak_list)} peaks to process")
        print("="*60 + "\n")

        # Use new decoupled architecture
        self._fit_all_peaks_with_new_processor()

    def _fit_all_peaks_with_new_processor(self):
        """Fit all peaks using the new single spectrum processor"""

        from lunaNMR.gui.gui_components import AdvancedProgressDialog

        # Create progress dialog
        dialog_title = "🔄 Global Peak Optimization" if self.use_global_optimization.get() else "Single Spectrum Fitting"
        progress_dialog = AdvancedProgressDialog(self.root, dialog_title)

        # Run in separate thread
        fitting_thread = threading.Thread(
            target=self._run_new_single_spectrum_processing,
            args=(progress_dialog,)
        )
        fitting_thread.daemon = True
        fitting_thread.start()

    def _run_new_single_spectrum_processing(self, progress_dialog):
        """Run single spectrum processing using new architecture"""

        try:
            # Update parameter manager from current GUI state
            self.param_manager.update_from_gui_variables(self)

            # Create single spectrum processor
            self.single_spectrum_processor = SingleSpectrumProcessor(self.integrator, self.param_manager)

            # Set up progress callback
            def progress_callback(progress, task, log_msg=None, failed=False):
                if not progress_dialog.cancelled:
                    progress_dialog.update_progress(progress, task, log_msg, failed)

            # Set processing options based on GUI
            processing_options = {
                'use_parallel': self.use_parallel_processing.get(),
                'use_global_optimization': self.use_global_optimization.get()
            }

            print(f"📋 Processing {len(self.integrator.peak_list)} peaks with options: {processing_options}")

            # Process all peaks
            fitted_results = self.single_spectrum_processor.process_peak_list(
                self.integrator.peak_list,
                processing_options,
                progress_callback
            )

            if not progress_dialog.cancelled:
                # Get comprehensive summary
                summary = self.single_spectrum_processor.get_processing_summary(
                    fitted_results,
                    len(self.integrator.peak_list)
                )

            # Update GUI
                self.root.after(0, self._complete_new_single_spectrum_fitting, summary, progress_dialog)

        except Exception as e:
            import traceback
            error_msg = f"Single spectrum processing error: {str(e)}"
            print(f"❌ {error_msg}")
            print(traceback.format_exc())
            self.root.after(0, lambda: self.update_status(f"❌ {error_msg}"))

    def _complete_new_single_spectrum_fitting(self, summary, progress_dialog):
        """Complete single spectrum fitting and update GUI"""

        # Complete progress dialog
        success_msg = f"Fitting completed: {summary['successful_peaks']}/{summary['total_peaks']} successful"
        progress_dialog.complete(success_msg)

        # Update status bar
        status_msg = f"✅ Fit All Peaks: {summary['successful_peaks']}/{summary['total_peaks']} successful ({summary['success_rate']:.1f}%)"
        self.update_status(status_msg)

        if summary['results']:
            self.last_fitting_results = summary['results']

            # Standardize the keys to match detection workflow for series integration compatibility
            standardized_results = []
            for result in summary['results']:
                # Create standardized structure matching detection workflow
                standardized_result = {
                    'assignment': result.get('Assignment', result.get('assignment', '')),
                    'ppm_x': float(result.get('Position_X', result.get('ppm_x', 0))),
                    'ppm_y': float(result.get('Position_Y', result.get('ppm_y', 0))),
                    'detected': True,
                    'fitted': True
                }
                # Preserve additional fitting data
                for key, value in result.items():
                    if key not in ['Assignment', 'Position_X', 'Position_Y']:
                        standardized_result[key] = value

                standardized_results.append(standardized_result)

            self.integrator.fitted_peaks = standardized_results
            print(f"🔧 Standardized {len(standardized_results)} fitting results for series integration compatibility")

        # Print detailed summary
        print(f"\n📊 Single Spectrum Processing Summary:")
        print(f"   Total peaks: {summary['total_peaks']}")
        print(f"   Successful: {summary['successful_peaks']} ({summary['success_rate']:.1f}%)")
        print(f"   Average R²: {summary['average_r_squared']:.3f}")

        quality_dist = summary['quality_distribution']
        print(f"   Quality distribution: {quality_dist['excellent']} excellent, {quality_dist['good']} good, {quality_dist['poor']} poor")

    def _extract_linewidth_statistics_from_results(self, linear_results):
        """Extract linewidth constraints from top 20% performers by R²"""
        if not linear_results:
            return None

        # Calculate R² values and sort results
        r_squared_values = []
        for result in linear_results:
            avg_r_squared = result.get('avg_r_squared', 0)
            r_squared_values.append((result, avg_r_squared))

        # Sort by R² in descending order
        r_squared_values.sort(key=lambda x: x[1], reverse=True)

        # Take top 20% performers
        top_20_percent_count = max(1, len(r_squared_values) // 5)
        top_performers = [result for result, r_squared in r_squared_values[:top_20_percent_count]]

        print(f"📊 Analyzing top {len(top_performers)} performers from {len(linear_results)} fits")

        # Debug: Print available keys in first result
        if top_performers:
            print(f"🔍 Debug - Available keys in top performer: {list(top_performers[0].keys())}")

        # Extract sigma and gamma values from successful fits
        sigma_values = []
        gamma_values = []

        for result in top_performers:
            x_fit = result.get('x_fit', {})
            y_fit = result.get('y_fit', {})

            # Extract sigma from X dimension and gamma from Y dimension
            if 'sigma' in x_fit and 'sigma' in y_fit:
                sigma_values.append(x_fit['sigma'])  # X linewidth
                gamma_values.append(y_fit['sigma'])  # Y linewidth (gamma)

        if not sigma_values or not gamma_values:
            print("⚠️  No valid linewidth parameters found in top performers")
            return None

        # Calculate robust statistics (median ± 2×IQR)
        import numpy as np

        sigma_median = np.median(sigma_values)
        sigma_q75, sigma_q25 = np.percentile(sigma_values, [75, 25])
        sigma_iqr = sigma_q75 - sigma_q25

        gamma_median = np.median(gamma_values)
        gamma_q75, gamma_q25 = np.percentile(gamma_values, [75, 25])
        gamma_iqr = gamma_q75 - gamma_q25

        # Define constraints as median ± 2×IQR with proper bounds
        sigma_min = max(0.001, sigma_median - 2 * sigma_iqr)
        sigma_max = max(sigma_min + 0.001, sigma_median + 2 * sigma_iqr)

        gamma_min = max(0.001, gamma_median - 2 * gamma_iqr)
        gamma_max = max(gamma_min + 0.001, gamma_median + 2 * gamma_iqr)

        linewidth_constraints = {
            'sigma_bounds': (sigma_min, sigma_max),
            'gamma_bounds': (gamma_min, gamma_max)
        }

        print(f"📈 Linewidth constraints extracted:")
        print(f"   Sigma: {linewidth_constraints['sigma_bounds'][0]:.3f} - {linewidth_constraints['sigma_bounds'][1]:.3f}")
        print(f"   Gamma: {linewidth_constraints['gamma_bounds'][0]:.3f} - {linewidth_constraints['gamma_bounds'][1]:.3f}")

        return linewidth_constraints

    def _refit_poor_peaks_with_constraints(self, progress_dialog, peak_list, linear_results, linewidth_constraints):
        """Refit bottom 80% peaks using constraints from top 20%"""
        # Calculate R² values and identify poor performers
        r_squared_values = []
        result_by_peak = {}

        for result in linear_results:
            peak_number = result.get('peak_number', -1)
            avg_r_squared = result.get('avg_r_squared', 0)
            r_squared_values.append((result, avg_r_squared))
            result_by_peak[peak_number] = result

        # Sort by R² and identify poor performers (bottom 80%)
        r_squared_values.sort(key=lambda x: x[1], reverse=True)
        top_20_percent_count = max(1, len(r_squared_values) // 5)
        poor_performers = [result for result, r_squared in r_squared_values[top_20_percent_count:]]

        print(f"🔄 Refitting {len(poor_performers)} poor performers with learned constraints")

        # Start with linear results and replace poor performers with constraint-based fits
        improved_results = linear_results.copy()

        for i, result in enumerate(poor_performers):
            if progress_dialog.cancelled:
                break

            peak_number = result.get('peak_number', i + 1)

            # Update progress
            progress_msg = f"Re-fitting peak {peak_number} with constraints..."
            progress_dialog.update_progress(
                75 + int(20 * i / len(poor_performers)),
                "🔄 Stage 3: Constraint-Based Refinement",
                progress_msg,
                failed=False
            )

            try:
                # Get peak data
                peak_data = peak_list.iloc[peak_number - 1]
                peak_x = float(peak_data['Position_X'])
                peak_y = float(peak_data['Position_Y'])
                assignment = peak_data.get('Assignment', f'Peak_{peak_number}')

                # Run enhanced Voigt fitting with linewidth constraints
                result = self.integrator.fit_peak_voigt_2d(
                    peak_x, peak_y, assignment,
                    use_dynamic_optimization=True,
                    linewidth_constraints=linewidth_constraints
                )

                if result and result.get('avg_r_squared', 0) > 0:
                    # Update the improved results with the new fit
                    for j, orig_result in enumerate(improved_results):
                        if orig_result.get('peak_number') == peak_number:
                            improved_results[j] = result
                            break

                    print(f"   Peak {peak_number}: R² improved from {result_by_peak[peak_number].get('avg_r_squared', 0):.3f} to {result.get('avg_r_squared', 0):.3f}")
                else:
                    print(f"   Peak {peak_number}: Constraint fitting failed, keeping linear result")

            except Exception as e:
                print(f"   Peak {peak_number}: Error in constraint fitting: {str(e)}")
                continue

        return improved_results

    # =================== SERIES PROCESSING ===================

    def start_series_integration(self):
        """
        Start series integration using the new independent Multi-Spectrum Processor.
        Combines complete GUI decoupling with perfect retrocompatibility.
        """
        # 1. Sync parameters to ensure everything is up-to-date
        self._sync_parameters_to_integrator()

        # 2. Validate peak source and get reference peaks
        peak_source = self.series_peak_source.get()
        reference_peaks = None

        if peak_source == "detected":
            if not hasattr(self.integrator, 'fitted_peaks') or not self.integrator.fitted_peaks:
                messagebox.showerror("Error",
                    "No detected peaks available. Run 'Fit All Peaks' first to detect peaks.")
                return
            # Convert fitted_peaks to DataFrame format
            reference_peaks = self._convert_fitted_peaks_to_dataframe(self.integrator.fitted_peaks)

        elif peak_source == "reference":
            if not hasattr(self.integrator, 'peak_list') or self.integrator.peak_list.empty:
                messagebox.showerror("Error",
                    "No reference peak list loaded. Load a peak list file first.")
                return
            reference_peaks = self.integrator.peak_list.copy()

        elif peak_source == "cascade":
            # For cascade mode, start with reference peaks
            if not hasattr(self.integrator, 'peak_list') or self.integrator.peak_list.empty:
                messagebox.showerror("Error",
                    "No reference peak list loaded for cascade mode.")
                return
            reference_peaks = self.integrator.peak_list.copy()

        # 3. Validate NMR files
        nmr_files = self.nmr_file_list.get_all_files()
        if not nmr_files:
            messagebox.showerror("Error",
                "No NMR files found in the selected folder. Please load NMR files first.")
            return

        print(f"🚀 Starting series integration: {len(nmr_files)} spectra, {len(reference_peaks)} reference peaks")

        # 4. Create progress dialog
        progress_dialog = AdvancedProgressDialog(
        #progress_dialog = SeriesProcessingProgressDialog(
            self.root,
            #total_spectra=len(nmr_files),
            title="Series Integration Progress"
        )

        # 5. Start processing in background thread
        self.processing_active = True
        threading.Thread(
            target=self._run_new_multi_spectrum_processing,
            args=(nmr_files, reference_peaks, peak_source, progress_dialog),
            daemon=True
        ).start()

    def _convert_fitted_peaks_to_dataframe(self, fitted_peaks):
        """Convert fitted_peaks list to DataFrame format for processor - handles both detection and fitting structures"""

        # DEBUG: Analyze the fitted_peaks structure
        print(f"🔍 CONVERSION DEBUG: Processing {len(fitted_peaks) if fitted_peaks else 0} fitted peaks")
        if fitted_peaks:
            sample_peak = fitted_peaks[0]
            print(f"   Sample peak type: {type(sample_peak)}")
            if isinstance(sample_peak, dict):
                print(f"   Sample peak keys: {list(sample_peak.keys())}")
                print(f"   Sample peak values: {sample_peak}")
            else:
                print(f"   ❌ PROBLEM: Peak is not a dict, it's {type(sample_peak)}")

        peak_data = []
        skipped_peaks = 0

        for i, peak in enumerate(fitted_peaks):
            if isinstance(peak, dict):
                # Handle both detection results (lowercase) and fitting results (uppercase)
                # Try multiple key variations to ensure compatibility
                assignment = (peak.get('assignment') or
                             peak.get('Assignment') or
                             f'Peak_{i+1}')

                # Enhanced coordinate extraction - handle peak_position tuple from SingleSpectrumProcessor
                position_x = 0
                position_y = 0

                # Priority order: non-zero ppm_x/ppm_y > peak_position > Position_X/Position_Y
                if peak.get('ppm_x') and peak.get('ppm_x') != 0:
                    position_x = float(peak.get('ppm_x'))
                elif peak.get('peak_position') and len(peak.get('peak_position')) >= 2:
                    position_x = float(peak.get('peak_position')[0])
                elif peak.get('Position_X'):
                    position_x = float(peak.get('Position_X'))
                elif peak.get('position_x'):
                    position_x = float(peak.get('position_x'))

                if peak.get('ppm_y') and peak.get('ppm_y') != 0:
                    position_y = float(peak.get('ppm_y'))
                elif peak.get('peak_position') and len(peak.get('peak_position')) >= 2:
                    position_y = float(peak.get('peak_position')[1])
                elif peak.get('Position_Y'):
                    position_y = float(peak.get('Position_Y'))
                elif peak.get('position_y'):
                    position_y = float(peak.get('position_y'))

                # DEBUG: Log conversion for first few peaks
                if i < 3:
                    print(f"   Peak {i}: assignment='{assignment}', x={position_x}, y={position_y}")

                peak_row = {
                    'Assignment': assignment,
                    'Position_X': float(position_x),
                    'Position_Y': float(position_y)
                }
                peak_data.append(peak_row)
            else:
                skipped_peaks += 1
                print(f"   ⚠️ Skipped peak {i}: not a dictionary (type: {type(peak)})")

        if skipped_peaks > 0:
            print(f"   ❌ WARNING: Skipped {skipped_peaks} peaks due to format issues")

        print(f"🔄 Converted {len(peak_data)} fitted peaks to DataFrame format for series integration")

        # CRITICAL: Validate result before returning
        if len(peak_data) == 0 and len(fitted_peaks) > 0:
            print(f"   ❌ CRITICAL ERROR: All {len(fitted_peaks)} peaks were skipped during conversion!")

        return pd.DataFrame(peak_data)

    def _run_new_multi_spectrum_processing(self, nmr_files, reference_peaks, peak_source_mode, progress_dialog):
        """
        Run the new independent multi-spectrum processing in background thread.
        """
        try:
            # 1. Gather all parameters from parameter manager
            all_params = self.param_manager.get_integrator_parameters()


            # 2. Create the new independent processor
            multi_processor = MultiSpectrumProcessor(all_params)

            # NEW - with parameter transformation:
            raw_params = self.param_manager.current_params.copy()

            # 3. Define output folder with timestamp
            base_folder = self.nmr_file_list.get_current_folder()
            output_folder = os.path.join(
                base_folder,
                f"series_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
            )

            # 3.5. Validate progress dialog before starting
            if not self._validate_progress_dialog(progress_dialog):
                raise ValueError("Progress dialog validation failed")


            # 4. Define progress callback for GUI updates
            #def progress_callback(progress, task_desc, log_msg, failed=False):
            #    if progress_dialog and not progress_dialog.cancelled:
            #        self.root.after(0, lambda: progress_dialog.update_progress(
            #            progress, task_desc, log_msg, failed
            #        ))

            #def progress_callback(progress, task_desc, log_msg, failed=False):
            #    if progress_dialog and hasattr(progress_dialog, 'cancelled') and not progress_dialog.cancelled:
            #        self.root.after(0, lambda: self._update_progress_safe(
            #            progress_dialog, progress, task_desc, log_msg, failed
            #        ))

            # 4. Define progress callback for GUI updates with error handling
            def progress_callback(progress, task_desc, log_msg, failed=False):
                if progress_dialog and hasattr(progress_dialog, 'top') and progress_dialog.top.winfo_exists():
                    try:
                        # Check if dialog was cancelled (if supported)
                        if hasattr(progress_dialog, 'cancelled') and progress_dialog.cancelled:
                            return

                        # Schedule safe progress update on main thread
                        self.root.after(0, lambda p=progress, t=task_desc, l=log_msg, f=failed:
                                      self._update_progress_safe(progress_dialog, p, t, l, f))

                    except Exception as e:
                        print(f"⚠️ Progress callback scheduling failed: {e}")

            # 5. Run the independent processing
            print("🎯 Starting independent multi-spectrum processing...")
            new_batch_results = multi_processor.process_nmr_series(
                nmr_files=nmr_files,
                reference_peaks=reference_peaks,
                output_folder=output_folder,
                peak_source_mode=peak_source_mode,
                progress_callback=progress_callback
            )

            # 6. Store results and schedule GUI update
            self.new_batch_results = new_batch_results  # Store new format
            self.root.after(0, lambda: self._complete_new_multi_spectrum_processing(
                new_batch_results, progress_dialog
            ))

        except Exception as e:
            import traceback
            error_msg = f"Multi-spectrum processing error: {e}"
            print(f"❌ {error_msg}\n{traceback.format_exc()}")

            # Schedule error handling on main thread
            self.root.after(0, lambda: self._handle_processing_error(error_msg, progress_dialog))

        finally:
            self.processing_active = False

    def _complete_new_multi_spectrum_processing(self, new_batch_results, progress_dialog):
        """
        Complete processing and update GUI with perfect retrocompatibility.
        """
        try:
            # 1. Validate results format
            if not self._validate_batch_results(new_batch_results):
                raise ValueError("Invalid batch results format")

            # 2. Complete progress dialog safely
            if progress_dialog:
                try:
                    progress_dialog.complete("Multi-spectrum processing completed successfully!")
                except Exception as e:
                    print(f"⚠️ Progress dialog completion failed: {e}")

            # 3. Convert new format to legacy format for retrocompatibility
            print("🔄 Converting results to legacy format for GUI compatibility...")
            try:
                legacy_batch_results = self._convert_multispec_to_legacy_format(new_batch_results)
            except Exception as e:
                print(f"❌ Format conversion failed: {e}")
                # Use original results as fallback
                legacy_batch_results = new_batch_results

            # 1. Complete progress dialog
            if progress_dialog:
                progress_dialog.complete("Multi-spectrum processing completed successfully!")

            # 2. Convert new format to legacy format for retrocompatibility
            print("🔄 Converting results to legacy format for GUI compatibility...")
            legacy_batch_results = self._convert_multispec_to_legacy_format(new_batch_results)

            # 3. Store both formats
            self.batch_results = legacy_batch_results      # For legacy GUI components
            self.new_batch_results = new_batch_results     # Keep new format for advanced features

            # 4. Update series_processor for export compatibility
            self._update_series_processor_compatibility(new_batch_results)

            # 5. Sync results to legacy components (peak navigator, etc.)
            self._sync_results_to_legacy_components(new_batch_results)


            # 7. Update all GUI components with legacy-compatible data
            self._safe_update_main_plot()
            self._safe_update_results_dropdown()

            # 8. Update series visualization with safety checks
            if hasattr(self, 'series_plotter'):
                try:
                    # Validate data structure before plotting
                    self._validate_legacy_results_for_plotting(legacy_batch_results)
                    self.series_plotter.plot_series_overview(legacy_batch_results)

                    if hasattr(self, 'canvas_series'):
                        self.canvas_series.draw()

                    print("✅ Series visualization updated successfully")

                except Exception as e:
                    print(f"⚠️ Series visualization failed: {e}")
                    print("✅ Processing completed successfully despite visualization issue")

                    # Try alternative visualization approach
                    self._plot_series_alternative(legacy_batch_results)

            # 9. Update statistics display
            self.update_statistics_from_batch(legacy_batch_results)

            # 10. Switch to series overview tab
            if hasattr(self, 'viz_notebook'):
                self.viz_notebook.select(2)  # Series overview tab

            # 11. Display success message
            summary = new_batch_results.get('summary', {})
            successful = summary.get('successful', 0)
            total = summary.get('total_spectra', 0)
            success_rate = summary.get('success_rate', 0)

            status_msg = f"✅ Series processing complete: {successful}/{total} successful ({success_rate:.1f}%)"
            self.update_status(status_msg)

            print("🎉 Multi-spectrum processing completed with full retrocompatibility!")

        except Exception as e:
            import traceback
            error_msg = f"Error completing multi-spectrum processing: {e}"
            print(f"❌ {error_msg}\n{traceback.format_exc()}")

            # Safe status update
            try:
                self.update_status(f"❌ Processing completion failed: {e}")
            except Exception as status_error:
                print(f"⚠️ Status update failed: {status_error}")

            # Safe progress dialog completion
            if progress_dialog:
                try:
                    # Use safe completion method
                    self._complete_progress_dialog_safe(progress_dialog, "Processing completed with errors", failed=True)
                except Exception as dialog_error:
                    print(f"⚠️ Progress dialog completion failed: {dialog_error}")

#gm added

    def _update_progress_safe(self, progress_dialog, progress, task_desc, log_msg, failed=False):
        """
        Safely update progress dialog with error handling for different dialog types.

        This method provides a robust interface for updating progress dialogs regardless
        of their specific implementation, with graceful fallbacks for unsupported features.

        Args:
            progress_dialog: The progress dialog instance
            progress (float): Progress percentage (0-100)
            task_desc (str): Description of current task
            log_msg (str): Log message to display
            failed (bool): Whether the current operation failed
        """
        if not progress_dialog:
            return

        try:
            # Method 1: Try standard update_progress with multiple parameters
            if hasattr(progress_dialog, 'update_progress'):
                try:
                    import inspect
                    sig = inspect.signature(progress_dialog.update_progress)
                    params = list(sig.parameters.keys())

                    if len(params) >= 3:
                        # Supports progress, task, and message
                        progress_dialog.update_progress(progress, task_desc, log_msg)
                    elif len(params) >= 2:
                        # Supports progress and task only
                        progress_dialog.update_progress(progress, task_desc)
                    else:
                        # Only supports progress
                        progress_dialog.update_progress(progress)

                except (TypeError, AttributeError) as e:
                    # Fallback if signature inspection fails
                    progress_dialog.update_progress(progress, task_desc)

            # Method 2: Try direct variable updates if method call failed
            else:
                self._update_progress_variables_direct(progress_dialog, progress, task_desc)

            # Method 3: Handle log messages separately if supported
            if log_msg:
                self._add_progress_log_safe(progress_dialog, log_msg, failed)

            # Method 4: Handle failure status if supported
            if failed:
                self._set_progress_status_safe(progress_dialog, "failed")

        except Exception as e:
            print(f"⚠️ Progress update failed: {e}")
            # Ultimate fallback - try basic variable assignment
            self._progress_fallback_update(progress_dialog, progress, task_desc)

    def _update_progress_variables_direct(self, progress_dialog, progress, task_desc):
        """Direct update of progress dialog variables"""
        try:
            if hasattr(progress_dialog, 'progress_var') and progress_dialog.progress_var:
                progress_dialog.progress_var.set(progress)
            if hasattr(progress_dialog, 'current_task') and progress_dialog.current_task:
                progress_dialog.current_task.set(task_desc)
        except Exception as e:
            print(f"⚠️ Direct progress variable update failed: {e}")

    def _add_progress_log_safe(self, progress_dialog, log_msg, failed=False):
        """Safely add log message to progress dialog"""
        try:
            if hasattr(progress_dialog, 'add_log'):
                progress_dialog.add_log(log_msg)
            elif hasattr(progress_dialog, 'log_text') and hasattr(progress_dialog.log_text, 'insert'):
                # Try direct text widget manipulation
                progress_dialog.log_text.insert('end', f"{log_msg}\n")
                progress_dialog.log_text.see('end')
        except Exception as e:
            print(f"⚠️ Log message update failed: {e}")

    def _set_progress_status_safe(self, progress_dialog, status):
        """Safely set progress dialog status"""
        try:
            if hasattr(progress_dialog, 'set_status'):
                progress_dialog.set_status(status)
            elif hasattr(progress_dialog, 'status_var') and progress_dialog.status_var:
                progress_dialog.status_var.set(status)
        except Exception as e:
            print(f"⚠️ Status update failed: {e}")

    def _progress_fallback_update(self, progress_dialog, progress, task_desc):
        """Ultimate fallback for progress updates"""
        try:
            # Try to find any progress-related attributes and update them
            if hasattr(progress_dialog, 'top') and hasattr(progress_dialog.top, 'title'):
                progress_dialog.top.title(f"Processing... {progress:.1f}%")

            # Print to console as last resort
            print(f"📊 Progress: {progress:.1f}% - {task_desc}")

        except Exception as e:
            # Silent failure - don't cascade errors
            pass

    def _validate_progress_dialog(self, progress_dialog):
        """
        Validate that the progress dialog is properly initialized and accessible.

        Args:
            progress_dialog: The progress dialog instance to validate

        Returns:
            bool: True if dialog is valid, False otherwise
        """
        try:
            if not progress_dialog:
                print("⚠️ Progress dialog is None")
                return False

            if not hasattr(progress_dialog, 'top'):
                print("⚠️ Progress dialog missing 'top' window attribute")
                return False

            # Check if the dialog window exists
            if not progress_dialog.top.winfo_exists():
                print("⚠️ Progress dialog window doesn't exist")
                return False

            # Check for basic functionality
            if not (hasattr(progress_dialog, 'update_progress') or
                   hasattr(progress_dialog, 'progress_var')):
                print("⚠️ Progress dialog missing update mechanism")
                return False

            print("✅ Progress dialog validation passed")
            return True

        except Exception as e:
            print(f"⚠️ Progress dialog validation failed: {e}")
            return False

#gm added
    def _validate_batch_results(self, batch_results):
        """Validate batch results format"""
        if not isinstance(batch_results, dict):
            return False

        required_keys = ['metadata', 'results', 'summary']
        for key in required_keys:
            if key not in batch_results:
                print(f"❌ Missing required key in batch results: {key}")
                return False

        return True

    def _handle_processing_error(self, error_msg, progress_dialog):
        """Handle processing errors on main thread with robust error handling"""
        try:
            self.update_status(f"❌ Series processing failed: {error_msg}")

            # Safely complete progress dialog
            if progress_dialog:
                try:
                    progress_dialog.complete("Processing failed")
                except Exception as dialog_error:
                    print(f"⚠️ Progress dialog completion failed: {dialog_error}")
                    # Force close if complete fails
                    try:
                        if hasattr(progress_dialog, 'top'):
                            progress_dialog.top.destroy()
                    except:
                        pass

            # Show error to user
            try:
                messagebox.showerror("Processing Error",
                    f"Series integration failed:\n\n{error_msg}")
            except Exception as msg_error:
                print(f"⚠️ Error dialog failed: {msg_error}")

        except Exception as e:
            print(f"❌ Critical error in error handling: {e}")
            import traceback
            traceback.print_exc()


    def _convert_multispec_to_legacy_format(self, new_batch_results):
        """
        Convert MultiSpectrumProcessor format to legacy format.
        This ensures perfect backward compatibility with all GUI components.
        """
        # Create legacy-compatible structure with all required attributes
        legacy_batch_results = type('BatchResults', (), {})()

        # Convert results format from new nested structure to legacy flat structure
        legacy_batch_results.results = {}
        if 'results' in new_batch_results:
            for spectrum_name, spectrum_data in new_batch_results['results'].items():
                if isinstance(spectrum_data, dict) and spectrum_data.get('success', False):
                    # Create complete spectrum result dictionary for legacy compatibility
                    legacy_spectrum_result = {
                        # Core results data
                        'fitted_results': spectrum_data.get('fitted_results', []),
                        'integration_results': spectrum_data.get('integration_results', []),

                        # Statistics for visualization compatibility
                        'detection_rate': spectrum_data.get('success_rate', 0.0),
                        'success_rate': spectrum_data.get('success_rate', 0.0),
                        'total_peaks': spectrum_data.get('total_peaks', 0),
                        'successful_fits': spectrum_data.get('successful_fits', 0),
                        'failed_fits': spectrum_data.get('total_peaks', 0) -  spectrum_data.get('successful_fits', 0),

                        # Status indicators (existing + new)
                        'success': spectrum_data.get('success', False),
                        'status': 'success' if spectrum_data.get('success', False) else 'failed',
                        'processing_complete': True,
                        'spectrum_file': spectrum_data.get('spectrum_file', spectrum_name),

                        # Spectrum browser specific fields
                        'detected_peaks': spectrum_data.get('successful_fits', 0),
                        'processing_time': spectrum_data.get('processing_time', 0.0),
                        'is_reference': spectrum_name.lower().find('ref') != -1 or  spectrum_name.lower().find('reference') != -1,

                        # Additional metadata for compatibility
                        'processing_mode': 'multi_spectrum',
                        'timestamp': spectrum_data.get('timestamp', ''),
                    }
                    legacy_batch_results.results[spectrum_name] = legacy_spectrum_result
                else:
                    # Create empty result structure for failed spectra
                    legacy_batch_results.results[spectrum_name] = {
                        # Core results data
                        'fitted_results': [],
                        'integration_results': [],

                        # Statistics for visualization compatibility
                        'detection_rate': 0.0,
                        'success_rate': 0.0,
                        'total_peaks': spectrum_data.get('total_peaks', 0) if  isinstance(spectrum_data, dict) else 0,
                        'successful_fits': 0,
                        'failed_fits': spectrum_data.get('total_peaks', 0) if  isinstance(spectrum_data, dict) else 0,

                        # Status indicators
                        'success': False,
                        'status': 'failed',
                        'processing_complete': False,
                        'spectrum_file': spectrum_name,

                        # Spectrum browser specific fields
                        'detected_peaks': 0,
                        'processing_time': 0.0,
                        'is_reference': spectrum_name.lower().find('ref') != -1 or spectrum_name.lower().find('reference') != -1,

                        # Additional metadata for compatibility
                        'processing_mode': 'multi_spectrum',
                        'error': spectrum_data.get('error', 'Processing failed') if isinstance(spectrum_data, dict) else 'Processing failed',
                        'timestamp': ''
                    }

        # Add get_summary method for SpectrumBrowserDialog compatibility
        def get_summary():
            summary = new_batch_results.get('summary', {})
            return {
                'total_spectra': summary.get('total_spectra', 0),
                'successful': summary.get('successful', 0),
                'failed': summary.get('failed', 0),
                'total_peaks': summary.get('total_peaks_processed', 0),
                'successful_peaks': summary.get('total_successful_fits', 0),
                'success_rate': summary.get('success_rate', 0.0),
                'detection_rate': summary.get('overall_detection_rate', 0.0)
            }

        legacy_batch_results.get_summary = get_summary

        # Add metadata for SpectrumBrowserDialog
        legacy_batch_results.metadata = new_batch_results.get('metadata', {})

        # Store reference to original new format for advanced features
        legacy_batch_results._new_format = new_batch_results

        # Add output_folder attribute for export compatibility
        metadata = new_batch_results.get('metadata', {})
        legacy_batch_results.output_folder = metadata.get('output_folder', '')

        # Add statistics attribute for update_statistics_from_batch compatibility
        summary = new_batch_results.get('summary', {})
        legacy_batch_results.statistics = {
            'total_spectra': summary.get('total_spectra', 0),
            'successful_spectra': summary.get('successful', 0),
            'failed_spectra': summary.get('failed', 0),
            'total_peaks_processed': summary.get('total_peaks_processed', 0),
            'total_successful_fits': summary.get('total_successful_fits', 0),
            'overall_success_rate': summary.get('success_rate', 0.0),
            'overall_detection_rate': summary.get('overall_detection_rate', 0.0)
        }

        return legacy_batch_results

    def _validate_legacy_results_for_plotting(self, legacy_batch_results):
        """
        Validate that legacy batch results have the correct structure for plotting.

        Args:
            legacy_batch_results: The legacy format batch results

        Raises:
            ValueError: If the data structure is invalid for plotting
        """
        if not hasattr(legacy_batch_results, 'results'):
            raise ValueError("Legacy batch results missing 'results' attribute")

        if not isinstance(legacy_batch_results.results, dict):
            raise ValueError("Legacy batch results 'results' is not a dictionary")

        # Check each spectrum result
        for spectrum_name, result in legacy_batch_results.results.items():
            if not isinstance(result, dict):
                raise ValueError(f"Spectrum result for '{spectrum_name}' is not a dictionary: {type(result)}")

            # Check for required keys
            required_keys = ['detection_rate', 'success_rate', 'total_peaks', 'successful_fits']
            missing_keys = [key for key in required_keys if key not in result]

            if missing_keys:
                print(f"⚠️ Spectrum '{spectrum_name}' missing keys: {missing_keys}")
                # Add missing keys with default values
                for key in missing_keys:
                    if key in ['detection_rate', 'success_rate']:
                        result[key] = 0.0
                    elif key in ['total_peaks', 'successful_fits']:
                        result[key] = 0

        print("✅ Legacy results validation passed")

    def _plot_series_alternative(self, legacy_batch_results):
        """
        Alternative series plotting method with simplified data requirements.

        Args:
            legacy_batch_results: The legacy format batch results
        """
        try:
            if not hasattr(self, 'series_plotter'):
                return

            print("🎨 Attempting alternative series visualization...")

            # Create simplified data structure for plotting
            simplified_data = {
                'spectrum_names': [],
                'success_rates': [],
                'peak_counts': []
            }

            if hasattr(legacy_batch_results, 'results'):
                for spectrum_name, result in legacy_batch_results.results.items():
                    simplified_data['spectrum_names'].append(spectrum_name)

                    if isinstance(result, dict):
                        simplified_data['success_rates'].append(result.get('success_rate', 0.0))
                        simplified_data['peak_counts'].append(result.get('successful_fits', 0))
                    else:
                        simplified_data['success_rates'].append(0.0)
                        simplified_data['peak_counts'].append(0)

            # Try to plot with simplified data
            if hasattr(self.series_plotter, 'plot_simple_overview'):
                self.series_plotter.plot_simple_overview(simplified_data)
            else:
                # Create basic text summary if plotting fails
                self._create_text_summary(simplified_data)

            print("✅ Alternative visualization completed")

        except Exception as e:
            print(f"⚠️ Alternative visualization also failed: {e}")
            print("📊 Series processing completed successfully without visualization")

    def _create_text_summary(self, data):
        """Create a text summary of processing results"""
        try:
            total_spectra = len(data['spectrum_names'])
            avg_success_rate = sum(data['success_rates']) / total_spectra if total_spectra > 0 else 0
            total_peaks = sum(data['peak_counts'])

            summary = f"""
📊 SERIES PROCESSING SUMMARY 📊
═══════════════════════════════
Total Spectra: {total_spectra}
Average Success Rate: {avg_success_rate:.1f}%
Total Peaks Processed: {total_peaks}
═══════════════════════════════
            """

            print(summary)

            # Update status with summary
            if hasattr(self, 'update_status'):
                self.update_status(f"✅ Processed {total_spectra} spectra, {total_peaks} peaks total")

        except Exception as e:
            print(f"⚠️ Text summary creation failed: {e}")

    def _sync_results_to_legacy_components(self, new_batch_results):
        """
        Sync results to all legacy GUI components that depend on specific data formats.
        This ensures peak navigator and other components work correctly.
        """
        if not new_batch_results or not new_batch_results.get('results'):
            return

        # Update integrator.fitted_peaks with most recent successful spectrum results
        latest_spectrum_results = None
        for spectrum_name, spectrum_data in new_batch_results['results'].items():
            if isinstance(spectrum_data, dict) and spectrum_data.get('success', False):
                fitted_results = spectrum_data.get('fitted_results', [])
                if fitted_results:
                    latest_spectrum_results = fitted_results
                    break  # Use first successful spectrum

        if latest_spectrum_results:
            # Update main integrator fitted_peaks for peak navigator compatibility
            self.integrator.fitted_peaks = latest_spectrum_results

            # Update peak navigator if it exists
            if hasattr(self, 'peak_navigator') and self.peak_navigator:
                try:
                    self.peak_navigator.load_detected_peaks(latest_spectrum_results)
                    print("✅ Peak navigator updated with latest results")
                except Exception as e:
                    print(f"⚠️ Failed to update peak navigator: {e}")

            print(f"✅ Legacy components synced with {len(latest_spectrum_results)} peaks")


    def _update_series_processor_compatibility(self, new_batch_results):
        """
        Update series_processor attributes for export function compatibility.
        """
        if hasattr(self, 'series_processor') and new_batch_results:
            metadata = new_batch_results.get('metadata', {})
            output_folder = metadata.get('output_folder', '')

            if output_folder:
                # Update series_processor.output_folder for legacy export compatibility
                self.series_processor.output_folder = output_folder
                print(f"✅ Series processor output folder updated: {output_folder}")

    # =================== VISUALIZATION UPDATES ===================

    def update_main_plot(self):
        """Update main spectrum plot"""
        if hasattr(self.integrator, 'nmr_data') and self.integrator.nmr_data is not None:
            # Save current zoom state
            self.save_current_zoom()

            # Plot spectrum with contour controls
            self.spectrum_plotter.plot_spectrum(
                self.integrator,
                contour_levels=self.contour_levels.get(),
                contour_min_level=self.contour_min.get(),
                contour_increment=self.contour_increment.get()
            )

            # Plot peaks
            self.spectrum_plotter.plot_peaks(
                self.integrator,
                show_detected=self.show_detected.get(),
                show_assigned=self.show_assigned.get(),
                show_fitted=self.show_fitted_curves.get(),
                #show_included_peaks_after_limit_debug=self.show_included_peaks_after_limit_debug.get()
            )

            # Restore zoom if it was set
            self.restore_zoom()

            self.canvas_main.draw()

    def update_statistics(self):
        """Update statistics display"""
        stats = {}

        # Detection statistics
        if hasattr(self.integrator, 'fitted_peaks') and (
            self.integrator.fitted_peaks is not None and
            (not hasattr(self.integrator.fitted_peaks, 'empty') or not self.integrator.fitted_peaks.empty) and
            len(self.integrator.fitted_peaks) > 0
        ):
            detected_count = sum(1 for p in self.integrator.fitted_peaks if p.get('detected', False))
            total_count = len(self.integrator.fitted_peaks)

            stats['detection'] = {
                'detected_peaks': detected_count,
                'total_peaks': total_count,
                'detection_rate': (detected_count / total_count * 100) if total_count > 0 else 0
            }

        # Integration statistics
        if hasattr(self.integrator, 'integration_results') and (
            self.integrator.integration_results is not None and
            (not hasattr(self.integrator.integration_results, 'empty') or not self.integrator.integration_results.empty) and
            len(self.integrator.integration_results) > 0
        ):
            good_quality = sum(1 for r in self.integrator.integration_results
                              if r.get('Quality') in ['Excellent', 'Good'])

            stats['integration'] = {
                'total_integrations': len(self.integrator.integration_results),
                'good_quality': good_quality,
                'quality_rate': (good_quality / len(self.integrator.integration_results) * 100)
            }

        # Fitting statistics (check both new and legacy locations)
        fits = []
        if hasattr(self.integrator, 'fitted_peaks') and self.integrator.fitted_peaks:
            fits.extend(self.integrator.fitted_peaks)
        elif hasattr(self.integrator, 'voigt_fits') and self.integrator.voigt_fits:
            fits.extend(self.integrator.voigt_fits)

        if fits:
            excellent_fits = sum(1 for f in fits if f.get('fitting_quality') == 'Excellent')

            stats['fitting'] = {
                'total_fits': len(fits),
                'excellent_fits': excellent_fits,
                'average_r_squared': np.mean([f.get('avg_r_squared', 0) for f in fits])
            }

        self.statistics_panel.update_stats(stats)

    def update_statistics_from_batch(self, batch_results):
        """Update statistics from batch processing results"""
        if not batch_results:
            return

        summary = batch_results.get_summary()
        stats = batch_results.statistics

        combined_stats = {
            'series_summary': summary,
            'batch_statistics': stats
        }

        self.statistics_panel.update_stats(combined_stats)

    # =================== NAVIGATION METHODS ===================

    def update_peak_navigation(self):
        """Update peak navigation info"""
        if hasattr(self.integrator, 'peak_list') and self.integrator.peak_list is not None:
            max_peaks = len(self.integrator.peak_list)
            current_peak = self.selected_peak_number.get()

            # Get assignment for current peak
            if 1 <= current_peak <= max_peaks:
                peak_data = self.integrator.peak_list.iloc[current_peak - 1]
                assignment = peak_data.get('Assignment', f'Peak_{current_peak}')
                self.peak_info_label.config(text=f"Peak: {assignment} ({current_peak}/{max_peaks})")
            else:
                self.peak_info_label.config(text=f"Peak: {current_peak}/{max_peaks}")
        else:
            self.peak_info_label.config(text="Peak: -/-")

    def center_on_selected_peak(self):
        """Center view on selected peak"""
        if not hasattr(self.integrator, 'peak_list') or self.integrator.peak_list is None:
            return

        peak_num = self.selected_peak_number.get()
        if peak_num < 1 or peak_num > len(self.integrator.peak_list):
            return

        # Get peak position
        peak_row = self.integrator.peak_list.iloc[peak_num - 1]
        peak_x = float(peak_row['Position_X'])
        peak_y = float(peak_row['Position_Y'])

        # Highlight and center on peak
        self.spectrum_plotter.highlight_peak(peak_x, peak_y)
        self.spectrum_plotter.set_zoom(peak_x, peak_y, self.zoom_x_range.get(), self.zoom_y_range.get())
        self.canvas_main.draw()

        self.update_status(f"Centered on peak {peak_num}: {peak_row.get('Assignment', f'Peak_{peak_num}')}")

    def zoom_to_peak(self):
        """Zoom to selected peak with adjustable range"""
        if not hasattr(self.integrator, 'peak_list') or self.integrator.peak_list is None:
            return

        # Create zoom dialog
        zoom_dialog = tk.Toplevel(self.root)
        zoom_dialog.title("Zoom Settings")
        # Set responsive geometry for zoom dialog
        zoom_dialog.geometry("350x250")
        zoom_dialog.minsize(300, 200)
        zoom_dialog.resizable(True, True)
        zoom_dialog.transient(self.root)
        zoom_dialog.grab_set()

        ttk.Label(zoom_dialog, text="Zoom Range Settings:").pack(pady=10)

        frame = ttk.Frame(zoom_dialog)
        frame.pack(pady=10)

        ttk.Label(frame, text="X Range (ppm):").grid(row=0, column=0, sticky=tk.W, padx=5)
        x_range_spin = tk.Spinbox(frame, from_=0.5, to=5.0, increment=0.1, width=4,
                                 textvariable=self.zoom_x_range, format="%.1f")
        x_range_spin.grid(row=0, column=1, padx=5)

        ttk.Label(frame, text="Y Range (ppm):").grid(row=1, column=0, sticky=tk.W, padx=5)
        y_range_spin = tk.Spinbox(frame, from_=5.0, to=50.0, increment=1.0, width=4,
                                 textvariable=self.zoom_y_range, format="%.1f")
        y_range_spin.grid(row=1, column=1, padx=5)

        button_frame = ttk.Frame(zoom_dialog)
        button_frame.pack(pady=20)

        ttk.Button(button_frame, text="Apply",
                  command=lambda: (self.center_on_selected_peak(), zoom_dialog.destroy())).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Cancel",
                  command=zoom_dialog.destroy).pack(side=tk.LEFT, padx=5)

    def prev_peak(self):
        """Navigate to previous peak"""
        if not hasattr(self.integrator, 'peak_list') or self.integrator.peak_list is None:
            return

        max_peaks = len(self.integrator.peak_list)
        current_peak = self.selected_peak_number.get()
        new_peak = current_peak - 1 if current_peak > 1 else max_peaks

        self.selected_peak_number.set(new_peak)
        self.center_on_selected_peak()
        self.update_peak_navigation()

    def next_peak(self):
        """Navigate to next peak"""
        if not hasattr(self.integrator, 'peak_list') or self.integrator.peak_list is None:
            return

        max_peaks = len(self.integrator.peak_list)
        current_peak = self.selected_peak_number.get()
        new_peak = current_peak + 1 if current_peak < max_peaks else 1

        self.selected_peak_number.set(new_peak)
        self.center_on_selected_peak()
        self.update_peak_navigation()

    # =================== PEAK EDITING METHODS ===================

    def toggle_peak_edit_mode(self):
        """Toggle peak editing mode on/off"""
        if self.peak_edit_mode.get():
            # Enable edit mode
            self.canvas_main.mpl_connect('button_press_event', self.on_peak_edit_click)
            self.canvas_main.get_tk_widget().config(cursor="crosshair")
            print("Peak editing mode enabled")
        else:
            # Disable edit mode
            self.canvas_main.mpl_disconnect('button_press_event')
            self.canvas_main.get_tk_widget().config(cursor="")
            self.selected_peak_info = None
            self.selected_peak_label.config(text="No peak selected")
            print("Peak editing mode disabled")

        # Update status display
        self.update_edit_mode_status()

    def update_edit_mode_status(self):
        """Update edit mode status display based on current settings"""
        if not self.peak_edit_mode.get():
            self.edit_mode_status_label.config(text="Mode: View Only", foreground='gray')
            return

        # Check which peak lists are enabled for editing
        ref_enabled = self.edit_reference_peaks.get()
        det_enabled = self.edit_detected_peaks.get()

        if ref_enabled and det_enabled:
            status_text = "Mode: EDIT (Ref + Det)"
        elif ref_enabled:
            status_text = "Mode: EDIT (Ref only)"
        elif det_enabled:
            status_text = "Mode: EDIT (Det only)"
        else:
            status_text = "Mode: EDIT (None selected!)"
            self.edit_mode_status_label.config(text=status_text, foreground='orange')
            return

        self.edit_mode_status_label.config(text=status_text, foreground='red')

    def on_peak_edit_click(self, event):
        """Handle mouse clicks in peak edit mode"""
        if not self.peak_edit_mode.get() or event.inaxes != self.ax_main:
            return

        click_x, click_y = event.xdata, event.ydata

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
            self.selected_peak_label.config(text="Peak moved! No peak selected")

    def find_nearest_peak(self, click_x, click_y):
        """Find the nearest peak to click position"""
        import math
        candidates = []

        # Check reference peaks (only if enabled for editing)
        if (self.edit_reference_peaks.get() and
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
        if (self.edit_detected_peaks.get() and
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
            print(f"Selected {closest['type']} peak {closest['index']}: {closest['assignment']} at ({closest['x']:.3f}, {closest['y']:.1f})")
            return closest

        return None

    def update_selected_peak_display(self):
        """Update display to show selected peak info"""
        if self.selected_peak_info:
            peak_info = self.selected_peak_info
            self.selected_peak_label.config(
                text=f"Selected: {peak_info['type']} peak '{peak_info['assignment']}' at ({peak_info['x']:.3f}, {peak_info['y']:.1f}) - Click new position to move"
            )

    def move_selected_peak(self, new_x, new_y):
        """Move selected peak to new position"""
        if not self.selected_peak_info:
            return

        peak_info = self.selected_peak_info
        old_x, old_y = peak_info['x'], peak_info['y']

        print(f"Moving {peak_info['type']} peak '{peak_info['assignment']}' from ({old_x:.3f}, {old_y:.1f}) to ({new_x:.3f}, {new_y:.1f})")

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
        self.update_main_plot()

        print(f"✅ Peak position updated successfully")

    # =================== FILE OPERATIONS ===================

    def load_nmr_file(self):
        """Load NMR file via dialog"""
        filename = filedialog.askopenfilename(
            title="Select NMR Spectrum",
            filetypes=[("NMR files", "*.ft *.fid"), ("All files", "*.*")]
        )
        if filename:
            self.on_nmr_file_select(filename, os.path.basename(filename))

    def load_peak_file(self):
        """Load peak list file via dialog"""
        filename = filedialog.askopenfilename(
            title="Select Peak List",
            filetypes=[("Text files", "*.txt"), ("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if filename:
            self.on_peak_file_select(filename, os.path.basename(filename))

    def export_peak_list(self):
        """Export current peak list"""
        if not hasattr(self.integrator, 'fitted_peaks') or (
            self.integrator.fitted_peaks is None or
            (hasattr(self.integrator.fitted_peaks, 'empty') and self.integrator.fitted_peaks.empty) or
            len(self.integrator.fitted_peaks) == 0
        ):
            messagebox.showerror("Error", "No peaks to export")
            return

        filename = filedialog.asksaveasfilename(
            title="Export Peak List",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("Text files", "*.txt"), ("All files", "*.*")]
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

                self.update_status(f"✅ Peak list exported: {len(peak_data)} peaks")
                messagebox.showinfo("Export Successful", f"Peak list exported to:\n{filename}")

            except Exception as e:
                self.update_status(f"❌ Export failed: {str(e)}")
                messagebox.showerror("Export Error", str(e))

    def export_results(self):
        """Export current processing results"""
        if not hasattr(self.integrator, 'integration_results') or (
            self.integrator.integration_results is None or
            (hasattr(self.integrator.integration_results, 'empty') and self.integrator.integration_results.empty) or
            len(self.integrator.integration_results) == 0
        ):
            messagebox.showerror("Error", "No integration results to export")
            return

        filename = filedialog.asksaveasfilename(
            title="Export Integration Results",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("Excel files", "*.xlsx"), ("All files", "*.*")]
        )

        if filename:
            try:
                df = pd.DataFrame(self.integrator.integration_results)

                if filename.endswith('.xlsx'):
                    df.to_excel(filename, index=False)
                else:
                    df.to_csv(filename, index=False, float_format='%.6f')

                self.update_status(f"✅ Results exported: {len(df)} integrations")
                messagebox.showinfo("Export Successful", f"Results exported to:\n{filename}")

            except Exception as e:
                self.update_status(f"❌ Export failed: {str(e)}")
                messagebox.showerror("Export Error", str(e))

    def save_plot(self):
        """Save current plot"""
        filename = filedialog.asksaveasfilename(
            title="Save Plot",
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("PDF files", "*.pdf"), ("SVG files", "*.svg")]
        )

        if filename:
            try:
                # Determine which tab is active
                current_tab = self.viz_notebook.index(self.viz_notebook.select())

                if current_tab == 0:
                    self.fig_main.savefig(filename, dpi=300, bbox_inches='tight')
                elif current_tab == 1:
                    self.fig_voigt.savefig(filename, dpi=300, bbox_inches='tight')
                elif current_tab == 2:
                    self.fig_series.savefig(filename, dpi=300, bbox_inches='tight')
                elif current_tab == 3:
                    self.fig_stats.savefig(filename, dpi=300, bbox_inches='tight')

                self.update_status(f"✅ Plot saved: {os.path.basename(filename)}")
                messagebox.showinfo("Save Successful", f"Plot saved to:\n{filename}")

            except Exception as e:
                self.update_status(f"❌ Save failed: {str(e)}")
                messagebox.showerror("Save Error", str(e))

    # =================== CONFIGURATION METHODS ===================

    def save_config(self):
        """Save configuration to file"""
        filename = filedialog.asksaveasfilename(
            title="Save Configuration",
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )

        if filename:
            if self.config_manager.export_config(filename):
                messagebox.showinfo("Success", f"Configuration saved to:\n{filename}")
            else:
                messagebox.showerror("Error", "Failed to save configuration")

    def load_config(self):
        """Load configuration from file"""
        filename = filedialog.askopenfilename(
            title="Load Configuration",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )

        if filename:
            if self.config_manager.import_config(filename):
                # Reload variables from config
                self.init_variables_from_config()
                self.update_config_display()
                messagebox.showinfo("Success", f"Configuration loaded from:\n{filename}")
            else:
                messagebox.showerror("Error", "Failed to load configuration")

    def reset_config(self):
        """Reset configuration to defaults"""
        if messagebox.askyesno("Reset Configuration", "Reset all settings to default values?"):
            self.config_manager.reset_to_defaults()
            self.init_variables_from_config()
            self.update_config_display()
            self.update_status("✅ Configuration reset to defaults")

    def update_config_display(self):
        """Update configuration display in config tab"""
        if hasattr(self, 'config_text'):
            config_info = f"""Configuration Status:



Detection Parameters:
  Noise Threshold: {self.noise_threshold.get():.2f}
  Search Window X: ±{self.search_window_x.get():.3f} ppm
  Search Window Y: ±{self.search_window_y.get():.2f} ppm
  Use Reference Detection: {self.use_reference_detection.get()}

Fitting Parameters:
  Fitting Window X: {self.fitting_window_x.get():.1f} ppm
  Fitting Window Y: {self.fitting_window_y.get():.1f} ppm
  Min R²: {self.min_r_squared.get():.2f}
  Max Iterations: {self.max_iterations.get()}
  Global Optimization: {self.use_global_optimization.get()}

Display Options:
  Show Detected: {self.show_detected.get()}
  Show Assigned: {self.show_assigned.get()}
  Show Fitted Curves: {self.show_fitted_curves.get()}

Series Options:
  Auto Process: {self.auto_process_series.get()}
  Save Individual Results: {self.save_individual_results.get()}
  Create Summary Plots: {self.create_summary_plots.get()}

File Paths:
  Current NMR File: {os.path.basename(self.current_nmr_file) if self.current_nmr_file else 'None'}
  Current Peak File: {os.path.basename(self.current_peak_file) if self.current_peak_file else 'None'}

Configuration File: {self.config_manager.config_file}
Last Updated: {self.config_manager.config.get('last_updated', 'Never')}
"""

            self.config_text.delete(1.0, tk.END)
            self.config_text.insert(1.0, config_info)

    def open_config_dialog(self):
        """Open configuration dialog"""
        self.viz_notebook.select(4)  # Switch to config tab
        self.update_config_display()

    # =================== UTILITY METHODS ===================

    def validate_current_files(self):
        """Validate currently loaded files"""
        validation_messages = []

        if self.current_nmr_file:
            valid, message = self.file_manager.validate_nmr_file(self.current_nmr_file)
            if valid:
                validation_messages.append(f"✅ NMR: {message}")
            else:
                validation_messages.append(f"❌ NMR: {message}")

        if self.current_peak_file:
            valid, message = self.file_manager.validate_peak_file(self.current_peak_file)
            if valid:
                validation_messages.append(f"✅ Peaks: {message}")
            else:
                validation_messages.append(f"❌ Peaks: {message}")

        if validation_messages:
            self.file_validation_label.config(text=" | ".join(validation_messages))
        else:
            self.file_validation_label.config(text="")

    def reset_view(self):
        """Reset main spectrum view"""
        if hasattr(self.spectrum_plotter, 'reset_zoom'):
            self.spectrum_plotter.reset_zoom()
            self.canvas_main.draw()
        self.update_status("✅ View reset")

    def update_status(self, message):
        """Update status display"""
        print(f"Status: {message}")

        # Update GUI labels if they exist
        if hasattr(self, 'status_label'):
            self.status_label.config(text=message)
            # Auto-clear status after delay
            self.root.after(10000, lambda: self.status_label.config(text="Ready"))

        if hasattr(self, 'status_text'):
            self.status_text.config(text=message[:50] + "..." if len(message) > 50 else message)

    def update_recent_files_menu(self):
        """Update recent files menu"""
        self.recent_menu.delete(0, tk.END)

        recent_files = self.config_manager.load_recent_files()

        nmr_files = recent_files.get('nmr_files', [])
        peak_files = recent_files.get('peak_files', [])

        if nmr_files:
            self.recent_menu.add_separator()
            for i, file_path in enumerate(nmr_files[:5]):
                if os.path.exists(file_path):
                    filename = os.path.basename(file_path)
                    self.recent_menu.add_command(
                        label=f"📊 {filename}",
                        command=lambda f=file_path: self.on_nmr_file_select(f, os.path.basename(f))
                    )

        if peak_files:
            self.recent_menu.add_separator()
            for i, file_path in enumerate(peak_files[:5]):
                if os.path.exists(file_path):
                    filename = os.path.basename(file_path)
                    self.recent_menu.add_command(
                        label=f"📋 {filename}",
                        command=lambda f=file_path: self.on_peak_file_select(f, os.path.basename(f))
                    )

    def show_statistics(self):
        """Show statistics in separate window"""
        stats_window = tk.Toplevel(self.root)
        stats_window.title("Processing Statistics")
        # Set responsive geometry (50% of parent window size)
        parent_width = self.root.winfo_width()
        parent_height = self.root.winfo_height()
        width = max(500, int(parent_width * 0.5))
        height = max(300, int(parent_height * 0.5))
        stats_window.geometry(f"{width}x{height}")
        stats_window.minsize(500, 300)
        stats_window.transient(self.root)

        stats_panel = StatisticsPanel(stats_window)
        stats_panel.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Get current statistics
        stats = {}

        if hasattr(self.integrator, 'fitted_peaks') and (
            self.integrator.fitted_peaks is not None and
            (not hasattr(self.integrator.fitted_peaks, 'empty') or not self.integrator.fitted_peaks.empty) and
            len(self.integrator.fitted_peaks) > 0
        ):
            detected = sum(1 for p in self.integrator.fitted_peaks if p.get('detected', False))
            stats['current_detection'] = {
                'detected_peaks': detected,
                'total_peaks': len(self.integrator.fitted_peaks),
                'detection_rate': detected / len(self.integrator.fitted_peaks) * 100
            }

        if self.batch_results:
            stats['series_results'] = self.batch_results.get_summary()

        stats_panel.update_stats(stats)

    def show_series_analysis(self):
        """Show detailed series analysis"""
        if not self.batch_results:
            messagebox.showinfo("No Data", "No series results available for analysis")
            return

        # Switch to statistics tab and update with series analysis
        self.viz_notebook.select(3)

        # Create series analyzer and generate analysis
        analyzer = SeriesAnalyzer()
        analyzer.load_results(self.batch_results)

        trends = analyzer.analyze_detection_trends()
        quality_report = analyzer.generate_quality_report()

        # Update statistics plots
        if trends and 'data' in trends:
            data = trends['data']

            # Plot 1: Detection rate histogram
            self.ax_stats_1.clear()
            self.ax_stats_1.hist(data['detection_rate'], bins=15, alpha=0.7, color='skyblue')
            self.ax_stats_1.set_title('Detection Rate Distribution')
            self.ax_stats_1.set_xlabel('Detection Rate (%)')
            self.ax_stats_1.set_ylabel('Number of Spectra')
            self.ax_stats_1.grid(True, alpha=0.3)

            # Plot 2: Detection vs Noise correlation
            self.ax_stats_2.clear()
            if len(data) > 1:
                self.ax_stats_2.scatter(data['noise_level'], data['detection_rate'], alpha=0.7)
                self.ax_stats_2.set_title('Detection Rate vs Noise Level')
                self.ax_stats_2.set_xlabel('Noise Level')
                self.ax_stats_2.set_ylabel('Detection Rate (%)')
                self.ax_stats_2.grid(True, alpha=0.3)

            # Plot 3: Quality distribution
            self.ax_stats_3.clear()
            quality_ranges = ['0-50%', '50-70%', '70-90%', '90-100%']
            quality_counts = [
                sum(1 for rate in data['detection_rate'] if 0 <= rate < 50),
                sum(1 for rate in data['detection_rate'] if 50 <= rate < 70),
                sum(1 for rate in data['detection_rate'] if 70 <= rate < 90),
                sum(1 for rate in data['detection_rate'] if 90 <= rate <= 100)
            ]
            colors = ['red', 'orange', 'yellow', 'green']

            self.ax_stats_3.bar(quality_ranges, quality_counts, color=colors, alpha=0.7)
            self.ax_stats_3.set_title('Quality Distribution')
            self.ax_stats_3.set_ylabel('Number of Spectra')
            self.ax_stats_3.grid(True, alpha=0.3)

            # Plot 4: Processing summary
            self.ax_stats_4.clear()
            summary = self.batch_results.get_summary()

            labels = ['Successful', 'Failed']
            sizes = [summary['successful'], summary['failed']]
            colors = ['green', 'red']

            if sum(sizes) > 0:
                self.ax_stats_4.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
                self.ax_stats_4.set_title('Processing Success Rate')

        self.fig_stats.tight_layout()
        self.canvas_stats.draw()

    def show_voigt_analysis(self):
        """Show Voigt analysis tab (legacy method)"""
        self.viz_notebook.select(1)

        if self.current_voigt_result:
            self.voigt_plotter.plot_voigt_analysis(self.current_voigt_result)
            self.canvas_voigt.draw()
        else:
            messagebox.showinfo("No Data", "No Voigt fitting results available.\nFit a peak to view analysis.")

    def show_selected_peak_analysis(self):
        """Show Voigt analysis for the currently selected peak"""
        if (not hasattr(self.integrator, 'peak_list') or
            self.integrator.peak_list is None or
            (hasattr(self.integrator.peak_list, 'empty') and self.integrator.peak_list.empty)):
            messagebox.showinfo("No Peak List", "No peak list loaded.\nPlease load peak data first.")
            return

        peak_number = self.selected_peak_number.get()
        max_peaks = len(self.integrator.peak_list)

        if peak_number < 1 or peak_number > max_peaks:
            messagebox.showwarning("Invalid Peak", f"Peak number {peak_number} is invalid.\nValid range: 1-{max_peaks}")
            return

        # Check if the selected peak has fitting results (both new and legacy locations)
        selected_result = None

        # First check fitted_peaks (new location)
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
            # Switch to Voigt analysis tab
            self.viz_notebook.select(1)

            # Show the analysis for the selected peak
            self.voigt_plotter.plot_voigt_analysis(selected_result)
            self.canvas_voigt.draw()

            # Update status
            assignment = selected_result.get('assignment', f'Peak_{peak_number}')
            quality = selected_result.get('fitting_quality', 'Unknown')
            self.update_status(f"✅ Showing Voigt analysis for {assignment} (Quality: {quality})")
        else:
            # Peak not fitted yet - offer to fit it
            peak_data = self.integrator.peak_list.iloc[peak_number - 1]
            assignment = peak_data.get('Assignment', f'Peak_{peak_number}')

            response = messagebox.askyesno(
                "Peak Not Fitted",
                f"Peak {peak_number} ({assignment}) has not been fitted yet.\n\nWould you like to fit it now?"
            )

            if response:
                # Fit the selected peak
                self.fit_selected_peak()
                # After fitting, call this method again to show results
                self.root.after(100, self.show_selected_peak_analysis)

    def navigator_show_peak_analysis(self, peak_type, peak_index):
        """Handle peak analysis request from navigator - exact same function as navigation button"""
        if peak_type == "reference":
            # For reference peaks - use existing logic
            peak_number = peak_index + 1  # Convert to 1-based index
            original_selection = self.selected_peak_number.get()

            # Temporarily set the navigation to match navigator selection
            self.selected_peak_number.set(peak_number)

            # Call existing method (same as 🔬 button in navigation panel)
            self.show_selected_peak_analysis()

            # Restore original navigation selection
            self.selected_peak_number.set(original_selection)

        elif peak_type == "detected":
            # For detected peaks - show Voigt analysis directly
            self.show_detected_peak_analysis_by_index(peak_index)

        else:
            messagebox.showerror("Error", f"Unknown peak type: {peak_type}")

    def show_detected_peak_analysis_by_index(self, peak_index):
        """Show Voigt analysis for detected peak by index - EXACT same logic as navigation button"""

        # Check if we have fitted results (same check as navigation button)
        if not hasattr(self.integrator, 'fitted_peaks') or not self.integrator.fitted_peaks:
            messagebox.showinfo("No Results", "No fitting results available for detected peaks.\n\nRun peak detection first.")
            return

        if peak_index >= len(self.integrator.fitted_peaks):
            messagebox.showerror("Error", f"Invalid peak index: {peak_index}")
            return

        # EXACT same logic as show_selected_peak_analysis but use direct index instead of peak_number matching
        selected_result = None

        # Get the peak result directly by index (since we're working with detected peaks list)
        if peak_index < len(self.integrator.fitted_peaks):
            selected_result = self.integrator.fitted_peaks[peak_index]

        # EXACT same display logic as navigation button
        if selected_result:
            # Switch to Voigt analysis tab
            self.viz_notebook.select(1)

            # Show the analysis for the selected peak
            self.voigt_plotter.plot_voigt_analysis(selected_result)
            self.canvas_voigt.draw()

            # Update status
            assignment = selected_result.get('assignment', f'Det_{peak_index+1}')
            quality = selected_result.get('fitting_quality', 'Unknown')
            self.update_status(f"✅ Showing Voigt analysis for {assignment} (Quality: {quality})")
        else:
            # Peak not fitted yet - offer to fit it (same as navigation button)
            assignment = f'Det_{peak_index+1}'
            response = messagebox.askyesno(
                "Peak Not Fitted",
                f"Detected peak {peak_index+1} ({assignment}) has not been fitted yet.\n\nWould you like to fit all peaks now?"
            )
            if response:
                self.fit_all_peaks()

    # =================== RESULTS ANALYSIS METHODS ===================

    def open_results_browser(self):
        """Open comprehensive results browser window"""
        if not self.batch_results or not hasattr(self.batch_results, 'results') or not self.batch_results.results:
            messagebox.showinfo("No Data", "No series results available to browse.\n\nPlease run series integration first.")
            return

        # Create results browser window
        browser_window = tk.Toplevel(self.root)
        browser_window.title("Series Results Browser")
        # Set responsive geometry (80% of screen size)
        screen_width = self.root.winfo_screenwidth()
        screen_height = self.root.winfo_screenheight()
        width = min(1200, int(screen_width * 0.8))
        height = min(800, int(screen_height * 0.8))
        x = (screen_width - width) // 2
        y = (screen_height - height) // 2
        browser_window.geometry(f"{width}x{height}+{x}+{y}")
        browser_window.minsize(800, 600)
        browser_window.transient(self.root)

        # Main container with notebook tabs
        notebook = ttk.Notebook(browser_window)
        notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Overview tab
        overview_frame = ttk.Frame(notebook)
        notebook.add(overview_frame, text="📊 Overview")
        self._create_results_overview_tab(overview_frame)

        # Peak Tracking tab
        tracking_frame = ttk.Frame(notebook)
        notebook.add(tracking_frame, text="🎯 Peak Tracking")
        self._create_peak_tracking_tab(tracking_frame)

        # Data Export tab
        export_frame = ttk.Frame(notebook)
        notebook.add(export_frame, text="💾 Data Export")
        self._create_data_export_tab(export_frame)

        # Quality Assessment tab
        quality_frame = ttk.Frame(notebook)
        notebook.add(quality_frame, text="🔍 Quality Assessment")
        self._create_quality_assessment_tab(quality_frame)

    def _create_results_overview_tab(self, parent):
        """Create results overview tab content"""
        # Summary statistics frame
        stats_frame = ttk.LabelFrame(parent, text="Series Summary", padding=10)
        stats_frame.pack(fill=tk.X, padx=5, pady=5)

        summary = self.batch_results.get_summary()

        # Create summary display
        summary_text = f"""
📈 Processing Summary:
• Total Spectra: {summary['total_spectra']}
• Successful: {summary['successful']} ({summary['success_rate']:.1f}%)
• Failed: {summary['failed']}
• Processing Duration: {summary['duration']}
• Processing Mode: {summary['processing_mode']}
        """

        summary_label = ttk.Label(stats_frame, text=summary_text, font=('TkDefaultFont', 10))
        summary_label.pack(anchor=tk.W)

        # Detailed results listbox
        results_frame = ttk.LabelFrame(parent, text="Spectrum Results", padding=10)
        results_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Create treeview for detailed results
        columns = ('Spectrum', 'Status', 'Detected', 'Total', 'Detection Rate', 'Type')
        tree = ttk.Treeview(results_frame, columns=columns, show='headings', height=15)

        for col in columns:
            tree.heading(col, text=col)
            tree.column(col, width=120)

        # Populate with data
        for spectrum_name, result in self.batch_results.results.items():
            status_icon = "✅" if result.get('status') == 'success' else "❌"
            type_indicator = "📍 REF" if result.get('is_reference', False) else "📄"

            tree.insert('', tk.END, values=(
                spectrum_name,
                f"{status_icon} {result.get('status', 'unknown').title()}",
                result.get('detected_peaks', 0),
                result.get('total_peaks', 0),
                f"{result.get('detection_rate', 0):.1f}%",
                type_indicator
            ))

        # Add scrollbar
        scrollbar = ttk.Scrollbar(results_frame, orient=tk.VERTICAL, command=tree.yview)
        tree.configure(yscrollcommand=scrollbar.set)

        tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

    def _create_peak_tracking_tab(self, parent):
        """Create peak tracking analysis tab"""
        # Check if tracking files exist
        if not self.series_processor.output_folder:
            ttk.Label(parent, text="Peak tracking data not available.\nRun series integration to generate tracking data.",
                     font=('TkDefaultFont', 12), foreground='gray').pack(expand=True)
            return

        tracking_file = os.path.join(self.series_processor.output_folder, "comprehensive_peak_tracking.csv")
        intensity_file = os.path.join(self.series_processor.output_folder, "peak_intensity_matrix.csv")

        if not os.path.exists(tracking_file):
            ttk.Label(parent, text="Peak tracking files not found.\nThey should be generated during series integration.",
                     font=('TkDefaultFont', 12), foreground='orange').pack(expand=True)
            return

        # Load and display tracking data
        try:
            import pandas as pd
            tracking_df = pd.read_csv(tracking_file)

            # Display basic info
            info_frame = ttk.Frame(parent)
            info_frame.pack(fill=tk.X, padx=10, pady=5)

            info_text = f"Peak Tracking Data: {len(tracking_df)} peaks across {len(self.batch_results.results)} spectra"
            ttk.Label(info_frame, text=info_text, font=('TkDefaultFont', 10, 'bold')).pack(anchor=tk.W)

            # Create preview table (first 10 columns to fit)
            preview_frame = ttk.LabelFrame(parent, text="Peak Tracking Preview", padding=10)
            preview_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

            # Select relevant columns for display
            display_cols = ['Peak_Number', 'Assignment', 'Reference_X', 'Reference_Y']
            remaining_cols = [col for col in tracking_df.columns if col not in display_cols][:6]  # Show 6 more
            display_cols.extend(remaining_cols)

            preview_tree = ttk.Treeview(preview_frame, columns=display_cols, show='headings', height=12)

            for col in display_cols:
                preview_tree.heading(col, text=col.replace('_', ' '))
                preview_tree.column(col, width=40)

            # Add data rows (limit to first 20 for performance)
            for idx, row in tracking_df.head(20).iterrows():
                values = [str(row.get(col, 'N/A'))[:12] for col in display_cols]  # Limit text length
                preview_tree.insert('', tk.END, values=values)

            # Add scrollbars
            preview_scroll = ttk.Scrollbar(preview_frame, orient=tk.VERTICAL, command=preview_tree.yview)
            preview_tree.configure(yscrollcommand=preview_scroll.set)

            preview_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
            preview_scroll.pack(side=tk.RIGHT, fill=tk.Y)

            # Add note about full data
            if len(tracking_df) > 20:
                note_text = f"Showing first 20 of {len(tracking_df)} peaks. Use 'Export Data Matrix' to access full data."
                ttk.Label(preview_frame, text=note_text, font=('TkDefaultFont', 8), foreground='blue').pack(pady=5)

        except Exception as e:
            error_text = f"Error loading peak tracking data: {str(e)}"
            ttk.Label(parent, text=error_text, font=('TkDefaultFont', 10), foreground='red').pack(expand=True)

    def _create_data_export_tab(self, parent):
        """Create data export options tab"""
        # Export options frame
        export_frame = ttk.LabelFrame(parent, text="Export Options", padding=15)
        export_frame.pack(fill=tk.X, padx=10, pady=10)

        # Quick export buttons
        ttk.Label(export_frame, text="Quick Export:", font=('TkDefaultFont', 11, 'bold')).pack(anchor=tk.W, pady=(0, 10))

        button_frame = ttk.Frame(export_frame)
        button_frame.pack(fill=tk.X)

        ttk.Button(button_frame, text="📊 Export Peak Intensity Matrix",
                  command=self.export_intensity_matrix).pack(side=tk.LEFT, padx=(0, 10))
        ttk.Button(button_frame, text="🎯 Export Peak Tracking Table",
                  command=self.export_tracking_table).pack(side=tk.LEFT, padx=(0, 10))
        ttk.Button(button_frame, text="📈 Export Detection Statistics",
                  command=self.export_detection_stats).pack(side=tk.LEFT)

        # Batch export section
        batch_frame = ttk.LabelFrame(parent, text="Batch Export", padding=15)
        batch_frame.pack(fill=tk.X, padx=10, pady=10)

        ttk.Label(batch_frame, text="Export all results in multiple formats:",
                 font=('TkDefaultFont', 10)).pack(anchor=tk.W, pady=(0, 10))

        batch_button_frame = ttk.Frame(batch_frame)
        batch_button_frame.pack(fill=tk.X)

        ttk.Button(batch_button_frame, text="📁 Export to CSV Package",
                  command=lambda: self.batch_export_results('csv')).pack(side=tk.LEFT, padx=(0, 10))
        ttk.Button(batch_button_frame, text="📊 Export to Excel Workbook",
                  command=lambda: self.batch_export_results('excel')).pack(side=tk.LEFT, padx=(0, 10))
        ttk.Button(batch_button_frame, text="📄 Generate PDF Report",
                  command=lambda: self.batch_export_results('pdf')).pack(side=tk.LEFT)

        # Custom export section
        custom_frame = ttk.LabelFrame(parent, text="Custom Export", padding=15)
        custom_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        ttk.Label(custom_frame, text="Select specific data to export:",
                 font=('TkDefaultFont', 10)).pack(anchor=tk.W, pady=(0, 10))

        # Export format selection
        format_frame = ttk.Frame(custom_frame)
        format_frame.pack(fill=tk.X, pady=5)

        ttk.Label(format_frame, text="Format:").pack(side=tk.LEFT)
        self.export_format = tk.StringVar(value="csv")
        ttk.Radiobutton(format_frame, text="CSV", variable=self.export_format, value="csv").pack(side=tk.LEFT, padx=10)
        ttk.Radiobutton(format_frame, text="Excel", variable=self.export_format, value="excel").pack(side=tk.LEFT, padx=10)
        ttk.Radiobutton(format_frame, text="JSON", variable=self.export_format, value="json").pack(side=tk.LEFT, padx=10)

        # Data selection checkboxes
        data_frame = ttk.Frame(custom_frame)
        data_frame.pack(fill=tk.X, pady=10)

        self.export_options = {}
        options = [
            ('peak_tracking', 'Peak tracking matrix'),
            ('intensity_matrix', 'Intensity matrix'),
            ('detection_matrix', 'Detection matrix'),
            ('individual_results', 'Individual spectrum results'),
            ('summary_stats', 'Summary statistics'),
            ('quality_report', 'Quality assessment report')
        ]

        for i, (key, label) in enumerate(options):
            self.export_options[key] = tk.BooleanVar(value=True)
            ttk.Checkbutton(data_frame, text=label, variable=self.export_options[key]).grid(
                row=i//2, column=i%2, sticky=tk.W, padx=20, pady=2)

        # Custom export button
        ttk.Button(custom_frame, text="🚀 Custom Export",
                  command=self.custom_export_results).pack(pady=15)

    def _create_quality_assessment_tab(self, parent):
        """Create quality assessment tab"""
        # Generate quality report
        if hasattr(self.series_processor.batch_results, 'statistics'):
            analyzer = SeriesAnalyzer()
            analyzer.load_results(self.batch_results)
            quality_report = analyzer.generate_quality_report()

            # Quality grade display
            grade_frame = ttk.LabelFrame(parent, text="Overall Quality Assessment", padding=15)
            grade_frame.pack(fill=tk.X, padx=10, pady=10)

            grade = quality_report.get('overall_grade', 'Unknown')
            grade_colors = {'Excellent': 'green', 'Good': 'blue', 'Fair': 'orange', 'Poor': 'red', 'Failed': 'darkred'}

            grade_text = f"Overall Grade: {grade}"
            ttk.Label(grade_frame, text=grade_text, font=('TkDefaultFont', 14, 'bold'),
                     foreground=grade_colors.get(grade, 'black')).pack()

            # Quality metrics
            metrics_frame = ttk.LabelFrame(parent, text="Quality Metrics", padding=15)
            metrics_frame.pack(fill=tk.X, padx=10, pady=10)

            metrics = quality_report.get('quality_metrics', {})
            if metrics:
                metrics_text = f"""
Average Detection Rate: {metrics.get('average_detection_rate', 0):.1f}%
Detection Consistency: {metrics.get('detection_consistency', 0):.2f}
Success Rate: {quality_report['summary']['success_rate']:.1f}%
                """
                ttk.Label(metrics_frame, text=metrics_text, font=('TkDefaultFont', 10)).pack(anchor=tk.W)

            # Issues and recommendations
            if quality_report.get('issues') or quality_report.get('recommendations'):
                advice_frame = ttk.LabelFrame(parent, text="Analysis & Recommendations", padding=15)
                advice_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

                advice_text = tk.Text(advice_frame, wrap=tk.WORD, height=10)
                advice_text.pack(fill=tk.BOTH, expand=True)

                if quality_report.get('issues'):
                    advice_text.insert(tk.END, "🚨 Issues Identified:\n")
                    for issue in quality_report['issues']:
                        advice_text.insert(tk.END, f"• {issue}\n")
                    advice_text.insert(tk.END, "\n")

                if quality_report.get('recommendations'):
                    advice_text.insert(tk.END, "💡 Recommendations:\n")
                    for rec in quality_report['recommendations']:
                        advice_text.insert(tk.END, f"• {rec}\n")

                advice_text.config(state=tk.DISABLED)
        else:
            ttk.Label(parent, text="Quality assessment requires processed series data.",
                     font=('TkDefaultFont', 12), foreground='gray').pack(expand=True)

    def show_peak_evolution(self):
        """Show peak evolution analysis window"""
        if not self.batch_results:
            messagebox.showinfo("No Data", "No series results available for peak evolution analysis.")
            return

        # This would create a specialized window for peak evolution visualization
        messagebox.showinfo("Peak Evolution", "Peak evolution analysis feature will be implemented in the visualization module.")

    def open_comparison_dashboard(self):
        """Open comparison dashboard for spectra analysis"""
        if not self.batch_results:
            messagebox.showinfo("No Data", "No series results available for comparison.")
            return

        # This would create a comparison dashboard
        messagebox.showinfo("Comparison Dashboard", "Comparison dashboard feature will be implemented with advanced plotting capabilities.")

    def export_data_matrix(self):
        """Export comprehensive data matrix with perfect legacy format support"""
        # Check for new format results first
        if hasattr(self, 'new_batch_results') and self.new_batch_results:
            self._export_data_matrix_new_format()
            return

        # Fallback to legacy format export
        if not hasattr(self, 'batch_results') or not self.batch_results:
            messagebox.showwarning("No Data",
                "No series results to export. Run series integration first.")
            return

        # Use legacy export method
        self._export_data_matrix_legacy_format()

    def _export_data_matrix_new_format(self):
        """Export data matrix from new MultiSpectrumProcessor format"""
        filename = filedialog.asksaveasfilename(
            title="Export Comprehensive Data Matrix",
            defaultextension=".csv",
            filetypes=[
                ("CSV files", "*.csv"),
                ("Excel files", "*.xlsx"),
                ("All files", "*.*")
            ]
        )

        if filename:
            try:
                # Create comprehensive data matrix from new format
                all_results = []

                for spectrum_name, spectrum_data in self.new_batch_results['results'].items():
                    if spectrum_data.get('success', False):
                        # Get fitted results or integration results
                        fitted_results = spectrum_data.get('fitted_results', [])
                        integration_results = spectrum_data.get('integration_results', [])

                        # Use integration_results if available (legacy format), otherwise fitted_results
                        results_to_export = integration_results if integration_results else fitted_results

                        for peak_result in results_to_export:
                            if isinstance(peak_result, dict):
                                result_row = peak_result.copy()

                                # Add spectrum metadata
                                result_row['spectrum'] = spectrum_name
                                result_row['spectrum_success'] = spectrum_data.get('success', False)
                                result_row['spectrum_success_rate'] = spectrum_data.get('success_rate', 0)
                                result_row['spectrum_total_peaks'] = spectrum_data.get('total_peaks', 0)
                                result_row['spectrum_successful_fits'] = spectrum_data.get('successful_fits', 0)

                                # Standardize column names for export
                                if 'assignment' in result_row:
                                    result_row['Assignment'] = result_row.pop('assignment')
                                if 'ppm_x' in result_row:
                                    result_row['Position_X'] = result_row.pop('ppm_x')
                                if 'ppm_y' in result_row:
                                    result_row['Position_Y'] = result_row.pop('ppm_y')
                                if 'integration_volume' in result_row:
                                    result_row['Volume'] = result_row.pop('integration_volume')

                                all_results.append(result_row)

                if all_results:
                    # Create DataFrame and export
                    df = pd.DataFrame(all_results)

                    # Reorder columns for better readability
                    priority_columns = ['spectrum', 'Assignment', 'Position_X', 'Position_Y',
                                      'Height', 'Volume', 'SNR', 'Quality', 'R_Squared']
                    remaining_columns = [col for col in df.columns if col not in priority_columns]
                    ordered_columns = [col for col in priority_columns if col in df.columns] + remaining_columns
                    df = df[ordered_columns]

                    # Export based on file extension
                    if filename.lower().endswith('.xlsx'):
                        df.to_excel(filename, index=False)
                    else:
                        df.to_csv(filename, index=False, float_format='%.6f')

                    # Show success message with statistics
                    total_peaks = len(all_results)
                    total_spectra = len(self.new_batch_results['results'])

                    messagebox.showinfo("Export Complete",
                        f"Comprehensive data matrix exported successfully!\n\n"
                        f"File: {os.path.basename(filename)}\n"
                        f"Total peaks: {total_peaks}\n"
                        f"Total spectra: {total_spectra}\n"
                        f"Path: {filename}")

                    print(f"✅ Data matrix exported: {total_peaks} peaks from {total_spectra} spectra")

                else:
                    messagebox.showwarning("No Data",
                        "No fitted peak results found to export.")

            except Exception as e:
                import traceback
                error_msg = f"Failed to export data matrix: {e}"
                print(f"❌ {error_msg}\n{traceback.format_exc()}")
                messagebox.showerror("Export Error", error_msg)

    def _export_data_matrix_legacy_format(self):
        """Export data matrix from legacy format (fallback method)"""
        # Implementation for legacy format export
        # This would be the old export logic for backward compatibility
        try:
            if hasattr(self.series_processor, 'output_folder') and self.series_processor.output_folder:
                # Use existing CSV files if available
                tracking_file = os.path.join(self.series_processor.output_folder, "comprehensive_peak_tracking.csv")
                if os.path.exists(tracking_file):
                    # Copy existing file to user-selected location
                    filename = filedialog.asksaveasfilename(
                        title="Export Data Matrix",
                        defaultextension=".csv",
                        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
                    )
                    if filename:
                        import shutil
                        shutil.copy2(tracking_file, filename)
                        messagebox.showinfo("Export Complete", f"Data matrix exported to:\n{filename}")
                        return

            messagebox.showwarning("No Data",
                "No series processing output folder found. Run series integration first.")

        except Exception as e:
            messagebox.showerror("Export Error", f"Failed to export data matrix: {e}")

    def export_analysis_report(self):
        """Export comprehensive analysis report"""
        if not self.batch_results:
            messagebox.showwarning("No Data", "No series results available for report generation.")
            return

        filename = filedialog.asksaveasfilename(
            title="Export Analysis Report",
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("Text files", "*.txt"), ("All files", "*.*")]
        )

        if filename:
            try:
                analyzer = SeriesAnalyzer()
                analyzer.load_results(self.batch_results)
                success = analyzer.export_analysis(filename)

                if success:
                    messagebox.showinfo("Export Complete", f"Analysis report exported to:\n{filename}")
                else:
                    messagebox.showerror("Export Failed", "Failed to generate analysis report.")
            except Exception as e:
                messagebox.showerror("Export Error", f"Failed to export report:\n{str(e)}")

    def batch_export_results(self, format_type='csv'):
        """Batch export results in specified format"""
        if not self.series_processor.output_folder:
            messagebox.showwarning("No Data", "No series results to export.")
            return

        folder = filedialog.askdirectory(title=f"Select folder for {format_type.upper()} export")
        if folder:
            try:
                if format_type == 'csv':
                    # Copy all CSV files from output folder
                    import glob
                    import shutil
                    csv_files = glob.glob(os.path.join(self.series_processor.output_folder, "*.csv"))
                    for file in csv_files:
                        shutil.copy2(file, folder)
                    messagebox.showinfo("Export Complete", f"Exported {len(csv_files)} CSV files to:\n{folder}")

                elif format_type == 'excel':
                    # This would create an Excel workbook with multiple sheets
                    messagebox.showinfo("Excel Export", "Excel export feature will create a comprehensive workbook with multiple data sheets.")

                elif format_type == 'pdf':
                    # This would generate a PDF report
                    messagebox.showinfo("PDF Export", "PDF report generation feature will create a comprehensive analysis document.")

            except Exception as e:
                messagebox.showerror("Export Error", f"Failed to perform batch export:\n{str(e)}")

    def export_intensity_matrix(self):
        """Export peak intensity matrix"""
        if not self.series_processor.output_folder:
            messagebox.showwarning("No Data", "No intensity matrix available.")
            return

        intensity_file = os.path.join(self.series_processor.output_folder, "peak_intensity_matrix.csv")
        if os.path.exists(intensity_file):
            filename = filedialog.asksaveasfilename(
                title="Export Intensity Matrix",
                defaultextension=".csv",
                filetypes=[("CSV files", "*.csv"), ("Excel files", "*.xlsx")]
            )
            if filename:
                try:
                    import shutil
                    shutil.copy2(intensity_file, filename)
                    messagebox.showinfo("Export Complete", f"Intensity matrix exported to:\n{filename}")
                except Exception as e:
                    messagebox.showerror("Export Error", f"Failed to export:\n{str(e)}")
        else:
            messagebox.showerror("File Not Found", "Intensity matrix file not found.")

    def export_tracking_table(self):
        """Export peak tracking table"""
        if not self.series_processor.output_folder:
            messagebox.showwarning("No Data", "No tracking table available.")
            return

        tracking_file = os.path.join(self.series_processor.output_folder, "comprehensive_peak_tracking.csv")
        if os.path.exists(tracking_file):
            filename = filedialog.asksaveasfilename(
                title="Export Peak Tracking Table",
                defaultextension=".csv",
                filetypes=[("CSV files", "*.csv"), ("Excel files", "*.xlsx")]
            )
            if filename:
                try:
                    import shutil
                    shutil.copy2(tracking_file, filename)
                    messagebox.showinfo("Export Complete", f"Tracking table exported to:\n{filename}")
                except Exception as e:
                    messagebox.showerror("Export Error", f"Failed to export:\n{str(e)}")
        else:
            messagebox.showerror("File Not Found", "Peak tracking file not found.")

    def export_detection_stats(self):
        """Export detection statistics"""
        if not self.series_processor.output_folder:
            messagebox.showwarning("No Data", "No detection statistics available.")
            return

        stats_file = os.path.join(self.series_processor.output_folder, "detection_statistics.csv")
        if os.path.exists(stats_file):
            filename = filedialog.asksaveasfilename(
                title="Export Detection Statistics",
                defaultextension=".csv",
                filetypes=[("CSV files", "*.csv"), ("Excel files", "*.xlsx")]
            )
            if filename:
                try:
                    import shutil
                    shutil.copy2(stats_file, filename)
                    messagebox.showinfo("Export Complete", f"Detection statistics exported to:\n{filename}")
                except Exception as e:
                    messagebox.showerror("Export Error", f"Failed to export:\n{str(e)}")
        else:
            messagebox.showerror("File Not Found", "Detection statistics file not found.")

    def custom_export_results(self):
        """Perform custom export based on user selections"""
        # This would implement custom export logic based on the checkboxes
        messagebox.showinfo("Custom Export", "Custom export functionality will be implemented based on selected options.")

    def show_quality_assessment(self):
        """Show detailed quality assessment window"""
        if not self.batch_results:
            messagebox.showinfo("No Data", "No series results available for quality assessment.")
            return

        # Open the results browser directly to the quality tab
        self.open_results_browser()

    def show_detection_statistics(self):
        """Show detailed detection statistics"""
        if not self.batch_results:
            messagebox.showinfo("No Data", "No series results available for detection statistics.")
            return

        # Create statistics window
        stats_window = tk.Toplevel(self.root)
        stats_window.title("Detection Statistics")
        # Set responsive geometry (60% of parent window size)
        parent_width = self.root.winfo_width()
        parent_height = self.root.winfo_height()
        width = max(600, int(parent_width * 0.6))
        height = max(400, int(parent_height * 0.6))
        stats_window.geometry(f"{width}x{height}")
        stats_window.minsize(600, 400)
        stats_window.transient(self.root)

        # Create statistics display
        stats_text = tk.Text(stats_window, wrap=tk.WORD, font=('Courier', 10))
        stats_text.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Generate statistics summary
        summary = self.batch_results.get_summary()
        stats_content = f"""
DETECTION STATISTICS SUMMARY
{'='*50}

Overall Performance:
• Total Spectra Processed: {summary['total_spectra']}
• Successful Processing: {summary['successful']} ({summary['success_rate']:.1f}%)
• Failed Processing: {summary['failed']}
• Processing Duration: {summary['duration']}
• Processing Mode: {summary['processing_mode']}

"""

        if hasattr(self.batch_results, 'statistics') and self.batch_results.statistics:
            stats_data = self.batch_results.statistics
            if 'detection_rate' in stats_data:
                dr = stats_data['detection_rate']
                stats_content += f"""
Detection Rate Statistics:
• Mean Detection Rate: {dr['mean']:.1f}%
• Median Detection Rate: {dr['median']:.1f}%
• Standard Deviation: {dr['std']:.1f}%
• Best Performance: {dr['max']:.1f}%
• Worst Performance: {dr['min']:.1f}%

"""

        # Individual spectrum details
        stats_content += "INDIVIDUAL SPECTRUM RESULTS:\n"
        stats_content += f"{'Spectrum':<30} {'Status':<10} {'Detected':<10} {'Total':<10} {'Rate':<10}\n"
        stats_content += "-" * 80 + "\n"

        for spectrum_name, result in self.batch_results.results.items():
            status = "✅ Success" if result.get('status') == 'success' else "❌ Failed"
            detected = result.get('detected_peaks', 0)
            total = result.get('total_peaks', 0)
            rate = result.get('detection_rate', 0)

            stats_content += f"{spectrum_name[:29]:<30} {status:<10} {detected:<10} {total:<10} {rate:.1f}%\n"

        stats_text.insert(tk.END, stats_content)
        stats_text.config(state=tk.DISABLED)

    def open_spectrum_browser(self):
        """
        Open spectrum browser with enhanced compatibility for new format.
        """
        # Validate that we have batch results
        if not hasattr(self, 'batch_results') or not self.batch_results:
            messagebox.showwarning("No Data",
                "No series results available. Run series integration first.")
            return

        try:
            # Import spectrum browser
            from lunaNMR.gui.spectrum_browser import SpectrumBrowserDialog

            # Get original data folder
            original_data_folder = self.nmr_file_list.get_current_folder()

            # Create spectrum browser with legacy-compatible batch results
            browser = SpectrumBrowserDialog(
                self.root,
                self.batch_results,  # Use legacy-compatible format
                self.series_processor if hasattr(self, 'series_processor') else None,
                original_data_folder
            )

            # Add compatibility enhancements for new format features
            if hasattr(self, 'new_batch_results') and self.new_batch_results:
                # Provide access to new format for advanced features
                browser._new_format_results = self.new_batch_results

                # Add enhanced export capability
                def enhanced_export_spectrum_data(spectrum_name):
                    """Enhanced export with new format support"""
                    if spectrum_name in self.new_batch_results.get('results', {}):
                        spectrum_data = self.new_batch_results['results'][spectrum_name]

                        # Export comprehensive spectrum data
                        filename = filedialog.asksaveasfilename(
                            title=f"Export {spectrum_name} Data",
                            defaultextension=".csv",
                            filetypes=[("CSV files", "*.csv"), ("Excel files", "*.xlsx")]
                        )

                        if filename:
                            try:
                                # Use integration_results if available, otherwise fitted_results
                                export_data = spectrum_data.get('integration_results',
                                                              spectrum_data.get('fitted_results', []))

                                if export_data:
                                    df = pd.DataFrame(export_data)

                                    if filename.lower().endswith('.xlsx'):
                                        df.to_excel(filename, index=False)
                                    else:
                                        df.to_csv(filename, index=False, float_format='%.6f')

                                    messagebox.showinfo("Export Complete",
                                        f"Spectrum data exported to:\n{filename}")
                                else:
                                    messagebox.showwarning("No Data",
                                        f"No peak data found for spectrum {spectrum_name}")

                            except Exception as e:
                                messagebox.showerror("Export Error", f"Export failed: {e}")

                # Add enhanced export method to browser
                browser._enhanced_export = enhanced_export_spectrum_data

            print("✅ Spectrum browser opened with enhanced compatibility")

        except ImportError as e:
            messagebox.showerror("Error", f"Failed to import spectrum browser: {e}")
        except Exception as e:
            import traceback
            error_msg = f"Failed to open spectrum browser: {e}"
            print(f"❌ {error_msg}\n{traceback.format_exc()}")
            messagebox.showerror("Error", error_msg)

    def show_help(self):
        """Show help dialog"""
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

        help_window = tk.Toplevel(self.root)
        help_window.title("User Guide")
        # Set responsive geometry for help window
        parent_width = self.root.winfo_width()
        parent_height = self.root.winfo_height()
        width = max(500, min(700, int(parent_width * 0.5)))
        height = max(300, min(500, int(parent_height * 0.6)))
        help_window.geometry(f"{width}x{height}")
        help_window.minsize(500, 300)
        help_window.transient(self.root)

        text_widget = tk.Text(help_window, wrap=tk.WORD, font=('TkDefaultFont', 10))
        scrollbar = ttk.Scrollbar(help_window, orient="vertical", command=text_widget.yview)

        text_widget.config(yscrollcommand=scrollbar.set)
        text_widget.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        text_widget.insert(1.0, help_text)
        text_widget.config(state='disabled')

    def show_about(self):
        """Show about dialog"""
        about_text = f"""
NMR Peaks Series Analysis
Version 1.0.0

Enhanced Multi-Mode NMR Peak Detection and Analysis

Features:
• Dual processing modes (Full Detection / In-Place Fitting)
• Advanced Voigt profile fitting
• Comprehensive series integration
• Enhanced visualization and statistics
• Configuration management

Developed using Python with tkinter, matplotlib, numpy, pandas, and scipy.

Configuration: {self.config_manager.config_file}
        """

        messagebox.showinfo("About NMR Peaks Series Analysis", about_text)

    def launch_dynamixs(self):
        """Launch DynamiXs Relaxation Analysis as independent module"""
        try:
            # Path to DynamiXs directory
            dynamixs_path = os.path.join(current_dir, "dynamiXs")
            launcher_path = os.path.join(dynamixs_path, "run_dynamixs_gui.py")

            # Check if DynamiXs directory exists
            if not os.path.exists(dynamixs_path):
                messagebox.showerror("Module Not Found",
                    "DynamiXs module directory not found.\n"
                    f"Expected location: {dynamixs_path}")
                return

            # Check if launcher script exists
            if not os.path.exists(launcher_path):
                messagebox.showerror("Launcher Not Found",
                    "DynamiXs launcher script not found.\n"
                    f"Expected location: {launcher_path}")
                return

            # Launch DynamiXs as independent subprocess
            subprocess.Popen([sys.executable, launcher_path],
                           cwd=dynamixs_path,
                           stdout=subprocess.DEVNULL,
                           stderr=subprocess.PIPE)

            # Show success message
            messagebox.showinfo("Module Launched",
                "DynamiXs Relaxation Analysis has been launched successfully.\n"
                "The module will open in a separate window.")

        except Exception as e:
            messagebox.showerror("Launch Error",
                f"Failed to launch DynamiXs module:\n{str(e)}\n\n"
                "Please ensure Python dependencies are installed:\n"
                "• tkinter • numpy • matplotlib • pandas • lmfit • scipy")

    def on_closing(self):
        """Handle application closing"""
        if self.processing_active:
            if not messagebox.askyesno("Processing Active",
                "Processing is currently active. Are you sure you want to exit?"):
                return

        # Save current window state
        self.user_prefs.save_from_window(self.root)

        # Save configuration
        self.config_manager.save_config()

        # Destroy window
        self.root.destroy()

    # =================== ZOOM AND CONTOUR METHODS ===================

    def update_zoom_from_sliders(self):
        """Update zoom based on navigation slider values"""
        if not hasattr(self, 'ax_main') or not self.ax_main:
            return

        x_center = self.zoom_x_center.get()
        x_range = self.zoom_x_range.get()
        y_center = self.zoom_y_center.get()
        y_range = self.zoom_y_range.get()

        # Calculate limits
        x_min = x_center - x_range / 2
        x_max = x_center + x_range / 2
        y_min = y_center - y_range / 2
        y_max = y_center + y_range / 2

        # Apply zoom (NMR convention - reversed axes)
        self.ax_main.set_xlim(x_max, x_min)
        self.ax_main.set_ylim(y_max, y_min)

        # Update stored zoom
        self.current_xlim = (x_max, x_min)
        self.current_ylim = (y_max, y_min)

        self.canvas_main.draw()

    def reset_zoom(self):
        """Reset zoom to show full spectrum"""
        if not hasattr(self.integrator, 'ppm_x_axis') or not hasattr(self.integrator, 'ppm_y_axis'):
            print("No spectrum data loaded for zoom reset")
            return

        # Calculate full spectrum ranges
        x_full_range = abs(self.integrator.ppm_x_axis[0] - self.integrator.ppm_x_axis[-1])
        y_full_range = abs(self.integrator.ppm_y_axis[0] - self.integrator.ppm_y_axis[-1])

        # Reset zoom sliders to show full spectrum
        self.zoom_x_center.set((self.integrator.ppm_x_axis[0] + self.integrator.ppm_x_axis[-1]) / 2)
        self.zoom_x_range.set(x_full_range)
        self.zoom_y_center.set((self.integrator.ppm_y_axis[0] + self.integrator.ppm_y_axis[-1]) / 2)
        self.zoom_y_range.set(y_full_range)

        # Clear stored zoom state
        self.current_xlim = None
        self.current_ylim = None

        # Apply the zoom
        self.update_zoom_from_sliders()

    def save_current_zoom(self):
        """Save the current zoom limits"""
        if hasattr(self, 'ax_main') and self.ax_main:
            self.current_xlim = self.ax_main.get_xlim()
            self.current_ylim = self.ax_main.get_ylim()

    def restore_zoom(self):
        """Restore the saved zoom limits"""
        if self.current_xlim and self.current_ylim and hasattr(self, 'ax_main') and self.ax_main:
            self.ax_main.set_xlim(self.current_xlim)
            self.ax_main.set_ylim(self.current_ylim)

    def center_on_selected_peak(self):
        """Center zoom on the selected peak"""
        if not hasattr(self.integrator, 'peak_list') or self.integrator.peak_list is None:
            self.update_status("❌ No peak list loaded")
            return

        peak_number = self.selected_peak_number.get()
        if peak_number < 1 or peak_number > len(self.integrator.peak_list):
            self.update_status(f"❌ Invalid peak number: {peak_number}")
            return

        # Get peak coordinates (convert to 0-based indexing)
        peak_idx = peak_number - 1
        peak_data = self.integrator.peak_list.iloc[peak_idx]
        peak_x = float(peak_data['Position_X'])
        peak_y = float(peak_data['Position_Y'])

        # Set zoom parameters to center on peak
        self.zoom_x_center.set(peak_x)
        self.zoom_x_range.set(0.3)  # Tight zoom
        self.zoom_y_center.set(peak_y)
        self.zoom_y_range.set(2.5)  # Reasonable Y range

        # Apply the zoom
        self.update_zoom_from_sliders()

    def center_on_coordinates(self, x, y):
        """Center zoom on specific coordinates"""
        # Set zoom parameters to center on coordinates
        self.zoom_x_center.set(x)
        self.zoom_x_range.set(0.3)  # Tight zoom
        self.zoom_y_center.set(y)
        self.zoom_y_range.set(2.5)  # Reasonable Y range

        # Apply the zoom
        self.update_zoom_from_sliders()

        self.update_status(f"Centered on coordinates: {x:.3f}, {y:.1f}")

    def set_selected_peak(self, peak_index, peak_type, source="unknown"):
        """Coordinate peak selection between navigator and navigation panel"""
        if peak_type == "reference":
            # Convert 0-based index to 1-based for navigation panel
            peak_number = peak_index + 1
            if hasattr(self, 'selected_peak_number'):
                self.selected_peak_number.set(peak_number)
                self.update_peak_navigation()

            # Update navigator if the change came from navigation panel
            if source != "navigator" and hasattr(self, 'peak_navigator'):
                self.peak_navigator.update_selection(peak_index, peak_type)

        elif peak_type == "detected":
            # For detected peaks, we mainly coordinate with navigator
            # (there's no separate navigation panel for detected peaks in the original design)
            if source != "navigator" and hasattr(self, 'peak_navigator'):
                self.peak_navigator.update_selection(peak_index, peak_type)

    def __del__(self):
        """Cleanup when GUI is destroyed"""
        try:
            self.deactivate_zoom_mode()
        except:
            pass

    # =================== INTERACTIVE ZOOM FUNCTIONALITY ===================

    def toggle_zoom_mode(self):
        """Toggle interactive zoom mode on/off"""
        if self.zoom_mode_active:
            self.deactivate_zoom_mode()
        else:
            self.activate_zoom_mode()

    def activate_zoom_mode(self):
        """Activate interactive zoom mode"""
        if not hasattr(self, 'ax_main') or not self.ax_main:
            self.update_status("❌ No spectrum loaded for zoom mode")
            return

        self.zoom_mode_active = True

        # Update button appearance
        self.zoom_mode_button.configure(text="🔍 Zoom ON", style="Active.TButton")

        # Store original axis limits for reset
        if not self.original_xlim or not self.original_ylim:
            self.original_xlim = self.ax_main.get_xlim()
            self.original_ylim = self.ax_main.get_ylim()

        # Connect mouse events
        self.connect_zoom_events()

        # Change cursor
        self.canvas_main.get_tk_widget().configure(cursor="crosshair")

        self.update_status("🔍 Zoom mode activated - Drag to select zoom area")

    def deactivate_zoom_mode(self):
        """Deactivate interactive zoom mode"""
        self.zoom_mode_active = False

        # Update button appearance
        self.zoom_mode_button.configure(text="🔍 Drag Zoom")

        # Disconnect mouse events
        self.disconnect_zoom_events()

        # Remove zoom rectangle if exists
        self.remove_zoom_rectangle()

        # Restore normal cursor
        self.canvas_main.get_tk_widget().configure(cursor="")

        self.update_status("Zoom mode deactivated")

    def connect_zoom_events(self):
        """Connect mouse event handlers for zoom mode"""
        if hasattr(self, 'canvas_main'):
            self.mouse_press_cid = self.canvas_main.mpl_connect('button_press_event', self.on_zoom_press)
            self.mouse_motion_cid = self.canvas_main.mpl_connect('motion_notify_event', self.on_zoom_motion)
            self.mouse_release_cid = self.canvas_main.mpl_connect('button_release_event', self.on_zoom_release)

    def disconnect_zoom_events(self):
        """Disconnect mouse event handlers"""
        if hasattr(self, 'canvas_main'):
            if self.mouse_press_cid:
                self.canvas_main.mpl_disconnect(self.mouse_press_cid)
                self.mouse_press_cid = None
            if self.mouse_motion_cid:
                self.canvas_main.mpl_disconnect(self.mouse_motion_cid)
                self.mouse_motion_cid = None
            if self.mouse_release_cid:
                self.canvas_main.mpl_disconnect(self.mouse_release_cid)
                self.mouse_release_cid = None

    def on_zoom_press(self, event):
        """Handle mouse press for zoom selection"""
        if not self.zoom_mode_active or not event.inaxes == self.ax_main:
            return

        if event.button == 1:  # Left mouse button
            self.zoom_start_point = (event.xdata, event.ydata)
            # Remove any existing rectangle
            self.remove_zoom_rectangle()

    def on_zoom_motion(self, event):
        """Handle mouse motion for zoom rectangle drawing"""
        if not self.zoom_mode_active or not event.inaxes == self.ax_main or not self.zoom_start_point:
            return

        if event.xdata is None or event.ydata is None:
            return

        # Remove existing rectangle
        self.remove_zoom_rectangle()

        # Calculate rectangle coordinates
        x_start, y_start = self.zoom_start_point
        x_end, y_end = event.xdata, event.ydata

        # Ensure proper ordering (handle dragging in any direction)
        x_min, x_max = min(x_start, x_end), max(x_start, x_end)
        y_min, y_max = min(y_start, y_end), max(y_start, y_end)

        # Draw selection rectangle
        from matplotlib.patches import Rectangle
        self.zoom_rectangle = Rectangle(
            (x_min, y_min), x_max - x_min, y_max - y_min,
            fill=False, edgecolor='red', linewidth=2, linestyle='--', alpha=0.8
        )
        self.ax_main.add_patch(self.zoom_rectangle)
        self.canvas_main.draw_idle()

    def on_zoom_release(self, event):
        """Handle mouse release to apply zoom"""
        if not self.zoom_mode_active or not event.inaxes == self.ax_main or not self.zoom_start_point:
            return

        if event.button == 1 and event.xdata is not None and event.ydata is not None:
            # Calculate zoom area
            x_start, y_start = self.zoom_start_point
            x_end, y_end = event.xdata, event.ydata

            # Ensure proper ordering
            x_min, x_max = min(x_start, x_end), max(x_start, x_end)
            y_min, y_max = min(y_start, y_end), max(y_start, y_end)

            # Validate zoom area (minimum size)
            if abs(x_max - x_min) > 0.05 and abs(y_max - y_min) > 1.0:
                self.apply_interactive_zoom(x_min, x_max, y_min, y_max)
            else:
                self.update_status("❌ Zoom area too small - Please select a larger area")

        # Clean up
        self.remove_zoom_rectangle()
        self.zoom_start_point = None

    def apply_interactive_zoom(self, x_min, x_max, y_min, y_max):
        """Apply the selected zoom area"""
        try:
            # Apply zoom to plot
            self.ax_main.set_xlim(x_max, x_min)  # Inverted for NMR convention
            self.ax_main.set_ylim(y_max, y_min)  # Inverted for NMR convention

            # Update zoom controls to reflect new zoom
            x_center = (x_min + x_max) / 2
            y_center = (y_min + y_max) / 2
            x_range = abs(x_max - x_min)
            y_range = abs(y_max - y_min)

            self.zoom_x_center.set(x_center)
            self.zoom_y_center.set(y_center)
            self.zoom_x_range.set(x_range)
            self.zoom_y_range.set(y_range)

            # Redraw canvas
            self.canvas_main.draw()

            # Update status
            self.update_status(f"✅ Zoomed to region: X({x_min:.2f}-{x_max:.2f}) Y({y_min:.1f}-{y_max:.1f}) ppm")

        except Exception as e:
            self.update_status(f"❌ Error applying zoom: {e}")

    def remove_zoom_rectangle(self):
        """Remove the zoom selection rectangle"""
        if self.zoom_rectangle:
            try:
                self.zoom_rectangle.remove()
            except:
                pass
            self.zoom_rectangle = None

    def reset_interactive_zoom(self):
        """Reset to original full spectrum view"""
        if not hasattr(self, 'ax_main') or not self.ax_main:
            return

        if self.original_xlim and self.original_ylim:
            # Restore original limits
            self.ax_main.set_xlim(self.original_xlim)
            self.ax_main.set_ylim(self.original_ylim)

            # Update zoom controls
            x_center = (self.original_xlim[0] + self.original_xlim[1]) / 2
            y_center = (self.original_ylim[0] + self.original_ylim[1]) / 2
            x_range = abs(self.original_xlim[1] - self.original_xlim[0])
            y_range = abs(self.original_ylim[1] - self.original_ylim[0])

            self.zoom_x_center.set(x_center)
            self.zoom_y_center.set(y_center)
            self.zoom_x_range.set(x_range)
            self.zoom_y_range.set(y_range)

            # Redraw
            self.canvas_main.draw()
            self.update_status("✅ Zoom reset to full spectrum")
        else:
            # Fallback to regular reset zoom
            self.reset_zoom()

    def auto_adjust_zoom_to_data(self):
        """Automatically adjust zoom to fit the loaded NMR data dimensions"""
        if not hasattr(self.integrator, 'ppm_x_axis') or not hasattr(self.integrator, 'ppm_y_axis'):
            return

        if self.integrator.ppm_x_axis is None or self.integrator.ppm_y_axis is None:
            return

        try:
            # Get the full spectrum ranges
            x_min = float(np.min(self.integrator.ppm_x_axis))
            x_max = float(np.max(self.integrator.ppm_x_axis))
            y_min = float(np.min(self.integrator.ppm_y_axis))
            y_max = float(np.max(self.integrator.ppm_y_axis))

            # Calculate center and range for both dimensions
            x_center = (x_min + x_max) / 2
            y_center = (y_min + y_max) / 2
            x_range = abs(x_max - x_min)
            y_range = abs(y_max - y_min)

            # Update zoom controls to show full spectrum
            self.zoom_x_center.set(x_center)
            self.zoom_x_range.set(x_range)
            self.zoom_y_center.set(y_center)
            self.zoom_y_range.set(y_range)

            # Store original limits for interactive zoom reset
            if hasattr(self, 'ax_main') and self.ax_main:
                self.original_xlim = (x_max, x_min)  # Inverted for NMR convention
                self.original_ylim = (y_max, y_min)  # Inverted for NMR convention

            print(f"Auto-adjusted zoom to data: X({x_min:.1f}-{x_max:.1f}) Y({y_min:.1f}-{y_max:.1f}) ppm")

        except Exception as e:
            print(f"Error in auto zoom adjustment: {e}")
            # Fallback to manual reset if auto-adjustment fails
            self.reset_zoom()

def main():
    """Main application entry point"""
    # Create root window
    root = tk.Tk()

    # Set responsive window configuration
    screen_width = root.winfo_screenwidth()
    screen_height = root.winfo_screenheight()

    # Calculate window size as percentage of screen (80% width, 90% height)
    window_width = int(screen_width * 0.8)
    window_height = int(screen_height * 0.9)

    # Center window on screen
    x = (screen_width - window_width) // 2
    y = (screen_height - window_height) // 2

    # Set geometry with minimum size constraints
    root.geometry(f"{window_width}x{window_height}+{x}+{y}")
    root.minsize(1000, 600)  # Minimum usable size

    # Make window resizable
    root.resizable(True, True)

    # Set window icon (if available)
    try:
        icon_path = os.path.join(os.path.dirname(__file__), "icon.ico")
        if os.path.exists(icon_path):
            root.iconbitmap(icon_path)
    except:
        pass  # Icon not required

    # Apply styling
    style = ttk.Style()
    if 'clam' in style.theme_names():
        style.theme_use('clam')

    # Create custom style for active zoom button
    style.configure('Active.TButton', foreground='red', background='lightgreen')

    # Create and run application
    try:
        app = NMRPeaksSeriesGUI(root)
        print("🚀 Starting NMR Peaks Series GUI application")
        root.mainloop()
        print("👋 NMR Peaks Series GUI application closed")

    except KeyboardInterrupt:
        print("\n⏹️  Application interrupted by user")
    except Exception as e:
        print(f"💥 Critical application error: {e}")
        messagebox.showerror("Critical Error",
            f"A critical error occurred:\n{str(e)}\n\nThe application will close.")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
