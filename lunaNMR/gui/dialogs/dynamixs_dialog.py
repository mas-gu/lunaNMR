# ABOUTME: Embedded DynamiXs dialog wrapper for lunaNMR integration
# ABOUTME: Directly imports existing DynamiXs Qt pages - NO code duplication

import os
import sys
from pathlib import Path
from typing import Optional, Dict, Any, List

from PySide6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QStackedWidget, QWidget, QFrame, QLabel
)
from PySide6.QtCore import Signal, Qt

from lunaNMR.gui.base.base_dialog import BaseDialog

# Add DynamiXs module to path for imports
_dynamixs_path = Path(__file__).parent.parent.parent.parent / "modules" / "dynamiXs_v2o0"
if str(_dynamixs_path) not in sys.path:
    sys.path.insert(0, str(_dynamixs_path))

# Direct imports from existing DynamiXs
from dynamiXs_gui import (
    T1T2FittingPage,
    SpectralDensityPage,
    PlottingPage,
    IntegratedAnalysisPage
)
from gui_components import (
    create_v_layout, create_h_layout,
    create_header_label, create_secondary_label, create_label,
    create_primary_button, create_secondary_button,
    get_font, open_directory_dialog
)
from constants import (
    SPACING_XS, SPACING_SM, SPACING_MD, SPACING_XL,
    PANEL_BG_COLOR, SECONDARY_TEXT, FONT_SIZE_SMALL
)


class DynamiXsDialog(BaseDialog):
    """Embedded DynamiXs dialog - wraps existing DynamiXs pages.

    This dialog provides the full DynamiXs NMR relaxation analysis workflow
    embedded within the lunaNMR application. It directly imports and uses
    the existing DynamiXs page classes for:
    - T1/T2 Fitting Analysis
    - Spectral Density Analysis
    - Plotting
    - Integrated Model-Free Analysis

    State Management:
        The dialog supports get_state()/set_state() and get_file_refs()/set_file_refs()
        methods for project save/load functionality. State is automatically saved
        to main_window when the dialog closes.
    """

    analysis_complete = Signal(dict)

    def __init__(self, parent=None, main_window=None):
        """Initialize the DynamiXs dialog.

        Args:
            parent: Parent widget
            main_window: Reference to main lunaNMR window (for state sync)
        """
        super().__init__(
            parent=parent,
            title="DynamiXs - NMR Relaxation Analysis Suite",
            default_size=(900, 700),
            min_size=(800, 600),
            modal=False
        )

        # Store reference to main window
        self.main_window = main_window

        # Current working directory
        self.current_dir = os.getcwd()

        # Reference to directory label for updates
        self.dir_label = None

        # Path to dynamiXs module
        self.dynamixs_path = _dynamixs_path

        # Setup UI
        self._setup_ui()

        # Restore state if available
        if main_window:
            if hasattr(main_window, 'dynamixs_state') and main_window.dynamixs_state:
                self.set_state(main_window.dynamixs_state)
            if hasattr(main_window, 'dynamixs_file_refs') and main_window.dynamixs_file_refs:
                self.set_file_refs(main_window.dynamixs_file_refs)

    def _setup_ui(self):
        """Setup the UI by embedding existing DynamiXs pages."""
        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        # Create stacked widget for page navigation
        self.stack = QStackedWidget()

        # Create main menu page
        self.main_menu_page = self._create_main_menu_page()
        self.stack.addWidget(self.main_menu_page)

        # Create all analysis pages (reusing existing classes)
        # Pass 'self' as main_window since we provide the same navigation interface
        self.t1t2_page = T1T2FittingPage(self)
        self.spectral_page = SpectralDensityPage(self)
        self.plotting_page = PlottingPage(self)
        self.integrated_page = IntegratedAnalysisPage(self)

        self.stack.addWidget(self.t1t2_page)
        self.stack.addWidget(self.spectral_page)
        self.stack.addWidget(self.plotting_page)
        self.stack.addWidget(self.integrated_page)

        layout.addWidget(self.stack)
        self.setLayout(layout)

        # Show main menu initially
        self.show_main_menu()

    def _create_main_menu_page(self) -> QWidget:
        """Create the main menu page with navigation buttons.

        This recreates the DynamiXs main menu layout for navigation
        between different analysis types.
        """
        page = QWidget()
        layout = create_v_layout(SPACING_MD, (SPACING_XL, SPACING_XL, SPACING_XL, SPACING_XL))
        page.setLayout(layout)

        # Title section
        title_frame = QFrame()
        title_layout = create_v_layout(SPACING_XS)
        title_frame.setLayout(title_layout)

        title_label = create_header_label("DynamiXs - NMR Relaxation Analysis Suite", level=1)
        title_layout.addWidget(title_label)

        subtitle_label = create_secondary_label("Choose your analysis type:")
        title_layout.addWidget(subtitle_label)

        layout.addWidget(title_frame)

        # Main options panel
        options_panel = QFrame()
        options_panel.setProperty("class", "panel")
        options_layout = create_v_layout(SPACING_SM, (SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD))
        options_panel.setLayout(options_layout)

        # Model Free Analysis button (primary)
        btn_integrated = create_primary_button(
            "Model Free Analysis",
            clicked=self.show_integrated_analysis,
            width=200
        )
        btn_integrated.setFont(get_font(18, bold=True))
        options_layout.addWidget(btn_integrated, alignment=Qt.AlignCenter)

        # Separator
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setStyleSheet(f"background-color: {PANEL_BG_COLOR}; max-height: 2px;")
        options_layout.addWidget(separator)

        # T1/T2 Fitting button
        btn_t1t2 = create_primary_button(
            "T1/T2 Fitting Analysis",
            clicked=self.show_t1t2_fitting,
            width=200
        )
        btn_t1t2.setFont(get_font(18, bold=True))
        options_layout.addWidget(btn_t1t2, alignment=Qt.AlignCenter)

        # Spectral Density button
        btn_spectral = create_primary_button(
            "Spectral Density Analysis",
            clicked=self.show_spectral_density,
            width=200
        )
        btn_spectral.setFont(get_font(18, bold=True))
        options_layout.addWidget(btn_spectral, alignment=Qt.AlignCenter)

        # Plot Data button
        btn_plot = create_secondary_button(
            "Plot Data",
            clicked=self.show_plotting,
            width=200
        )
        btn_plot.setFont(get_font(18, bold=True))
        options_layout.addWidget(btn_plot, alignment=Qt.AlignCenter)

        layout.addWidget(options_panel, stretch=1)

        # Working directory info
        dir_frame = QFrame()
        dir_frame.setProperty("class", "card")
        dir_layout = create_h_layout(SPACING_SM, (SPACING_MD, SPACING_SM, SPACING_MD, SPACING_SM))
        dir_frame.setLayout(dir_layout)

        self.dir_label = create_label(
            f"Working Directory: {self.current_dir}",
            font_size=FONT_SIZE_SMALL
        )
        self.dir_label.setStyleSheet(f"color: {SECONDARY_TEXT};")
        dir_layout.addWidget(self.dir_label, stretch=1)

        change_dir_btn = create_secondary_button(
            "Change Directory",
            clicked=self.change_working_directory,
            width=140
        )
        dir_layout.addWidget(change_dir_btn)

        layout.addWidget(dir_frame)

        return page

    # -------------------------------------------------------------------------
    # Navigation Methods (same interface as DynamiXsMainWindow)
    # -------------------------------------------------------------------------

    def show_main_menu(self):
        """Show the main menu page."""
        self.stack.setCurrentWidget(self.main_menu_page)

    def show_t1t2_fitting(self):
        """Show T1/T2 fitting page."""
        self.stack.setCurrentWidget(self.t1t2_page)

    def show_spectral_density(self):
        """Show spectral density analysis page."""
        self.stack.setCurrentWidget(self.spectral_page)

    def show_plotting(self):
        """Show plotting page."""
        self.stack.setCurrentWidget(self.plotting_page)

    def show_integrated_analysis(self):
        """Show integrated analysis page."""
        self.stack.setCurrentWidget(self.integrated_page)

    def change_working_directory(self):
        """Allow user to change working directory."""
        new_dir = open_directory_dialog(
            self,
            title="Select Working Directory",
            directory=self.current_dir
        )
        if new_dir:
            self.current_dir = new_dir
            os.chdir(new_dir)
            if self.dir_label:
                self.dir_label.setText(f"Working Directory: {self.current_dir}")

    # -------------------------------------------------------------------------
    # State Management (for project save/load)
    # -------------------------------------------------------------------------

    def get_state(self) -> Dict[str, Any]:
        """Collect state from all pages for project save.

        Returns:
            Dictionary containing state from all pages
        """
        return {
            'active_page': self.stack.currentIndex(),
            'current_dir': self.current_dir,
            't1t2': self._get_t1t2_state(),
            'spectral': self._get_spectral_state(),
            'integrated': self._get_integrated_state(),
        }

    def set_state(self, state: Dict[str, Any]):
        """Restore state to all pages from project load.

        Args:
            state: Dictionary containing state for all pages
        """
        if 'current_dir' in state:
            self.current_dir = state['current_dir']
            if self.dir_label:
                self.dir_label.setText(f"Working Directory: {self.current_dir}")

        if 't1t2' in state:
            self._set_t1t2_state(state['t1t2'])
        if 'spectral' in state:
            self._set_spectral_state(state['spectral'])
        if 'integrated' in state:
            self._set_integrated_state(state['integrated'])
        if 'active_page' in state:
            self.stack.setCurrentIndex(state['active_page'])

    def get_file_refs(self) -> Dict[str, Any]:
        """Collect file references from all pages.

        Returns:
            Dictionary containing file paths from all pages
        """
        return {
            't1t2': self._get_t1t2_file_refs(),
            'spectral': self._get_spectral_file_refs(),
            'integrated': self._get_integrated_file_refs(),
        }

    def set_file_refs(self, refs: Dict[str, Any]):
        """Restore file references to all pages.

        Args:
            refs: Dictionary containing file paths for all pages
        """
        if 't1t2' in refs:
            self._set_t1t2_file_refs(refs['t1t2'])
        if 'spectral' in refs:
            self._set_spectral_file_refs(refs['spectral'])
        if 'integrated' in refs:
            self._set_integrated_file_refs(refs['integrated'])

    # -------------------------------------------------------------------------
    # T1/T2 Page State Helpers
    # -------------------------------------------------------------------------

    def _get_t1t2_state(self) -> Dict[str, Any]:
        """Extract state from T1/T2 fitting page using session state API."""
        page = self.t1t2_page
        if hasattr(page, 'get_session_state'):
            return page.get_session_state()
        return {}

    def _set_t1t2_state(self, state: Dict[str, Any]):
        """Restore state to T1/T2 fitting page using session state API."""
        page = self.t1t2_page
        if hasattr(page, 'restore_session_state'):
            page.restore_session_state(state)

    def _get_t1t2_file_refs(self) -> Dict[str, Any]:
        """Extract file references from T1/T2 fitting page."""
        page = self.t1t2_page
        refs = {}

        # Field 1 files
        if hasattr(page, 'field1_t1_file') and page.field1_t1_file:
            refs['field1_t1_file'] = page.field1_t1_file
        if hasattr(page, 'field1_t2_file') and page.field1_t2_file:
            refs['field1_t2_file'] = page.field1_t2_file
        if hasattr(page, 'field1_t1rho_file') and page.field1_t1rho_file:
            refs['field1_t1rho_file'] = page.field1_t1rho_file

        # Field 2 files
        if hasattr(page, 'field2_t1_file') and page.field2_t1_file:
            refs['field2_t1_file'] = page.field2_t1_file
        if hasattr(page, 'field2_t2_file') and page.field2_t2_file:
            refs['field2_t2_file'] = page.field2_t2_file
        if hasattr(page, 'field2_t1rho_file') and page.field2_t1rho_file:
            refs['field2_t1rho_file'] = page.field2_t1rho_file

        # T1rho specific files (legacy)
        if hasattr(page, 'peak_list_file') and page.peak_list_file:
            refs['peak_list_file'] = page.peak_list_file

        return refs

    def _set_t1t2_file_refs(self, refs: Dict[str, Any]):
        """Restore file references to T1/T2 fitting page."""
        page = self.t1t2_page

        # Field 1 files
        if 'field1_t1_file' in refs:
            page.field1_t1_file = refs['field1_t1_file']
            if hasattr(page, 'field1_t1_display'):
                page.field1_t1_display.setText(os.path.basename(refs['field1_t1_file']))
        if 'field1_t2_file' in refs:
            page.field1_t2_file = refs['field1_t2_file']
            if hasattr(page, 'field1_t2_display'):
                page.field1_t2_display.setText(os.path.basename(refs['field1_t2_file']))
        if 'field1_t1rho_file' in refs:
            page.field1_t1rho_file = refs['field1_t1rho_file']
            if hasattr(page, 'field1_t1rho_display'):
                page.field1_t1rho_display.setText(os.path.basename(refs['field1_t1rho_file']))

        # Field 2 files
        if 'field2_t1_file' in refs:
            page.field2_t1_file = refs['field2_t1_file']
            if hasattr(page, 'field2_t1_display'):
                page.field2_t1_display.setText(os.path.basename(refs['field2_t1_file']))
        if 'field2_t2_file' in refs:
            page.field2_t2_file = refs['field2_t2_file']
            if hasattr(page, 'field2_t2_display'):
                page.field2_t2_display.setText(os.path.basename(refs['field2_t2_file']))
        if 'field2_t1rho_file' in refs:
            page.field2_t1rho_file = refs['field2_t1rho_file']
            if hasattr(page, 'field2_t1rho_display'):
                page.field2_t1rho_display.setText(os.path.basename(refs['field2_t1rho_file']))

        # Legacy T1rho files
        if 'peak_list_file' in refs:
            page.peak_list_file = refs['peak_list_file']
            if hasattr(page, 'peak_list_display'):
                page.peak_list_display.setText(os.path.basename(refs['peak_list_file']))

    # -------------------------------------------------------------------------
    # Spectral Density Page State Helpers
    # -------------------------------------------------------------------------

    def _get_spectral_state(self) -> Dict[str, Any]:
        """Extract state from Spectral Density page using session state API."""
        page = self.spectral_page
        if hasattr(page, 'get_session_state'):
            return page.get_session_state()
        return {}

    def _set_spectral_state(self, state: Dict[str, Any]):
        """Restore state to Spectral Density page using session state API."""
        page = self.spectral_page
        if hasattr(page, 'restore_session_state'):
            page.restore_session_state(state)

    def _get_spectral_file_refs(self) -> Dict[str, Any]:
        """Extract file references from Spectral Density page."""
        page = self.spectral_page
        refs = {}

        if hasattr(page, 'input_file1') and page.input_file1:
            refs['input_file1'] = page.input_file1
        if hasattr(page, 'input_file2') and page.input_file2:
            refs['input_file2'] = page.input_file2

        return refs

    def _set_spectral_file_refs(self, refs: Dict[str, Any]):
        """Restore file references to Spectral Density page."""
        page = self.spectral_page

        if 'input_file1' in refs:
            page.input_file1 = refs['input_file1']
            if hasattr(page, 'file1_display'):
                import os
                page.file1_display.setText(os.path.basename(refs['input_file1']))
        if 'input_file2' in refs:
            page.input_file2 = refs['input_file2']
            if hasattr(page, 'file2_display'):
                import os
                page.file2_display.setText(os.path.basename(refs['input_file2']))

    # -------------------------------------------------------------------------
    # Integrated Analysis Page State Helpers
    # -------------------------------------------------------------------------

    def _get_integrated_state(self) -> Dict[str, Any]:
        """Extract state from Integrated Analysis page using session state API."""
        page = self.integrated_page
        if hasattr(page, 'get_session_state'):
            return page.get_session_state()
        return {}

    def _set_integrated_state(self, state: Dict[str, Any]):
        """Restore state to Integrated Analysis page using session state API."""
        page = self.integrated_page
        if hasattr(page, 'restore_session_state'):
            page.restore_session_state(state)

    def _get_integrated_file_refs(self) -> Dict[str, Any]:
        """Extract file references from Integrated Analysis page."""
        page = self.integrated_page
        refs = {}

        # Field 1 files
        if hasattr(page, 'field1_t1_file') and page.field1_t1_file:
            refs['field1_t1'] = page.field1_t1_file
        if hasattr(page, 'field1_t2_file') and page.field1_t2_file:
            refs['field1_t2'] = page.field1_t2_file
        if hasattr(page, 'field1_noe_sat_file') and page.field1_noe_sat_file:
            refs['field1_noe_sat'] = page.field1_noe_sat_file
        if hasattr(page, 'field1_noe_unsat_file') and page.field1_noe_unsat_file:
            refs['field1_noe_unsat'] = page.field1_noe_unsat_file

        # Field 2 files
        if hasattr(page, 'field2_t1_file') and page.field2_t1_file:
            refs['field2_t1'] = page.field2_t1_file
        if hasattr(page, 'field2_t2_file') and page.field2_t2_file:
            refs['field2_t2'] = page.field2_t2_file
        if hasattr(page, 'field2_noe_sat_file') and page.field2_noe_sat_file:
            refs['field2_noe_sat'] = page.field2_noe_sat_file
        if hasattr(page, 'field2_noe_unsat_file') and page.field2_noe_unsat_file:
            refs['field2_noe_unsat'] = page.field2_noe_unsat_file

        return refs

    def _set_integrated_file_refs(self, refs: Dict[str, Any]):
        """Restore file references to Integrated Analysis page."""
        page = self.integrated_page
        import os

        # Field 1 files
        if 'field1_t1' in refs:
            page.field1_t1_file = refs['field1_t1']
            if hasattr(page, 'field1_t1_display'):
                page.field1_t1_display.setText(os.path.basename(refs['field1_t1']))
        if 'field1_t2' in refs:
            page.field1_t2_file = refs['field1_t2']
            if hasattr(page, 'field1_t2_display'):
                page.field1_t2_display.setText(os.path.basename(refs['field1_t2']))
        if 'field1_noe_sat' in refs:
            page.field1_noe_sat_file = refs['field1_noe_sat']
            if hasattr(page, 'field1_noe_sat_display'):
                page.field1_noe_sat_display.setText(os.path.basename(refs['field1_noe_sat']))
        if 'field1_noe_unsat' in refs:
            page.field1_noe_unsat_file = refs['field1_noe_unsat']
            if hasattr(page, 'field1_noe_unsat_display'):
                page.field1_noe_unsat_display.setText(os.path.basename(refs['field1_noe_unsat']))

        # Field 2 files
        if 'field2_t1' in refs:
            page.field2_t1_file = refs['field2_t1']
            if hasattr(page, 'field2_t1_display'):
                page.field2_t1_display.setText(os.path.basename(refs['field2_t1']))
        if 'field2_t2' in refs:
            page.field2_t2_file = refs['field2_t2']
            if hasattr(page, 'field2_t2_display'):
                page.field2_t2_display.setText(os.path.basename(refs['field2_t2']))
        if 'field2_noe_sat' in refs:
            page.field2_noe_sat_file = refs['field2_noe_sat']
            if hasattr(page, 'field2_noe_sat_display'):
                page.field2_noe_sat_display.setText(os.path.basename(refs['field2_noe_sat']))
        if 'field2_noe_unsat' in refs:
            page.field2_noe_unsat_file = refs['field2_noe_unsat']
            if hasattr(page, 'field2_noe_unsat_display'):
                page.field2_noe_unsat_display.setText(os.path.basename(refs['field2_noe_unsat']))

    # -------------------------------------------------------------------------
    # Result Files Management
    # -------------------------------------------------------------------------

    def get_output_directories(self) -> Dict[str, str]:
        """Get output directories that may contain result files.

        Returns:
            Dict mapping page name to output directory path
        """
        dirs = {}

        # T1T2 page - output_dir attribute or input file directory
        if hasattr(self.t1t2_page, 'output_dir') and self.t1t2_page.output_dir:
            dirs['t1t2'] = self.t1t2_page.output_dir
        elif hasattr(self.t1t2_page, 'input_file') and self.t1t2_page.input_file:
            import os
            dirs['t1t2'] = os.path.dirname(self.t1t2_page.input_file)

        # Integrated page - output_dir attribute
        if hasattr(self.integrated_page, 'output_dir') and self.integrated_page.output_dir:
            dirs['integrated'] = self.integrated_page.output_dir

        return dirs

    def get_result_files(self) -> Dict[str, List[str]]:
        """Get list of result files from output directories.

        Returns:
            Dict mapping page name to list of result file paths (*.csv, *.txt)
        """
        import os
        import glob

        result_files = {}
        output_dirs = self.get_output_directories()

        for page_name, output_dir in output_dirs.items():
            if output_dir and os.path.isdir(output_dir):
                files = []
                # Collect CSV and TXT result files
                files.extend(glob.glob(os.path.join(output_dir, "*.csv")))
                files.extend(glob.glob(os.path.join(output_dir, "*.txt")))
                # Filter to only include actual result files (not input files)
                result_patterns = [
                    '_fit_results', '_basic', '_detailed', '_results',
                    '_T2_from_T1rho', '_spectral_density', '_series_tidy',
                    '_data', '_intensities', 'hetNOE_ratios'
                ]
                filtered = [f for f in files
                           if any(p in os.path.basename(f) for p in result_patterns)]
                if filtered:
                    result_files[page_name] = filtered

        return result_files

    def set_output_directory(self, page_name: str, output_dir: str):
        """Set output directory for a page after project load.

        Args:
            page_name: 't1t2', 'spectral', or 'integrated'
            output_dir: New output directory path
        """
        import os

        if page_name == 't1t2' and hasattr(self.t1t2_page, 'output_dir'):
            self.t1t2_page.output_dir = output_dir
        elif page_name == 'integrated':
            if hasattr(self.integrated_page, 'output_dir'):
                self.integrated_page.output_dir = output_dir
            if hasattr(self.integrated_page, 'outdir_display'):
                self.integrated_page.outdir_display.setText(os.path.basename(output_dir))

    # -------------------------------------------------------------------------
    # Dialog Event Handlers
    # -------------------------------------------------------------------------

    def closeEvent(self, event):
        """Save state when dialog closes."""
        if self.main_window:
            self.main_window.dynamixs_state = self.get_state()
            self.main_window.dynamixs_file_refs = self.get_file_refs()
        super().closeEvent(event)
