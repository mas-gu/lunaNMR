# ABOUTME: Embedded DynamiXs dialog wrapper for lunaNMR integration
# ABOUTME: Directly imports existing DynamiXs Qt pages - NO code duplication

import os
import sys
from pathlib import Path
from typing import Optional, Dict, Any, List

from PySide6.QtWidgets import (
    QVBoxLayout, QStackedWidget, QWidget, QFrame, QInputDialog
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

        # Current analysis name (required before entering analysis pages)
        self.current_analysis_name = None

        # Analysis metadata (source_series, analysis_type, created_at, etc.)
        self._analysis_metadata = {}

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

        from methyl_t2_page import MethylT2FittingPage
        self.methyl_t2_page = MethylT2FittingPage(self)

        self.stack.addWidget(self.t1t2_page)
        self.stack.addWidget(self.spectral_page)
        self.stack.addWidget(self.plotting_page)
        self.stack.addWidget(self.integrated_page)
        self.stack.addWidget(self.methyl_t2_page)

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

        # T2 Methyl Fitting button (bi-exponential)
        btn_methyl_t2 = create_primary_button(
            "T2 Methyl Fitting Analysis",
            clicked=self.show_methyl_t2_fitting,
            width=200
        )
        btn_methyl_t2.setFont(get_font(18, bold=True))
        options_layout.addWidget(btn_methyl_t2, alignment=Qt.AlignCenter)

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
        """Show T1/T2 fitting page.

        Prompts for analysis name if not already set.
        """
        if not self.current_analysis_name:
            name = self._prompt_analysis_name("T1/T2 Fitting")
            if name is None:
                return  # User cancelled
            self.current_analysis_name = name

        self.stack.setCurrentWidget(self.t1t2_page)

    def show_methyl_t2_fitting(self):
        """Show methyl bi-exp T2 fitting page."""
        self.stack.setCurrentWidget(self.methyl_t2_page)

    def show_spectral_density(self):
        """Show spectral density analysis page."""
        self.stack.setCurrentWidget(self.spectral_page)

    def show_plotting(self):
        """Show plotting page."""
        self.stack.setCurrentWidget(self.plotting_page)

    def show_integrated_analysis(self):
        """Show integrated analysis page.

        Prompts for analysis name if not already set.
        """
        if not self.current_analysis_name:
            name = self._prompt_analysis_name("Model-Free Analysis")
            if name is None:
                return  # User cancelled
            self.current_analysis_name = name

        self.stack.setCurrentWidget(self.integrated_page)

    def _prompt_analysis_name(self, analysis_type: str) -> Optional[str]:
        """Prompt user for analysis name before starting analysis.

        Args:
            analysis_type: Type of analysis (for dialog title)

        Returns:
            Analysis name if provided, None if cancelled.
        """
        # Generate suggested name
        suggested_name = "Analysis_1"

        # Check for existing analyses to avoid duplicates
        existing_analyses = set()
        if self.main_window and hasattr(self.main_window, 'saved_dynamixs'):
            existing_analyses = set(self.main_window.saved_dynamixs.keys())

        if suggested_name in existing_analyses:
            # Append number to make unique
            base_name = "Analysis"
            counter = 2
            while f"{base_name}_{counter}" in existing_analyses:
                counter += 1
            suggested_name = f"{base_name}_{counter}"

        # Show input dialog
        name, ok = QInputDialog.getText(
            self,
            "Analysis Name",
            f"Enter a name for this {analysis_type}:",
            text=suggested_name
        )

        if not ok or not name.strip():
            return None

        # Validate name (no special characters that could cause file issues)
        name = name.strip()
        invalid_chars = '<>:"/\\|?*'
        for char in invalid_chars:
            name = name.replace(char, '_')

        return name

    def clear_analysis_name(self):
        """Clear current analysis name (e.g., when starting fresh)."""
        self.current_analysis_name = None
        self._analysis_metadata = {}

    def get_analysis_metadata(self) -> Dict[str, Any]:
        """Get metadata for current analysis.

        Returns:
            Dict with analysis metadata including:
            - analysis_name: Name of this analysis
            - analysis_type: 't1t2' or 'model_free'
            - source_series: Name of series that provided input data (if known)
            - created_at: ISO timestamp when analysis was started
        """
        from datetime import datetime

        metadata = getattr(self, '_analysis_metadata', {}).copy()

        # Always include current analysis name
        metadata['analysis_name'] = self.current_analysis_name

        # Ensure created_at is set
        if 'created_at' not in metadata and self.current_analysis_name:
            metadata['created_at'] = datetime.now().isoformat()

        # Determine analysis type from current page
        current_widget = self.stack.currentWidget()
        if current_widget == self.t1t2_page:
            metadata['analysis_type'] = 't1t2'
        elif current_widget == self.integrated_page:
            metadata['analysis_type'] = 'model_free'
        elif current_widget == self.spectral_page:
            metadata['analysis_type'] = 'spectral_density'

        return metadata

    def set_analysis_metadata(self, metadata: Dict[str, Any]):
        """Set metadata for current analysis.

        Args:
            metadata: Dict with analysis metadata (source_series, etc.)
        """
        self._analysis_metadata = metadata.copy()
        if 'analysis_name' in metadata:
            self.current_analysis_name = metadata['analysis_name']

    def set_source_series(self, series_name: str):
        """Set the source series for this analysis.

        Called when DynamiXs is launched from series results.

        Args:
            series_name: Name of the series that provided input data
        """
        if not hasattr(self, '_analysis_metadata'):
            self._analysis_metadata = {}
        self._analysis_metadata['source_series'] = series_name

    def detect_source_series_from_path(self, file_path: str) -> Optional[str]:
        """Detect source series from a file path.

        Checks if the file is within a known series folder structure:
        - <project>/.lunaNMR/series_results/<series_name>/
        - Or if the path contains a series name from saved_series

        Args:
            file_path: Path to a CSV file loaded in DynamiXs

        Returns:
            Series name if detected, None otherwise
        """
        if not file_path:
            return None

        file_path = os.path.normpath(file_path)

        # Check for series_results folder pattern
        if 'series_results' in file_path:
            parts = file_path.split(os.sep)
            try:
                series_idx = parts.index('series_results')
                if series_idx + 1 < len(parts):
                    series_name = parts[series_idx + 1]
                    return series_name
            except ValueError:
                pass

        # Check against saved series names if available
        if self.main_window and hasattr(self.main_window, 'saved_series'):
            for series_name in self.main_window.saved_series.keys():
                # Check if series name appears in path
                if series_name in file_path:
                    return series_name

        return None

    def try_auto_detect_source_series(self):
        """Try to auto-detect source series from loaded file paths.

        Checks file references in the current state for series patterns.
        Sets source_series in metadata if detected.
        """
        file_refs = self.get_file_refs()

        for page_name, refs in file_refs.items():
            if not refs:
                continue

            for field_name, file_path in refs.items():
                if file_path:
                    series_name = self.detect_source_series_from_path(file_path)
                    if series_name:
                        self.set_source_series(series_name)
                        return series_name

        return None

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
            'analysis_name': self.current_analysis_name,
            'analysis_metadata': self.get_analysis_metadata(),
            't1t2': self._get_t1t2_state(),
            'spectral': self._get_spectral_state(),
            'integrated': self._get_integrated_state(),
            'methyl_t2': self._get_methyl_t2_state(),
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

        if 'analysis_name' in state:
            self.current_analysis_name = state['analysis_name']

        if 'analysis_metadata' in state:
            self.set_analysis_metadata(state['analysis_metadata'])

        if 't1t2' in state:
            self._set_t1t2_state(state['t1t2'])
        if 'spectral' in state:
            self._set_spectral_state(state['spectral'])
        if 'integrated' in state:
            self._set_integrated_state(state['integrated'])
        if 'methyl_t2' in state:
            self._set_methyl_t2_state(state['methyl_t2'])
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
            'methyl_t2': self._get_methyl_t2_file_refs(),
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
        if 'methyl_t2' in refs:
            self._set_methyl_t2_file_refs(refs['methyl_t2'])

    # -------------------------------------------------------------------------
    # Named analyses (project save/reopen, parity with the Kd module)
    # -------------------------------------------------------------------------

    _ANALYSIS_TYPE_LABELS = {
        't1t2': 'T1T2', 'methyl_t2': 'methylT2',
        'spectral': 'spectralDensity', 'integrated': 'modelFree',
    }

    def _active_page_key(self) -> Optional[str]:
        """Key of the page the user is currently on (a result page), or None."""
        w = self.stack.currentWidget()
        return {self.t1t2_page: 't1t2', self.spectral_page: 'spectral',
                self.integrated_page: 'integrated',
                self.methyl_t2_page: 'methyl_t2'}.get(w)

    @staticmethod
    def _page_has_results(state: Dict[str, Any], page: Optional[str]) -> bool:
        """Whether the given page's state holds a completed run worth saving."""
        ps = (state.get(page) or {}) if page else {}
        if page == 't1t2':
            return bool(ps.get('fitted_experiments') or ps.get('output_dir'))
        if page == 'methyl_t2':
            return bool(ps.get('last_results_file'))
        if page in ('spectral', 'integrated'):
            sr = ps.get('session_results') or {}
            return bool(sr.get('analysis_complete') or ps.get('output_dir'))
        return False

    def _analysis_save_name(self, state: Dict[str, Any]) -> Optional[str]:
        """Name for the active run: <source_series>_<type> (e.g. HSPA1A_T1T2), or None
        if the active page has no completed results to save."""
        page = self._active_page_key()
        if not self._page_has_results(state, page):
            return None
        meta = state.get('analysis_metadata') or {}
        series = str(meta.get('source_series') or '').strip()
        type_label = self._ANALYSIS_TYPE_LABELS.get(page, page)
        base = f"{series}_{type_label}" if series else (self.current_analysis_name or type_label)
        return base.strip('_ ') or type_label

    def ensure_current_saved(self):
        """Auto-capture the active page's current run into the project's saved DynamiXs
        analyses (upsert by <series>_<type>). No-op if the active page has no results."""
        state = self.get_state()
        name = self._analysis_save_name(state)
        if not name:
            return None
        mw = self.main_window
        if getattr(mw, 'dynamixs_analyses', None) is None:
            mw.dynamixs_analyses = {}
        mw.dynamixs_analyses[name] = {'state': state, 'file_refs': self.get_file_refs()}
        return name

    def open_analysis(self, entry: Dict[str, Any], name: Optional[str] = None):
        """Reopen a saved DynamiXs analysis: restore its state (which switches to its
        page and repoints result paths) and its file references."""
        if not isinstance(entry, dict):
            return
        self.set_state(entry.get('state') or {})
        refs = entry.get('file_refs')
        if refs:
            self.set_file_refs(refs)
        if name:
            self.current_analysis_name = name

    # -------------------------------------------------------------------------
    # Methyl T2 Page State Helpers
    # -------------------------------------------------------------------------

    def _get_methyl_t2_state(self) -> Dict[str, Any]:
        """Extract state from Methyl T2 fitting page using session state API."""
        page = self.methyl_t2_page
        if hasattr(page, 'get_session_state'):
            return page.get_session_state()
        return {}

    def _set_methyl_t2_state(self, state: Dict[str, Any]):
        """Restore state to Methyl T2 fitting page using session state API."""
        page = self.methyl_t2_page
        if hasattr(page, 'restore_session_state'):
            page.restore_session_state(state)

    def _get_methyl_t2_file_refs(self) -> Dict[str, Any]:
        """Extract file references from Methyl T2 fitting page."""
        page = self.methyl_t2_page
        refs = {}
        if getattr(page, 'input_file', None):
            refs['input_file'] = page.input_file
        return refs

    def _set_methyl_t2_file_refs(self, refs: Dict[str, Any]):
        """Restore file references to Methyl T2 fitting page."""
        page = self.methyl_t2_page
        if 'input_file' in refs:
            page.input_file = refs['input_file']
            if hasattr(page, 'file_drop'):
                page.file_drop.setText(os.path.basename(refs['input_file']))

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

    # -------------------------------------------------------------------------
    # Inspect Peak Integration (FitViewer → MultiSpectrumViewer)
    # -------------------------------------------------------------------------

    def connect_fit_viewer_signals(self, fit_viewer):
        """Connect FitViewer's inspect_peak_requested signal to handler.

        Called after a FitViewer is opened to enable the Inspect Peak flow.

        Args:
            fit_viewer: The FitViewer instance to connect signals from
        """
        if fit_viewer and hasattr(fit_viewer, 'inspect_peak_requested'):
            fit_viewer.inspect_peak_requested.connect(self._on_inspect_peak_requested)

    def _on_inspect_peak_requested(self, residue_id: str, series_name: str, all_fit_data: dict = None):
        """Handle Inspect Peak request from FitViewer.

        Opens the MultiSpectrumViewerDialog with the selected residue,
        using the specified series name directly.

        Args:
            residue_id: Residue ID from FitViewer (e.g., "142")
            series_name: Name of the NMR series to inspect (e.g., "T1_asyn_series")
            all_fit_data: Dict mapping residue ID to exponential fit data from DynamiXs.
                E.g., {'142': {fit_data}, '143': {fit_data}, ...}
                Each fit_data contains: time_points, intensities, fit_curve, t_value, t_error, etc.
        """
        from PySide6.QtWidgets import QMessageBox

        if not series_name:
            QMessageBox.warning(
                self,
                "No Series",
                "No series name provided.\n\n"
                "Please assign a series using drag-and-drop in the Project Series panel."
            )
            return

        # Get series data from main_window
        if not self.main_window or not hasattr(self.main_window, 'saved_series'):
            QMessageBox.warning(
                self,
                "No Series Data",
                "Cannot open MultiSpectrum Viewer: No series data available."
            )
            return

        series_data = self.main_window.saved_series.get(series_name)
        if not series_data:
            QMessageBox.warning(
                self,
                "Series Not Found",
                f"Cannot find series '{series_name}' in saved series.\n\n"
                "Available series: " + ", ".join(self.main_window.saved_series.keys())
            )
            return

        # Open MultiSpectrumViewerDialog with the series data and all fit data
        self._open_multi_spectrum_viewer(series_data, residue_id, all_fit_data)

    def _get_source_series_for_experiment(self, experiment_type: str) -> Optional[str]:
        """Get source series name for a specific experiment type.

        Looks up the series from:
        1. T1T2FittingPage's source_series_map (from drag-and-drop)
        2. File paths used for the fit (detects series from path)

        Args:
            experiment_type: "t1", "t2", or "t1rho"

        Returns:
            Series name or None if not found
        """
        # Try T1T2FittingPage's source_series_map first (from drag-and-drop)
        if hasattr(self, 't1t2_page') and hasattr(self.t1t2_page, 'source_series_map'):
            source_map = self.t1t2_page.source_series_map
            for field in ['field1', 'field2']:
                key = f"{field}_{experiment_type}"
                if key in source_map:
                    return source_map[key]

        # Try detecting from file paths used in T1T2FittingPage
        if hasattr(self, 't1t2_page'):
            for field in ['field1', 'field2']:
                file_attr = f"{field}_{experiment_type}_file"
                file_path = getattr(self.t1t2_page, file_attr, None)
                if file_path:
                    series_name = self.detect_source_series_from_path(file_path)
                    if series_name:
                        return series_name

        # Try IntegratedAnalysisPage (Model Free)
        if hasattr(self, 'integrated_page'):
            # Check source_series_map
            if hasattr(self.integrated_page, 'source_series_map'):
                source_map = self.integrated_page.source_series_map
                for field in ['field1', 'field2']:
                    key = f"{field}_{experiment_type}"
                    if key in source_map:
                        return source_map[key]

            # Try detecting from file paths
            for field in ['field1', 'field2']:
                file_attr = f"{field}_{experiment_type}_file"
                file_path = getattr(self.integrated_page, file_attr, None)
                if file_path:
                    series_name = self.detect_source_series_from_path(file_path)
                    if series_name:
                        return series_name

        return None

    def _open_multi_spectrum_viewer(self, series_data, initial_assignment: str = None,
                                      exponential_fit_data: dict = None):
        """Open MultiSpectrumViewerDialog with series data.

        Args:
            series_data: Series integration results (BatchResults object)
            initial_assignment: Optional assignment to auto-select in Peak Mode
            exponential_fit_data: Dict mapping residue ID to exponential fit data from DynamiXs.
                E.g., {'142': {fit_data}, '143': {fit_data}, ...}
        """
        try:
            from lunaNMR.gui.dialogs.multi_spectrum_viewer_dialog import MultiSpectrumViewerDialog

            # Convert BatchResults to list format expected by MultiSpectrumViewerDialog
            all_results = self._convert_batch_results_to_list(series_data)

            if not all_results:
                from PySide6.QtWidgets import QMessageBox
                QMessageBox.warning(
                    self,
                    "No Results",
                    "The series data contains no spectrum results."
                )
                return

            # Get file_manager from main_window if available
            file_manager = None
            if self.main_window and hasattr(self.main_window, 'file_manager'):
                file_manager = self.main_window.file_manager

            # Create and show the viewer
            self.multi_spectrum_viewer = MultiSpectrumViewerDialog(
                parent=self,
                all_results=all_results,
                file_manager=file_manager,
                initial_assignment=initial_assignment,
                exponential_fit_data=exponential_fit_data
            )
            self.multi_spectrum_viewer.show()

        except Exception as e:
            from PySide6.QtWidgets import QMessageBox
            QMessageBox.critical(
                self,
                "Error",
                f"Failed to open MultiSpectrum Viewer:\n{str(e)}"
            )

    def _convert_batch_results_to_list(self, batch_results) -> List[Dict[str, Any]]:
        """Convert BatchResults object to list format for MultiSpectrumViewerDialog.

        Args:
            batch_results: BatchResults object with .results dict

        Returns:
            List of result dicts with spectrum_name, spectrum_file, fitted_peaks
        """
        import os

        all_results = []

        # Get data_folder for constructing full paths (same as main_window.open_multi_spectrum_viewer)
        data_folder = None
        if hasattr(batch_results, 'metadata'):
            data_folder = batch_results.metadata.get('data_folder')
        elif isinstance(batch_results, dict):
            metadata = batch_results.get('metadata', {})
            data_folder = metadata.get('data_folder')

        # Fallback to main_window.current_nmr_folder if not in metadata
        if not data_folder and self.main_window:
            data_folder = getattr(self.main_window, 'current_nmr_folder', None)

        # BatchResults has .results attribute which is a dict keyed by spectrum_name
        if hasattr(batch_results, 'results') and batch_results.results:
            for spectrum_name, result_data in batch_results.results.items():
                if not isinstance(result_data, dict):
                    continue

                # Skip failed spectra
                if not result_data.get('success', True):
                    continue

                # Extract fitted peaks (may be 'fitted_results' or 'integration_results')
                fitted_peaks = result_data.get('fitted_results', [])
                if not fitted_peaks:
                    fitted_peaks = result_data.get('integration_results', [])

                # Construct full path for spectrum_file (same logic as main_window)
                spectrum_file = result_data.get('spectrum_file', '')
                if not spectrum_file or not os.path.isabs(spectrum_file):
                    if data_folder:
                        spectrum_file = os.path.join(data_folder, spectrum_name)

                result_entry = {
                    'spectrum_name': spectrum_name,
                    'spectrum_file': spectrum_file,
                    'fitted_peaks': fitted_peaks,
                    'integration_results': fitted_peaks,
                }
                all_results.append(result_entry)

        return all_results

    def _normalize_peak_data(self, peak: Dict) -> Dict:
        """Normalize peak data to ensure numeric values are floats, not lists.

        After JSON deserialization, numpy scalars may become single-element lists.
        This method converts them back to proper floats.

        Args:
            peak: Peak data dict

        Returns:
            Normalized peak dict with proper float values
        """
        if not isinstance(peak, dict):
            return peak

        normalized = {}
        # Keys that should be scalar floats
        scalar_keys = {
            # Position keys (various naming conventions)
            'pos_f1', 'pos_f2', 'position_x', 'position_y',
            'center_x', 'center_y', 'peak_x', 'peak_y',
            'ppm_x', 'ppm_y', 'Position_X', 'Position_Y',
            # Linewidth keys
            'lw_gau_f1', 'lw_gau_f2', 'lw_lor_f1', 'lw_lor_f2',
            'sigma_x', 'sigma_y', 'gamma_x', 'gamma_y',
            'gaussian_fwhm_x', 'gaussian_fwhm_y',
            'lorentzian_fwhm_x', 'lorentzian_fwhm_y',
            # Intensity/amplitude keys
            'intensity', 'amplitude', 'volume', 'height',
            # Quality metrics
            'r_squared', 'avg_r_squared', 'r_squared_x', 'r_squared_y',
            'r_squared_local',
            # Other scalars
            'baseline', 'noise_level', 'window_size',
        }

        for key, value in peak.items():
            if key in scalar_keys:
                normalized[key] = self._to_float(value)
            elif isinstance(value, dict):
                # Recursively normalize nested dicts (like 'all_peaks', 'region_2d')
                normalized[key] = self._normalize_peak_data(value)
            elif isinstance(value, list) and len(value) > 0 and isinstance(value[0], dict):
                # Normalize list of peak dicts (like 'all_peaks')
                normalized[key] = [self._normalize_peak_data(p) for p in value]
            else:
                normalized[key] = value

        return normalized

    def _to_float(self, value) -> float:
        """Convert a value to float, handling lists and numpy types.

        Args:
            value: A float, int, list, or numpy scalar

        Returns:
            Float value
        """
        if value is None:
            return 0.0
        if isinstance(value, (list, tuple)):
            if len(value) == 0:
                return 0.0
            return float(value[0])
        try:
            return float(value)
        except (TypeError, ValueError):
            return 0.0
