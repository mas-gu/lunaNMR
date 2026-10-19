#!/usr/bin/env python3
"""
DynamiXs v2.0 GUI - A comprehensive interface for NMR relaxation data analysis

This GUI provides access to:
- T1/T2 fitting analysis
- Spectral density function analysis
- Data plotting and comparison tools
- Data formatting utilities

Built with PySide6 (Qt) following LunaNMR v1.0 Style Guide
for cross-platform compatibility (macOS, Windows, Linux)
"""

import os
import sys
import glob
import time
from pathlib import Path
from typing import Optional, List, Dict, Any

from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QPushButton, QLineEdit, QTextEdit, QPlainTextEdit,
    QFrame, QGroupBox, QScrollArea, QStackedWidget,
    QCheckBox, QRadioButton, QButtonGroup, QComboBox,
    QSpinBox, QDoubleSpinBox, QProgressBar,
    QFileDialog, QMessageBox, QSizePolicy, QListWidget, QListWidgetItem,
    QAbstractItemView
)
from PySide6.QtCore import Qt, Signal, Slot, QThread, QTimer, QMimeData
from PySide6.QtGui import QFont, QDragEnterEvent, QDropEvent

# Import design constants
from constants import (
    BG_COLOR, PANEL_BG_COLOR, FRAME_BG_COLOR,
    PRIMARY_TEXT, SECONDARY_TEXT,
    SUCCESS_GREEN,
    SPACING_XS, SPACING_SM, SPACING_MD, SPACING_LG, SPACING_XL,
    FONT_SIZE_LARGE_HEADER, FONT_SIZE_SECTION_LABEL, FONT_SIZE_BODY, FONT_SIZE_SMALL,
    WINDOW_MIN_WIDTH, WINDOW_MIN_HEIGHT
)

# Import GUI components
from gui_components import (
    load_stylesheet, get_font, get_monospace_font,
    CollapsibleGroupBox, StyledGroupBox,
    create_primary_button, create_secondary_button, create_destructive_button,
    create_success_button, create_icon_button,
    create_label, create_header_label, create_secondary_label,
    create_line_edit, create_text_edit, create_plain_text_edit,
    create_spin_box, create_double_spin_box, create_combo_box,
    create_check_box, create_radio_button,
    create_h_layout, create_v_layout, create_grid_layout, create_scroll_area,
    show_info, show_warning, show_error, ask_yes_no,
    open_file_dialog, open_files_dialog, save_file_dialog, open_directory_dialog
)

# Import workers for background tasks
from workers import (
    T1T2FittingWorker, T1T2FittingParams,
    SpectralDensityWorker, SpectralDensityParams,
    IntegratedAnalysisWorker, IntegratedAnalysisParams
)
from dynamiXs_T1_T2.delay_units import (
    TIME_UNITS_DEFAULT, resolve_time_units, series_sidecar
)
from dynamiXs_T1_T2.series_inputs import tidy_table_refusal

# The delay units a page can state for its series, and the factor turning ms into each:
# the initial-guess fields are in ms, the fitters work in the series' own units.
DELAY_UNITS = ('ms', 's', 'us')
MS_TO_UNITS = {'ms': 1.0, 's': 1e-3, 'us': 1e3}

# The intensity a page fits, its series column, and the series matrix that holds it.
OBSERVABLES = (('Height', 'height', 'peak_intensity_matrix.csv'),
               ('Volume', 'volume', 'peak_volume_matrix.csv'))


# =============================================================================
# LUNАНMR IMPORT HELPER
# =============================================================================

def _ensure_lunaNMR_path():
    """Ensure lunaNMR package is importable by adding its root to sys.path."""
    # Path: dynamiXs_v2o0 -> modules -> lunaNMR_v1o0 (root)
    lunaNMR_root = Path(__file__).parent.parent.parent
    if str(lunaNMR_root) not in sys.path:
        sys.path.insert(0, str(lunaNMR_root))


# =============================================================================
# MAIN APPLICATION WINDOW
# =============================================================================

class DynamiXsMainWindow(QMainWindow):
    """Main application window with navigation between analysis pages."""

    def __init__(self):
        super().__init__()

        self.setWindowTitle("DynamiXs v2.0 - NMR Relaxation Analysis Suite")
        self.setMinimumSize(WINDOW_MIN_WIDTH, WINDOW_MIN_HEIGHT)
        self.resize(900, 700)

        # Get the dynamiXs directory path
        self.dynamixs_path = Path(__file__).parent

        # Current working directory for file operations
        self.current_dir = os.getcwd()

        # Create central widget and stacked layout
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.main_layout = QVBoxLayout(self.central_widget)
        self.main_layout.setContentsMargins(0, 0, 0, 0)
        self.main_layout.setSpacing(0)

        # Stacked widget for page navigation
        self.stack = QStackedWidget()
        self.main_layout.addWidget(self.stack)

        # Create pages
        self._create_pages()

        # Show main menu
        self.show_main_menu()

    def _create_pages(self):
        """Create all application pages."""
        # Main menu page
        self.main_menu_page = self._create_main_menu_page()
        self.stack.addWidget(self.main_menu_page)

        # T1/T2 Fitting page
        self.t1t2_page = T1T2FittingPage(self)
        self.stack.addWidget(self.t1t2_page)

        # T2 Methyl bi-exp Fitting page
        from methyl_t2_page import MethylT2FittingPage
        self.methyl_t2_page = MethylT2FittingPage(self)
        self.stack.addWidget(self.methyl_t2_page)

        # Spectral Density page
        self.spectral_page = SpectralDensityPage(self)
        self.stack.addWidget(self.spectral_page)

        # Integrated Analysis page
        self.integrated_page = IntegratedAnalysisPage(self)
        self.stack.addWidget(self.integrated_page)

    def _create_main_menu_page(self) -> QWidget:
        """Create the main menu page."""
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

        layout.addWidget(options_panel, stretch=1)

        # Working directory info
        dir_frame = QFrame()
        dir_frame.setProperty("class", "card")
        dir_layout = create_h_layout(SPACING_SM, (SPACING_MD, SPACING_SM, SPACING_MD, SPACING_SM))
        dir_frame.setLayout(dir_layout)

        self.dir_label = create_label(f"Working Directory: {self.current_dir}",
                                      font_size=FONT_SIZE_SMALL)
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
    # Navigation Methods
    # -------------------------------------------------------------------------

    def show_main_menu(self):
        """Show the main menu page."""
        self.stack.setCurrentWidget(self.main_menu_page)

    def show_t1t2_fitting(self):
        """Show T1/T2 fitting page."""
        self.stack.setCurrentWidget(self.t1t2_page)

    def show_methyl_t2_fitting(self):
        """Show methyl bi-exp T2 fitting page."""
        self.stack.setCurrentWidget(self.methyl_t2_page)

    def show_spectral_density(self):
        """Show spectral density analysis page."""
        self.stack.setCurrentWidget(self.spectral_page)

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
            self.dir_label.setText(f"Working Directory: {self.current_dir}")

    def get_state(self) -> dict:
        """Get complete DynamiXs state for project saving.

        Collects state from all pages including field data, session results,
        and source_series_map for series-based data.
        """
        state = {
            't1t2': {},
            'spectral_density': {},
            'model_free': {},
        }

        # T1/T2 Fitting page state
        if hasattr(self, 't1t2_page'):
            page = self.t1t2_page
            state['t1t2'] = {
                'session_results': getattr(page, 'session_results', {}),
                'fitted_experiments': getattr(page, 'fitted_experiments', []),
                'output_dir': getattr(page, 'output_dir', None),
                'field2_enabled': getattr(page, 'field2_enabled', False),
                'source_series_map': getattr(page, 'source_series_map', {}),
                # Field frequencies
                'field1_freq': page.field1_freq_spin.value() if hasattr(page, 'field1_freq_spin') else 600.0,
                'field2_freq': page.field2_freq_spin.value() if hasattr(page, 'field2_freq_spin') else 700.0,
            }

        # Spectral Density page state
        if hasattr(self, 'spectral_page'):
            page = self.spectral_page
            state['spectral_density'] = {
                'session_results': getattr(page, 'session_results', {}),
                'output_dir': getattr(page, 'output_dir', None),
                'source_series_map': getattr(page, 'source_series_map', {}),
            }

        # Model Free / Integrated Analysis page state
        if hasattr(self, 'integrated_page'):
            page = self.integrated_page
            state['model_free'] = {
                'session_results': getattr(page, 'session_results', {}),
                'output_dir': getattr(page, 'output_dir', None),
                'source_series_map': getattr(page, 'source_series_map', {}),
            }

        return state

    def get_file_refs(self) -> dict:
        """Get file references from all pages for project saving.

        Returns dict of all file paths used by DynamiXs pages.
        """
        file_refs = {
            't1t2': {},
            'spectral_density': {},
            'model_free': {},
        }

        # T1/T2 Fitting page file refs
        if hasattr(self, 't1t2_page'):
            page = self.t1t2_page
            file_refs['t1t2'] = {
                'field1_t1_file': getattr(page, 'field1_t1_file', None),
                'field1_t2_file': getattr(page, 'field1_t2_file', None),
                'field1_t1rho_file': getattr(page, 'field1_t1rho_file', None),
                'field2_t1_file': getattr(page, 'field2_t1_file', None),
                'field2_t2_file': getattr(page, 'field2_t2_file', None),
                'field2_t1rho_file': getattr(page, 'field2_t1rho_file', None),
                'peak_list_file': getattr(page, 'peak_list_file', None),
            }

        # Spectral Density page file refs
        if hasattr(self, 'spectral_page'):
            page = self.spectral_page
            file_refs['spectral_density'] = {
                'field1_t1_file': getattr(page, 'field1_t1_file', None),
                'field1_t2_file': getattr(page, 'field1_t2_file', None),
                'field1_noe_sat_file': getattr(page, 'field1_noe_sat_file', None),
                'field1_noe_unsat_file': getattr(page, 'field1_noe_unsat_file', None),
                'field2_t1_file': getattr(page, 'field2_t1_file', None),
                'field2_t2_file': getattr(page, 'field2_t2_file', None),
                'field2_noe_sat_file': getattr(page, 'field2_noe_sat_file', None),
                'field2_noe_unsat_file': getattr(page, 'field2_noe_unsat_file', None),
            }

        return file_refs

    def restore_state(self, state: dict, file_refs: dict):
        """Restore DynamiXs state from project load.

        Args:
            state: State dict from get_state()
            file_refs: File refs dict from get_file_refs()
        """
        # Restore T1/T2 Fitting page
        if hasattr(self, 't1t2_page') and 't1t2' in state:
            page = self.t1t2_page
            t1t2_state = state['t1t2']
            t1t2_refs = file_refs.get('t1t2', {})

            # Restore session results
            if 'session_results' in t1t2_state:
                page.session_results = t1t2_state['session_results']
            if 'fitted_experiments' in t1t2_state:
                page.fitted_experiments = t1t2_state['fitted_experiments']
            if 'output_dir' in t1t2_state:
                page.output_dir = t1t2_state['output_dir']
            if 'field2_enabled' in t1t2_state:
                page.field2_enabled = t1t2_state['field2_enabled']
            if 'source_series_map' in t1t2_state:
                page.source_series_map = t1t2_state['source_series_map']

            # Restore frequencies
            if hasattr(page, 'field1_freq_spin') and 'field1_freq' in t1t2_state:
                page.field1_freq_spin.setValue(t1t2_state['field1_freq'])
            if hasattr(page, 'field2_freq_spin') and 'field2_freq' in t1t2_state:
                page.field2_freq_spin.setValue(t1t2_state['field2_freq'])

            # Restore file references and update displays
            for attr in ['field1_t1_file', 'field1_t2_file', 'field1_t1rho_file',
                        'field2_t1_file', 'field2_t2_file', 'field2_t1rho_file',
                        'peak_list_file']:
                if attr in t1t2_refs and t1t2_refs[attr]:
                    setattr(page, attr, t1t2_refs[attr])

            # Update displays for field data
            page._update_field_displays()

        # Restore Spectral Density page
        if hasattr(self, 'spectral_page') and 'spectral_density' in state:
            page = self.spectral_page
            sd_state = state['spectral_density']
            sd_refs = file_refs.get('spectral_density', {})

            if 'session_results' in sd_state:
                page.session_results = sd_state['session_results']
            if 'output_dir' in sd_state:
                page.output_dir = sd_state['output_dir']
            if 'source_series_map' in sd_state:
                page.source_series_map = sd_state['source_series_map']

            # Restore file references
            for attr in ['field1_t1_file', 'field1_t2_file', 'field1_noe_sat_file',
                        'field1_noe_unsat_file', 'field2_t1_file', 'field2_t2_file',
                        'field2_noe_sat_file', 'field2_noe_unsat_file']:
                if attr in sd_refs and sd_refs[attr]:
                    setattr(page, attr, sd_refs[attr])


# =============================================================================
# BASE PAGE CLASS
# =============================================================================

def resolve_series_csv(batch_results):
    """Path to a series' series_analysis_tidy.csv for the Kd input.

    A series opened from a project reads the copy bundled beside its
    batch_results.json, never its recorded run folder, which may hold a later run; a
    series run in this session reads its own run's. If none exists on disk, returns
    "" for a series opened from a project, else the recorded csv_path for display.
    """
    _ensure_lunaNMR_path()
    from lunaNMR.utils.project_manager import series_csv_path
    path = series_csv_path(batch_results, "series_analysis_tidy.csv")
    if path is not None:
        return str(path)
    if getattr(batch_results, 'bundled_folder', None):
        return ""
    return (getattr(batch_results, 'metadata', None) or {}).get('csv_path') or ""


class BasePage(QWidget):
    """Base class for analysis pages with common header and back button."""

    def __init__(self, main_window: DynamiXsMainWindow, title: str):
        super().__init__()
        self.main_window = main_window
        self.title = title

        self._setup_base_layout()

    def _setup_base_layout(self):
        """Setup the base layout with header."""
        self.page_layout = create_v_layout(SPACING_MD, (SPACING_XL, SPACING_XL, SPACING_XL, SPACING_XL))
        self.setLayout(self.page_layout)

        # Header
        header = QFrame()
        header_layout = create_h_layout(SPACING_SM)
        header.setLayout(header_layout)

        title_label = create_header_label(self.title, level=1)
        header_layout.addWidget(title_label)
        header_layout.addStretch()

        back_btn = create_secondary_button(
            "← Back to Main",
            clicked=self.main_window.show_main_menu,
            width=140
        )
        header_layout.addWidget(back_btn)

        self.page_layout.addWidget(header)

        # Content area (subclasses add to this)
        self.content_frame = QFrame()
        self.content_frame.setProperty("class", "panel")
        self.content_layout = create_v_layout(SPACING_SM, (SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD))
        self.content_frame.setLayout(self.content_layout)
        self.page_layout.addWidget(self.content_frame, stretch=1)

    # ---------- the LunaNMR project's series ----------

    def _analysis_running(self, log) -> bool:
        """True while this page's analysis thread still runs, said through `log`.

        A second start would replace `self.worker`, the only reference to the running
        QThread; Python then collects it mid-run and Qt aborts the application.
        """
        worker = getattr(self, 'worker', None)
        if worker is not None and worker.isRunning():
            log("An analysis is already running on this page: wait for it to finish.")
            return True
        return False

    def _lunanmr_main_window(self):
        """The real LunaNMR main window. In the app the page's main_window is the
        hosting dialog, which holds the main window as .main_window; in isolation the
        main window is passed directly."""
        mw = self.main_window
        return getattr(mw, 'main_window', None) or mw

    def _get_available_series(self, include=None) -> List[Dict[str, Any]]:
        """The project's series as {'name', 'csv_path'} dicts, each tidy CSV resolved by
        resolve_series_csv. `include(batch_results)` keeps only the series it accepts;
        by default, the filter the page's series panel was built with."""
        include = include or getattr(self, '_series_include', None)
        main = self._lunanmr_main_window()
        saved_series = getattr(main, 'saved_series', {}) or {}
        return [{'name': name, 'csv_path': resolve_series_csv(br)}
                for name, br in saved_series.items() if include is None or include(br)]

    def _series_intensity_matrix(self, series_name):
        """A series' matrix of the intensity the page fits on, or None: for a series
        opened from a project, the copy bundled beside it; for a series run in this
        session, the one beside its outputs."""
        _ensure_lunaNMR_path()
        from lunaNMR.utils.project_manager import series_csv_path
        main = self._lunanmr_main_window()
        batch = (getattr(main, 'saved_series', {}) or {}).get(series_name)
        path = series_csv_path(batch, self._observable_matrix_name())
        return str(path) if path is not None else None

    # ---------- what is fitted: height or volume ----------

    def _make_observable_combo(self, label="Fit on:", redrop=True):
        """The row choosing the intensity the page fits, Height (the default, as the CLI
        reads peak_intensity_matrix.csv) or Volume. With `redrop`, a series already
        dropped takes its other matrix when the choice changes."""
        self.observable_combo = create_combo_box([o[0] for o in OBSERVABLES], width=110)
        if redrop:
            self.observable_combo.currentTextChanged.connect(self._redrop_series)
        return self._make_field_row(label, self.observable_combo, width=110)

    def _observable(self) -> str:
        """The series column fitted: 'height' or 'volume'."""
        combo = getattr(self, 'observable_combo', None)
        text = combo.currentText() if combo is not None else OBSERVABLES[0][0]
        return next(o[1] for o in OBSERVABLES if o[0] == text)

    def _observable_matrix_name(self) -> str:
        """The series matrix holding the intensity the page fits."""
        return next(o[2] for o in OBSERVABLES if o[1] == self._observable())

    def _set_observable(self, column):
        """Select the observable whose series column is `column`, if it is one."""
        label = next((o[0] for o in OBSERVABLES if o[1] == column), None)
        if label and getattr(self, 'observable_combo', None) is not None:
            self.observable_combo.setCurrentText(label)

    def _redrop_series(self, *_):
        """Give every slot a series was dropped on that series' matrix of the chosen
        intensity."""
        for field_name, series_name in list(getattr(self, 'source_series_map', {}).items()):
            self._on_series_dropped(field_name, series_name, "")

    # ---------- delay units ----------

    def _make_units_combo(self):
        """The page's one delay-units selector (ms, s, us; ms by default, the CLI's
        --time-units default): the units the series' delay headers are written in."""
        self.units_combo = create_combo_box(list(DELAY_UNITS), width=70)
        self.units_combo.setToolTip(
            "Units of the series' delay values. A dropped series' series_metadata.json "
            "sets it; units contradicting that record, or headers that state their own, "
            "are refused.")
        return self.units_combo

    def _take_recorded_units(self, path, log):
        """Set the selector to the delay units recorded in the series_metadata.json beside
        `path`, and say so through `log`."""
        sidecar = series_sidecar(path) if path else None
        recorded = (sidecar or {}).get('value_units')
        if recorded in DELAY_UNITS:
            self.units_combo.setCurrentText(recorded)
            log(f"Delay units: {recorded}, as recorded in {sidecar['path']}")

    def _checked_units(self, paths) -> Optional[str]:
        """The selected delay units, or None when they contradict the series_metadata.json
        beside one of `paths` or that file's unit-stating delay headers: the refusal is
        shown, as the CLI's is."""
        units = self.units_combo.currentText() or TIME_UNITS_DEFAULT
        try:
            for path in paths:
                if path:
                    resolve_time_units(units, series_sidecar(path), "Delay units", path)
        except ValueError as refusal:
            show_error(self, "Delay units", str(refusal))
            return None
        return units

    def _add_series_panel(self, layout, hint, empty_text, include=None, max_height=120):
        """Add the draggable list of the project's series, with its hint, an empty-state
        line and a refresh button, to `layout`. `include` filters as _get_available_series."""
        self._series_include = include
        hint_label = create_secondary_label(hint)
        hint_label.setWordWrap(True)
        layout.addWidget(hint_label)
        self.no_series_label = create_secondary_label(empty_text)
        self.no_series_label.setStyleSheet(f"color: {SECONDARY_TEXT}; font-style: italic; padding: 8px;")
        self.no_series_label.setWordWrap(True)
        layout.addWidget(self.no_series_label)
        self.series_list_widget = DraggableSeriesList()
        self.series_list_widget.setMaximumHeight(max_height)
        self.series_list_widget.setVisible(False)  # Hidden until populated
        layout.addWidget(self.series_list_widget)
        refresh_row = QFrame()
        refresh_layout = create_h_layout(SPACING_SM)
        refresh_row.setLayout(refresh_layout)
        refresh_layout.addWidget(create_secondary_button(
            "Refresh Series List", clicked=self._populate_series_list, width=150))
        refresh_layout.addStretch()
        layout.addWidget(refresh_row)
        QTimer.singleShot(100, self._populate_series_list)

    def _populate_series_list(self):
        """Fill the draggable series list; show the empty-state line when there are none."""
        if not hasattr(self, 'series_list_widget'):
            return
        self.series_list_widget.clear()
        series = self._get_available_series()
        self.no_series_label.setVisible(not series)
        self.series_list_widget.setVisible(bool(series))
        for s in series:
            item = QListWidgetItem(s['name'])
            item.setData(Qt.UserRole, s['name'])
            item.setData(Qt.UserRole + 1, s['csv_path'])
            if s['csv_path']:
                item.setToolTip(f"CSV: {s['csv_path']}")
            self.series_list_widget.addItem(item)

    # ---------- layout helpers shared by the pages ----------

    def _make_section(self, title: str) -> QFrame:
        """A titled card the page fills through its layout()."""
        frame = QFrame()
        frame.setProperty("class", "card")
        frame.setStyleSheet(
            f"QFrame[class='card'] {{ background-color: {FRAME_BG_COLOR}; "
            f"border-radius: 8px; }}"
        )
        body = create_v_layout(SPACING_XS, (SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM))
        frame.setLayout(body)
        body.addWidget(create_header_label(title, level=3))
        return frame

    def _make_field_row(self, label_text: str, widget, width: int = 180):
        """A label and its control on one row, the control at most `width` wide."""
        row = QFrame()
        row_layout = create_h_layout(SPACING_SM)
        row.setLayout(row_layout)
        row_layout.addWidget(create_label(label_text), stretch=1)
        widget.setMaximumWidth(width)
        row_layout.addWidget(widget)
        return row

    def _make_float_spin(self, attr_name, lo, hi, default, decimals=2, step=1.0):
        """A QDoubleSpinBox, kept on the page as `attr_name`."""
        spin = QDoubleSpinBox()
        spin.setRange(lo, hi)
        spin.setDecimals(decimals)
        spin.setSingleStep(step)
        spin.setValue(default)
        setattr(self, attr_name, spin)
        return spin

    def _log(self, msg: str):
        """Append to the page's log."""
        self.log_text.appendPlainText(str(msg))

    def _page_log(self, msg: str):
        """Append to whichever log the page shows its runs in."""
        for name in ('log_text', 'results_text', 'progress_text'):
            widget = getattr(self, name, None)
            if widget is not None:
                widget.appendPlainText(str(msg))
                return

    def _refuses_tidy_table(self, path) -> bool:
        """Whether `path` is a series' tidy table, which the fitters cannot read as a
        relaxation matrix; the refusal, naming the matrix to use, goes to the page's log."""
        refusal = tidy_table_refusal(path, self._observable_matrix_name())
        if refusal:
            self._page_log(refusal)
        return bool(refusal)

    # ---------- saving the run to the project ----------

    def _make_save_run_button(self):
        """The "Save run to project" button: the hosting dialog names the page's run and
        keeps it in the project session. Disabled where no project hosts the page."""
        host = self.main_window
        self.save_run_btn = create_secondary_button(
            "Save run to project", width=170,
            clicked=lambda: host.save_run_to_project(self))
        self.save_run_btn.setEnabled(hasattr(host, 'save_run_to_project'))
        return self.save_run_btn


# =============================================================================
# DRAG-DROP WIDGETS FOR SERIES PICKER
# =============================================================================

class DraggableSeriesList(QListWidget):
    """List widget with drag support for series items."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setDragEnabled(True)
        self.setDragDropMode(QAbstractItemView.DragOnly)
        self.setDefaultDropAction(Qt.CopyAction)
        self.setSelectionMode(QAbstractItemView.SingleSelection)
        self.setStyleSheet(f"""
            QListWidget {{
                background-color: {FRAME_BG_COLOR};
                border: 1px solid {PANEL_BG_COLOR};
                border-radius: 4px;
                padding: 4px;
            }}
            QListWidget::item {{
                padding: 6px 8px;
                border-radius: 3px;
            }}
            QListWidget::item:selected {{
                background-color: {PANEL_BG_COLOR};
            }}
            QListWidget::item:hover {{
                background-color: {PANEL_BG_COLOR};
            }}
        """)

    def mimeData(self, items):
        """Build the drag payload the drop target reads. Overriding this (rather
        than startDrag) lets the base view own the drag lifecycle, so an item
        stays draggable after being dropped once."""
        mime_data = QMimeData()
        if items:
            item = items[0]
            series_name = item.data(Qt.UserRole)
            csv_path = item.data(Qt.UserRole + 1) or ""
            mime_data.setText(f"series:{series_name}:{csv_path}")
        return mime_data


class DropTargetLabel(QLabel):
    """Label that accepts dropped series items."""

    series_dropped = Signal(str, str, str)  # field_name, series_name, csv_path

    def __init__(self, field_name: str, text: str = "No file selected", parent=None):
        super().__init__(text, parent)
        self.field_name = field_name
        self.setAcceptDrops(True)
        self._update_style(False)

    def _update_style(self, is_hover: bool):
        """Update style based on hover state."""
        if is_hover:
            self.setStyleSheet(f"""
                color: {PRIMARY_TEXT};
                background-color: {SUCCESS_GREEN}40;
                padding: 4px 8px;
                border-radius: 4px;
                border: 2px dashed {SUCCESS_GREEN};
            """)
        else:
            self.setStyleSheet(f"""
                color: {SECONDARY_TEXT};
                background-color: {FRAME_BG_COLOR};
                padding: 4px 8px;
                border-radius: 4px;
            """)

    def dragEnterEvent(self, event: QDragEnterEvent):
        """Accept drag if it contains series data."""
        if event.mimeData().hasText():
            text = event.mimeData().text()
            if text.startswith("series:"):
                event.acceptProposedAction()
                self._update_style(True)
                return
        event.ignore()

    def dragLeaveEvent(self, event):
        """Reset style when drag leaves."""
        self._update_style(False)

    def dropEvent(self, event: QDropEvent):
        """Handle dropped series."""
        self._update_style(False)
        text = event.mimeData().text()
        if text.startswith("series:"):
            parts = text.split(":", 2)
            if len(parts) >= 2:
                series_name = parts[1]
                csv_path = parts[2] if len(parts) > 2 else ""
                self.series_dropped.emit(self.field_name, series_name, csv_path)
                event.acceptProposedAction()
                return
        event.ignore()


# =============================================================================
# T1/T2 FITTING PAGE
# =============================================================================

# Residue names are identifiers, read as text: a list named 3.0, 4.0, ... otherwise loads
# as floats, and the T1rho -> T2 conversion, matching names as strings, fails on them.
RESIDUE_NAME_DTYPES = {'Residue': str, 'residue': str, 'Assignment': str, 'res': str}


class T1T2FittingPage(BasePage):
    """T1/T2 exponential decay fitting analysis page with dual-field support."""

    def __init__(self, main_window: DynamiXsMainWindow):
        super().__init__(main_window, "T1/T2 Fitting Analysis")
        self.worker = None
        self._init_session_state()
        self._init_file_vars()
        self._setup_content()

    def _init_session_state(self):
        """Initialize session state for persistent results across fits."""
        # Session results store fit outputs without overwriting
        # Each entry: {output_dir, json_folder, results_file, n_fitted, frequency_mhz}
        self.session_results = {
            'field1': {
                'frequency_mhz': 600.0,
                't1': None,
                't2': None,
                't1rho': None
            },
            'field2': None  # Enabled when user activates dual-field
        }
        # Track which experiments have been successfully fitted
        self.fitted_experiments = set()  # e.g., {"field1_t1", "field1_t2", "field2_t1"}

    def _init_file_vars(self):
        """Initialize file path variables for dual-field support."""
        # Legacy single-file input (for backward compatibility)
        self.input_file = None

        # Field 1 files
        self.field1_t1_file = None
        self.field1_t2_file = None
        self.field1_t1rho_file = None

        # Field 2 files
        self.field2_t1_file = None
        self.field2_t2_file = None
        self.field2_t1rho_file = None

        # Output settings
        self.output_dir = None

        # Peak list for the T1rho -> T2 conversion
        self.peak_list_file = None

        # Field 2 enabled state
        self.field2_enabled = False

    def _create_file_row(self, parent_layout, label_text: str, field_name: str, enabled: bool = True) -> DropTargetLabel:
        """Create a file selection row with drop-enabled display for data input."""
        row = QFrame()
        row_layout = create_h_layout(SPACING_SM)
        row.setLayout(row_layout)

        label = create_label(label_text)
        label.setFixedWidth(100)
        row_layout.addWidget(label)

        # Use DropTargetLabel for drag-and-drop from series picker
        display = DropTargetLabel(field_name, "No file selected")
        display.setFont(get_font(FONT_SIZE_SMALL))
        display.setEnabled(enabled)
        display.series_dropped.connect(self._on_series_dropped)
        row_layout.addWidget(display, stretch=1)

        browse_btn = create_secondary_button("Browse", width=80, enabled=enabled)
        browse_btn.clicked.connect(lambda: self._browse_data_file(field_name))
        row_layout.addWidget(browse_btn)

        # Store button reference for enable/disable
        setattr(self, f"{field_name}_browse_btn", browse_btn)

        parent_layout.addWidget(row)
        return display

    def _browse_data_file(self, field_name: str):
        """Browse for a data file and store path."""
        filename = open_file_dialog(
            self,
            title=f"Select {field_name.replace('_', ' ').title()} File",
            directory=self.main_window.current_dir,
            filter="CSV files (*.csv);;All files (*.*)"
        )
        if filename and not self._refuses_tidy_table(filename):
            self._assign_data_file(field_name, filename)

    def _assign_data_file(self, field_name: str, path: str):
        """Make `path` the slot's data. A series dropped there before no longer is: a
        fit reads a dropped series from memory first."""
        setattr(self, f"{field_name}_file", path)
        getattr(self, 'source_series_map', {}).pop(field_name, None)
        display = getattr(self, f"{field_name}_display", None)
        if display:
            display.setText(os.path.basename(path))
            display.setToolTip(path)

    # =========================================================================
    # SERIES PICKER METHODS
    # =========================================================================

    def _get_fitting_dataframe_from_series(self, series_name: str):
        """Get fitting DataFrame from BatchResults in memory.

        Accesses lunaNMR's saved_series to get BatchResults and calls
        to_fitting_dataframe() to generate the pivot format needed for T1/T2 fitting.

        Args:
            series_name: Name of the series to retrieve

        Returns:
            DataFrame in fitting format (Assignment + delay columns), or None if not found
        """
        # Access lunaNMR main window through DynamiXsDialog
        lunaNMR_main = None
        if hasattr(self.main_window, 'main_window'):
            lunaNMR_main = self.main_window.main_window

        if not lunaNMR_main:
            return None

        # Get saved_series from lunaNMR main window
        saved_series = getattr(lunaNMR_main, 'saved_series', {}) or {}

        if series_name not in saved_series:
            return None

        batch_results = saved_series[series_name]

        # Check if BatchResults has the to_fitting_dataframe method
        if not hasattr(batch_results, 'to_fitting_dataframe'):
            return None

        try:
            fitting_df = batch_results.to_fitting_dataframe(value_column=self._observable())
            if fitting_df.empty:
                return None
            return fitting_df
        except Exception as e:
            print(f"Error generating fitting DataFrame from series '{series_name}': {e}")
            return None

    def _on_series_dropped(self, field_name: str, series_name: str, csv_path: str):
        """Handle a series being dropped on a data field.

        Args:
            field_name: Target field (e.g., 'field1_t1', 'field2_t1rho')
            series_name: Name of the dropped series
            csv_path: Path to series_analysis_tidy.csv
        """
        # Set the file path
        if csv_path:
            setattr(self, f"{field_name}_file", csv_path)
        else:
            # CSV path not found - warn user
            show_warning(
                self,
                "CSV Not Found",
                f"Could not find series_analysis_tidy.csv for series '{series_name}'.\n\n"
                "The series was saved but the CSV output file could not be located.\n"
                "Please use the Browse button to select the CSV file manually."
            )

        # Update the display
        display = getattr(self, f"{field_name}_display", None)
        if display:
            if csv_path:
                display.setText(f"📊 {series_name}")
                display.setToolTip(f"Series: {series_name}\nCSV: {csv_path}")
            else:
                display.setText(f"⚠️ {series_name} (CSV not found)")
                display.setToolTip(f"Series: {series_name}\nCSV: Not found - use Browse to select")

        # Store source_series metadata for later use in Inspect Peak
        if not hasattr(self, 'source_series_map'):
            self.source_series_map = {}
        self.source_series_map[field_name] = series_name

        # Log to results area
        if hasattr(self, 'results_text'):
            self.results_text.appendPlainText(
                f"Assigned series '{series_name}' to {field_name.replace('_', ' ').title()}"
            )
            self._take_recorded_units(csv_path, self.results_text.appendPlainText)

    def _update_field_displays(self):
        """Update field display widgets after state restore.

        Called after loading a project to show the assigned files/series
        in the UI displays.
        """
        source_series_map = getattr(self, 'source_series_map', {})

        # Field mapping: attribute name -> display widget attribute
        field_mappings = [
            ('field1_t1', 'field1_t1_display'),
            ('field1_t2', 'field1_t2_display'),
            ('field1_t1rho', 'field1_t1rho_display'),
            ('field2_t1', 'field2_t1_display'),
            ('field2_t2', 'field2_t2_display'),
            ('field2_t1rho', 'field2_t1rho_display'),
        ]

        for field_key, display_attr in field_mappings:
            file_attr = f"{field_key}_file"
            file_path = getattr(self, file_attr, None)
            display = getattr(self, display_attr, None)

            if display is None:
                continue

            # Check if this field was set from a series
            series_name = source_series_map.get(field_key)
            if series_name:
                display.setText(f"📊 {series_name}")
                display.setToolTip(f"Series: {series_name}\nFile: {file_path or 'In memory'}")
            elif file_path:
                # Show filename only
                display.setText(os.path.basename(file_path))
                display.setToolTip(file_path)
            else:
                display.setText("No file selected")
                display.setToolTip("")

        # Update output directory display
        if hasattr(self, 'outdir_display') and self.output_dir:
            self.outdir_display.setText(os.path.basename(self.output_dir))
            self.outdir_display.setToolTip(self.output_dir)

    @staticmethod
    def _create_carrier_spin():
        """A carrier input that starts unset: its minimum reads 'required'."""
        spin = create_double_spin_box(0.0, 0.0, 300.0, 3, " ppm", width=110)
        spin.setSpecialValueText("required")
        return spin

    def _toggle_field2(self):
        """Enable or disable Field 2 data section."""
        self.field2_enabled = not self.field2_enabled

        if self.field2_enabled:
            self.enable_field2_btn.setText("Disable")
            # Initialize field2 session structure
            if self.session_results.get('field2') is None:
                self.session_results['field2'] = {
                    'frequency_mhz': self.field2_freq_spin.value(),
                    't1': None,
                    't2': None,
                    't1rho': None
                }
        else:
            self.enable_field2_btn.setText("Enable")

        # Enable/disable Field 2 widgets
        self.field2_freq_spin.setEnabled(self.field2_enabled)
        self.field2_auto_btn.setEnabled(self.field2_enabled)
        for name in ("theta", "omega1", "carrier"):
            getattr(self, f"field2_{name}_spin").setEnabled(self.field2_enabled)

        for exp in ['t1', 't2', 't1rho']:
            display = getattr(self, f"field2_{exp}_display", None)
            btn = getattr(self, f"field2_{exp}_browse_btn", None)
            if display:
                display.setEnabled(self.field2_enabled)
            if btn:
                btn.setEnabled(self.field2_enabled)

        # Update analysis dropdown to show/hide field2 options
        self._update_analysis_options()

    def _auto_load_folder(self, field: str):
        """Auto-detect and load T1, T2, and T1rho files from a folder."""
        folder = open_directory_dialog(
            self,
            f"Select {field.title()} Data Folder",
            self.main_window.current_dir
        )
        if not folder:
            return

        # Pattern matching for relaxation data files
        patterns = {
            't1': ['*T1*.csv', '*t1*.csv', '*R1*.csv'],
            't2': ['*T2*.csv', '*t2*.csv', '*R2*.csv'],
            't1rho': ['*T1rho*.csv', '*t1rho*.csv', '*R1rho*.csv', '*T1r*.csv']
        }

        found_files = {}
        for exp, file_patterns in patterns.items():
            for pattern in file_patterns:
                matches = glob.glob(os.path.join(folder, pattern))
                if matches:
                    found_files[exp] = matches[0]
                    break

        # Update file displays
        for exp, filepath in found_files.items():
            self._assign_data_file(f"{field}_{exp}", filepath)

        if found_files:
            self.results_text.appendPlainText(
                f"Auto-loaded from {os.path.basename(folder)}: {', '.join(found_files.keys())}"
            )
        else:
            self.results_text.appendPlainText(
                f"No matching files found in {os.path.basename(folder)}"
            )

    def _update_analysis_options(self):
        """Update analysis dropdown based on field2 enabled state."""
        current_selection = self.analysis_combo.currentText()

        # Clear and rebuild options
        self.analysis_combo.clear()

        # Field 1 options (always available)
        self.analysis_combo.addItems([
            "Field 1 - T1",
            "Field 1 - T2",
            "Field 1 - T1ρ",
            "Field 1 - All",
        ])

        # Field 2 options (only when enabled)
        if self.field2_enabled:
            self.analysis_combo.addItems([
                "Field 2 - T1",
                "Field 2 - T2",
                "Field 2 - T1ρ",
                "Field 2 - All",
                "Both Fields - All",
            ])

        # Try to restore previous selection
        index = self.analysis_combo.findText(current_selection)
        if index >= 0:
            self.analysis_combo.setCurrentIndex(index)

    def _on_analysis_selection_changed(self, selection: str):
        """Update the results prefix based on analysis selection."""
        if not hasattr(self, 'prefix_edit'):
            return

        # Extract experiment type from selection
        if "T1ρ" in selection:
            prefix = "T1rho"
        elif "T2" in selection:
            prefix = "T2"
        elif "T1" in selection:
            prefix = "T1"
        elif "All" in selection:
            prefix = "all"
        else:
            prefix = "analysis"

        self.prefix_edit.setText(prefix)

    def _setup_content(self):
        """Setup the page content."""
        # Replace default content frame with scroll area (like IntegratedAnalysisPage)
        self.page_layout.removeWidget(self.content_frame)
        self.content_frame.deleteLater()

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)

        scroll_content = QWidget()
        scroll_content.setProperty("class", "panel")
        scroll_layout = create_v_layout(SPACING_MD, (SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD))
        scroll_content.setLayout(scroll_layout)

        # --- Project Series Section (drag series to data fields) ---
        series_group = CollapsibleGroupBox("Project Series", collapsible=True, initially_collapsed=False)
        series_content = series_group.content_layout

        self._add_series_panel(
            series_content, "Drag a series to assign it to a T1/T2/T1ρ data field below",
            "No project series available. Use Browse buttons or Auto-Load from Folder.")

        scroll_layout.addWidget(series_group)

        # --- Field 1 Data Section (Required) ---
        field1_group = CollapsibleGroupBox("Field 1 Data (Required)")
        field1_content = field1_group.content_layout

        # Field 1 frequency in header
        freq1_row = QFrame()
        freq1_layout = create_h_layout(SPACING_XS)
        freq1_row.setLayout(freq1_layout)
        self.field1_freq_spin = create_double_spin_box(600.0, 100.0, 1200.0, 1, " MHz", width=100)
        freq1_layout.addWidget(self.field1_freq_spin)
        field1_group.add_header_widget(freq1_row)

        # Auto-load button for Field 1
        auto_load_row1 = QFrame()
        auto_load_layout1 = create_h_layout(SPACING_SM)
        auto_load_row1.setLayout(auto_load_layout1)
        self.field1_auto_btn = create_primary_button(
            "Auto-Load from Folder",
            clicked=lambda: self._auto_load_folder("field1"),
            width=200
        )
        auto_load_layout1.addWidget(self.field1_auto_btn)
        auto_load_layout1.addWidget(create_secondary_label("Automatically detect T1, T2, and T1ρ files"))
        auto_load_layout1.addStretch()
        field1_content.addWidget(auto_load_row1)

        # File rows for Field 1
        self.field1_t1_display = self._create_file_row(field1_content, "T1 Data:", "field1_t1")
        self.field1_t2_display = self._create_file_row(field1_content, "T2 Data:", "field1_t2")
        self.field1_t1rho_display = self._create_file_row(field1_content, "T1ρ Data:", "field1_t1rho")

        scroll_layout.addWidget(field1_group)

        # --- Field 2 Data Section (Optional) ---
        field2_group = CollapsibleGroupBox("Field 2 Data (Optional - for dual-field)")
        field2_content = field2_group.content_layout

        # Field 2 frequency and enable button in header
        freq2_row = QFrame()
        freq2_layout = create_h_layout(SPACING_SM)
        freq2_row.setLayout(freq2_layout)
        self.field2_freq_spin = create_double_spin_box(700.0, 100.0, 1200.0, 1, " MHz", width=100)
        self.field2_freq_spin.setEnabled(False)
        freq2_layout.addWidget(self.field2_freq_spin)
        self.enable_field2_btn = create_secondary_button("Enable", clicked=self._toggle_field2, width=80)
        freq2_layout.addWidget(self.enable_field2_btn)
        field2_group.add_header_widget(freq2_row)

        # Auto-load button for Field 2
        auto_load_row2 = QFrame()
        auto_load_layout2 = create_h_layout(SPACING_SM)
        auto_load_row2.setLayout(auto_load_layout2)
        self.field2_auto_btn = create_primary_button(
            "Auto-Load from Folder",
            clicked=lambda: self._auto_load_folder("field2"),
            width=200,
            enabled=False
        )
        auto_load_layout2.addWidget(self.field2_auto_btn)
        auto_load_layout2.addWidget(create_secondary_label("Automatically detect T1, T2, and T1ρ files"))
        auto_load_layout2.addStretch()
        field2_content.addWidget(auto_load_row2)

        # File rows for Field 2 (disabled by default)
        self.field2_t1_display = self._create_file_row(field2_content, "T1 Data:", "field2_t1", enabled=False)
        self.field2_t2_display = self._create_file_row(field2_content, "T2 Data:", "field2_t2", enabled=False)
        self.field2_t1rho_display = self._create_file_row(field2_content, "T1ρ Data:", "field2_t1rho", enabled=False)

        scroll_layout.addWidget(field2_group)

        # --- T1ρ → T2 Conversion Section ---
        self.t1rho_conversion_group = CollapsibleGroupBox("T1ρ → T2 Conversion", collapsible=True, initially_collapsed=True)
        t1rho_conv_content = self.t1rho_conversion_group.content_layout

        # Requirements info
        req_label = create_secondary_label("Requires: T1 fit results + T1ρ fit results")
        t1rho_conv_content.addWidget(req_label)

        # Status indicators for required fits
        status_row = QFrame()
        status_layout = create_h_layout(SPACING_MD)
        status_row.setLayout(status_layout)
        status_layout.addWidget(create_label("Available:"))
        self.conv_t1_status = create_label("T1: ○", font_size=FONT_SIZE_SMALL)
        self.conv_t1_status.setStyleSheet(f"color: {SECONDARY_TEXT};")
        status_layout.addWidget(self.conv_t1_status)
        self.conv_t1rho_status = create_label("T1ρ: ○", font_size=FONT_SIZE_SMALL)
        self.conv_t1rho_status.setStyleSheet(f"color: {SECONDARY_TEXT};")
        status_layout.addWidget(self.conv_t1rho_status)
        status_layout.addStretch()
        t1rho_conv_content.addWidget(status_row)

        # Acquisition parameters, PER FIELD: each field is its own acquisition, and none
        # of these can be read off the spectra or the filenames. The carrier has no
        # default - a guessed carrier tilts every residue by the wrong angle.
        acq_grid = QGridLayout()
        acq_grid.setHorizontalSpacing(SPACING_MD)
        acq_grid.addWidget(create_label("Field 1"), 0, 1)
        acq_grid.addWidget(create_label("Field 2"), 0, 2)
        rows = (("Tilt angle θ (cnst28):", "theta", lambda: create_double_spin_box(
                    45.0, 1.0, 89.0, 1, " °", width=110)),
                ("Spin-lock ω₁ (cnst27):", "omega1", lambda: create_double_spin_box(
                    2000.0, 100.0, 10000.0, 3, " Hz", width=110)),
                ("¹⁵N carrier (O3):", "carrier", self._create_carrier_spin))
        for row, (label, name, make) in enumerate(rows, start=1):
            acq_grid.addWidget(create_label(label), row, 0)
            for column, field in enumerate(("field1", "field2"), start=1):
                spin = make()
                setattr(self, f"{field}_{name}_spin", spin)
                acq_grid.addWidget(spin, row, column)
        acq_grid.setColumnStretch(3, 1)
        acq_frame = QFrame()
        acq_frame.setLayout(acq_grid)
        t1rho_conv_content.addWidget(acq_frame)
        for name in ("theta", "omega1", "carrier"):
            getattr(self, f"field2_{name}_spin").setEnabled(False)

        # Peak list file for T1rho calculation
        peak_list_row = QFrame()
        peak_list_layout = create_h_layout(SPACING_SM)
        peak_list_row.setLayout(peak_list_layout)
        peak_list_layout.addWidget(create_label("Peak List:"))
        self.peak_list_display = create_label("No file selected", font_size=FONT_SIZE_SMALL)
        self.peak_list_display.setStyleSheet(f"color: {SECONDARY_TEXT}; background-color: {FRAME_BG_COLOR}; padding: 4px 8px; border-radius: 4px;")
        peak_list_layout.addWidget(self.peak_list_display, stretch=1)
        peak_list_browse_btn = create_secondary_button("Browse", clicked=self._browse_peak_list_file, width=80)
        peak_list_layout.addWidget(peak_list_browse_btn)
        t1rho_conv_content.addWidget(peak_list_row)

        # Calculate button
        calc_btn_row = QFrame()
        calc_btn_layout = create_h_layout(SPACING_SM)
        calc_btn_row.setLayout(calc_btn_layout)
        self.calc_t2_btn = create_primary_button(
            "Calculate T2 from T1ρ",
            clicked=self._calculate_t2_from_t1rho,
            width=200
        )
        calc_btn_layout.addWidget(self.calc_t2_btn)
        calc_btn_layout.addStretch()
        t1rho_conv_content.addWidget(calc_btn_row)

        scroll_layout.addWidget(self.t1rho_conversion_group)

        # --- Analysis Selection ---
        analysis_frame = QFrame()
        analysis_layout = create_h_layout(SPACING_SM)
        analysis_frame.setLayout(analysis_layout)

        analysis_layout.addWidget(create_label("Analysis to Run:"))
        self.analysis_combo = create_combo_box([
            "Field 1 - T1",
            "Field 1 - T2",
            "Field 1 - T1ρ",
            "Field 1 - All",
        ], width=200)
        self.analysis_combo.currentTextChanged.connect(self._on_analysis_selection_changed)
        analysis_layout.addWidget(self.analysis_combo)
        analysis_layout.addStretch()

        scroll_layout.addWidget(analysis_frame)

        # --- Analysis Status ---
        status_group = StyledGroupBox("Analysis Status")
        status_layout = create_h_layout(SPACING_MD, (SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM))
        status_group.setLayout(status_layout)

        # Field 1 status
        f1_status = QFrame()
        f1_layout = create_h_layout(SPACING_XS)
        f1_status.setLayout(f1_layout)
        f1_layout.addWidget(create_label("Field 1:"))
        self.field1_t1_indicator = create_label("T1 ○")
        self.field1_t2_indicator = create_label("T2 ○")
        self.field1_t1rho_indicator = create_label("T1ρ ○")
        f1_layout.addWidget(self.field1_t1_indicator)
        f1_layout.addWidget(self.field1_t2_indicator)
        f1_layout.addWidget(self.field1_t1rho_indicator)
        status_layout.addWidget(f1_status)

        # Field 2 status
        f2_status = QFrame()
        f2_layout = create_h_layout(SPACING_XS)
        f2_status.setLayout(f2_layout)
        f2_layout.addWidget(create_label("Field 2:"))
        self.field2_t1_indicator = create_label("T1 ○")
        self.field2_t2_indicator = create_label("T2 ○")
        self.field2_t1rho_indicator = create_label("T1ρ ○")
        f2_layout.addWidget(self.field2_t1_indicator)
        f2_layout.addWidget(self.field2_t2_indicator)
        f2_layout.addWidget(self.field2_t1rho_indicator)
        status_layout.addWidget(f2_status)

        status_layout.addStretch()
        scroll_layout.addWidget(status_group)

        # Results prefix
        prefix_frame = QFrame()
        prefix_layout = create_h_layout(SPACING_SM)
        prefix_frame.setLayout(prefix_layout)

        prefix_layout.addWidget(create_label("Results Prefix:"))
        self.prefix_edit = create_line_edit("T1", width=300)
        prefix_layout.addWidget(self.prefix_edit)
        prefix_layout.addStretch()

        scroll_layout.addWidget(prefix_frame)

        # Output Configuration
        output_group = StyledGroupBox("Output Configuration")
        output_layout = create_v_layout(SPACING_XS, (SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD))
        output_group.setLayout(output_layout)

        # Output directory row
        outdir_row = QFrame()
        outdir_layout = create_h_layout(SPACING_SM)
        outdir_row.setLayout(outdir_layout)
        outdir_layout.addWidget(create_label("Output Directory:"))
        self.outdir_display = create_label("(same as input file)", font_size=FONT_SIZE_SMALL)
        self.outdir_display.setStyleSheet(f"color: {SECONDARY_TEXT}; background-color: {FRAME_BG_COLOR}; padding: 4px 8px; border-radius: 4px;")
        outdir_layout.addWidget(self.outdir_display, stretch=1)
        outdir_btn = create_secondary_button("Browse", clicked=self._browse_output_dir, width=80)
        outdir_layout.addWidget(outdir_btn)
        output_layout.addWidget(outdir_row)

        # JSON folder row
        json_row = QFrame()
        json_layout = create_h_layout(SPACING_SM)
        json_row.setLayout(json_layout)
        json_layout.addWidget(create_label("JSON Data Folder:"))
        self.json_display = create_label("json", font_size=FONT_SIZE_SMALL)
        self.json_display.setStyleSheet(f"color: {SECONDARY_TEXT}; background-color: {FRAME_BG_COLOR}; padding: 4px 8px; border-radius: 4px;")
        json_layout.addWidget(self.json_display, stretch=1)
        json_btn = create_secondary_button("Browse", clicked=self._browse_json_folder, width=80)
        json_layout.addWidget(json_btn)
        output_layout.addWidget(json_row)

        scroll_layout.addWidget(output_group)

        # Advanced parameters
        advanced_group = StyledGroupBox("Advanced Parameters")
        advanced_layout = create_v_layout(SPACING_XS, (SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD))
        advanced_group.setLayout(advanced_layout)

        # Initial Amplitude
        amp_row = QFrame()
        amp_layout = create_h_layout(SPACING_SM)
        amp_row.setLayout(amp_layout)
        amp_layout.addWidget(create_label("Initial Amplitude:"))
        self.amp_spin = create_double_spin_box(5.0, 0.001, 1000.0, 3, width=100)
        amp_layout.addWidget(self.amp_spin)
        amp_layout.addStretch()
        advanced_layout.addWidget(amp_row)

        # Initial Time Constant
        time_row = QFrame()
        time_layout = create_h_layout(SPACING_SM)
        time_row.setLayout(time_layout)
        time_layout.addWidget(create_label("Initial Time Constant:"))
        self.time_spin = create_double_spin_box(100.0, 1.0, 10000.0, 1, " ms", width=100)
        time_layout.addWidget(self.time_spin)
        time_layout.addWidget(create_label("Delay units:"))
        time_layout.addWidget(self._make_units_combo())
        time_layout.addStretch()
        advanced_layout.addWidget(time_row)
        advanced_layout.addWidget(self._make_observable_combo(redrop=False))

        # Error estimation method
        error_row = QFrame()
        error_layout = create_h_layout(SPACING_SM)
        error_row.setLayout(error_layout)
        error_layout.addWidget(create_label("Error Method:"))
        self.error_method_combo = create_combo_box(
            ["Analytical (fast)", "Bootstrap"],
            width=150
        )
        self.error_method_combo.currentIndexChanged.connect(self._on_error_method_changed)
        error_layout.addWidget(self.error_method_combo)
        error_layout.addStretch()
        advanced_layout.addWidget(error_row)

        # Bootstrap iterations (only visible when Bootstrap selected)
        self.boot_row = QFrame()
        boot_layout = create_h_layout(SPACING_SM)
        self.boot_row.setLayout(boot_layout)
        boot_layout.addWidget(create_label("Bootstrap Iterations:"))
        self.boot_spin = create_spin_box(1000, 100, 10000, width=100)
        boot_layout.addWidget(self.boot_spin)
        boot_layout.addStretch()
        self.boot_row.hide()  # Hidden by default (analytical is default)
        advanced_layout.addWidget(self.boot_row)

        scroll_layout.addWidget(advanced_group)

        # Start button
        start_btn = create_primary_button("Start Analysis", clicked=self._start_analysis, width=200)
        scroll_layout.addWidget(start_btn, alignment=Qt.AlignCenter)

        # Results area
        results_group = StyledGroupBox("Analysis Results")
        results_layout = create_v_layout(SPACING_SM, (SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD))
        results_group.setLayout(results_layout)

        self.results_text = create_plain_text_edit(read_only=True)
        self.results_text.setMinimumHeight(150)
        results_layout.addWidget(self.results_text)

        # Progress bar for fitting progress
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(True)
        self.progress_bar.setFormat("%p% - Fitting residues...")
        self.progress_bar.hide()  # Hidden until analysis starts
        results_layout.addWidget(self.progress_bar)

        # Button row for results actions
        btn_row = QFrame()
        btn_row_layout = create_h_layout(SPACING_SM)
        btn_row.setLayout(btn_row_layout)

        view_results_btn = create_secondary_button("View Results", clicked=self._display_results, width=120)
        btn_row_layout.addWidget(view_results_btn)

        view_fits_btn = create_secondary_button("View Fits", clicked=self._open_fit_viewer, width=100)
        btn_row_layout.addWidget(view_fits_btn)
        btn_row_layout.addWidget(self._make_save_run_button())

        btn_row_layout.addStretch()
        results_layout.addWidget(btn_row)

        scroll_layout.addWidget(results_group)

        # Finalize scroll area
        scroll.setWidget(scroll_content)
        self.page_layout.addWidget(scroll, stretch=1)

    def _on_error_method_changed(self, index: int):
        """Show/hide bootstrap iterations based on error method selection."""
        if index == 1:  # Bootstrap selected
            self.boot_row.show()
        else:  # Analytical selected
            self.boot_row.hide()

    def _browse_input_file(self):
        """Browse for input CSV file."""
        filename = open_file_dialog(
            self,
            title="Select Input CSV File",
            directory=self.main_window.current_dir,
            filter="CSV files (*.csv);;All files (*.*)"
        )
        if filename:
            self.input_file = filename
            self.file_display.setText(os.path.basename(filename))

    def _start_analysis(self):
        """Start the T1/T2 analysis based on selected analysis type."""
        if self._analysis_running(self.results_text.appendPlainText):
            return
        # Require output directory to be set before fitting
        if not self.output_dir:
            show_warning(
                self,
                "Output Directory Required",
                "Please select an output directory before running analysis.\n"
                "All fit results and JSON files will be saved there."
            )

            # Open directory selection dialog
            initial_dir = self.main_window.current_dir
            folder = open_directory_dialog(self, "Select Output Directory", initial_dir)
            if folder:
                self.output_dir = folder
                self.outdir_display.setText(os.path.basename(folder))
            else:
                # User cancelled - abort
                return

        self.results_text.clear()

        # Parse selected analysis from dropdown
        selection = self.analysis_combo.currentText()
        self.results_text.appendPlainText(f"Starting analysis: {selection}")

        # Parse field and experiment type from selection
        if selection.startswith("Field 1"):
            field_name = "field1"
            field_freq = self.field1_freq_spin.value()
        else:
            field_name = "field2"
            field_freq = self.field2_freq_spin.value()

        # Handle "All" option
        if "All" in selection:
            self._run_all_analyses(field_name)
            return

        # Determine experiment type and get data source
        if "T1ρ" in selection:
            exp_type = "T1rho"
            data_field_key = f"{field_name}_t1rho"
        elif "T2" in selection:
            exp_type = "T2"
            data_field_key = f"{field_name}_t2"
        else:  # T1
            exp_type = "T1"
            data_field_key = f"{field_name}_t1"

        # Try to get data from memory (series) first, then fall back to file path
        input_df = None
        input_file = None
        series_name = None

        # Check if data was set from a series (memory-based)
        source_series_map = getattr(self, 'source_series_map', {})
        if data_field_key in source_series_map:
            series_name = source_series_map[data_field_key]
            # Try to get BatchResults from lunaNMR's saved_series
            input_df = self._get_fitting_dataframe_from_series(series_name)
            if input_df is not None:
                self.results_text.appendPlainText(f"Using data from series '{series_name}' in memory")

        # Fall back to file path if no memory data
        if input_df is None:
            input_file = getattr(self, f"{data_field_key}_file", None)

        if input_df is None and not input_file:
            show_error(self, "Error", f"Please select a {exp_type} data file for {field_name.replace('field', 'Field ')}.")
            return

        # The slot's file - a dropped series' tidy CSV - sits beside the series' record
        time_units = self._checked_units([getattr(self, f"{data_field_key}_file", None)])
        if time_units is None:
            return

        # Use the required output directory
        output_dir = self.output_dir

        # Build JSON folder path
        json_folder_name = self.json_display.text() if hasattr(self, 'json_display') else "json"
        json_folder = os.path.join(output_dir, json_folder_name)

        # Determine error method from dropdown
        error_method = "bootstrap" if self.error_method_combo.currentIndex() == 1 else "analytical"

        # Every file of the run starts with its prefix, the CLI's <prefix>_<field> naming:
        # runs sharing an output folder keep their own, and a project bundles each
        # analysis with only its own.
        run_prefix = self._run_prefix()

        # Create parameters (input_df takes priority over input_file)
        params = T1T2FittingParams(
            input_file=input_file or "",
            output_dir=output_dir,
            experiment_type=exp_type,
            results_prefix=f"{run_prefix}_{field_name}_{exp_type}",
            multicore=True,
            initial_amplitude=self.amp_spin.value(),
            initial_decay=self.time_spin.value() * MS_TO_UNITS[time_units],
            time_units=time_units,
            n_bootstrap=self.boot_spin.value(),
            error_method=error_method,
            json_folder=json_folder,
            field_name=field_name,
            run_prefix=run_prefix,
            field_freq=field_freq,
            input_df=input_df,
            series_name=series_name or ""
        )

        # Show progress bar and reset
        self.progress_bar.setValue(0)
        self.progress_bar.setFormat(f"%p% - Fitting {exp_type} for {field_name}...")
        self.progress_bar.show()

        # Create and start worker
        self.worker = T1T2FittingWorker(params)
        self.worker.progress.connect(self._on_progress)
        self.worker.progress_value.connect(self._on_progress_value)
        self.worker.finished.connect(self._on_finished)
        self.worker.error.connect(self._on_error)
        self.worker.start()

    def _run_prefix(self) -> str:
        """The Results Prefix every file of the next run starts with."""
        return self.prefix_edit.text().strip() or "analysis"

    def _session_entries(self):
        """Every stored fit of this page's session, across fields and experiments."""
        for field_data in self.session_results.values():
            for exp_type in ('t1', 't2', 't1rho'):
                if field_data and field_data.get(exp_type):
                    yield field_data[exp_type]

    def _file_prefixes(self) -> list:
        """What the files of each run stored in the session start with."""
        return sorted({e['file_prefix'] for e in self._session_entries()
                       if e.get('file_prefix')})

    def _run_all_analyses(self, field_name: str):
        """Run all analyses (T1, T2) for the specified field."""
        self.results_text.appendPlainText(f"Running all analyses for {field_name}...")

        # Get available files for this field
        t1_file = getattr(self, f"{field_name}_t1_file", None)
        t2_file = getattr(self, f"{field_name}_t2_file", None)

        analyses_to_run = []
        if t1_file:
            analyses_to_run.append(("T1", t1_file))
        if t2_file:
            analyses_to_run.append(("T2", t2_file))

        if not analyses_to_run:
            show_error(self, "Error", f"No data files selected for {field_name.replace('field', 'Field ')}.")
            return

        # For now, just run the first one (batch processing will be added later)
        exp_type, input_file = analyses_to_run[0]
        self.results_text.appendPlainText(f"Starting {exp_type} analysis...")

        # Update combo to match and call regular analysis
        index = 0 if exp_type == "T1" else 1
        if field_name == "field1":
            self.analysis_combo.setCurrentIndex(index)
        self._start_analysis()

    def _calculate_t2_from_t1rho(self):
        """Calculate T2 from the T1 and T1ρ fit results of every field that has both.

        Each field converts on its own tilt angle, spin-lock field, carrier and
        spectrometer frequency: they are separate acquisitions. The carrier must be
        entered - it cannot be read off the spectra, and a guessed one tilts every residue
        by the wrong angle. Nothing is written unless every field to convert has one.
        """
        fields = []
        # Field 2 only while it is enabled: disabled, it keeps its results but its
        # parameters cannot be edited, so it must neither convert nor ask for a carrier.
        for field in ('field1', 'field2') if self.field2_enabled else ('field1',):
            data = self.session_results.get(field) or {}
            t1, t1rho = data.get('t1') or {}, data.get('t1rho') or {}
            if t1.get('results_file') and t1rho.get('results_file'):
                fields.append((field, t1['results_file'], t1rho['results_file']))
        if not fields:
            show_error(self, "Error", "No field has both T1 and T1ρ fit results. "
                                      "Please run T1 and T1ρ fitting first.")
            return
        if not self.peak_list_file:
            show_error(self, "Error", "Please select a peak list file for T1ρ → T2 calculation.")
            return
        unset = [field.replace('field', 'Field ') for field, _, _ in fields
                 if getattr(self, f'{field}_carrier_spin').value() <= 0]
        if unset:
            show_error(self, "Carrier required",
                       f"Enter the ¹⁵N carrier (O3, in ppm) for {', '.join(unset)}. It is an "
                       f"acquisition parameter: it cannot be read from the spectra, and a "
                       f"guessed carrier gives every residue the wrong tilt angle.")
            return

        self.results_text.clear()
        self.results_text.appendPlainText("Starting T1ρ → T2 Conversion...")
        for field, t1_file, t1rho_file in fields:
            try:
                self._convert_t1rho_field(field, t1_file, t1rho_file)
            except Exception as e:
                self.results_text.appendPlainText(
                    f"\nError during T1ρ → T2 conversion ({field}):\n{str(e)}")
                import traceback
                self.results_text.appendPlainText(traceback.format_exc())
        self._update_status_indicators()

    def _convert_t1rho_field(self, field, t1_file, t1rho_file):
        """Convert one field's T1 + T1ρ fits to T2, and hold it as that field's T2."""
        _ensure_lunaNMR_path()
        from lunaNMR.utils.t1rho_calculator import calculate_T2_from_T1rho, write_t2_fit_results

        theta = getattr(self, f'{field}_theta_spin').value()
        omega1 = getattr(self, f'{field}_omega1_spin').value()
        carrier = getattr(self, f'{field}_carrier_spin').value()
        spec_freq_mhz = getattr(self, f'{field}_freq_spin').value()
        label = field.replace('field', 'Field ')
        self.results_text.appendPlainText(
            f"\n{label}: T1 {os.path.basename(t1_file)}, T1ρ {os.path.basename(t1rho_file)}")
        self.results_text.appendPlainText(
            f"  θ = {theta}°, ω₁ = {omega1} Hz, carrier = {carrier:.3f} ppm, "
            f"{spec_freq_mhz} MHz")

        r1_data = self._normalize_relaxation_columns(self._load_relaxation_results(t1_file), 'R1')
        r1rho_data = self._normalize_relaxation_columns(
            self._load_relaxation_results(t1rho_file), 'R1rho')
        result = calculate_T2_from_T1rho(
            r1_data=r1_data,
            r1rho_data=r1rho_data,
            theta_deg=theta,
            omega1_hz=omega1,
            carrier_ppm=carrier,
            spec_freq_mhz=spec_freq_mhz,
            peak_list=self._load_peak_list(self.peak_list_file)
        )

        file_prefix = f"{self._run_prefix()}_{field}_T2_from_T1rho"
        output_file = os.path.join(self.output_dir, f"{file_prefix}.csv")
        result.to_csv(output_file, index=False)

        # Also create a fit_results.txt format for compatibility with View Results
        results_txt_file = os.path.join(self.output_dir, f"{file_prefix}_fit_results.txt")
        write_t2_fit_results(result, results_txt_file)

        # The derived T2 becomes the field's T2 unless a measured one is there: a fitted
        # T2 is the more direct number, and `project build` keeps it the same way.
        current = self.session_results[field].get('t2')
        if current and current.get('derived_from') != 't1rho':
            self.results_text.appendPlainText(
                f"  {label} kept the measured T2; the derived one is in "
                f"{os.path.basename(output_file)}")
        else:
            self.session_results[field]['t2'] = {
                'output_dir': self.output_dir,
                'json_folder': None,  # No JSON for derived T2
                'file_prefix': file_prefix,
                'results_file': results_txt_file,
                'n_fitted': len(result),
                'frequency_mhz': spec_freq_mhz,
                'derived_from': 't1rho'  # Mark as derived
            }
            self.fitted_experiments.add(f'{field}_t2')
            self.results_text.appendPlainText(f"  T2 results now available in 'View Results'")

        self.results_text.appendPlainText(f"  Residues processed: {len(result)}")
        self.results_text.appendPlainText(
            f"  Mean T2: {result['T2'].mean():.1f} ± {result['T2'].std():.1f} ms")
        self.results_text.appendPlainText(f"  Results saved to: {output_file}")

    def _load_relaxation_results(self, file_path: str) -> 'pd.DataFrame':
        """Load relaxation results file, auto-detecting delimiter.

        Supports:
        - Tab-delimited files (fit results .txt from dynamiXs)
        - Comma-delimited CSV files
        - Semicolon-delimited files
        """
        import pandas as pd

        # Try to auto-detect delimiter by reading first line
        with open(file_path, 'r') as f:
            first_line = f.readline()

        if '\t' in first_line:
            return pd.read_csv(file_path, dtype=RESIDUE_NAME_DTYPES, sep='\t')
        elif ';' in first_line:
            return pd.read_csv(file_path, dtype=RESIDUE_NAME_DTYPES, sep=';')
        else:
            return pd.read_csv(file_path, dtype=RESIDUE_NAME_DTYPES)

    def _normalize_relaxation_columns(self, df: 'pd.DataFrame', rate_type: str) -> 'pd.DataFrame':
        """Normalize column names in relaxation data DataFrame.

        Handles both rate data (R1, R1rho) and time data (T1, T1rho) by
        converting times to rates if necessary. Times are assumed to be in ms.
        """
        import numpy as np
        df = df.copy()

        # Common column name variations
        residue_cols = ['residue', 'Residue', 'res', 'assignment', 'Assignment']

        if rate_type == 'R1':
            rate_cols = ['R1', 'r1', 'R1_value', '1/T1']
            time_cols = ['T1', 't1', 'T1_value']  # Time data to convert
            err_cols = ['R1_err', 'R1_error', 'r1_err', 'R1_std']
            time_err_cols = ['T1_err', 't1_err', 'T1_error']
        else:  # R1rho
            rate_cols = ['R1rho', 'R1ρ', 'r1rho', 'R1rho_value', '1/T1rho']
            time_cols = ['T1rho', 'T1ρ', 't1rho', 'T1rho_value']
            err_cols = ['R1rho_err', 'R1ρ_err', 'R1rho_error', 'r1rho_err']
            time_err_cols = ['T1rho_err', 'T1ρ_err', 't1rho_err']

        # Normalize residue column
        for col in residue_cols:
            if col in df.columns:
                df = df.rename(columns={col: 'residue'})
                break

        # Check for rate data first
        found_rate = False
        for col in rate_cols:
            if col in df.columns:
                df = df.rename(columns={col: rate_type})
                found_rate = True
                break

        # If no rate data found, try to convert from time data
        if not found_rate:
            for col in time_cols:
                if col in df.columns:
                    # Convert T (ms) to R (s^-1): R = 1000/T
                    df[rate_type] = 1000.0 / df[col]
                    # Also convert error: dR/dT = -1000/T^2, so R_err = T_err * 1000/T^2
                    for err_col in time_err_cols:
                        if err_col in df.columns:
                            df[f'{rate_type}_err'] = df[err_col] * 1000.0 / (df[col] ** 2)
                            break
                    found_rate = True
                    break

        # Check if we found the required data
        if not found_rate:
            available_cols = list(df.columns)
            if rate_type == 'R1':
                expected = "R1, T1, or similar"
            else:
                expected = "R1rho, T1rho, or similar"
            raise ValueError(
                f"Could not find {rate_type} data in file.\n"
                f"Expected columns: {expected}\n"
                f"Available columns: {available_cols}\n"
                f"Make sure you selected the correct file type."
            )

        # Normalize error column (if not already done from time conversion)
        if f'{rate_type}_err' not in df.columns:
            for col in err_cols:
                if col in df.columns:
                    df = df.rename(columns={col: f'{rate_type}_err'})
                    break

        # Add default error if not present
        if f'{rate_type}_err' not in df.columns:
            df[f'{rate_type}_err'] = 0.0

        # Handle NaN/inf from division
        df = df.replace([np.inf, -np.inf], np.nan)
        df = df.dropna(subset=[rate_type])

        return df

    def _load_peak_list(self, file_path: str) -> 'pd.DataFrame':
        """Load peak list file supporting multiple formats.

        Supports:
        - CSV files (.csv) - comma or semicolon delimited
        - Text files (.txt) - tab or space delimited
        - SPARKY peaks (.peaks) - SPARKY format
        """
        import pandas as pd

        ext = os.path.splitext(file_path)[1].lower()

        try:
            if ext == '.peaks':
                # SPARKY format: typically space/tab delimited with header
                # Format: Assignment w1 w2 [other columns...]
                df = pd.read_csv(file_path, dtype=RESIDUE_NAME_DTYPES, delim_whitespace=True, comment='#')
                # Normalize SPARKY column names
                col_map = {
                    'Assignment': 'residue',
                    'w1': 'N_ppm',
                    'w2': 'H_ppm',
                }
                df = df.rename(columns={k: v for k, v in col_map.items() if k in df.columns})

            elif ext == '.txt':
                # Auto-detect delimiter by checking first line
                with open(file_path, 'r') as f:
                    first_line = f.readline()

                if ',' in first_line:
                    # Comma-delimited (lunaNMR format)
                    df = pd.read_csv(file_path, dtype=RESIDUE_NAME_DTYPES, skipinitialspace=True)
                elif '\t' in first_line:
                    # Tab-delimited
                    df = pd.read_csv(file_path, dtype=RESIDUE_NAME_DTYPES, sep='\t', skipinitialspace=True)
                else:
                    # Whitespace-delimited
                    df = pd.read_csv(file_path, dtype=RESIDUE_NAME_DTYPES, delim_whitespace=True)

            else:
                # CSV - try comma first, then semicolon
                try:
                    df = pd.read_csv(file_path, dtype=RESIDUE_NAME_DTYPES, skipinitialspace=True)
                    if len(df.columns) == 1 and ';' in open(file_path).readline():
                        df = pd.read_csv(file_path, dtype=RESIDUE_NAME_DTYPES, sep=';', skipinitialspace=True)
                except:
                    df = pd.read_csv(file_path, dtype=RESIDUE_NAME_DTYPES, sep=';', skipinitialspace=True)

            # Strip whitespace from column names
            df.columns = df.columns.str.strip()

            # Normalize common column name variations
            col_map = {
                # Residue/Assignment columns
                'Assignment': 'residue',
                'Residue': 'residue',
                'res': 'residue',
                # 15N chemical shift columns
                '15N': 'N_ppm',
                'N': 'N_ppm',
                'w1': 'N_ppm',
                'Position_Y': 'N_ppm',  # lunaNMR format (Y is always 15N)
                # 1H chemical shift columns
                '1H': 'H_ppm',
                'H': 'H_ppm',
                'w2': 'H_ppm',
                'Position_X': 'H_ppm',  # lunaNMR format (X is 1H)
            }
            df = df.rename(columns={k: v for k, v in col_map.items() if k in df.columns})

            return df

        except Exception as e:
            raise RuntimeError(f"Failed to load peak list: {e}")

    @Slot(str)
    def _on_progress(self, message: str):
        """Handle progress updates."""
        self.results_text.appendPlainText(message)

    @Slot(int)
    def _on_progress_value(self, value: int):
        """Handle progress bar value updates."""
        self.progress_bar.setValue(value)

    @Slot(object)
    def _on_finished(self, results: dict):
        """Handle analysis completion and store results in session."""
        self.progress_bar.setValue(100)
        n_fitted = results.get('n_fitted') or 0
        for line in self._fit_report(results):
            self.results_text.appendPlainText(line)
        if not n_fitted:
            # Nothing to view or to feed model-free: the session keeps what it had.
            self.progress_bar.setFormat("Nothing fitted")
            self._release_worker()
            return
        self.progress_bar.setFormat("Complete!")

        # Store results in session (persistent across fits)
        field_name = results.get('field_name', 'field1')
        exp_type = results.get('experiment_type', 'T1').lower()  # t1, t2, or t1rho

        # Create result entry
        result_entry = {
            'output_dir': results.get('output_dir'),
            'json_folder': results.get('json_folder'),
            'json_file': results.get('json_file'),
            'file_prefix': results.get('file_prefix'),
            'results_file': results.get('results_file'),
            'n_fitted': results.get('n_fitted'),
            'frequency_mhz': results.get('field_freq', 600.0)
        }

        # Store in session_results
        if field_name not in self.session_results or self.session_results[field_name] is None:
            self.session_results[field_name] = {
                'frequency_mhz': result_entry['frequency_mhz'],
                't1': None,
                't2': None,
                't1rho': None
            }

        self.session_results[field_name][exp_type] = result_entry
        self.session_results[field_name]['frequency_mhz'] = result_entry['frequency_mhz']

        # Track fitted experiments
        session_key = f"{field_name}_{exp_type}"
        self.fitted_experiments.add(session_key)

        # Update status indicators (if they exist)
        self._update_status_indicators()

        self.results_text.appendPlainText(f"\nSession updated: {session_key} stored")
        self.results_text.appendPlainText(f"Fitted experiments: {sorted(self.fitted_experiments)}")
        self._release_worker()

    @staticmethod
    def _fit_report(results: dict) -> list:
        """The page's account of a finished fit: what the CLI's t1t2 summary reports."""
        n_fitted = results.get('n_fitted') or 0
        n_total = results.get('n_total', n_fitted)
        exp = results.get('experiment_type', 'T1')
        if n_fitted:
            lines = [f"\n{exp} analysis complete: {n_fitted}/{n_total} residues fitted"]
        else:
            lines = [f"\nNo residue fitted (0/{n_total}): nothing is stored for this "
                     f"experiment"]
        failed = results.get('failed_residues') or []
        if failed and n_fitted:
            lines.append(f"Not fitted ({len(failed)}): {', '.join(map(str, failed))}")
        dropped = results.get('dropped_columns') or []
        if dropped:
            lines.append(f"WARNING: {len(dropped)} column(s) with no parseable delay were "
                         f"left out of the fit: {', '.join(map(str, dropped))}")
        marginal = results.get('n_window_marginal')
        if marginal:
            lines.append(f"NOTE: {marginal} residues have t_max/T < 2 - the delay window is "
                         f"too short to pin the time constant down for those")
        lines.append(f"Results saved to: {results.get('results_file', 'N/A')}")
        return lines

    def _release_worker(self):
        """Clean up the worker thread to prevent a crash on exit."""
        if hasattr(self, 'worker') and self.worker is not None:
            self.worker.wait()
            self.worker.deleteLater()
            self.worker = None

    @Slot(str)
    def _on_error(self, error_msg: str):
        """Handle analysis error."""
        self.progress_bar.hide()
        self.results_text.appendPlainText(f"\nError during analysis:\n{error_msg}")
        self._release_worker()

    def _update_status_indicators(self):
        """Update visual status indicators for fitted experiments."""
        # Map experiment keys to indicator widgets
        indicator_map = {
            'field1_t1': getattr(self, 'field1_t1_indicator', None),
            'field1_t2': getattr(self, 'field1_t2_indicator', None),
            'field1_t1rho': getattr(self, 'field1_t1rho_indicator', None),
            'field2_t1': getattr(self, 'field2_t1_indicator', None),
            'field2_t2': getattr(self, 'field2_t2_indicator', None),
            'field2_t1rho': getattr(self, 'field2_t1rho_indicator', None),
        }

        for key, indicator in indicator_map.items():
            if indicator is None:
                continue

            # Get the experiment type label (T1, T2, T1ρ)
            exp_label = key.split('_')[-1].upper()
            if exp_label == 'T1RHO':
                exp_label = 'T1ρ'

            if key in self.fitted_experiments:
                indicator.setText(f"{exp_label} ✓")
                indicator.setStyleSheet("color: #4CAF50; font-weight: bold;")  # Green
            else:
                indicator.setText(f"{exp_label} ○")
                indicator.setStyleSheet("color: #9E9E9E;")  # Gray

        # Update T1ρ → T2 conversion status indicators
        if hasattr(self, 'conv_t1_status'):
            if 'field1_t1' in self.fitted_experiments:
                self.conv_t1_status.setText("T1: ✓")
                self.conv_t1_status.setStyleSheet("color: #4CAF50; font-weight: bold;")
            else:
                self.conv_t1_status.setText("T1: ○")
                self.conv_t1_status.setStyleSheet(f"color: {SECONDARY_TEXT};")

        if hasattr(self, 'conv_t1rho_status'):
            if 'field1_t1rho' in self.fitted_experiments:
                self.conv_t1rho_status.setText("T1ρ: ✓")
                self.conv_t1rho_status.setStyleSheet("color: #4CAF50; font-weight: bold;")
            else:
                self.conv_t1rho_status.setText("T1ρ: ○")
                self.conv_t1rho_status.setStyleSheet(f"color: {SECONDARY_TEXT};")

    def get_session_state(self) -> dict:
        """Get current session state for project save.

        Returns serializable dict with all UI settings and session results.
        """
        return {
            # Session results (fit outputs)
            'session_results': self.session_results,
            'fitted_experiments': list(self.fitted_experiments),

            # Field 1 files
            'field1_t1_file': self.field1_t1_file,
            'field1_t2_file': self.field1_t2_file,
            'field1_t1rho_file': self.field1_t1rho_file,

            # Field 2 files
            'field2_enabled': self.field2_enabled,
            'field2_t1_file': self.field2_t1_file,
            'field2_t2_file': self.field2_t2_file,
            'field2_t1rho_file': self.field2_t1rho_file,

            # Source series mapping (for memory-based data from lunaNMR)
            'source_series_map': getattr(self, 'source_series_map', {}),

            # UI settings
            'field1_freq': self.field1_freq_spin.value(),
            'field2_freq': self.field2_freq_spin.value(),
            'analysis_selection': self.analysis_combo.currentText(),
            'prefix': self.prefix_edit.text(),
            'file_prefixes': self._file_prefixes(),
            'output_dir': self.output_dir,

            # T1rho parameters
            'theta_angle': self.field1_theta_spin.value(),
            'omega1': self.field1_omega1_spin.value(),
            'field1_carrier': self.field1_carrier_spin.value(),
            'field2_theta_angle': self.field2_theta_spin.value(),
            'field2_omega1': self.field2_omega1_spin.value(),
            'field2_carrier': self.field2_carrier_spin.value(),
            'peak_list_file': self.peak_list_file,

            # Advanced parameters
            'initial_amplitude': self.amp_spin.value(),
            'initial_decay': self.time_spin.value(),
            'error_method': self.error_method_combo.currentIndex(),
            'bootstrap_iter': self.boot_spin.value(),
            'delay_units': self.units_combo.currentText(),
            'fit_on': self._observable(),
        }

    def restore_session_state(self, state: dict):
        """Restore session state from project load.

        Args:
            state: Dict from get_session_state()
        """
        if not state:
            return

        # Restore session results
        self.session_results = state.get('session_results', self.session_results)
        self.fitted_experiments = set(state.get('fitted_experiments', []))

        # Restore Field 1 files
        self.field1_t1_file = state.get('field1_t1_file')
        self.field1_t2_file = state.get('field1_t2_file')
        self.field1_t1rho_file = state.get('field1_t1rho_file')

        # Restore Field 2 settings
        self.field2_t1_file = state.get('field2_t1_file')
        self.field2_t2_file = state.get('field2_t2_file')
        self.field2_t1rho_file = state.get('field2_t1rho_file')

        # Restore source series mapping (for memory-based data)
        self.source_series_map = state.get('source_series_map', {})

        # Enable Field 2 if it was enabled
        if state.get('field2_enabled', False) and not self.field2_enabled:
            self._toggle_field2()

        # Restore UI settings
        self.field1_freq_spin.setValue(state.get('field1_freq', 600.0))
        self.field2_freq_spin.setValue(state.get('field2_freq', 700.0))

        # Restore analysis selection
        analysis_text = state.get('analysis_selection', '')
        if analysis_text:
            idx = self.analysis_combo.findText(analysis_text)
            if idx >= 0:
                self.analysis_combo.setCurrentIndex(idx)

        self.prefix_edit.setText(state.get('prefix', 'T1_analysis'))
        self.output_dir = state.get('output_dir')
        if self.output_dir:
            self.outdir_display.setText(os.path.basename(self.output_dir))

        # Restore T1rho parameters
        self.field1_theta_spin.setValue(state.get('theta_angle', 45.0))
        self.field1_omega1_spin.setValue(state.get('omega1', 2000.0))
        self.field1_carrier_spin.setValue(state.get('field1_carrier', 0.0))
        self.field2_theta_spin.setValue(state.get('field2_theta_angle', 45.0))
        self.field2_omega1_spin.setValue(state.get('field2_omega1', 2000.0))
        self.field2_carrier_spin.setValue(state.get('field2_carrier', 0.0))
        self.peak_list_file = state.get('peak_list_file')
        if self.peak_list_file:
            self.peak_list_display.setText(os.path.basename(self.peak_list_file))

        # Restore advanced parameters
        self.amp_spin.setValue(state.get('initial_amplitude', 5.0))
        self.time_spin.setValue(state.get('initial_decay', 100.0))
        self.error_method_combo.setCurrentIndex(state.get('error_method', 0))
        self.boot_spin.setValue(state.get('bootstrap_iter', 1000))
        self.units_combo.setCurrentText(state.get('delay_units', TIME_UNITS_DEFAULT))
        self._set_observable(state.get('fit_on'))

        # Update status indicators
        self._update_status_indicators()

        # Update field displays (handles both series and file paths)
        self._update_field_displays()

    def _browse_peak_list_file(self):
        """Browse for peak list file."""
        filename = open_file_dialog(
            self,
            title="Select Peak List File",
            directory=self.main_window.current_dir,
            filter="Peak Lists (*.csv *.txt *.peaks);;CSV files (*.csv);;Text files (*.txt);;SPARKY peaks (*.peaks);;All files (*.*)"
        )
        if filename:
            self.peak_list_file = filename
            self.peak_list_display.setText(os.path.basename(filename))

    def _display_results(self):
        """Open the T1T2 results viewer window with session results."""
        try:
            from visualization.t1t2_results_viewer import T1T2ResultsViewer

            # Build session_results dict for the viewer
            # Contains all fitted experiments with their result files
            viewer_session = {}
            for field_name, field_data in self.session_results.items():
                if field_data is None:
                    continue
                for exp_type in ['t1', 't2', 't1rho']:
                    if field_data.get(exp_type):
                        result_entry = field_data[exp_type]
                        key = f"{field_name}_{exp_type}"
                        viewer_session[key] = {
                            'results_file': result_entry.get('results_file'),
                            'frequency_mhz': result_entry.get('frequency_mhz'),
                            'experiment_type': exp_type.upper().replace('T1RHO', 'T1rho')
                        }

            if viewer_session:
                # Multi-experiment viewer with stacked subplots
                self.results_viewer = T1T2ResultsViewer(
                    parent=None,
                    session_results=viewer_session
                )
                self.results_viewer.show()
            else:
                # No session data - show message
                self.results_text.appendPlainText("No fitted results to display. Run fitting first.")
        except ImportError as e:
            self.results_text.appendPlainText(f"Could not open Results Viewer: {e}")
        except Exception as e:
            self.results_text.appendPlainText(f"Error opening Results Viewer: {e}")

    def _browse_output_dir(self):
        """Browse for output directory."""
        initial_dir = self.output_dir if self.output_dir else self.main_window.current_dir
        folder = open_directory_dialog(self, "Select Output Directory", initial_dir)
        if folder:
            self.output_dir = folder
            self.outdir_display.setText(os.path.basename(folder))

    def _browse_json_folder(self):
        """Browse for JSON folder."""
        initial_dir = self.output_dir if self.output_dir else self.main_window.current_dir
        folder = open_directory_dialog(self, "Select JSON Data Folder", initial_dir)
        if folder:
            self.json_display.setText(os.path.basename(folder))

    def _open_fit_viewer(self):
        """Open the fit viewer window with session data."""
        try:
            from visualization.fit_viewer import FitViewer

            # The session's own fit JSONs: a folder shared with other runs holds theirs too
            json_files = sorted({e['json_file'] for e in self._session_entries()
                                 if e.get('json_file') and os.path.exists(e['json_file'])})
            if not json_files:
                self.results_text.appendPlainText(
                    "The session records no fit JSON on disk: nothing to view.")
                return

            # Get source series map for Inspect Peak functionality
            series_map = getattr(self, 'source_series_map', {})

            # Store as instance variable to prevent garbage collection
            self.fit_viewer = FitViewer(
                parent=None,  # Separate window
                json_files=json_files,
                source_series_map=series_map
            )

            # Connect FitViewer's Inspect Peak signal to DynamiXsDialog handler
            # (main_window is DynamiXsDialog when embedded in lunaNMR)
            if hasattr(self.main_window, 'connect_fit_viewer_signals'):
                self.main_window.connect_fit_viewer_signals(self.fit_viewer)

            self.fit_viewer.show()
        except ImportError as e:
            self.results_text.appendPlainText(f"Could not open Fit Viewer: {e}")
        except Exception as e:
            self.results_text.appendPlainText(f"Error opening Fit Viewer: {e}")


# =============================================================================
# SPECTRAL DENSITY PAGE
# =============================================================================

class SpectralDensityPage(BasePage):
    """Spectral density analysis page."""

    def __init__(self, main_window: DynamiXsMainWindow):
        super().__init__(main_window, "Spectral Density Analysis")
        self.input_file1 = None
        self.input_file2 = None
        self.worker = None
        self.output_dir = None  # Set when analysis runs

        # Session state tracking
        self.session_results = {
            'analysis_complete': False,
            'results_prefix': None,
            'output_dir': None,
            'results_file': None,  # Path to detailed results CSV
        }

        self._setup_content()

    def _setup_content(self):
        """Setup the page content."""
        # Create scroll area for content
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)

        scroll_content = QWidget()
        scroll_layout = create_v_layout(SPACING_SM)
        scroll_content.setLayout(scroll_layout)

        # Analysis type selection
        analysis_frame = QFrame()
        analysis_layout = create_h_layout(SPACING_SM)
        analysis_frame.setLayout(analysis_layout)

        analysis_label = create_label("Analysis Type:")
        analysis_layout.addWidget(analysis_label)

        self.analysis_group = QButtonGroup(self)
        self.single_field_radio = create_radio_button("Single Field", checked=True)
        self.dual_field_radio = create_radio_button("Dual Field")
        self.analysis_group.addButton(self.single_field_radio, 1)
        self.analysis_group.addButton(self.dual_field_radio, 2)
        self.single_field_radio.toggled.connect(self._update_analysis_type)
        analysis_layout.addWidget(self.single_field_radio)
        analysis_layout.addWidget(self.dual_field_radio)
        analysis_layout.addStretch()

        scroll_layout.addWidget(analysis_frame)

        # Spectral density method
        method_frame = QFrame()
        method_layout = create_h_layout(SPACING_SM)
        method_frame.setLayout(method_layout)

        method_label = create_label("Spectral Density Method:")
        method_layout.addWidget(method_label)

        self.method_group = QButtonGroup(self)
        self.jwh_radio = create_radio_button("J(ωH)", checked=True)
        self.j087_radio = create_radio_button("J(0.87ωH)")
        self.method_group.addButton(self.jwh_radio, 1)
        self.method_group.addButton(self.j087_radio, 2)
        method_layout.addWidget(self.jwh_radio)
        method_layout.addWidget(self.j087_radio)
        method_layout.addStretch()

        scroll_layout.addWidget(method_frame)

        # Field parameters
        field_group = StyledGroupBox("Field Parameters")
        field_layout = create_v_layout(SPACING_XS, (SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD))
        field_group.setLayout(field_layout)

        # Field 1
        field1_row = QFrame()
        field1_layout = create_h_layout(SPACING_SM)
        field1_row.setLayout(field1_layout)
        field1_layout.addWidget(create_label("Field 1 Frequency (MHz):"))
        self.field1_spin = create_double_spin_box(600.0, 100.0, 1200.0, 1, width=100)
        field1_layout.addWidget(self.field1_spin)
        field1_layout.addStretch()
        field_layout.addWidget(field1_row)

        # Field 2
        field2_row = QFrame()
        field2_layout = create_h_layout(SPACING_SM)
        field2_row.setLayout(field2_layout)
        field2_layout.addWidget(create_label("Field 2 Frequency (MHz):"))
        self.field2_spin = create_double_spin_box(700.0, 100.0, 1200.0, 1, width=100)
        self.field2_spin.setEnabled(False)
        field2_layout.addWidget(self.field2_spin)
        field2_layout.addStretch()
        field_layout.addWidget(field2_row)

        scroll_layout.addWidget(field_group)

        # File import section
        files_group = StyledGroupBox("Data Import")
        files_layout = create_v_layout(SPACING_SM, (SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD))
        files_group.setLayout(files_layout)

        # File 1
        file1_label = create_label("Data File 1:")
        files_layout.addWidget(file1_label)

        file1_row = QFrame()
        file1_row_layout = create_h_layout(SPACING_SM)
        file1_row.setLayout(file1_row_layout)
        self.file1_display = create_label("No file selected", font_size=FONT_SIZE_SMALL)
        self.file1_display.setStyleSheet(f"color: {SECONDARY_TEXT};")
        file1_row_layout.addWidget(self.file1_display, stretch=1)
        browse1_btn = create_secondary_button("Browse", clicked=self._browse_file1, width=80)
        file1_row_layout.addWidget(browse1_btn)
        files_layout.addWidget(file1_row)

        # File 2
        file2_label = create_label("Data File 2:")
        files_layout.addWidget(file2_label)

        file2_row = QFrame()
        file2_row_layout = create_h_layout(SPACING_SM)
        file2_row.setLayout(file2_row_layout)
        self.file2_display = create_label("Not required for single field", font_size=FONT_SIZE_SMALL)
        self.file2_display.setStyleSheet(f"color: {SECONDARY_TEXT};")
        file2_row_layout.addWidget(self.file2_display, stretch=1)
        self.browse2_btn = create_secondary_button("Browse", clicked=self._browse_file2, width=80)
        self.browse2_btn.setEnabled(False)
        file2_row_layout.addWidget(self.browse2_btn)
        files_layout.addWidget(file2_row)

        scroll_layout.addWidget(files_group)

        # Advanced parameters
        advanced_group = StyledGroupBox("Advanced Parameters")
        advanced_layout = create_v_layout(SPACING_XS, (SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD))
        advanced_group.setLayout(advanced_layout)

        # N-H Bond Length
        rnh_row = QFrame()
        rnh_layout = create_h_layout(SPACING_SM)
        rnh_row.setLayout(rnh_layout)
        rnh_layout.addWidget(create_label("N-H Bond Length (Å):"))
        self.rnh_spin = create_double_spin_box(1.015, 0.9, 1.2, 3, width=100)
        rnh_layout.addWidget(self.rnh_spin)
        rnh_layout.addStretch()
        advanced_layout.addWidget(rnh_row)

        # 15N CSA
        csa_row = QFrame()
        csa_layout = create_h_layout(SPACING_SM)
        csa_row.setLayout(csa_layout)
        csa_layout.addWidget(create_label("15N CSA (ppm):"))
        self.csa_spin = create_double_spin_box(-160.0, -200.0, -100.0, 1, width=100)
        csa_layout.addWidget(self.csa_spin)
        csa_layout.addStretch()
        advanced_layout.addWidget(csa_row)

        # Monte Carlo
        mc_row = QFrame()
        mc_layout = create_h_layout(SPACING_SM)
        mc_row.setLayout(mc_layout)
        self.mc_check = create_check_box("Use Monte Carlo Errors")
        mc_layout.addWidget(self.mc_check)
        mc_layout.addStretch()
        advanced_layout.addWidget(mc_row)

        mc_samples_row = QFrame()
        mc_samples_layout = create_h_layout(SPACING_SM)
        mc_samples_row.setLayout(mc_samples_layout)
        mc_samples_layout.addWidget(create_label("Monte Carlo Samples:"))
        self.mc_spin = create_spin_box(25, 1, 1000, width=100)
        mc_samples_layout.addWidget(self.mc_spin)
        mc_samples_layout.addStretch()
        advanced_layout.addWidget(mc_samples_row)

        scroll_layout.addWidget(advanced_group)

        # Output prefix
        prefix_frame = QFrame()
        prefix_layout = create_v_layout(SPACING_XS)
        prefix_frame.setLayout(prefix_layout)
        prefix_label = create_label("Output Prefix:")
        prefix_layout.addWidget(prefix_label)
        self.prefix_edit = create_line_edit("spectral_density_analysis", width=400)
        prefix_layout.addWidget(self.prefix_edit)
        scroll_layout.addWidget(prefix_frame)

        # Start button
        start_btn = create_primary_button("Start Analysis", clicked=self._start_analysis, width=200)
        scroll_layout.addWidget(start_btn, alignment=Qt.AlignCenter)
        scroll_layout.addWidget(self._make_save_run_button(), alignment=Qt.AlignCenter)

        # Results area
        results_group = StyledGroupBox("Analysis Results")
        results_layout = create_v_layout(SPACING_SM, (SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD))
        results_group.setLayout(results_layout)

        self.results_text = create_plain_text_edit(read_only=True)
        self.results_text.setMinimumHeight(150)
        results_layout.addWidget(self.results_text)

        scroll_layout.addWidget(results_group)

        scroll.setWidget(scroll_content)
        self.content_layout.addWidget(scroll)

    def _update_analysis_type(self):
        """Update UI based on analysis type selection."""
        dual = self.dual_field_radio.isChecked()
        self.field2_spin.setEnabled(dual)
        self.browse2_btn.setEnabled(dual)
        if not dual:
            self.file2_display.setText("Not required for single field")

    def _browse_file1(self):
        """Browse for file 1."""
        filename = open_file_dialog(
            self,
            title="Select Data File 1",
            directory=self.main_window.current_dir,
            filter="CSV files (*.csv);;All files (*.*)"
        )
        if filename:
            self.input_file1 = filename
            self.file1_display.setText(os.path.basename(filename))

    def _browse_file2(self):
        """Browse for file 2."""
        filename = open_file_dialog(
            self,
            title="Select Data File 2",
            directory=self.main_window.current_dir,
            filter="CSV files (*.csv);;All files (*.*)"
        )
        if filename:
            self.input_file2 = filename
            self.file2_display.setText(os.path.basename(filename))

    def _start_analysis(self):
        """Start the spectral density analysis."""
        if self._analysis_running(self.results_text.appendPlainText):
            return
        if not self.input_file1:
            show_error(self, "Error", "Please select at least one input file.")
            return

        # Check for dual-field requirements
        is_dual = self.dual_field_radio.isChecked()
        if is_dual and not self.input_file2:
            show_error(self, "Error", "Dual-field analysis requires two input files.\nPlease select File 2.")
            return

        self.results_text.clear()
        self.results_text.appendPlainText("Starting analysis...")

        # Store output settings for session tracking
        self.output_dir = self.main_window.current_dir
        results_prefix = self.prefix_edit.text()

        # Create parameters
        params = SpectralDensityParams(
            input_file=self.input_file1,
            input_file2=self.input_file2 if self.input_file2 else "",
            output_dir=self.output_dir,
            results_prefix=results_prefix,
            field1_freq=self.field1_spin.value(),
            field2_freq=self.field2_spin.value(),
            dual_field=is_dual,
            use_087=self.j087_radio.isChecked(),
            r_nh=self.rnh_spin.value(),
            csa=self.csa_spin.value(),
            monte_carlo=self.mc_check.isChecked(),
            n_samples=self.mc_spin.value()
        )

        # Create and start worker
        self.worker = SpectralDensityWorker(params)
        self.worker.progress.connect(self._on_progress)
        self.worker.finished.connect(self._on_finished)
        self.worker.error.connect(self._on_error)
        self.worker.start()

    @Slot(str)
    def _on_progress(self, message: str):
        """Handle progress updates."""
        self.results_text.appendPlainText(message)

    @Slot(object)
    def _on_finished(self, results: dict):
        """Handle analysis completion."""
        self.results_text.appendPlainText("\nAnalysis completed successfully!")

        # Store session results for project save/load
        results_prefix = self.prefix_edit.text()
        results_file = os.path.join(self.output_dir, f"{results_prefix}_spectral_density_basic.csv")

        self.session_results = {
            'analysis_complete': True,
            'results_prefix': results_prefix,
            'output_dir': self.output_dir,
            'results_file': results_file if os.path.exists(results_file) else None,
        }

    @Slot(str)
    def _on_error(self, error_msg: str):
        """Handle analysis error."""
        self.results_text.appendPlainText(f"\nError during analysis:\n{error_msg}")

    def get_session_state(self) -> dict:
        """Get current session state for project save.

        Returns serializable dict with all UI settings and session results.
        """
        return {
            # Session results
            'session_results': self.session_results,

            # Input files
            'input_file1': self.input_file1,
            'input_file2': self.input_file2,
            'output_dir': self.output_dir,

            # Analysis type
            'mode': 'dual_field' if self.dual_field_radio.isChecked() else 'single_field',
            'use_087': self.j087_radio.isChecked(),

            # Field frequencies
            'field1_freq': self.field1_spin.value(),
            'field2_freq': self.field2_spin.value(),

            # Physical parameters
            'rnh_angstrom': self.rnh_spin.value(),
            'csa_ppm': self.csa_spin.value(),

            # Analysis options
            'monte_carlo': self.mc_check.isChecked(),
            'n_samples': self.mc_spin.value(),

            # Output prefix, and what the files of the run stored here start with
            'prefix': self.prefix_edit.text(),
            'file_prefixes': ([self.session_results['results_prefix']]
                              if (self.session_results or {}).get('results_prefix') else []),
        }

    def restore_session_state(self, state: dict):
        """Restore session state from project load.

        Args:
            state: Dict from get_session_state()
        """
        if not state:
            return

        # Restore session results
        self.session_results = state.get('session_results', self.session_results)

        # Restore input files
        self.input_file1 = state.get('input_file1')
        self.input_file2 = state.get('input_file2')
        self.output_dir = state.get('output_dir')

        # Update file displays
        if self.input_file1:
            self.file1_display.setText(os.path.basename(self.input_file1))
        if self.input_file2:
            self.file2_display.setText(os.path.basename(self.input_file2))

        # Restore analysis type
        mode = state.get('mode', 'single_field')
        if mode == 'dual_field':
            self.dual_field_radio.setChecked(True)
        else:
            self.single_field_radio.setChecked(True)

        # Restore J(ωH) vs J(0.87ωH) method
        if state.get('use_087', False):
            self.j087_radio.setChecked(True)
        else:
            self.jwh_radio.setChecked(True)

        # Restore field frequencies
        self.field1_spin.setValue(state.get('field1_freq', 600.0))
        self.field2_spin.setValue(state.get('field2_freq', 700.0))

        # Restore physical parameters
        self.rnh_spin.setValue(state.get('rnh_angstrom', 1.015))
        self.csa_spin.setValue(state.get('csa_ppm', -160.0))

        # Restore analysis options
        self.mc_check.setChecked(state.get('monte_carlo', False))
        self.mc_spin.setValue(state.get('n_samples', 50))

        # Restore output prefix
        self.prefix_edit.setText(state.get('prefix', 'spectral_density_analysis'))

        # Update UI based on analysis type
        self._update_analysis_type()


# =============================================================================
# INTEGRATED ANALYSIS PAGE
# =============================================================================

class IntegratedAnalysisPage(BasePage):
    """Integrated Model Free Analysis - automated workflow from raw data."""

    def __init__(self, main_window: DynamiXsMainWindow):
        super().__init__(main_window, "Model Free Analysis")
        self.worker = None
        self._analysis_start_time = None
        self._elapsed_timer = QTimer()
        self._elapsed_timer.timeout.connect(self._update_elapsed_time)
        self._init_file_vars()
        self._setup_content()

    def _init_file_vars(self):
        """Initialize file path variables."""
        # Field 1 files
        self.field1_t1_file = None
        self.field1_t2_file = None
        self.field1_noe_sat_file = None
        self.field1_noe_unsat_file = None

        # Field 2 files
        self.field2_t1_file = None
        self.field2_t2_file = None
        self.field2_noe_sat_file = None
        self.field2_noe_unsat_file = None

        # Output settings - default to current working directory
        self.output_dir = self.main_window.current_dir
        self.json_folder = None

        # Session state tracking
        self.session_results = {
            'analysis_complete': False,
            'results_prefix': None,
            'output_dir': None,
            'json_folder': None,
            'results_file': None,  # Path to spectral density results
        }

    def _setup_content(self):
        """Setup the page content."""
        # Replace default content frame with scroll area
        self.page_layout.removeWidget(self.content_frame)
        self.content_frame.deleteLater()

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)

        scroll_content = QWidget()
        scroll_content.setProperty("class", "panel")
        scroll_layout = create_v_layout(SPACING_MD, (SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD))
        scroll_content.setLayout(scroll_layout)

        # Subtitle
        # Note: "From NMR Spectra" option removed - use LunaNMR Fit Series with
        # "Extract delays" checkbox to produce CSVs, then load them here.
        subtitle = create_secondary_label("Automated workflow from NMR relaxation data to model-free parameters")
        scroll_layout.addWidget(subtitle)

        # --- Project Series Section (drag series to data fields) ---
        series_group = CollapsibleGroupBox("Project Series", collapsible=True, initially_collapsed=False)
        series_content = series_group.content_layout

        self._add_series_panel(
            series_content, "Drag a series to assign it to a T1/T2/hetNOE data field below",
            "No project series available. Use Browse buttons or Auto-Load from Folder.")

        scroll_layout.addWidget(series_group)

        # --- CSV Input Section ---
        csv_input_layout = create_v_layout(SPACING_MD)

        # Field 1 Data Section
        field1_group = CollapsibleGroupBox("Field 1 Data (Required)")
        field1_content = field1_group.content_layout

        # Field 1 frequency in header
        freq1_row = QFrame()
        freq1_layout = create_h_layout(SPACING_XS)
        freq1_row.setLayout(freq1_layout)
        self.field1_freq_spin = create_double_spin_box(600.0, 100.0, 1200.0, 1, " MHz", width=100)
        freq1_layout.addWidget(self.field1_freq_spin)
        field1_group.add_header_widget(freq1_row)

        # Auto-load button
        auto_load_row = QFrame()
        auto_load_layout = create_h_layout(SPACING_SM)
        auto_load_row.setLayout(auto_load_layout)
        self.field1_auto_btn = create_primary_button(
            "Auto-Load from Folder",
            clicked=lambda: self._auto_load_folder("field1"),
            width=200
        )
        auto_load_layout.addWidget(self.field1_auto_btn)
        auto_load_layout.addWidget(create_secondary_label("Automatically detect and load T1, T2, and hetNOE files"))
        auto_load_layout.addStretch()
        field1_content.addWidget(auto_load_row)

        # Manual file rows for Field 1
        self.field1_t1_display = self._create_file_row(field1_content, "T1 Data:", "field1_t1")
        self.field1_t2_display = self._create_file_row(field1_content, "T2 Data:", "field1_t2")
        self.field1_noe_sat_display = self._create_file_row(field1_content, "hetNOE Saturated:", "field1_noe_sat")
        self.field1_noe_unsat_display = self._create_file_row(field1_content, "hetNOE Unsaturated:", "field1_noe_unsat")

        csv_input_layout.addWidget(field1_group)

        # Field 2 Data Section
        field2_group = CollapsibleGroupBox("Field 2 Data (Optional - for dual-field)")
        field2_content = field2_group.content_layout

        # Field 2 frequency and enable button in header
        freq2_row = QFrame()
        freq2_layout = create_h_layout(SPACING_SM)
        freq2_row.setLayout(freq2_layout)
        self.field2_freq_spin = create_double_spin_box(700.0, 100.0, 1200.0, 1, " MHz", width=100)
        self.field2_freq_spin.setEnabled(False)
        freq2_layout.addWidget(self.field2_freq_spin)
        self.enable_field2_btn = create_secondary_button("Enable", clicked=self._toggle_field2, width=80)
        freq2_layout.addWidget(self.enable_field2_btn)
        field2_group.add_header_widget(freq2_row)

        # Auto-load button for Field 2
        auto_load_row2 = QFrame()
        auto_load_layout2 = create_h_layout(SPACING_SM)
        auto_load_row2.setLayout(auto_load_layout2)
        self.field2_auto_btn = create_primary_button(
            "Auto-Load from Folder",
            clicked=lambda: self._auto_load_folder("field2"),
            width=200,
            enabled=False
        )
        auto_load_layout2.addWidget(self.field2_auto_btn)
        auto_load_layout2.addWidget(create_secondary_label("Automatically detect and load T1, T2, and hetNOE files"))
        auto_load_layout2.addStretch()
        field2_content.addWidget(auto_load_row2)

        # Manual file rows for Field 2
        self.field2_t1_display = self._create_file_row(field2_content, "T1 Data:", "field2_t1", enabled=False)
        self.field2_t2_display = self._create_file_row(field2_content, "T2 Data:", "field2_t2", enabled=False)
        self.field2_noe_sat_display = self._create_file_row(field2_content, "hetNOE Saturated:", "field2_noe_sat", enabled=False)
        self.field2_noe_unsat_display = self._create_file_row(field2_content, "hetNOE Unsaturated:", "field2_noe_unsat", enabled=False)

        csv_input_layout.addWidget(field2_group)

        # Add CSV input layout to main scroll layout
        csv_container = QFrame()
        csv_container.setLayout(csv_input_layout)
        scroll_layout.addWidget(csv_container)

        # Analysis Parameters (collapsible)
        params_group = CollapsibleGroupBox("Analysis Parameters", collapsible=True, initially_collapsed=True)
        params_content = params_group.content_layout

        # Analysis method
        method_label = create_label("Analysis Method:", font_size=FONT_SIZE_SECTION_LABEL, bold=True)
        params_content.addWidget(method_label)

        self.method_group = QButtonGroup(self)
        methods = [
            ("Single-field J(ωH)", "single_jwh"),
            ("Single-field J(0.87ωH)", "single_087"),
            ("Dual-field J(ωH) [requires Field 2]", "dual_jwh"),
            ("Dual-field J(0.87ωH) [requires Field 2]", "dual_087")
        ]

        self.method_radios = []
        for i, (label, value) in enumerate(methods):
            radio = create_radio_button(label, checked=(value == "dual_087"))
            self.method_group.addButton(radio, i)
            params_content.addWidget(radio)
            self.method_radios.append((radio, value))

        # Physical parameters
        phys_label = create_label("Physical Parameters:", font_size=FONT_SIZE_SECTION_LABEL, bold=True)
        params_content.addWidget(phys_label)

        phys_grid = QFrame()
        phys_layout = create_grid_layout(SPACING_SM, SPACING_XS)
        phys_grid.setLayout(phys_layout)

        phys_layout.addWidget(create_label("N-H bond length:"), 0, 0)
        self.rnh_spin = create_double_spin_box(1.015, 0.9, 1.2, 3, width=80)
        phys_layout.addWidget(self.rnh_spin, 0, 1)
        phys_layout.addWidget(create_label("Å"), 0, 2)

        phys_layout.addWidget(create_label("15N CSA:"), 1, 0)
        self.csa_spin = create_double_spin_box(-160.0, -200.0, -100.0, 1, width=80)
        phys_layout.addWidget(self.csa_spin, 1, 1)
        phys_layout.addWidget(create_label("ppm"), 1, 2)

        params_content.addWidget(phys_grid)

        # Fitting parameters
        fit_label = create_label("Fitting Parameters:", font_size=FONT_SIZE_SECTION_LABEL, bold=True)
        params_content.addWidget(fit_label)

        fit_grid = QFrame()
        fit_layout = create_grid_layout(SPACING_SM, SPACING_XS)
        fit_grid.setLayout(fit_layout)

        # Initial amplitude guess
        fit_layout.addWidget(create_label("Initial amplitude:"), 0, 0)
        self.init_amp_spin = create_double_spin_box(5.0, 0.001, 1000.0, 3, width=80)
        fit_layout.addWidget(self.init_amp_spin, 0, 1)

        # Initial T1/T2 time constant guess
        fit_layout.addWidget(create_label("Initial T1 (ms):"), 1, 0)
        self.init_t1_spin = create_double_spin_box(800.0, 1.0, 10000.0, 1, width=80)
        fit_layout.addWidget(self.init_t1_spin, 1, 1)
        # One for all the T1 and T2 series; a series not in ms is fitted from a
        # data-driven initial instead of these ms guesses, as `dynamixs modelfree` does.
        fit_layout.addWidget(create_label("Delay units:"), 1, 2)
        fit_layout.addWidget(self._make_units_combo(), 1, 3)

        fit_layout.addWidget(create_label("Initial T2 (ms):"), 2, 0)
        self.init_t2_spin = create_double_spin_box(100.0, 1.0, 5000.0, 1, width=80)
        fit_layout.addWidget(self.init_t2_spin, 2, 1)

        # Error estimation method
        fit_layout.addWidget(create_label("Error method:"), 3, 0)
        self.mf_error_method_combo = create_combo_box(
            ["Analytical (fast)", "Bootstrap"],
            width=130
        )
        self.mf_error_method_combo.currentIndexChanged.connect(self._on_mf_error_method_changed)
        fit_layout.addWidget(self.mf_error_method_combo, 3, 1)

        # Bootstrap iterations (row for bootstrap settings)
        self.mf_boot_label = create_label("Bootstrap:")
        fit_layout.addWidget(self.mf_boot_label, 4, 0)
        self.t1_boot_spin = create_spin_box(1000, 100, 10000, width=80)
        fit_layout.addWidget(self.t1_boot_spin, 4, 1)
        self.mf_boot_iter_label = create_label("iterations")
        fit_layout.addWidget(self.mf_boot_iter_label, 4, 2)
        # Hide bootstrap row by default (analytical is default)
        self.mf_boot_label.hide()
        self.t1_boot_spin.hide()
        self.mf_boot_iter_label.hide()

        # Monte Carlo iterations
        fit_layout.addWidget(create_label("Monte Carlo:"), 5, 0)
        self.mc_spin = create_spin_box(25, 1, 1000, width=80)
        fit_layout.addWidget(self.mc_spin, 5, 1)
        fit_layout.addWidget(create_label("iterations"), 5, 2)

        fit_layout.addWidget(self._make_observable_combo(), 6, 0, 1, 3)

        params_content.addWidget(fit_grid)

        scroll_layout.addWidget(params_group)

        # Output configuration
        output_group = StyledGroupBox("Output Configuration")
        output_layout = create_v_layout(SPACING_SM, (SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD))
        output_group.setLayout(output_layout)

        # Output directory
        outdir_row = QFrame()
        outdir_layout = create_h_layout(SPACING_SM)
        outdir_row.setLayout(outdir_layout)
        outdir_layout.addWidget(create_label("Output Directory:"))
        self.outdir_display = create_label(os.path.basename(self.output_dir), font_size=FONT_SIZE_SMALL)
        self.outdir_display.setStyleSheet(f"color: {SECONDARY_TEXT}; background-color: {FRAME_BG_COLOR}; padding: 4px 8px; border-radius: 4px;")
        outdir_layout.addWidget(self.outdir_display, stretch=1)
        outdir_btn = create_secondary_button("Browse", clicked=self._browse_output_dir, width=80)
        outdir_layout.addWidget(outdir_btn)
        output_layout.addWidget(outdir_row)

        # Output prefix
        prefix_row = QFrame()
        prefix_layout = create_h_layout(SPACING_SM)
        prefix_row.setLayout(prefix_layout)
        prefix_layout.addWidget(create_label("File Prefix:"))
        self.prefix_edit = create_line_edit("integrated_analysis", width=250)
        prefix_layout.addWidget(self.prefix_edit)
        prefix_layout.addStretch()
        output_layout.addWidget(prefix_row)

        # JSON folder
        json_row = QFrame()
        json_layout = create_h_layout(SPACING_SM)
        json_row.setLayout(json_layout)
        json_layout.addWidget(create_label("JSON Data Folder:"))
        self.json_display = create_label("json", font_size=FONT_SIZE_SMALL)
        self.json_display.setStyleSheet(f"color: {SECONDARY_TEXT}; background-color: {FRAME_BG_COLOR}; padding: 4px 8px; border-radius: 4px;")
        json_layout.addWidget(self.json_display, stretch=1)
        json_btn = create_secondary_button("Browse", clicked=self._browse_json_folder, width=80)
        json_layout.addWidget(json_btn)
        output_layout.addWidget(json_row)

        scroll_layout.addWidget(output_group)

        # Action buttons
        btn_frame = QFrame()
        btn_layout = create_h_layout(SPACING_SM)
        btn_frame.setLayout(btn_layout)

        start_btn = create_primary_button("Start Integrated Analysis", clicked=self._start_analysis, width=200)
        btn_layout.addWidget(start_btn)

        reset_btn = create_secondary_button("Reset", clicked=self._reset_form, width=80)
        btn_layout.addWidget(reset_btn)

        view_results_btn = create_secondary_button("View Results", clicked=self._open_results_viewer, width=120)
        btn_layout.addWidget(view_results_btn)

        view_fits_btn = create_secondary_button("View Fits", clicked=self._open_fit_viewer, width=100)
        btn_layout.addWidget(view_fits_btn)
        btn_layout.addWidget(self._make_save_run_button())

        btn_layout.addStretch()
        scroll_layout.addWidget(btn_frame)

        # Progress log
        log_group = StyledGroupBox("Progress Log")
        log_layout = create_v_layout(SPACING_SM, (SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD))
        log_group.setLayout(log_layout)

        # Timer display
        timer_row = QFrame()
        timer_layout = create_h_layout(SPACING_SM)
        timer_row.setLayout(timer_layout)
        timer_layout.addWidget(create_label("Elapsed Time:"))
        self.elapsed_time_label = create_label("--:--")
        self.elapsed_time_label.setStyleSheet(f"font-weight: bold; font-size: {FONT_SIZE_LARGE_HEADER}pt; color: {PRIMARY_TEXT};")
        timer_layout.addWidget(self.elapsed_time_label)
        timer_layout.addStretch()
        log_layout.addWidget(timer_row)

        self.progress_text = create_plain_text_edit(read_only=True)
        self.progress_text.setMinimumHeight(200)
        log_layout.addWidget(self.progress_text)

        scroll_layout.addWidget(log_group)

        scroll.setWidget(scroll_content)
        self.page_layout.addWidget(scroll, stretch=1)

    def _create_file_row(self, parent_layout, label_text: str, field_name: str, enabled: bool = True) -> DropTargetLabel:
        """Create a file selection row with drop-enabled display."""
        row = QFrame()
        row_layout = create_h_layout(SPACING_SM)
        row.setLayout(row_layout)

        label = create_label(label_text)
        label.setFixedWidth(150)
        row_layout.addWidget(label)

        # Use DropTargetLabel for drag-and-drop from series picker
        display = DropTargetLabel(field_name, "No file selected")
        display.setFont(get_font(FONT_SIZE_SMALL))
        display.setEnabled(enabled)
        display.series_dropped.connect(self._on_series_dropped)
        row_layout.addWidget(display, stretch=1)

        browse_btn = create_secondary_button(
            "Browse",
            clicked=lambda: self._browse_file(field_name, display),
            width=80,
            enabled=enabled
        )
        row_layout.addWidget(browse_btn)

        # Store browse button for enabling/disabling
        setattr(self, f"{field_name}_browse_btn", browse_btn)

        parent_layout.addWidget(row)
        return display

    def _browse_file(self, field_name: str, display_label: DropTargetLabel):
        """Browse for a file."""
        filename = open_file_dialog(
            self,
            title=f"Select {field_name.replace('_', ' ').title()} File",
            directory=self.main_window.current_dir,
            filter="CSV files (*.csv);;All files (*.*)"
        )
        if not filename:
            return
        # A T1 or T2 slot reads a relaxation matrix; a hetNOE slot reads one plane's file
        if 'noe' not in field_name and self._refuses_tidy_table(filename):
            return
        setattr(self, f"{field_name}_file", filename)
        display_label.setText(os.path.basename(filename))

    # =========================================================================
    # SERIES PICKER METHODS
    # =========================================================================

    def _on_series_dropped(self, field_name: str, series_name: str, csv_path: str):
        """Handle a series being dropped on a data field.

        A T1 or T2 slot takes the series' intensity matrix, the file `dynamixs modelfree`
        reads. A hetNOE slot takes no series: it reads one plane's intensities, and a
        hetNOE series holds both planes in one matrix.

        Args:
            field_name: Target field (e.g., 'field1_t1', 'field2_noe_sat')
            series_name: Name of the dropped series
            csv_path: Path to series_analysis_tidy.csv
        """
        label = field_name.replace('_', ' ').title()
        if 'noe' in field_name:
            self.progress_text.appendPlainText(
                f"{label} takes one plane's intensities (residue,intensity[,error]), not a "
                f"series: compute the hetNOE planes from series '{series_name}' and browse "
                f"for the file.")
            return
        matrix = self._series_intensity_matrix(series_name)
        display = getattr(self, f"{field_name}_display", None)
        if not matrix:
            self.progress_text.appendPlainText(
                f"Series '{series_name}' has no {self._observable_matrix_name()} beside its outputs: "
                f"browse for the {label} file.")
            if display:
                display.setText(f"⚠️ {series_name} (no intensity matrix)")
            return
        setattr(self, f"{field_name}_file", matrix)
        self._take_recorded_units(matrix, self.progress_text.appendPlainText)
        if display:
            display.setText(f"📊 {series_name}")
            display.setToolTip(f"Series: {series_name}\nMatrix: {matrix}")

        # Store source_series metadata for later use in Inspect Peak
        if not hasattr(self, 'source_series_map'):
            self.source_series_map = {}
        self.source_series_map[field_name] = series_name

        self.progress_text.appendPlainText(f"Assigned series '{series_name}' to {label}")

    def _on_mf_error_method_changed(self, index: int):
        """Show/hide bootstrap iterations based on error method selection."""
        if index == 1:  # Bootstrap selected
            self.mf_boot_label.show()
            self.t1_boot_spin.show()
            self.mf_boot_iter_label.show()
        else:  # Analytical selected
            self.mf_boot_label.hide()
            self.t1_boot_spin.hide()
            self.mf_boot_iter_label.hide()

    def _toggle_field2(self):
        """Toggle Field 2 enabled state."""
        currently_enabled = self.field2_freq_spin.isEnabled()
        new_state = not currently_enabled

        self.field2_freq_spin.setEnabled(new_state)
        self.field2_auto_btn.setEnabled(new_state)
        self.enable_field2_btn.setText("Disable" if new_state else "Enable")

        # Enable/disable browse buttons
        for suffix in ['t1', 't2', 'noe_sat', 'noe_unsat']:
            btn = getattr(self, f"field2_{suffix}_browse_btn", None)
            if btn:
                btn.setEnabled(new_state)

    def _auto_load_folder(self, field_prefix: str):
        """Auto-detect and load files from a folder."""
        folder_path = open_directory_dialog(
            self,
            title=f"Select Folder for {field_prefix.upper()} Data",
            directory=self.main_window.current_dir
        )

        if not folder_path:
            return

        csv_files = glob.glob(os.path.join(folder_path, "*.csv"))
        if not csv_files:
            show_warning(self, "No Files", "No CSV files found in the selected folder.")
            return

        # Detection logic
        t1_file = t2_file = noe_sat_file = noe_unsat_file = None

        for f in csv_files:
            name = os.path.basename(f).upper()
            if 'T1' in name and 'T2' not in name:
                t1_file = f
            elif 'T2' in name:
                t2_file = f
            elif 'NOSAT' in name or 'UNSATURATED' in name:
                noe_unsat_file = f
            elif 'SAT' in name or 'SATURATED' in name:
                noe_sat_file = f

        # Set files and update displays
        files_found = []
        if t1_file:
            setattr(self, f"{field_prefix}_t1_file", t1_file)
            getattr(self, f"{field_prefix}_t1_display").setText(os.path.basename(t1_file))
            files_found.append("T1")
        if t2_file:
            setattr(self, f"{field_prefix}_t2_file", t2_file)
            getattr(self, f"{field_prefix}_t2_display").setText(os.path.basename(t2_file))
            files_found.append("T2")
        if noe_sat_file:
            setattr(self, f"{field_prefix}_noe_sat_file", noe_sat_file)
            getattr(self, f"{field_prefix}_noe_sat_display").setText(os.path.basename(noe_sat_file))
            files_found.append("hetNOE Sat")
        if noe_unsat_file:
            setattr(self, f"{field_prefix}_noe_unsat_file", noe_unsat_file)
            getattr(self, f"{field_prefix}_noe_unsat_display").setText(os.path.basename(noe_unsat_file))
            files_found.append("hetNOE Unsat")

        if files_found:
            # Change button to green
            btn = self.field1_auto_btn if field_prefix == "field1" else self.field2_auto_btn
            btn.setStyleSheet(f"background-color: {SUCCESS_GREEN};")
            self.progress_text.appendPlainText(f"Loaded {', '.join(files_found)} from {folder_path}")
        else:
            show_warning(self, "No Files Detected",
                        "Could not auto-detect file types.\n"
                        "Please load files manually or ensure filenames contain T1, T2, SAT, or NOSAT.")

    def _browse_output_dir(self):
        """Browse for output directory."""
        folder = open_directory_dialog(self, "Select Output Directory", self.output_dir)
        if folder:
            self.output_dir = folder
            self.outdir_display.setText(os.path.basename(folder))

    def _browse_json_folder(self):
        """Browse for JSON folder."""
        folder = open_directory_dialog(self, "Select JSON Data Folder", self.output_dir)
        if folder:
            self.json_display.setText(os.path.basename(folder))

    def _reset_form(self):
        """Reset the form to default values."""
        # Reset file variables
        for prefix in ['field1', 'field2']:
            for suffix in ['t1', 't2', 'noe_sat', 'noe_unsat']:
                setattr(self, f"{prefix}_{suffix}_file", None)
                display = getattr(self, f"{prefix}_{suffix}_display", None)
                if display:
                    display.setText("No file selected")

        # Reset button styles
        self.field1_auto_btn.setStyleSheet("")
        self.field2_auto_btn.setStyleSheet("")

        # Clear progress log
        self.progress_text.clear()
        self.progress_text.appendPlainText("Form reset.")

    def _start_analysis(self):
        """Start the integrated analysis."""
        if self._analysis_running(self.progress_text.appendPlainText):
            return
        # Validate required files
        if not all([self.field1_t1_file, self.field1_t2_file,
                   self.field1_noe_sat_file, self.field1_noe_unsat_file]):
            show_error(self, "Missing Files", "Please load all Field 1 data files.")
            return

        # Require output directory to be set before fitting
        if not self.output_dir or not os.path.isdir(self.output_dir):
            show_warning(
                self,
                "Output Directory Required",
                "Please select an output directory before running analysis.\n"
                "All fit results and JSON files will be saved there."
            )
            # Open directory selection dialog
            initial_dir = self.main_window.current_dir
            folder = open_directory_dialog(self, "Select Output Directory", initial_dir)
            if folder:
                self.output_dir = folder
                self.outdir_display.setText(os.path.basename(folder))
            else:
                # User cancelled - abort
                return

        # Get selected method
        method = "dual_087"
        for radio, value in self.method_radios:
            if radio.isChecked():
                method = value
                break

        # Check if dual-field method requires Field 2 data
        if method.startswith("dual") and not self.field2_freq_spin.isEnabled():
            show_error(self, "Field 2 Required",
                      "Dual-field analysis requires Field 2 data. Please enable and load Field 2 files.")
            return

        time_units = self._checked_units([self.field1_t1_file, self.field1_t2_file,
                                          self.field2_t1_file, self.field2_t2_file])
        if time_units is None:
            return

        self.progress_text.clear()
        self.progress_text.appendPlainText("Starting integrated analysis...")

        # Start timer
        self._analysis_start_time = time.time()
        self.elapsed_time_label.setText("00:00")
        self._elapsed_timer.start(1000)  # Update every second

        # Get error method from combo box
        error_method = "bootstrap" if self.mf_error_method_combo.currentIndex() == 1 else "analytical"

        # Create parameters
        params = IntegratedAnalysisParams(
            field1_t1_file=self.field1_t1_file,
            field1_t2_file=self.field1_t2_file,
            field1_noe_sat_file=self.field1_noe_sat_file,
            field1_noe_unsat_file=self.field1_noe_unsat_file,
            field1_freq=self.field1_freq_spin.value(),
            field2_t1_file=self.field2_t1_file or "",
            field2_t2_file=self.field2_t2_file or "",
            field2_noe_sat_file=self.field2_noe_sat_file or "",
            field2_noe_unsat_file=self.field2_noe_unsat_file or "",
            field2_freq=self.field2_freq_spin.value(),
            dual_field=self.field2_freq_spin.isEnabled(),
            analysis_method=method,
            r_nh=self.rnh_spin.value(),
            csa=self.csa_spin.value(),
            initial_amplitude=self.init_amp_spin.value(),
            initial_t1=self.init_t1_spin.value(),
            initial_t2=self.init_t2_spin.value(),
            n_bootstrap=self.t1_boot_spin.value(),
            error_method=error_method,
            n_monte_carlo=self.mc_spin.value(),
            output_dir=self.output_dir,
            results_prefix=self.prefix_edit.text(),
            json_folder=os.path.join(self.output_dir, self.json_display.text()),
            time_units=time_units
        )

        # Create and start worker
        self.worker = IntegratedAnalysisWorker(params)
        self.worker.progress.connect(self._on_progress)
        self.worker.finished.connect(self._on_finished)
        self.worker.error.connect(self._on_error)
        self.worker.start()

    @Slot(str)
    def _on_progress(self, message: str):
        """Handle progress updates."""
        self.progress_text.appendPlainText(message)

    @Slot(object)
    def _on_finished(self, results: dict):
        """Handle analysis completion."""
        # Stop timer and calculate final elapsed time
        self._elapsed_timer.stop()
        if self._analysis_start_time:
            elapsed = time.time() - self._analysis_start_time
            minutes, seconds = divmod(int(elapsed), 60)
            hours, minutes = divmod(minutes, 60)
            if hours > 0:
                time_str = f"{hours:02d}:{minutes:02d}:{seconds:02d}"
            else:
                time_str = f"{minutes:02d}:{seconds:02d}"
            self.elapsed_time_label.setText(time_str)
            self.progress_text.appendPlainText(f"\nTotal analysis time: {time_str}")

        self.progress_text.appendPlainText("\n" + "="*50)
        self.progress_text.appendPlainText("INTEGRATED ANALYSIS COMPLETE!")
        self.progress_text.appendPlainText("="*50)
        for line in self._fit_report(results):
            self.progress_text.appendPlainText(line)

        # Store session results for project save/load
        results_prefix = self.prefix_edit.text()
        json_folder = os.path.join(self.output_dir, self.json_display.text())
        results_file = os.path.join(self.output_dir, f"{results_prefix}_spectral_density_detailed.csv")

        self.json_folder = json_folder
        self.session_results = {
            'analysis_complete': True,
            'results_prefix': results_prefix,
            'output_dir': self.output_dir,
            'json_folder': json_folder,
            # The run's own fit JSONs (<prefix>_<field>_<T1|T2>): the folder may hold
            # another run's too.
            'json_files': sorted(glob.glob(os.path.join(
                glob.escape(json_folder), f"{glob.escape(results_prefix)}_field[12]_T[12]_fit_data.json"))),
            'results_file': results_file if os.path.exists(results_file) else None,
        }

        # Auto-open viewers
        self._open_results_viewer()
        self._open_fit_viewer()

    @staticmethod
    def _fit_report(results: dict) -> list:
        """The page's account of a finished model-free run: what the CLI's modelfree
        summary reports."""
        capped = results.get('at_iteration_cap_residues') or []
        lines = [f"Model-free complete: {results.get('n_successful')}/"
                 f"{results.get('n_residues')} residues",
                 f"{results.get('n_at_iteration_cap', 0)} fits stopped at the iteration cap "
                 f"without converging" + (f": {', '.join(capped)}" if capped else "")]
        for field, excluded in (results.get('excluded_residues') or {}).items():
            if excluded:
                lines.append(f"{field.replace('field', 'Field ')}: {len(excluded)} residues "
                             f"excluded at the merge (not in all datasets): "
                             f"{', '.join(map(str, excluded))}")
        flagged = {k: v for k, v in (results.get('validation_warnings') or {}).items() if v}
        if flagged:
            lines.append("Validation warnings (residues per check, all fields): " + ", ".join(
                f"{k.replace('_', ' ')} {v}" for k, v in flagged.items()))
        return lines

    @Slot(str)
    def _on_error(self, error_msg: str):
        """Handle analysis error."""
        self._elapsed_timer.stop()  # Stop timer on error
        self.progress_text.appendPlainText(f"\nError during analysis:\n{error_msg}")

    def _update_elapsed_time(self):
        """Update the elapsed time display."""
        if self._analysis_start_time:
            elapsed = time.time() - self._analysis_start_time
            minutes, seconds = divmod(int(elapsed), 60)
            hours, minutes = divmod(minutes, 60)
            if hours > 0:
                time_str = f"{hours:02d}:{minutes:02d}:{seconds:02d}"
            else:
                time_str = f"{minutes:02d}:{seconds:02d}"
            self.elapsed_time_label.setText(time_str)

    def _open_results_viewer(self):
        """Open the results viewer window."""
        try:
            from visualization.results_viewer import ResultsViewer

            # Build results file path
            prefix = self.prefix_edit.text() if hasattr(self, 'prefix_edit') else "analysis"
            results_file = os.path.join(self.output_dir, f"{prefix}_spectral_density_detailed.csv")

            # Get field frequencies
            field1_freq = self.field1_freq_spin.value() if hasattr(self, 'field1_freq_spin') else 600.0
            field2_freq = self.field2_freq_spin.value() if hasattr(self, 'field2_freq_spin') else None
            is_dual_field = hasattr(self, 'field2_enabled') and self.field2_enabled

            # Store as instance variable to prevent garbage collection
            self.results_viewer = ResultsViewer(
                parent=None,  # Separate window
                results_file=results_file if os.path.exists(results_file) else None,
                field1_freq=field1_freq,
                field2_freq=field2_freq,
                is_dual_field=is_dual_field
            )
            self.results_viewer.show()
        except ImportError as e:
            self.progress_text.appendPlainText(f"Could not open Results Viewer: {e}")
        except Exception as e:
            self.progress_text.appendPlainText(f"Error opening Results Viewer: {e}")

    def _open_fit_viewer(self):
        """Open the fit viewer window."""
        try:
            from visualization.fit_viewer import FitViewer

            # The run's own fits, where it recorded them and they are still there
            json_files = [f for f in (self.session_results or {}).get('json_files') or []
                          if os.path.exists(f)]
            if not json_files:
                self.progress_text.appendPlainText(
                    "The run records no fit JSON on disk: nothing to view.")
                return

            # Get source series map for Inspect Peak functionality
            series_map = getattr(self, 'source_series_map', {})

            # Store as instance variable to prevent garbage collection
            self.fit_viewer = FitViewer(
                parent=None,  # Separate window
                json_files=json_files,
                source_series_map=series_map
            )

            # Connect FitViewer's Inspect Peak signal to DynamiXsDialog handler
            # (main_window is DynamiXsDialog when embedded in lunaNMR)
            if hasattr(self.main_window, 'connect_fit_viewer_signals'):
                self.main_window.connect_fit_viewer_signals(self.fit_viewer)

            self.fit_viewer.show()
        except ImportError as e:
            self.progress_text.appendPlainText(f"Could not open Fit Viewer: {e}")
        except Exception as e:
            self.progress_text.appendPlainText(f"Error opening Fit Viewer: {e}")

    def get_session_state(self) -> dict:
        """Get current session state for project save.

        Returns serializable dict with all UI settings and session results.
        """
        # Get selected analysis method
        analysis_method = "dual_087"
        for radio, value in self.method_radios:
            if radio.isChecked():
                analysis_method = value
                break

        return {
            # Session results
            'session_results': self.session_results,

            # Field 1 files
            'field1_t1_file': self.field1_t1_file,
            'field1_t2_file': self.field1_t2_file,
            'field1_noe_sat_file': self.field1_noe_sat_file,
            'field1_noe_unsat_file': self.field1_noe_unsat_file,

            # Field 2 files
            'field2_t1_file': self.field2_t1_file,
            'field2_t2_file': self.field2_t2_file,
            'field2_noe_sat_file': self.field2_noe_sat_file,
            'field2_noe_unsat_file': self.field2_noe_unsat_file,

            # Field configuration
            'field1_freq': self.field1_freq_spin.value(),
            'field2_freq': self.field2_freq_spin.value(),
            'enable_dual_field': self.field2_freq_spin.isEnabled(),

            # Analysis method
            'analysis_method': analysis_method,

            # Physical parameters
            'rnh_angstrom': self.rnh_spin.value(),
            'csa_ppm': self.csa_spin.value(),

            # Fitting parameters
            'initial_amplitude': self.init_amp_spin.value(),
            't1_initial_time': self.init_t1_spin.value(),
            't2_initial_time': self.init_t2_spin.value(),
            'error_method': self.mf_error_method_combo.currentText(),
            'bootstrap_iterations': self.t1_boot_spin.value(),
            'monte_carlo_iterations': self.mc_spin.value(),
            'delay_units': self.units_combo.currentText(),
            'fit_on': self._observable(),

            # Output configuration
            'output_dir': self.output_dir,
            'prefix': self.prefix_edit.text(),
            'file_prefixes': ([self.session_results['results_prefix']]
                              if (self.session_results or {}).get('results_prefix') else []),
            'json_folder_name': self.json_display.text(),
        }

    def restore_session_state(self, state: dict):
        """Restore session state from project load.

        Args:
            state: Dict from get_session_state()
        """
        if not state:
            return

        # Restore session results
        self.session_results = state.get('session_results', self.session_results)

        # Restore Field 1 files
        self.field1_t1_file = state.get('field1_t1_file')
        self.field1_t2_file = state.get('field1_t2_file')
        self.field1_noe_sat_file = state.get('field1_noe_sat_file')
        self.field1_noe_unsat_file = state.get('field1_noe_unsat_file')

        # Update Field 1 displays
        if self.field1_t1_file:
            self.field1_t1_display.setText(os.path.basename(self.field1_t1_file))
        if self.field1_t2_file:
            self.field1_t2_display.setText(os.path.basename(self.field1_t2_file))
        if self.field1_noe_sat_file:
            self.field1_noe_sat_display.setText(os.path.basename(self.field1_noe_sat_file))
        if self.field1_noe_unsat_file:
            self.field1_noe_unsat_display.setText(os.path.basename(self.field1_noe_unsat_file))

        # Restore Field 2 files
        self.field2_t1_file = state.get('field2_t1_file')
        self.field2_t2_file = state.get('field2_t2_file')
        self.field2_noe_sat_file = state.get('field2_noe_sat_file')
        self.field2_noe_unsat_file = state.get('field2_noe_unsat_file')

        # Enable Field 2 if it was enabled
        if state.get('enable_dual_field', False):
            if not self.field2_freq_spin.isEnabled():
                self._toggle_field2()

        # Update Field 2 displays
        if self.field2_t1_file:
            self.field2_t1_display.setText(os.path.basename(self.field2_t1_file))
        if self.field2_t2_file:
            self.field2_t2_display.setText(os.path.basename(self.field2_t2_file))
        if self.field2_noe_sat_file:
            self.field2_noe_sat_display.setText(os.path.basename(self.field2_noe_sat_file))
        if self.field2_noe_unsat_file:
            self.field2_noe_unsat_display.setText(os.path.basename(self.field2_noe_unsat_file))

        # Restore field frequencies
        self.field1_freq_spin.setValue(state.get('field1_freq', 600.0))
        self.field2_freq_spin.setValue(state.get('field2_freq', 700.0))

        # Restore analysis method
        analysis_method = state.get('analysis_method', 'dual_087')
        for radio, value in self.method_radios:
            if value == analysis_method:
                radio.setChecked(True)
                break

        # Restore physical parameters
        self.rnh_spin.setValue(state.get('rnh_angstrom', 1.015))
        self.csa_spin.setValue(state.get('csa_ppm', -160.0))

        # Restore fitting parameters
        self.init_amp_spin.setValue(state.get('initial_amplitude', 5.0))
        self.init_t1_spin.setValue(state.get('t1_initial_time', 800.0))
        self.init_t2_spin.setValue(state.get('t2_initial_time', 100.0))

        error_method = state.get('error_method', 'Analytical (fast)')
        idx = self.mf_error_method_combo.findText(error_method)
        if idx >= 0:
            self.mf_error_method_combo.setCurrentIndex(idx)

        self.t1_boot_spin.setValue(state.get('bootstrap_iterations', 1000))
        self.mc_spin.setValue(state.get('monte_carlo_iterations', 50))
        self.units_combo.setCurrentText(state.get('delay_units', TIME_UNITS_DEFAULT))
        self.observable_combo.blockSignals(True)   # restoring, not choosing: no re-drop
        self._set_observable(state.get('fit_on'))
        self.observable_combo.blockSignals(False)

        # Restore output configuration
        self.output_dir = state.get('output_dir', self.output_dir)
        self.outdir_display.setText(os.path.basename(self.output_dir))
        self.prefix_edit.setText(state.get('prefix', 'integrated_analysis'))
        self.json_display.setText(state.get('json_folder_name', 'json'))
