# ABOUTME: Qt dialog for configuring and running series integration on multiple NMR spectra
# ABOUTME: Port of v0.9 series integration workflow to PySide6 with progress tracking

import os
import glob
import re
import logging
from datetime import datetime, timedelta
from typing import List, Dict, Any

import pandas as pd

from PySide6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QRadioButton,
    QGroupBox, QListWidget, QListWidgetItem, QProgressBar,
    QTextEdit, QButtonGroup, QMessageBox, QFileDialog, QCheckBox, QFrame,
    QComboBox
)
from PySide6.QtCore import Qt, Signal, QThread, QObject

from lunaNMR.gui.base.base_dialog import BaseDialog
from lunaNMR.gui.styles.design_system import (
    SPACING_SM, SPACING_MD, SPACING_LG,
    FONT_SIZE_BODY, FONT_SIZE_SECTION_LABEL, FONT_SIZE_SMALL,
    PRIMARY_TEXT, SECONDARY_TEXT,
    SECONDARY_BUTTON_BG, SECONDARY_BUTTON_HOVER, SECONDARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BORDER, BUTTON_CORNER_RADIUS, BUTTON_HEIGHT_DIALOG,
    FRAME_BG_COLOR,
    SUCCESS_GREEN, DESTRUCTIVE_BUTTON_BG, DESTRUCTIVE_BUTTON_HOVER,
    DESTRUCTIVE_BUTTON_TEXT, INFO_BLUE
)

logger = logging.getLogger(__name__)


def natural_sort_key(s: str) -> list:
    """Generate a sort key for natural (human) sorting.

    Treats embedded numbers as integers for proper numeric ordering.
    Example: ["8ms", "54ms", "102ms"] instead of ["102ms", "54ms", "8ms"]

    Args:
        s: String to generate sort key for

    Returns:
        List of strings and integers for comparison
    """
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(r'(\d+)', str(s))]


class ProcessingWorker(QObject):
    """Worker for background series processing."""

    progress = Signal(int, int, str)  # current, total, message
    finished = Signal(object)  # batch_results
    error = Signal(str)

    def __init__(self, processor, nmr_files, reference_peaks, peak_source_mode, voigt_params,
                 extract_delays: bool = False, series_mode: str = "time",
                 pre_detected_peaks=None, sn_from_gui_locked=None,
                 peak_list_contour_min: float = None, lock_detection_threshold: bool = True):
        super().__init__()
        self.processor = processor
        self.nmr_files = nmr_files
        self.reference_peaks = reference_peaks
        self.peak_source_mode = peak_source_mode
        self.voigt_params = voigt_params
        self.extract_delays = extract_delays
        self.series_mode = series_mode
        self.pre_detected_peaks = pre_detected_peaks
        self.sn_from_gui_locked = sn_from_gui_locked
        self.peak_list_contour_min = peak_list_contour_min
        self.lock_detection_threshold = lock_detection_threshold
        self._cancelled = False
        self._paused = False

    def cancel(self):
        self._cancelled = True
        if self.processor:
            self.processor.processing_active = False

    def set_paused(self, paused: bool):
        """Set pause state."""
        self._paused = paused
        if self.processor:
            self.processor.paused = paused

    def run(self):
        """Run the series processing."""
        try:
            from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor
            from lunaNMR.utils.output_manager import OutputManager
            import os
            from datetime import datetime

            # Create processor with parameters
            processor = MultiSpectrumProcessor(self.voigt_params)
            self.processor = processor

            # Set up progress callback that works with MultiSpectrumProcessor format
            def progress_callback(progress_pct, task_desc, log_msg, failed=False):
                if self._cancelled:
                    return False  # Signal to stop
                # Convert to (current, total, message) format
                total = len(self.nmr_files) if self.nmr_files else 100
                current = int(progress_pct * total / 100)
                self.progress.emit(current, total, log_msg or task_desc)
                return True

            # Set OutputManager callback for this thread so log_* functions route to GUI
            OutputManager.set_callback(progress_callback)

            processor.processing_active = True

            # Create output folder
            if self.nmr_files:
                base_folder = os.path.dirname(self.nmr_files[0])
                output_folder = os.path.join(
                    base_folder,
                    f"series_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
                )
            else:
                output_folder = None

            # Process series using correct method name
            batch_results = processor.process_nmr_series(
                nmr_files=self.nmr_files,
                reference_peaks=self.reference_peaks,
                output_folder=output_folder,
                peak_source_mode=self.peak_source_mode,
                progress_callback=progress_callback,
                extract_delays=self.extract_delays,
                series_mode=self.series_mode,
                pre_detected_peaks=self.pre_detected_peaks,
                sn_from_gui_locked=self.sn_from_gui_locked,
                peak_list_contour_min=self.peak_list_contour_min,
                lock_detection_threshold=self.lock_detection_threshold
            )

            if not self._cancelled:
                self.finished.emit(batch_results)

        except Exception as e:
            import traceback
            logger.error(f"Series processing error: {e}\n{traceback.format_exc()}")
            self.error.emit(str(e))
        finally:
            # Clear OutputManager callback for this thread
            from lunaNMR.utils.output_manager import OutputManager
            OutputManager.set_callback(None)


class SeriesIntegrationDialog(BaseDialog):
    """Dialog for configuring and running series integration.

    This dialog allows users to:
        - Select peak source (detected, reference, cascade)
        - View available NMR files
        - Monitor processing progress
        - View results summary

    Based on v0.9 start_series_integration() (main_gui.py:4452-4600)

    Example:
        ```python
        dialog = SeriesIntegrationDialog(parent, main_window)
        if dialog.exec():
            batch_results = dialog.batch_results
        ```
    """

    processing_complete = Signal(object)  # Emits BatchResults when done

    def __init__(self, parent=None, main_window=None):
        """Initialize the series integration dialog.

        Args:
            parent: Parent widget
            main_window: Reference to MainWindow for accessing integrator and parameters
        """
        super().__init__(
            parent=parent,
            title="Series Integration - lunaNMR v1.0",
            default_size=(700, 600),
            modal=True
        )

        self.main_window = main_window
        self.batch_results = None
        self.worker = None
        self.thread = None
        self.nmr_files = []
        self.current_series_name = None  # Name for current series integration

        # Progress tracking
        self.start_time = None  # Set when processing starts
        self.completed_tasks = 0
        self.failed_tasks = 0
        self.total_spectra = 0
        self.current_spectrum_idx = 0
        self.current_spectrum_name = ""
        self.first_spectrum_time = None  # Time to complete first spectrum (for ETA)
        self.last_spectrum_start = None  # When current spectrum started

        # State flags
        self._is_paused = False

        # Build UI
        self.setup_ui()

        # Populate file list
        self._populate_file_list()

        # Center on parent
        if parent:
            self.center_on_parent()
        else:
            self.center_on_screen()

        logger.debug("SeriesIntegrationDialog initialized")

    def setup_ui(self):
        """Setup the dialog user interface."""
        layout = QVBoxLayout()
        layout.setSpacing(SPACING_MD)
        layout.setContentsMargins(SPACING_LG, SPACING_LG, SPACING_LG, SPACING_LG)

        # Title
        title_label = QLabel("Series Integration")
        title_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_SECTION_LABEL}px;
                font-weight: bold;
                color: {PRIMARY_TEXT};
            }}
        """)
        layout.addWidget(title_label)

        # Peak Source section
        peak_source_group = self._create_peak_source_section()
        layout.addWidget(peak_source_group)

        # File list section
        file_list_group = self._create_file_list_section()
        layout.addWidget(file_list_group, stretch=1)

        # Progress section
        progress_group = self._create_progress_section()
        layout.addWidget(progress_group)

        # Button row
        button_layout = self._create_button_row()
        layout.addLayout(button_layout)

        self.setLayout(layout)

    def _create_peak_source_section(self) -> QGroupBox:
        """Create peak source selection section."""
        group = QGroupBox("🎯 Peak Source")
        group.setStyleSheet(self._get_group_style())

        layout = QVBoxLayout()
        layout.setSpacing(SPACING_SM)

        self.peak_source_group = QButtonGroup(self)

        # Cascade mode option (default)
        self.cascade_radio = QRadioButton("Cascade Mode (refine positions across series)")
        self.cascade_radio.setToolTip("Start with detected/reference peaks, refine positions for each spectrum")
        self.peak_source_group.addButton(self.cascade_radio)
        layout.addWidget(self.cascade_radio)

        # Detected peaks option
        self.detected_radio = QRadioButton("Use Detected Peaks (from 'Fit All Peaks')")
        self.detected_radio.setToolTip("Use peaks detected and fitted in the current spectrum")
        self.peak_source_group.addButton(self.detected_radio)
        layout.addWidget(self.detected_radio)

        # Independent mode option
        self.independent_radio = QRadioButton("Independent Mode (full detect + fit per spectrum)")
        self.independent_radio.setToolTip(
            "Runs the exact same process as GUI 'Fit Spectrum' for each spectrum:\n"
            "- Uses original reference peaks for detection (no position drift)\n"
            "- Full PASS 1 learning + adaptive optimization\n"
            "- Complete fitting pipeline per spectrum\n"
            "Best for maximum accuracy when spectra may differ significantly."
        )
        self.peak_source_group.addButton(self.independent_radio)
        layout.addWidget(self.independent_radio)

        # Default selection: Cascade mode
        self.cascade_radio.setChecked(True)

        # Connect signals for terminal feedback
        self.cascade_radio.toggled.connect(self._on_peak_source_cascade)
        self.detected_radio.toggled.connect(self._on_peak_source_detected)
        self.independent_radio.toggled.connect(self._on_peak_source_independent)

        group.setLayout(layout)
        return group

    def _on_peak_source_detected(self, checked: bool):
        """Handle detected peaks radio button toggle."""
        if checked:
            print("Peak source: DETECTED - Using peaks from 'Fit All Peaks'")

    def _on_peak_source_cascade(self, checked: bool):
        """Handle cascade mode radio button toggle."""
        if checked:
            print("Peak source: CASCADE - Refining positions across series")

    def _on_peak_source_independent(self, checked: bool):
        """Handle independent mode radio button toggle."""
        if checked:
            print("Peak source: INDEPENDENT - Full detect + fit for each spectrum")

    def _create_file_list_section(self) -> QGroupBox:
        """Create NMR file list section."""
        group = QGroupBox("📁 NMR Files to Process")
        group.setStyleSheet(self._get_group_style())

        layout = QVBoxLayout()

        # Info label
        self.file_info_label = QLabel("0 files found")
        self.file_info_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}px; color: {PRIMARY_TEXT};")
        layout.addWidget(self.file_info_label)

        # File list
        self.file_list = QListWidget()
        self.file_list.setStyleSheet("""
            QListWidget {
                background-color: white;
                border: 1px solid #C7C7CC;
                border-radius: 6px;
            }
            QListWidget::item {
                padding: 2px;
            }
        """)
        layout.addWidget(self.file_list)

        # Bottom row with X-axis extraction mode and Browse button
        bottom_layout = QHBoxLayout()

        # X-axis value extracted from filenames: off / time series / titration
        xaxis_label = QLabel("X-axis from filenames:")
        xaxis_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}px; color: {PRIMARY_TEXT};")
        bottom_layout.addWidget(xaxis_label)

        self.xaxis_mode_combo = QComboBox()
        # (label, mode) — index 0 is off; modes match DelayExtractor(mode=...)
        for text in ("Off (use filenames)", "Time series (ms/s)", "Titration point"):
            self.xaxis_mode_combo.addItem(text)
        self.xaxis_mode_combo.setToolTip(
            "Off: use spectrum filenames as output column headers.\n"
            "Time series: read _50ms / _1s delays and use as headers (DynamiXs-ready).\n"
            "Titration point: read the _1o0 / _0o5 suffix (o means a decimal point)\n"
            "and plot intensities directly against the titration value."
        )
        self.xaxis_mode_combo.currentIndexChanged.connect(self._on_xaxis_mode_changed)
        bottom_layout.addWidget(self.xaxis_mode_combo)

        bottom_layout.addStretch()

        browse_btn = QPushButton("Change Folder...")
        browse_btn.clicked.connect(self._browse_folder)
        bottom_layout.addWidget(browse_btn)

        layout.addLayout(bottom_layout)

        # Delay info label (shown when extract delays is enabled)
        self.delay_info_label = QLabel("")
        self.delay_info_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}px; color: {SECONDARY_TEXT}; font-style: italic;")
        self.delay_info_label.setVisible(False)
        layout.addWidget(self.delay_info_label)

        group.setLayout(layout)
        return group

    def _create_progress_section(self) -> QGroupBox:
        """Create progress display section (Voigt Fitting dialog style)."""
        group = QGroupBox("📊 Progress")
        group.setStyleSheet(self._get_group_style())

        layout = QVBoxLayout()

        # Progress bar with built-in percentage text (styled via main.qss)
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setFormat("%p%")  # Show percentage inside bar
        self.progress_bar.setTextVisible(True)
        self.progress_bar.setMinimumHeight(24)
        layout.addWidget(self.progress_bar)

        # Current Task frame (like Voigt Fitting)
        task_frame = QFrame()
        task_frame.setStyleSheet("""
            QFrame {
                background-color: #F5F5F5;
                border: 1px solid #C7C7CC;
                border-radius: 6px;
            }
        """)
        task_layout = QVBoxLayout(task_frame)
        task_layout.setContentsMargins(10, 6, 10, 6)

        task_header = QLabel("Current Task:")
        task_header.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}px; color: {SECONDARY_TEXT}; font-weight: bold; border: none; background: transparent;")
        task_layout.addWidget(task_header)

        self.status_label = QLabel("Ready to start")
        self.status_label.setStyleSheet(f"font-size: {FONT_SIZE_BODY}px; color: {PRIMARY_TEXT}; border: none; background: transparent;")
        task_layout.addWidget(self.status_label)

        layout.addWidget(task_frame)

        # Statistics frame with icons (like Voigt Fitting)
        stats_frame = QFrame()
        stats_frame.setStyleSheet("""
            QFrame {
                background-color: #F5F5F5;
                border: 1px solid #C7C7CC;
                border-radius: 6px;
            }
        """)
        stats_layout = QVBoxLayout(stats_frame)
        stats_layout.setContentsMargins(10, 6, 10, 6)

        stats_header = QLabel("📈 Statistics")
        stats_header.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}px; color: {SECONDARY_TEXT}; font-weight: bold; border: none; background: transparent;")
        stats_layout.addWidget(stats_header)

        self.stats_label = QLabel("🔄 Elapsed: -- | ⏱ ETA: --\n✓ Spectra: 0/0")
        self.stats_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}px; color: {PRIMARY_TEXT}; border: none; background: transparent;")
        stats_layout.addWidget(self.stats_label)

        layout.addWidget(stats_frame)

        # Log area with header (like Voigt Fitting "Detailed Log")
        log_header = QLabel("📋 Detailed Log")
        log_header.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}px; color: {SECONDARY_TEXT}; font-weight: bold; margin-top: 4px; border: none; background: transparent;")
        layout.addWidget(log_header)

        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setMaximumHeight(100)
        self.log_text.setStyleSheet("""
            QTextEdit {
                font-family: monospace;
                font-size: 10px;
                background-color: #F5F5F5;
                border: 1px solid #C7C7CC;
                border-radius: 6px;
            }
        """)
        layout.addWidget(self.log_text)

        group.setLayout(layout)
        return group

    def _create_button_row(self) -> QHBoxLayout:
        """Create button row with Voigt Fitting-style colored buttons."""
        layout = QHBoxLayout()

        # Start button - Green with icon (matching Voigt Fitting style)
        self.start_btn = QPushButton("▶ Start Integration")
        self.start_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        self.start_btn.setStyleSheet(self._get_green_button_style())
        self.start_btn.clicked.connect(self._start_processing)
        layout.addWidget(self.start_btn)

        # Pause button - Red with icon (matching Voigt Fitting style)
        self.pause_btn = QPushButton("⏸ Pause")
        self.pause_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        self.pause_btn.setStyleSheet(self._get_red_button_style())
        self.pause_btn.clicked.connect(self._on_pause_clicked)
        self.pause_btn.setEnabled(False)
        layout.addWidget(self.pause_btn)

        # Cancel button - Red with icon (matching Voigt Fitting style)
        self.cancel_btn = QPushButton("✕ Cancel")
        self.cancel_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        self.cancel_btn.setStyleSheet(self._get_red_button_style())
        self.cancel_btn.clicked.connect(self._cancel_processing)
        self.cancel_btn.setEnabled(False)
        layout.addWidget(self.cancel_btn)

        layout.addStretch()

        # Results viewer buttons (enabled after processing completes)
        # Browse Results button hidden but code kept for future use
        self.browse_results_btn = QPushButton("📋 Browse Results")
        self.browse_results_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        self.browse_results_btn.setStyleSheet(self._get_secondary_button_style())
        self.browse_results_btn.setToolTip("Browse individual spectrum results")
        self.browse_results_btn.clicked.connect(self._open_spectrum_browser)
        self.browse_results_btn.setEnabled(False)
        self.browse_results_btn.setVisible(False)  # Hidden but functional

        # Overlay Viewer button - Blue with icon
        self.overlay_viewer_btn = QPushButton("📊 Overlay Viewer")
        self.overlay_viewer_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        self.overlay_viewer_btn.setStyleSheet(self._get_blue_button_style())
        self.overlay_viewer_btn.setToolTip("View multiple spectra overlaid")
        self.overlay_viewer_btn.clicked.connect(self._open_overlay_viewer)
        self.overlay_viewer_btn.setEnabled(False)
        layout.addWidget(self.overlay_viewer_btn)

        # Close button - Green with icon (matching Voigt Fitting style)
        self.close_btn = QPushButton("✓ Close")
        self.close_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        self.close_btn.setStyleSheet(self._get_green_button_style())
        self.close_btn.clicked.connect(self.close)
        layout.addWidget(self.close_btn)

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

    def _get_green_button_style(self) -> str:
        """Get green success button style (matching Voigt Fitting dialog)."""
        return f"""
            QPushButton {{
                background-color: {SUCCESS_GREEN};
                color: white;
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px {SPACING_MD}px;
                font-size: {FONT_SIZE_BODY}px;
                font-weight: bold;
            }}
            QPushButton:hover {{
                background-color: #2DB84E;
            }}
            QPushButton:disabled {{
                background-color: #A8E6B5;
            }}
        """

    def _get_red_button_style(self) -> str:
        """Get red control button style (matching Voigt Fitting dialog)."""
        return f"""
            QPushButton {{
                background-color: {DESTRUCTIVE_BUTTON_BG};
                color: {DESTRUCTIVE_BUTTON_TEXT};
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px {SPACING_MD}px;
                font-size: {FONT_SIZE_BODY}px;
                font-weight: bold;
            }}
            QPushButton:hover {{
                background-color: {DESTRUCTIVE_BUTTON_HOVER};
            }}
            QPushButton:disabled {{
                background-color: #F5B5B2;
            }}
        """

    def _get_blue_button_style(self) -> str:
        """Get blue info button style (matching Voigt Fitting dialog)."""
        return f"""
            QPushButton {{
                background-color: {INFO_BLUE};
                color: white;
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px {SPACING_MD}px;
                font-size: {FONT_SIZE_BODY}px;
                font-weight: bold;
            }}
            QPushButton:hover {{
                background-color: #4A8DD4;
            }}
            QPushButton:disabled {{
                background-color: #B5D4F5;
            }}
        """

    def _populate_file_list(self):
        """Populate the NMR file list from current folder."""
        self.file_list.clear()
        self.nmr_files = []

        if not self.main_window:
            return

        nmr_folder = getattr(self.main_window, 'current_nmr_folder', None)
        if not nmr_folder or not os.path.exists(nmr_folder):
            self.file_info_label.setText("No NMR folder selected")
            return

        # Get NMR files (unsorted - will be sorted by _refresh_file_list_display)
        self.nmr_files = self._get_nmr_files(nmr_folder)

        # Update info label
        self.file_info_label.setText(f"{len(self.nmr_files)} files found in {os.path.basename(nmr_folder)}")

        # Refresh display with correct ordering (by delay if checkbox enabled)
        self._refresh_file_list_display()

    def _get_nmr_files(self, folder: str) -> List[str]:
        """Get all NMR files from folder.

        Based on v0.9 get_all_nmr_files() (main_gui.py:2813-2830)
        """
        files = []
        for ext in ["ft", "ft2", "fid"]:
            pattern = os.path.join(folder, f"*.{ext}")
            files.extend(glob.glob(pattern))

        return sorted(files, key=natural_sort_key)

    def _xaxis_enabled(self) -> bool:
        """True when an x-axis value (time or titration) is extracted from filenames."""
        return self.xaxis_mode_combo.currentIndex() > 0

    def _series_mode(self) -> str:
        """Selected filename-value mode: 'time' or 'titration' (default 'time')."""
        return "titration" if self.xaxis_mode_combo.currentIndex() == 2 else "time"

    def _on_xaxis_mode_changed(self, index: int):
        """Handle x-axis extraction mode change (off / time / titration)."""
        if not self._xaxis_enabled():
            self.delay_info_label.setVisible(False)

        # Re-populate file list to show correct order
        self._refresh_file_list_display()

    def _update_delay_info(self):
        """Update the value-extraction info label based on current file list."""
        if not self._xaxis_enabled():
            self.delay_info_label.setVisible(False)
            return

        titration = self._series_mode() == "titration"
        noun = "titration points" if titration else "delays"
        pattern_hint = "_1o0 or _0o5" if titration else "_50ms or _1s"
        unit = "" if titration else "ms"

        if not self.nmr_files:
            self.delay_info_label.setText(f"No files to extract {noun} from")
            self.delay_info_label.setVisible(True)
            return

        from lunaNMR.utils.delay_extractor import DelayExtractor
        extractor = DelayExtractor(mode=self._series_mode())

        files_with_values = extractor.sort_files_by_delay(self.nmr_files)

        if not files_with_values:
            self.delay_info_label.setText(f"⚠ No {noun} found in filenames (expected pattern: {pattern_hint})")
            self.delay_info_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}px; color: #FF6B6B; font-style: italic;")
        else:
            values = [d for _, d in files_with_values]
            missing = len(self.nmr_files) - len(files_with_values)

            info_text = f"✓ Found {len(values)} {noun}: {min(values):g}{unit} - {max(values):g}{unit}"
            if missing > 0:
                info_text += f" ({missing} files without {noun})"

            self.delay_info_label.setText(info_text)
            self.delay_info_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}px; color: {SECONDARY_TEXT}; font-style: italic;")

        self.delay_info_label.setVisible(True)

    def _refresh_file_list_display(self):
        """Refresh file list display with correct ordering.

        When x-axis extraction is enabled, files are sorted by their extracted
        value (delay or titration point). Otherwise, files are sorted alphabetically.
        """
        if not self.nmr_files:
            return

        self.file_list.clear()

        if self._xaxis_enabled():
            titration = self._series_mode() == "titration"
            unit = "" if titration else "ms"

            from lunaNMR.utils.delay_extractor import DelayExtractor
            extractor = DelayExtractor(mode=self._series_mode())

            # Get files sorted by extracted value with sequence numbers
            files_with_sequence = extractor.sort_files_with_sequence(self.nmr_files)

            # Update nmr_files to be in sorted order
            self.nmr_files = [f for f, _, _ in files_with_sequence]

            # Populate list with the extracted value
            for filepath, value, sequence in files_with_sequence:
                basename = os.path.basename(filepath)
                col_name = extractor.get_column_name(value, sequence)
                display_text = f"{basename}  →  {col_name}{unit}"
                item = QListWidgetItem(display_text)
                item.setData(Qt.UserRole, filepath)
                self.file_list.addItem(item)

            # Update delay info label
            self._update_delay_info()
        else:
            # Sort alphabetically (natural sort)
            self.nmr_files = sorted(self.nmr_files, key=natural_sort_key)

            # Populate list without delay info
            for filepath in self.nmr_files:
                item = QListWidgetItem(os.path.basename(filepath))
                item.setData(Qt.UserRole, filepath)
                self.file_list.addItem(item)

    def _browse_folder(self):
        """Browse for NMR folder."""
        folder = QFileDialog.getExistingDirectory(
            self,
            "Select NMR Data Folder",
            getattr(self.main_window, 'current_nmr_folder', '') if self.main_window else ''
        )

        if folder:
            if self.main_window:
                self.main_window.current_nmr_folder = folder
            self._populate_file_list()

    def _prompt_series_name(self) -> str:
        """Prompt user for series name before starting integration.

        Returns:
            Series name if provided, None if cancelled.
        """
        from PySide6.QtWidgets import QInputDialog

        # Generate suggested name from folder
        suggested_name = "Series_1"
        if self.main_window and hasattr(self.main_window, 'current_nmr_folder'):
            folder = self.main_window.current_nmr_folder
            if folder:
                suggested_name = os.path.basename(folder)

        # Check for existing series to avoid duplicates
        existing_series = getattr(self.main_window, 'saved_series', {}) if self.main_window else {}
        if suggested_name in existing_series:
            # Append number to make unique
            base_name = suggested_name
            counter = 2
            while f"{base_name}_{counter}" in existing_series:
                counter += 1
            suggested_name = f"{base_name}_{counter}"

        # Show input dialog
        name, ok = QInputDialog.getText(
            self,
            "Series Name",
            "Enter a name for this series integration:",
            text=suggested_name
        )

        if not ok or not name.strip():
            return None

        # Validate name (no special characters that could cause file issues)
        name = name.strip()
        invalid_chars = '<>:"/\\|?*'
        for char in invalid_chars:
            name = name.replace(char, '_')

        # Check if name already exists
        if name in existing_series:
            result = QMessageBox.question(
                self,
                "Series Exists",
                f"A series named '{name}' already exists.\n\n"
                "Do you want to overwrite it?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No
            )
            if result != QMessageBox.Yes:
                return None

        return name

    def _start_processing(self):
        """Start series processing."""
        # Validate
        if not self.nmr_files:
            QMessageBox.warning(self, "No Files", "No NMR files found. Select a folder first.")
            return

        if not self.main_window:
            QMessageBox.warning(self, "Error", "Main window reference not available.")
            return

        # Prompt for series name
        series_name = self._prompt_series_name()
        if series_name is None:
            return  # User cancelled

        self.current_series_name = series_name

        # Get peak source
        peak_source_mode = "detected"
        reference_peaks = None

        if self.detected_radio.isChecked():
            peak_source_mode = "detected"
            if not hasattr(self.main_window.integrator, 'fitted_peaks') or not self.main_window.integrator.fitted_peaks:
                QMessageBox.warning(
                    self, "No Peaks",
                    "No detected peaks available. Run 'Fit All Peaks' first."
                )
                return
            reference_peaks = self._convert_fitted_peaks_to_dataframe(
                self.main_window.integrator.fitted_peaks
            )

        elif self.cascade_radio.isChecked():
            peak_source_mode = "cascade"
            if not hasattr(self.main_window.integrator, 'peak_list') or self.main_window.integrator.peak_list.empty:
                QMessageBox.warning(
                    self, "No Peak List",
                    "No reference peak list loaded for cascade mode."
                )
                return
            reference_peaks = self.main_window.integrator.peak_list.copy()

        elif self.independent_radio.isChecked():
            peak_source_mode = "independent"
            # Use fitted_peaks (includes dummies) - user cleans up peak list in GUI before series
            if not hasattr(self.main_window.integrator, 'fitted_peaks') or not self.main_window.integrator.fitted_peaks:
                QMessageBox.warning(
                    self, "No Peaks",
                    "No detected peaks available for independent mode. Run 'Fit All Peaks' first."
                )
                return
            reference_peaks = self._convert_fitted_peaks_to_dataframe(
                self.main_window.integrator.fitted_peaks
            )
            self._log(f"Independent mode: Using {len(reference_peaks)} GUI-detected peaks (includes dummies)")

        # For cascade mode, check if user already ran detection (fitted_peaks available)
        pre_detected_peaks = None
        if peak_source_mode == "cascade":
            if hasattr(self.main_window.integrator, 'fitted_peaks') and self.main_window.integrator.fitted_peaks:
                pre_detected_peaks = self.main_window.integrator.fitted_peaks
                self._log(f"Using {len(pre_detected_peaks)} pre-detected peaks from GUI")
            else:
                self._log("No pre-detected peaks - will run detection on first spectrum")

        # Check if from-GUI mode is active and lock threshold for series
        sn_from_gui_locked = None
        if hasattr(self.main_window, 'sn_from_gui_checkbox') and self.main_window.sn_from_gui_checkbox.isChecked():
            if hasattr(self.main_window.integrator, 'sn_absolute_threshold'):
                sn_from_gui_locked = {
                    'absolute_threshold': self.main_window.integrator.sn_absolute_threshold,
                    'contour_min': self.main_window.contour_min_spin.value()
                }
                self._log(f"Series integration: Locking threshold at {sn_from_gui_locked['absolute_threshold']:.2e} "
                         f"(contour_min={sn_from_gui_locked['contour_min']:.3f})")

        # Always capture contour_min for Peak List mode detection threshold
        # Each spectrum computes: peak_list_contour_threshold = max_intensity * contour_min * 1.05
        peak_list_contour_min = self.main_window.contour_min_spin.value()
        lock_detection_threshold = getattr(self.main_window, 'lock_detection_threshold', True)
        if lock_detection_threshold:
            self._log(f"Peak List mode: Using contour_min={peak_list_contour_min:.5f} (locked from spectrum 1)")
        else:
            self._log(f"Peak List mode: Using contour_min={peak_list_contour_min:.5f} (per-spectrum threshold)")

        # Get parameters
        voigt_params = self._get_voigt_params()

        # Get x-axis extraction options (off / time / titration)
        extract_delays = self._xaxis_enabled()
        series_mode = self._series_mode()

        # Update UI
        self.start_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)
        self.pause_btn.setEnabled(True)
        self.pause_btn.setText("⏸ Pause")
        self.progress_bar.setValue(0)
        self.status_label.setText("Starting...")
        self.log_text.clear()

        # Initialize timing and counters
        self.start_time = datetime.now()
        self.completed_tasks = 0
        self.failed_tasks = 0
        self.total_spectra = len(self.nmr_files)
        self.current_spectrum_idx = 0
        self.current_spectrum_name = ""
        self.first_spectrum_time = None
        self.last_spectrum_start = datetime.now()
        self._is_paused = False
        self.stats_label.setText("Starting...")

        if extract_delays:
            delay_msg = f" ({series_mode} x-axis enabled)"
        else:
            delay_msg = ""
        self._log(f"Starting series integration: {len(self.nmr_files)} spectra, {len(reference_peaks)} peaks{delay_msg}")

        # Create worker and thread
        self.thread = QThread()
        self.worker = ProcessingWorker(
            None, self.nmr_files, reference_peaks, peak_source_mode, voigt_params,
            extract_delays=extract_delays,
            series_mode=series_mode,
            pre_detected_peaks=pre_detected_peaks,
            sn_from_gui_locked=sn_from_gui_locked,
            peak_list_contour_min=peak_list_contour_min,
            lock_detection_threshold=lock_detection_threshold
        )
        self.worker.moveToThread(self.thread)

        # Connect signals
        self.thread.started.connect(self.worker.run)
        self.worker.progress.connect(self._on_progress)
        self.worker.finished.connect(self._on_finished)
        self.worker.error.connect(self._on_error)
        self.worker.finished.connect(self.thread.quit)
        self.worker.error.connect(self.thread.quit)

        # Start
        self.thread.start()

    def _cancel_processing(self):
        """Cancel the current processing."""
        if self.worker:
            self.worker.cancel()
        self._log("Cancelling...")
        self.status_label.setText("Cancelling...")

    def _on_pause_clicked(self):
        """Toggle pause state."""
        self._is_paused = not self._is_paused

        if self._is_paused:
            self.pause_btn.setText("▶ Resume")
            self.status_label.setText("Paused...")
            if self.worker:
                self.worker.set_paused(True)
        else:
            self.pause_btn.setText("⏸ Pause")
            if self.worker:
                self.worker.set_paused(False)

    def _on_progress(self, current: int, total: int, message: str):
        """Handle progress update with elapsed time and ETA."""
        import re

        # Parse spectrum index from message (format: "Processing spectrum N/M" or "Loading filename")
        spectrum_match = re.search(r'spectrum\s+(\d+)/(\d+)', message or '', re.IGNORECASE)
        if spectrum_match:
            spectrum_idx = int(spectrum_match.group(1))
            total_spectra = int(spectrum_match.group(2))
        else:
            # Fallback to passed values
            spectrum_idx = current + 1  # current is 0-based from callback
            total_spectra = total if total > 0 else len(self.nmr_files)

        # Track spectrum changes for ETA calculation
        if spectrum_idx != self.current_spectrum_idx:
            now = datetime.now()

            # If moving to spectrum 2, we just finished spectrum 1 - record first spectrum time
            if spectrum_idx == 2 and self.last_spectrum_start:
                self.first_spectrum_time = (now - self.last_spectrum_start).total_seconds()

            # Update tracking
            self.current_spectrum_idx = spectrum_idx
            self.last_spectrum_start = now
            self.total_spectra = total_spectra

        # Calculate progress: completed spectra / total
        # spectrum_idx is 1-based and represents current (in-progress), so completed = idx - 1
        completed = spectrum_idx - 1
        if total_spectra > 0:
            percent = int((completed / total_spectra) * 100)
            self.progress_bar.setValue(percent)

        # Extract spectrum name from message (e.g., "Loading T1_50ms.ft")
        loading_match = re.search(r'Loading\s+(\S+)', message or '')
        if loading_match:
            self.current_spectrum_name = loading_match.group(1)

        # Update status - simple display of current spectrum
        status = f"Processing: {self.current_spectrum_name}" if self.current_spectrum_name else f"Spectrum {spectrum_idx}/{total_spectra}"
        self.status_label.setText(status)

        # Update statistics with ETA
        self._update_statistics(spectrum_idx, total_spectra)

        # Log significant messages only (not every progress tick)
        if message and 'Loading' in message:
            self._log(message)

    def _on_finished(self, batch_results):
        """Handle processing completion."""
        self.batch_results = batch_results
        self.start_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)
        self.pause_btn.setEnabled(False)
        self.progress_bar.setValue(100)

        # Calculate total time
        if self.start_time:
            total_time = datetime.now() - self.start_time
            total_time_str = str(total_time).split('.')[0]
        else:
            total_time_str = "N/A"

        if batch_results:
            # Handle both dict format and BatchResults object
            if hasattr(batch_results, 'get_summary'):
                summary = batch_results.get_summary()
            elif isinstance(batch_results, dict):
                summary = batch_results.get('summary', {})
            else:
                summary = {}

            # Update final stats
            self.stats_label.setText(
                f"Completed in: {total_time_str}\n"
                f"Successful: {summary.get('successful', 0)} | Failed: {summary.get('failed', 0)}"
            )

            self.status_label.setText(
                f"Complete: {summary.get('successful', 0)}/{summary.get('total_spectra', 0)} spectra processed"
            )
            self._log("Series integration complete!")
            self._log(f"  Total: {summary.get('total_spectra', 0)} spectra")
            self._log(f"  Successful: {summary.get('successful', 0)}")
            self._log(f"  Failed: {summary.get('failed', 0)}")
            self._log(f"  Duration: {summary.get('duration', 'N/A')}")

            # Store results in main window
            if self.main_window:
                # Store as current batch_results (for backward compatibility)
                self.main_window.batch_results = batch_results

                # Also store in saved_series dict with the series name
                if not hasattr(self.main_window, 'saved_series'):
                    self.main_window.saved_series = {}

                series_name = self.current_series_name or "Unnamed_Series"

                # Add metadata to batch_results
                if hasattr(batch_results, 'metadata'):
                    batch_results.metadata['series_name'] = series_name
                    batch_results.metadata['spectrum_count'] = len(batch_results.results) if hasattr(batch_results, 'results') else 0

                self.main_window.saved_series[series_name] = batch_results
                self._log(f"Series saved as: '{series_name}'")

            self.processing_complete.emit(batch_results)

            # Enable results viewer buttons
            self.browse_results_btn.setEnabled(True)
            self.overlay_viewer_btn.setEnabled(True)

            QMessageBox.information(
                self,
                "Complete",
                f"Series integration complete!\n\n"
                f"Processed: {summary.get('total_spectra', 0)} spectra\n"
                f"Successful: {summary.get('successful', 0)}\n"
                f"Failed: {summary.get('failed', 0)}\n"
                f"Duration: {summary.get('duration', 'N/A')}"
            )
        else:
            self.status_label.setText("Processing completed (no results)")
            self.stats_label.setText(f"Completed in: {total_time_str}\nNo results available")
            self._log("Warning: No results returned")

    def _on_error(self, error_msg: str):
        """Handle processing error."""
        self.start_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)
        self.pause_btn.setEnabled(False)
        self.status_label.setText(f"Error: {error_msg}")
        self._log(f"ERROR: {error_msg}", is_error=True)

        QMessageBox.critical(
            self,
            "Processing Error",
            f"Series integration failed:\n\n{error_msg}"
        )

    def _open_spectrum_browser(self):
        """Open spectrum browser dialog to browse individual results."""
        from lunaNMR.gui.dialogs import SpectrumBrowserDialog

        if not self.batch_results:
            QMessageBox.warning(self, "No Data",
                "No series results available.")
            return

        # Get series processor and data folder
        series_processor = getattr(self.main_window, 'series_processor', None) if self.main_window else None
        original_data_folder = getattr(self.main_window, 'current_nmr_folder', None) if self.main_window else None

        dialog = SpectrumBrowserDialog(
            parent=self,
            batch_results=self.batch_results,
            series_processor=series_processor,
            original_data_folder=original_data_folder
        )
        dialog.show()

        self._log("Opened spectrum browser")

    def _open_overlay_viewer(self):
        """Open multi-spectrum overlay viewer."""
        import os
        from lunaNMR.gui.dialogs import MultiSpectrumViewerDialog

        if not self.batch_results:
            QMessageBox.warning(self, "No Data",
                "No series results available.")
            return

        # Get data folder from batch_results metadata (v0.9 main_gui.py:7551 pattern)
        # This is the folder containing the NMR files
        data_folder = None
        if hasattr(self.batch_results, 'metadata'):
            data_folder = self.batch_results.metadata.get('data_folder')
        elif isinstance(self.batch_results, dict):
            metadata = self.batch_results.get('metadata', {})
            data_folder = metadata.get('data_folder')

        # Fallback to main_window.current_nmr_folder if not in metadata
        if not data_folder:
            data_folder = getattr(self.main_window, 'current_nmr_folder', None) if self.main_window else None

        self._log(f"Data folder for spectrum paths: {data_folder}")

        # Get all results as list
        all_results = []
        if hasattr(self.batch_results, 'results'):
            for name, result in self.batch_results.results.items():
                result_copy = result.copy()
                result_copy['spectrum_name'] = name

                # CRITICAL: Construct full path for spectrum_file (v0.9 main_gui.py:7547-7553)
                # result may contain 'spectrum_file' but it's just a filename, not full path
                if 'spectrum_file' not in result_copy or not os.path.isabs(result_copy.get('spectrum_file', '')):
                    if data_folder:
                        result_copy['spectrum_file'] = os.path.join(data_folder, name)
                        self._log(f"Constructed full path for {name}: {result_copy['spectrum_file']}")
                    else:
                        self._log(f"⚠️ WARNING: No data folder - cannot construct path for {name}")

                all_results.append(result_copy)

        # Sort results: by extracted value if enabled, otherwise by natural sort
        if self._xaxis_enabled():
            # Sort by extracted value, handling duplicates (e.g., "50", "50_2", "100")
            from lunaNMR.utils.delay_extractor import DelayExtractor
            extractor = DelayExtractor(mode=self._series_mode())

            def delay_sort_key(result):
                name = result.get('spectrum_name', '')
                # Parse value from column name format (e.g., "50", "50_2")
                parsed = extractor.parse_column_name(str(name))
                if parsed:
                    value, sequence = parsed
                    return (value, sequence, str(name))
                # Try extracting from filename pattern (mode-aware)
                value = extractor.extract_value(str(name))
                if value is not None:
                    return (value, 1, str(name))
                # Put unrecognized at the end
                return (float('inf'), 0, str(name))

            all_results.sort(key=delay_sort_key)
        else:
            # Sort by spectrum name using natural sort (8ms < 54ms < 102ms, not alphabetical)
            all_results.sort(key=lambda x: natural_sort_key(x.get('spectrum_name', '')))

        # Get file manager from main window
        file_manager = getattr(self.main_window, 'file_manager', None) if self.main_window else None

        dialog = MultiSpectrumViewerDialog(
            parent=self,
            all_results=all_results,
            file_manager=file_manager
        )
        dialog.show()

        self._log(f"Opened multi-spectrum overlay viewer with {len(all_results)} spectra")

    def _log(self, message: str, is_error: bool = False):
        """Add timestamped message to log with icons (Voigt Fitting style)."""
        timestamp = datetime.now().strftime('%H:%M:%S')
        icon = "✕" if is_error else "✓"
        log_entry = f"[{timestamp}] {icon}  {message}"

        self.log_text.append(log_entry)

        # Auto-scroll to bottom
        self.log_text.verticalScrollBar().setValue(
            self.log_text.verticalScrollBar().maximum()
        )

    def _update_statistics(self, current_spectrum: int, total_spectra: int):
        """Update elapsed time and ETA display."""
        if self.start_time is None:
            return

        elapsed = datetime.now() - self.start_time
        elapsed_str = str(elapsed).split('.')[0]  # Remove microseconds

        # Calculate ETA based on first spectrum time
        # ETA = first_spectrum_time * 4 * remaining_spectra (rough estimate)
        completed = current_spectrum - 1  # current is 1-based, in-progress
        remaining = total_spectra - completed

        if self.first_spectrum_time and remaining > 0:
            # Use first spectrum time * 4 as estimate per spectrum
            eta_seconds = self.first_spectrum_time * 4 * remaining
            eta_str = str(timedelta(seconds=int(eta_seconds)))
        elif completed > 0:
            # Fallback: use elapsed / completed * remaining
            avg_time = elapsed.total_seconds() / completed
            eta_seconds = avg_time * remaining
            eta_str = str(timedelta(seconds=int(eta_seconds)))
        else:
            eta_str = "Calculating..."

        # Update stats label with icons (Voigt Fitting style)
        stats_text = (
            f"🔄 Elapsed: {elapsed_str} | ⏱ ETA: {eta_str}\n"
            f"✓ Spectra: {current_spectrum}/{total_spectra}"
        )
        self.stats_label.setText(stats_text)

    def _convert_fitted_peaks_to_dataframe(self, fitted_peaks: List[Dict]) -> pd.DataFrame:
        """Convert fitted peaks to DataFrame format.

        Based on v0.9 _convert_fitted_peaks_to_dataframe() (main_gui.py:4512-4587)
        """
        peak_data = []

        for i, peak in enumerate(fitted_peaks):
            if isinstance(peak, dict):
                assignment = (peak.get('assignment') or
                             peak.get('Assignment') or
                             f'Peak_{i+1}')

                # Get positions - try multiple key variants
                position_x = (peak.get('peak_x') or
                             peak.get('ppm_x') or
                             peak.get('Position_X') or
                             peak.get('position_x') or
                             peak.get('pos_f2') or 0)

                position_y = (peak.get('peak_y') or
                             peak.get('ppm_y') or
                             peak.get('Position_Y') or
                             peak.get('position_y') or
                             peak.get('pos_f1') or 0)

                peak_data.append({
                    'Assignment': assignment,
                    'Position_X': float(position_x),
                    'Position_Y': float(position_y)
                })

        return pd.DataFrame(peak_data)

    def _get_voigt_params(self) -> Dict[str, Any]:
        """Get Voigt fitting parameters from main window.

        Returns parameters in the nested structure expected by MultiSpectrumProcessor.
        Based on v0.9 param_manager.get_effective_parameters() structure.
        """
        if not self.main_window:
            return {}

        params = {
            # Detection parameters (CRITICAL - was missing, causing 'search_window_x' KeyError)
            'detection_params': {
                'search_window_x': getattr(self.main_window, 'search_window_x', 0.08),
                'search_window_y': getattr(self.main_window, 'search_window_y', 0.8),
                'noise_threshold': getattr(self.main_window, 'noise_threshold', 3.0),
            },
            # GUI parameters for fitting control
            'gui_params': {
                'fix_positions': getattr(self.main_window, 'fix_positions', False),
                'fix_linewidths': getattr(self.main_window, 'fix_linewidths', False),
                'use_parallel_processing': getattr(self.main_window, 'use_parallel_processing', False),
                'use_centroid_refinement': getattr(self.main_window, 'use_centroid_refinement', True),
                'centroid_window_x_ppm': getattr(self.main_window, 'centroid_window_x_ppm', 0.02),
                'centroid_window_y_ppm': getattr(self.main_window, 'centroid_window_y_ppm', 1.0),
                'centroid_noise_multiplier': getattr(self.main_window, 'centroid_noise_multiplier', 2.0),
                'use_ps2d_multi_peak': getattr(self.main_window, 'use_ps2d_multi_peak', True),
                'use_ps2d_linewidth_reuse': getattr(self.main_window, 'use_ps2d_linewidth_reuse', False),
                # ML training data collection
                'collect_training_data': self._get_ml_collect_enabled(),
                # Custom linewidth parameters (only if enabled)
                'lw_lorentz_1h': getattr(self.main_window, 'lw_lorentz_1h', None) if getattr(self.main_window, 'use_custom_linewidths', False) else None,
                'lw_gauss_1h': getattr(self.main_window, 'lw_gauss_1h', None) if getattr(self.main_window, 'use_custom_linewidths', False) else None,
                'lw_lorentz_15n': getattr(self.main_window, 'lw_lorentz_15n', None) if getattr(self.main_window, 'use_custom_linewidths', False) else None,
                'lw_gauss_15n': getattr(self.main_window, 'lw_gauss_15n', None) if getattr(self.main_window, 'use_custom_linewidths', False) else None,
                # Peak detection parameters
                'height_threshold': getattr(self.main_window, 'height_threshold', 0.1),
                'distance_factor': getattr(self.main_window, 'distance_factor', 2.0),
                'prominence_threshold': getattr(self.main_window, 'prominence_threshold', 0.05),
                'smoothing_sigma': getattr(self.main_window, 'smoothing_sigma', 1.0),
                'max_peaks_fit': getattr(self.main_window, 'max_peaks_fit', 50),
                'max_optimization_iterations': getattr(self.main_window, 'max_optimization_iterations', 50),
            },
            # Fitting parameters
            'fitting_params': {
                'min_r_squared': getattr(self.main_window, 'min_r_squared', 0.5),
                'max_iterations': getattr(self.main_window, 'max_iterations', 1000),
                'fitting_window_x': getattr(self.main_window, 'fitting_window_x', 0.2),
                'fitting_window_y': getattr(self.main_window, 'fitting_window_y', 2.0),
            },
            # Processing options
            'processing_options': {
                'use_parallel_processing': getattr(self.main_window, 'use_parallel_processing', False),
                'use_global_optimization': getattr(self.main_window, 'use_global_optimization', False),
                'enable_cascade_drift_limit': getattr(self.main_window, 'enable_cascade_drift_limit', True),
                'rerun_adaptive_per_spectrum': getattr(self.main_window, 'rerun_adaptive_per_spectrum', False),
                'lock_cluster_assignments': getattr(self.main_window, 'lock_cluster_assignments', False),
                'use_original_reference_for_detection': False,
            }
        }

        # CRITICAL: Pass clusters from reference fit to ensure series uses same clusters
        # The clusters are stored in enhanced_fitter.series_params after "Fit All Peaks"
        reference_clusters = None
        try:
            if (hasattr(self.main_window, 'integrator') and
                hasattr(self.main_window.integrator, 'enhanced_fitter') and
                hasattr(self.main_window.integrator.enhanced_fitter, 'series_params') and
                self.main_window.integrator.enhanced_fitter.series_params):
                reference_clusters = self.main_window.integrator.enhanced_fitter.series_params.get(
                    'locked_clusters_by_assignment'
                )
                if reference_clusters:
                    self._log(f"Passing {len(reference_clusters)} clusters from reference fit to series")
        except Exception as e:
            self._log(f"Warning: Could not get reference clusters: {e}")

        if reference_clusters:
            params['reference_clusters_by_assignment'] = reference_clusters

        # Independent mode: Force full pipeline (rerun adaptive + original reference)
        if hasattr(self, 'independent_radio') and self.independent_radio.isChecked():
            params['processing_options']['rerun_adaptive_per_spectrum'] = True
            params['processing_options']['use_original_reference_for_detection'] = True

        return params

    def _get_ml_collect_enabled(self) -> bool:
        """Check if ML training data collection is enabled via ML Learning Center config."""
        if not self.main_window:
            return False

        try:
            if hasattr(self.main_window, 'config_manager'):
                config = self.main_window.config_manager.config
                ml_config = config.get('ml_learning', {})
                return ml_config.get('collect_training_data', False)
        except Exception:
            pass

        return False
