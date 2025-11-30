# ABOUTME: Qt dialog for configuring and running series integration on multiple NMR spectra
# ABOUTME: Port of v0.9 series integration workflow to PySide6 with progress tracking

import os
import glob
import re
import logging
import threading
from typing import Optional, List, Dict, Any

import pandas as pd

from PySide6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QRadioButton,
    QGroupBox, QListWidget, QListWidgetItem, QProgressBar,
    QTextEdit, QButtonGroup, QMessageBox, QFileDialog
)
from PySide6.QtCore import Qt, Signal, QThread, QObject

from lunaNMR.gui.base.base_dialog import BaseDialog
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


class ProcessingWorker(QObject):
    """Worker for background series processing."""

    progress = Signal(int, int, str)  # current, total, message
    finished = Signal(object)  # batch_results
    error = Signal(str)

    def __init__(self, processor, nmr_files, reference_peaks, peak_source_mode, voigt_params):
        super().__init__()
        self.processor = processor
        self.nmr_files = nmr_files
        self.reference_peaks = reference_peaks
        self.peak_source_mode = peak_source_mode
        self.voigt_params = voigt_params
        self._cancelled = False

    def cancel(self):
        self._cancelled = True
        if self.processor:
            self.processor.processing_active = False

    def run(self):
        """Run the series processing."""
        try:
            from lunaNMR.processors.multi_spectrum_processor import MultiSpectrumProcessor
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
                progress_callback=progress_callback
            )

            if not self._cancelled:
                self.finished.emit(batch_results)

        except Exception as e:
            import traceback
            logger.error(f"Series processing error: {e}\n{traceback.format_exc()}")
            self.error.emit(str(e))


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
        group = QGroupBox("Peak Source")
        group.setStyleSheet(self._get_group_style())

        layout = QVBoxLayout()
        layout.setSpacing(SPACING_SM)

        self.peak_source_group = QButtonGroup(self)

        # Detected peaks option
        self.detected_radio = QRadioButton("Use Detected Peaks (from 'Fit All Peaks')")
        self.detected_radio.setToolTip("Use peaks detected and fitted in the current spectrum")
        self.peak_source_group.addButton(self.detected_radio)
        layout.addWidget(self.detected_radio)

        # Reference peaks option
        self.reference_radio = QRadioButton("Use Reference Peak List (from loaded file)")
        self.reference_radio.setToolTip("Use peak positions from the loaded peak list file")
        self.peak_source_group.addButton(self.reference_radio)
        layout.addWidget(self.reference_radio)

        # Cascade mode option
        self.cascade_radio = QRadioButton("Cascade Mode (refine positions across series)")
        self.cascade_radio.setToolTip("Start with reference peaks, refine positions for each spectrum")
        self.peak_source_group.addButton(self.cascade_radio)
        layout.addWidget(self.cascade_radio)

        # Default selection
        self.detected_radio.setChecked(True)

        # Connect signals for terminal feedback
        self.detected_radio.toggled.connect(self._on_peak_source_detected)
        self.reference_radio.toggled.connect(self._on_peak_source_reference)
        self.cascade_radio.toggled.connect(self._on_peak_source_cascade)

        group.setLayout(layout)
        return group

    def _on_peak_source_detected(self, checked: bool):
        """Handle detected peaks radio button toggle."""
        if checked:
            print("🔘 Peak source: DETECTED - Using peaks from 'Fit All Peaks'")

    def _on_peak_source_reference(self, checked: bool):
        """Handle reference peaks radio button toggle."""
        if checked:
            print("🔘 Peak source: REFERENCE - Using peaks from loaded peak list file")

    def _on_peak_source_cascade(self, checked: bool):
        """Handle cascade mode radio button toggle."""
        if checked:
            print("🔘 Peak source: CASCADE - Refining positions across series")

    def _create_file_list_section(self) -> QGroupBox:
        """Create NMR file list section."""
        group = QGroupBox("NMR Files to Process")
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

        # Browse button
        browse_layout = QHBoxLayout()
        browse_layout.addStretch()
        browse_btn = QPushButton("Change Folder...")
        browse_btn.clicked.connect(self._browse_folder)
        browse_layout.addWidget(browse_btn)
        layout.addLayout(browse_layout)

        group.setLayout(layout)
        return group

    def _create_progress_section(self) -> QGroupBox:
        """Create progress display section."""
        group = QGroupBox("Progress")
        group.setStyleSheet(self._get_group_style())

        layout = QVBoxLayout()

        # Progress bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        layout.addWidget(self.progress_bar)

        # Status label
        self.status_label = QLabel("Ready to start")
        self.status_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}px; color: {PRIMARY_TEXT};")
        layout.addWidget(self.status_label)

        # Log area
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setMaximumHeight(100)
        self.log_text.setStyleSheet("""
            QTextEdit {
                font-family: monospace;
                font-size: 10px;
                background-color: #F5F5F5;
                border: 1px solid #C7C7CC;
                border-radius: 4px;
            }
        """)
        layout.addWidget(self.log_text)

        group.setLayout(layout)
        return group

    def _create_button_row(self) -> QHBoxLayout:
        """Create button row."""
        layout = QHBoxLayout()

        # Start button
        self.start_btn = QPushButton("Start Integration")
        self.start_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        self.start_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {PRIMARY_BUTTON_BG};
                color: {PRIMARY_BUTTON_TEXT};
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px {SPACING_MD}px;
                font-size: {FONT_SIZE_BODY}px;
                font-weight: bold;
            }}
            QPushButton:hover {{
                background-color: {PRIMARY_BUTTON_HOVER};
            }}
            QPushButton:disabled {{
                background-color: #C7C7CC;
            }}
        """)
        self.start_btn.clicked.connect(self._start_processing)
        layout.addWidget(self.start_btn)

        # Cancel button
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        self.cancel_btn.setStyleSheet(self._get_secondary_button_style())
        self.cancel_btn.clicked.connect(self._cancel_processing)
        self.cancel_btn.setEnabled(False)
        layout.addWidget(self.cancel_btn)

        layout.addStretch()

        # Results viewer buttons (enabled after processing completes)
        # Browse Results button hidden but code kept for future use
        self.browse_results_btn = QPushButton("Browse Results")
        self.browse_results_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        self.browse_results_btn.setStyleSheet(self._get_secondary_button_style())
        self.browse_results_btn.setToolTip("Browse individual spectrum results")
        self.browse_results_btn.clicked.connect(self._open_spectrum_browser)
        self.browse_results_btn.setEnabled(False)
        self.browse_results_btn.setVisible(False)  # Hidden but functional

        self.overlay_viewer_btn = QPushButton("Overlay Viewer")
        self.overlay_viewer_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        self.overlay_viewer_btn.setStyleSheet(self._get_secondary_button_style())
        self.overlay_viewer_btn.setToolTip("View multiple spectra overlaid")
        self.overlay_viewer_btn.clicked.connect(self._open_overlay_viewer)
        self.overlay_viewer_btn.setEnabled(False)
        layout.addWidget(self.overlay_viewer_btn)

        # Close button
        self.close_btn = QPushButton("Close")
        self.close_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        self.close_btn.setStyleSheet(self._get_secondary_button_style())
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

        # Get NMR files
        self.nmr_files = self._get_nmr_files(nmr_folder)

        # Populate list
        for f in self.nmr_files:
            item = QListWidgetItem(os.path.basename(f))
            item.setData(Qt.UserRole, f)
            self.file_list.addItem(item)

        self.file_info_label.setText(f"{len(self.nmr_files)} files found in {os.path.basename(nmr_folder)}")

    def _get_nmr_files(self, folder: str) -> List[str]:
        """Get all NMR files from folder.

        Based on v0.9 get_all_nmr_files() (main_gui.py:2813-2830)
        """
        files = []
        for ext in ["ft", "ft2", "fid"]:
            pattern = os.path.join(folder, f"*.{ext}")
            files.extend(glob.glob(pattern))

        # Natural sort
        def natural_sort_key(s):
            return [int(text) if text.isdigit() else text.lower()
                    for text in re.split(r'(\d+)', s)]

        return sorted(files, key=natural_sort_key)

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

    def _start_processing(self):
        """Start series processing."""
        # Validate
        if not self.nmr_files:
            QMessageBox.warning(self, "No Files", "No NMR files found. Select a folder first.")
            return

        if not self.main_window:
            QMessageBox.warning(self, "Error", "Main window reference not available.")
            return

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

        elif self.reference_radio.isChecked():
            peak_source_mode = "reference"
            if not hasattr(self.main_window.integrator, 'peak_list') or self.main_window.integrator.peak_list.empty:
                QMessageBox.warning(
                    self, "No Peak List",
                    "No reference peak list loaded. Load a peak list file first."
                )
                return
            reference_peaks = self.main_window.integrator.peak_list.copy()

        elif self.cascade_radio.isChecked():
            peak_source_mode = "cascade"
            if not hasattr(self.main_window.integrator, 'peak_list') or self.main_window.integrator.peak_list.empty:
                QMessageBox.warning(
                    self, "No Peak List",
                    "No reference peak list loaded for cascade mode."
                )
                return
            reference_peaks = self.main_window.integrator.peak_list.copy()

        # Get parameters
        voigt_params = self._get_voigt_params()

        # Update UI
        self.start_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)
        self.progress_bar.setValue(0)
        self.status_label.setText("Starting...")
        self.log_text.clear()
        self._log(f"Starting series integration: {len(self.nmr_files)} spectra, {len(reference_peaks)} peaks")

        # Create worker and thread
        self.thread = QThread()
        self.worker = ProcessingWorker(
            None, self.nmr_files, reference_peaks, peak_source_mode, voigt_params
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

    def _on_progress(self, current: int, total: int, message: str):
        """Handle progress update."""
        if total > 0:
            percent = int(current / total * 100)
            self.progress_bar.setValue(percent)

        status = f"Processing spectrum {current}/{total}"
        if message:
            status += f": {message}"
        self.status_label.setText(status)

        if message:
            self._log(message)

    def _on_finished(self, batch_results):
        """Handle processing completion."""
        self.batch_results = batch_results
        self.start_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)
        self.progress_bar.setValue(100)

        if batch_results:
            # Handle both dict format and BatchResults object
            if hasattr(batch_results, 'get_summary'):
                summary = batch_results.get_summary()
            elif isinstance(batch_results, dict):
                summary = batch_results.get('summary', {})
            else:
                summary = {}

            self.status_label.setText(
                f"Complete: {summary.get('successful', 0)}/{summary.get('total_spectra', 0)} spectra processed"
            )
            self._log(f"Series integration complete!")
            self._log(f"  Total: {summary.get('total_spectra', 0)} spectra")
            self._log(f"  Successful: {summary.get('successful', 0)}")
            self._log(f"  Failed: {summary.get('failed', 0)}")
            self._log(f"  Duration: {summary.get('duration', 'N/A')}")

            # Store results in main window
            if self.main_window:
                self.main_window.batch_results = batch_results

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
            self._log("Warning: No results returned")

    def _on_error(self, error_msg: str):
        """Handle processing error."""
        self.start_btn.setEnabled(True)
        self.cancel_btn.setEnabled(False)
        self.status_label.setText(f"Error: {error_msg}")
        self._log(f"ERROR: {error_msg}")

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

        # Sort by spectrum name for consistent ordering
        all_results.sort(key=lambda x: x.get('spectrum_name', ''))

        # Get file manager from main window
        file_manager = getattr(self.main_window, 'file_manager', None) if self.main_window else None

        dialog = MultiSpectrumViewerDialog(
            parent=self,
            all_results=all_results,
            file_manager=file_manager
        )
        dialog.show()

        self._log(f"Opened multi-spectrum overlay viewer with {len(all_results)} spectra")

    def _log(self, message: str):
        """Add message to log."""
        self.log_text.append(message)

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
            }
        }

        return params
