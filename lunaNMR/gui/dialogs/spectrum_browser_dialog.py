# ABOUTME: Qt-based spectrum browser dialog for browsing individual spectra from series results
# ABOUTME: Port of Tkinter spectrum_browser.py to PySide6 for v1.0 Qt interface

import logging
from typing import Dict

from PySide6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QLineEdit, QComboBox, QTreeWidget, QTreeWidgetItem,
    QHeaderView, QMessageBox
)
from PySide6.QtCore import Qt, Signal

from lunaNMR.gui.base.base_dialog import BaseDialog
from lunaNMR.gui.styles.design_system import (
    SPACING_SM, SPACING_MD, SPACING_LG,
    FONT_SIZE_BODY, FONT_SIZE_SECTION_LABEL, PRIMARY_TEXT, PRIMARY_BUTTON_BG, PRIMARY_BUTTON_HOVER, PRIMARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BG, SECONDARY_BUTTON_HOVER, SECONDARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BORDER, BUTTON_CORNER_RADIUS, BUTTON_HEIGHT_DIALOG,
    FRAME_BG_COLOR
)

logger = logging.getLogger(__name__)


class SpectrumBrowserDialog(BaseDialog):
    """Dialog for browsing and selecting individual spectra from series results.

    This is a Qt port of the Tkinter SpectrumBrowserDialog from spectrum_browser.py.

    Features:
        - List of all spectra in series with status and quality info
        - Search and filter by name, quality, status
        - Sort by any column
        - Double-click to open detailed spectrum viewer
        - Quality indicators with color coding

    Based on v0.9 SpectrumBrowserDialog (spectrum_browser.py:116-478)

    Example:
        ```python
        dialog = SpectrumBrowserDialog(parent, batch_results, series_processor)
        if dialog.exec():
            selected = dialog.selected_spectrum
        ```
    """

    spectrum_selected = Signal(str)  # Emits spectrum name when selected

    def __init__(self, parent=None, batch_results=None, series_processor=None,
                 original_data_folder=None):
        """Initialize the spectrum browser dialog.

        Args:
            parent: Parent widget
            batch_results: BatchResults object containing series processing results
            series_processor: MultiSpectrumProcessor instance
            original_data_folder: Path to original data folder
        """
        super().__init__(
            parent=parent,
            title="Browse Individual Spectra - Quality Control",
            default_size=(1000, 700),
            modal=False  # Non-modal to allow interaction with main window
        )

        self.batch_results = batch_results
        self.series_processor = series_processor
        self.original_data_folder = original_data_folder
        self.selected_spectrum = None

        # Sort state
        self.sort_column = 1  # Default sort by spectrum name
        self.sort_ascending = True

        # Build UI
        self.setup_ui()

        # Populate list
        self.populate_spectrum_list()

        # Center on parent
        if parent:
            self.center_on_parent()
        else:
            self.center_on_screen()

        logger.debug("SpectrumBrowserDialog initialized")

    def setup_ui(self):
        """Setup the dialog user interface."""
        layout = QVBoxLayout()
        layout.setSpacing(SPACING_MD)
        layout.setContentsMargins(SPACING_LG, SPACING_LG, SPACING_LG, SPACING_LG)

        # Title and info row
        title_layout = QHBoxLayout()

        title_label = QLabel("Individual Spectrum Browser")
        title_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_SECTION_LABEL}px;
                font-weight: bold;
                color: {PRIMARY_TEXT};
            }}
        """)
        title_layout.addWidget(title_label)

        title_layout.addStretch()

        # Summary info
        if self.batch_results and hasattr(self.batch_results, 'get_summary'):
            summary = self.batch_results.get_summary()
            info_text = f"Series: {summary.get('total_spectra', 0)} spectra, "
            info_text += f"{summary.get('successful', 0)} successful "
            info_text += f"({summary.get('success_rate', 0):.1f}%)"
            info_label = QLabel(info_text)
            info_label.setStyleSheet(f"""
                QLabel {{
                    font-size: {FONT_SIZE_BODY}px;
                    color: #007AFF;
                }}
            """)
            title_layout.addWidget(info_label)

        layout.addLayout(title_layout)

        # Filter section
        filter_group = self.create_filter_section()
        layout.addWidget(filter_group)

        # Spectrum list section
        list_group = self.create_list_section()
        layout.addWidget(list_group, stretch=1)

        # Button row
        button_layout = self.create_button_row()
        layout.addLayout(button_layout)

        self.setLayout(layout)

    def create_filter_section(self) -> QGroupBox:
        """Create the filter and search section.

        Returns:
            QGroupBox containing filter controls
        """
        group = QGroupBox("Filter & Search")
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

        # Search box
        layout.addWidget(QLabel("Search:"))
        self.search_edit = QLineEdit()
        self.search_edit.setPlaceholderText("Filter by name...")
        self.search_edit.setMaximumWidth(200)
        self.search_edit.textChanged.connect(self.filter_spectra)
        layout.addWidget(self.search_edit)

        layout.addSpacing(SPACING_MD)

        # Quality filter
        layout.addWidget(QLabel("Quality:"))
        self.quality_filter = QComboBox()
        self.quality_filter.addItems(["All", "Excellent", "Good", "Fair", "Poor", "Failed"])
        self.quality_filter.setMaximumWidth(120)
        self.quality_filter.currentTextChanged.connect(self.filter_spectra)
        layout.addWidget(self.quality_filter)

        layout.addSpacing(SPACING_MD)

        # Status filter
        layout.addWidget(QLabel("Status:"))
        self.status_filter = QComboBox()
        self.status_filter.addItems(["All", "Success", "Failed"])
        self.status_filter.setMaximumWidth(100)
        self.status_filter.currentTextChanged.connect(self.filter_spectra)
        layout.addWidget(self.status_filter)

        layout.addStretch()

        group.setLayout(layout)
        return group

    def create_list_section(self) -> QGroupBox:
        """Create the spectrum list section.

        Returns:
            QGroupBox containing spectrum tree widget
        """
        group = QGroupBox("Spectrum List")
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

        layout = QVBoxLayout()
        layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)

        # Create tree widget
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels([
            "Status", "Spectrum", "Type", "Detected", "Total",
            "Detection Rate", "Quality", "Processing Time"
        ])

        # Configure columns
        header = self.tree.header()
        header.setSectionResizeMode(0, QHeaderView.ResizeToContents)  # Status
        header.setSectionResizeMode(1, QHeaderView.Stretch)  # Spectrum name
        header.setSectionResizeMode(2, QHeaderView.ResizeToContents)  # Type
        header.setSectionResizeMode(3, QHeaderView.ResizeToContents)  # Detected
        header.setSectionResizeMode(4, QHeaderView.ResizeToContents)  # Total
        header.setSectionResizeMode(5, QHeaderView.ResizeToContents)  # Detection Rate
        header.setSectionResizeMode(6, QHeaderView.ResizeToContents)  # Quality
        header.setSectionResizeMode(7, QHeaderView.ResizeToContents)  # Processing Time

        # Enable sorting
        self.tree.setSortingEnabled(True)
        self.tree.sortByColumn(1, Qt.AscendingOrder)

        # Connect signals
        self.tree.itemDoubleClicked.connect(self.on_spectrum_double_click)
        self.tree.itemClicked.connect(self.on_spectrum_select)

        # Style
        self.tree.setStyleSheet("""
            QTreeWidget {
                background-color: white;
                border: 1px solid #C7C7CC;
                border-radius: 6px;
            }
            QTreeWidget::item {
                padding: 4px;
            }
            QTreeWidget::item:selected {
                background-color: #007AFF;
                color: white;
            }
        """)

        layout.addWidget(self.tree)

        group.setLayout(layout)
        return group

    def create_button_row(self) -> QHBoxLayout:
        """Create the button row.

        Returns:
            QHBoxLayout containing action buttons
        """
        layout = QHBoxLayout()
        layout.setSpacing(SPACING_SM)

        # Open Spectrum Viewer button (primary)
        open_viewer_btn = QPushButton("Open Spectrum Viewer")
        open_viewer_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        open_viewer_btn.setStyleSheet(f"""
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
        open_viewer_btn.clicked.connect(self.open_spectrum_viewer)
        layout.addWidget(open_viewer_btn)

        # Quick Analysis button
        quick_analysis_btn = QPushButton("Quick Analysis")
        quick_analysis_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        quick_analysis_btn.setStyleSheet(f"""
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
        quick_analysis_btn.clicked.connect(self.show_quick_analysis)
        layout.addWidget(quick_analysis_btn)

        # Refresh button
        refresh_btn = QPushButton("Refresh List")
        refresh_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        refresh_btn.setStyleSheet(f"""
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
        refresh_btn.clicked.connect(self.refresh_list)
        layout.addWidget(refresh_btn)

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

    def populate_spectrum_list(self):
        """Populate the spectrum list with data from batch results.

        Based on v0.9 populate_spectrum_list() (spectrum_browser.py:242-299)
        """
        self.tree.clear()

        if not self.batch_results or not hasattr(self.batch_results, 'results'):
            logger.warning("No batch results to display")
            return

        # Add spectrum data
        for spectrum_name, result in self.batch_results.results.items():
            # Calculate quality score
            quality_score = self.calculate_quality_score(result)
            quality_text = self.get_quality_text(quality_score)

            # Status icon
            status_icon = "✓" if result.get('status') == 'success' else "✗"

            # Spectrum type
            spectrum_type = "REF" if result.get('is_reference', False) else "STD"

            # Processing time
            proc_time = result.get('processing_time', 0)
            proc_time_str = f"{proc_time:.2f}s" if proc_time else "N/A"

            # Detection rate
            detection_rate = result.get('detection_rate', 0)

            # Create item
            item = QTreeWidgetItem([
                status_icon,
                spectrum_name,
                spectrum_type,
                str(result.get('detected_peaks', 0)),
                str(result.get('total_peaks', 0)),
                f"{detection_rate:.1f}%",
                quality_text,
                proc_time_str
            ])

            # Store result data for later access
            item.setData(0, Qt.UserRole, result)
            item.setData(1, Qt.UserRole, spectrum_name)

            # Color coding based on quality and status
            if result.get('status') == 'failed':
                item.setBackground(0, Qt.red)
                for col in range(self.tree.columnCount()):
                    item.setBackground(col, Qt.lightGray)
            elif quality_score >= 90:
                for col in range(self.tree.columnCount()):
                    item.setBackground(col, Qt.green)
            elif quality_score >= 70:
                for col in range(self.tree.columnCount()):
                    item.setBackground(col, Qt.cyan)
            elif quality_score >= 50:
                for col in range(self.tree.columnCount()):
                    item.setBackground(col, Qt.yellow)

            self.tree.addTopLevelItem(item)

        logger.debug(f"Populated spectrum list with {self.tree.topLevelItemCount()} items")

    def calculate_quality_score(self, result: Dict) -> float:
        """Calculate quality score for a spectrum result.

        Based on v0.9 calculate_quality_score() (spectrum_browser.py:301-325)

        Args:
            result: Result dictionary for a spectrum

        Returns:
            Quality score from 0-100
        """
        if result.get('status') == 'failed':
            return 0

        # Weighted scoring based on different metrics
        detection_rate = result.get('detection_rate', 0)
        avg_r_squared = result.get('avg_r_squared', 0)
        avg_snr = result.get('avg_snr', 0)

        # Detection rate: 40% weight
        detection_score = detection_rate * 0.4

        # R-squared: 40% weight (scale 0-1 to 0-100)
        r_squared_score = avg_r_squared * 100 * 0.4

        # SNR: 20% weight (cap at 10 for 100% score)
        snr_score = min(avg_snr / 10, 1.0) * 100 * 0.2

        total_score = detection_score + r_squared_score + snr_score
        return min(total_score, 100)

    def get_quality_text(self, score: float) -> str:
        """Convert quality score to text label.

        Args:
            score: Quality score (0-100)

        Returns:
            Quality label string
        """
        if score >= 90:
            return "Excellent"
        elif score >= 70:
            return "Good"
        elif score >= 50:
            return "Fair"
        elif score > 0:
            return "Poor"
        else:
            return "Failed"

    def filter_spectra(self):
        """Filter spectrum list based on search and filter criteria.

        Based on v0.9 filter_spectra() (spectrum_browser.py:327-380)
        """
        search_text = self.search_edit.text().lower()
        quality_filter = self.quality_filter.currentText()
        status_filter = self.status_filter.currentText()

        # Show/hide items based on filters
        for i in range(self.tree.topLevelItemCount()):
            item = self.tree.topLevelItem(i)
            spectrum_name = item.text(1).lower()
            quality = item.text(6)
            status = "Success" if item.text(0) == "✓" else "Failed"

            # Check all filters
            show = True

            # Search filter
            if search_text and search_text not in spectrum_name:
                show = False

            # Quality filter
            if quality_filter != "All" and quality != quality_filter:
                show = False

            # Status filter
            if status_filter != "All" and status != status_filter:
                show = False

            item.setHidden(not show)

    def on_spectrum_select(self, item: QTreeWidgetItem, column: int):
        """Handle spectrum selection.

        Args:
            item: Selected tree item
            column: Column that was clicked
        """
        spectrum_name = item.data(1, Qt.UserRole)
        self.selected_spectrum = spectrum_name
        logger.debug(f"Selected spectrum: {spectrum_name}")

    def on_spectrum_double_click(self, item: QTreeWidgetItem, column: int):
        """Handle double-click on spectrum - open viewer.

        Args:
            item: Double-clicked tree item
            column: Column that was clicked
        """
        self.selected_spectrum = item.data(1, Qt.UserRole)
        self.open_spectrum_viewer()

    def open_spectrum_viewer(self):
        """Open detailed spectrum viewer for selected spectrum."""
        from .spectrum_viewer_dialog import SpectrumViewerDialog

        if not self.selected_spectrum:
            QMessageBox.warning(
                self,
                "No Selection",
                "Please select a spectrum from the list first."
            )
            return

        # Get result data
        result_data = self.batch_results.results.get(self.selected_spectrum)
        if not result_data:
            QMessageBox.warning(
                self,
                "Data Not Found",
                f"Could not find data for spectrum: {self.selected_spectrum}"
            )
            return

        # Get spectrum file path
        spectrum_file_path = result_data.get('spectrum_file', '')
        if not spectrum_file_path and self.original_data_folder:
            # Try to construct path from folder and name
            import os
            spectrum_file_path = os.path.join(self.original_data_folder, self.selected_spectrum)

        # Open SpectrumViewerDialog
        viewer = SpectrumViewerDialog(
            parent=self,
            spectrum_name=self.selected_spectrum,
            spectrum_file_path=spectrum_file_path,
            result_data=result_data,
            series_processor=self.series_processor
        )
        viewer.show()

        # Emit signal
        self.spectrum_selected.emit(self.selected_spectrum)

    def show_quick_analysis(self):
        """Show quick analysis summary for selected spectrum.

        Based on v0.9 show_quick_analysis() (spectrum_browser.py:430-478)
        """
        if not self.selected_spectrum:
            QMessageBox.warning(
                self,
                "No Selection",
                "Please select a spectrum from the list first."
            )
            return

        # Get result data
        result_data = self.batch_results.results.get(self.selected_spectrum)
        if not result_data:
            return

        # Build analysis summary
        summary = f"Quick Analysis: {self.selected_spectrum}\n"
        summary += "=" * 50 + "\n\n"

        summary += f"Status: {result_data.get('status', 'Unknown')}\n"
        summary += f"Processing Time: {result_data.get('processing_time', 0):.2f}s\n\n"

        summary += "Detection Statistics:\n"
        summary += f"  - Detected Peaks: {result_data.get('detected_peaks', 0)}\n"
        summary += f"  - Total Peaks: {result_data.get('total_peaks', 0)}\n"
        summary += f"  - Detection Rate: {result_data.get('detection_rate', 0):.1f}%\n\n"

        summary += "Quality Metrics:\n"
        summary += f"  - Average R²: {result_data.get('avg_r_squared', 0):.4f}\n"
        summary += f"  - Average SNR: {result_data.get('avg_snr', 0):.2f}\n"

        quality_score = self.calculate_quality_score(result_data)
        summary += f"  - Quality Score: {quality_score:.1f}/100 ({self.get_quality_text(quality_score)})\n"

        QMessageBox.information(self, "Quick Analysis", summary)

    def refresh_list(self):
        """Refresh the spectrum list."""
        self.populate_spectrum_list()
        self.filter_spectra()
        logger.debug("Spectrum list refreshed")
