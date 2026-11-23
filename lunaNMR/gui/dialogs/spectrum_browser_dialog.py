# ABOUTME: Qt-based spectrum browser dialog for browsing individual spectra from series results
# ABOUTME: Port of Tkinter spectrum_browser.py to PySide6 for v1.0 Qt interface

import logging
from typing import Dict

from PySide6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QLineEdit, QComboBox, QTreeWidget, QTreeWidgetItem,
    QHeaderView, QMessageBox
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor

from lunaNMR.gui.base.base_dialog import BaseDialog
from lunaNMR.gui.dialogs.spectrum_files import series_data_folder, series_spectrum_file
from lunaNMR.gui.dialogs.spectrum_viewer_dialog import SpectrumViewerDialog
from lunaNMR.processors.multi_spectrum_processor import spectrum_summary
from lunaNMR.gui.styles.design_system import (
    SPACING_SM, SPACING_MD, SPACING_LG,
    FONT_SIZE_BODY, FONT_SIZE_SECTION_LABEL, PRIMARY_TEXT, PRIMARY_BUTTON_BG, PRIMARY_BUTTON_HOVER, PRIMARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BG, SECONDARY_BUTTON_HOVER, SECONDARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BORDER, BUTTON_CORNER_RADIUS, BUTTON_HEIGHT_DIALOG,
    INFO_BLUE, SUCCESS_GREEN, WARNING_ORANGE, ERROR_RED,
    FRAME_BG_COLOR
)

logger = logging.getLogger(__name__)

COLUMNS = ["Status", "Spectrum", "Detected", "Total", "Detection Rate",
           "Mean R²", "Mean S/N", "Spectrum score"]
STATUS, NAME, DETECTED, TOTAL, DETECTION_RATE, MEAN_R_SQUARED, MEAN_SNR, SCORE = range(len(COLUMNS))

# The row tint of each band of the per-spectrum score; a Poor spectrum's row is left plain.
SCORE_BAND_COLORS = {
    "Excellent": SUCCESS_GREEN,
    "Good": INFO_BLUE,
    "Fair": WARNING_ORANGE,
    "Failed": ERROR_RED,
}


class SpectrumRow(QTreeWidgetItem):
    """A row of the spectrum list. A column holding a number is sorted by that number,
    kept in the cell's user data, and not by the text it is shown as."""

    def __lt__(self, other):
        column = self.treeWidget().sortColumn()
        mine, theirs = self.data(column, Qt.UserRole), other.data(column, Qt.UserRole)
        if isinstance(mine, (int, float)) and isinstance(theirs, (int, float)):
            return mine < theirs
        # The text is compared here: QTreeWidgetItem's own comparison, called from
        # Python, dispatches back to this method and never returns.
        return self.text(column) < other.text(column)


class SpectrumBrowserDialog(BaseDialog):
    """Dialog for browsing and selecting individual spectra from series results.

    This is a Qt port of the Tkinter SpectrumBrowserDialog from spectrum_browser.py.

    Features:
        - List of all spectra in series with status, detection, mean R², mean S/N
          and the per-spectrum score computed from those three
        - Search and filter by name, spectrum score, status
        - Sort by any column
        - Double-click to open detailed spectrum viewer
        - Rows tinted by the band of their spectrum score

    The spectrum score rates a spectrum of a series. It is not a fit's quality label,
    which comes from fit_quality_label and is read from each fit.

    Based on v0.9 SpectrumBrowserDialog (spectrum_browser.py:116-478)

    Example:
        ```python
        dialog = SpectrumBrowserDialog(parent, batch_results)
        dialog.show()
        ```
    """

    def __init__(self, parent=None, batch_results=None, original_data_folder=None):
        """Initialize the spectrum browser dialog.

        Args:
            parent: Parent widget
            batch_results: BatchResults object containing series processing results
            original_data_folder: Folder of the spectra, for a series that records none
        """
        super().__init__(
            parent=parent,
            title="Browse Individual Spectra - Quality Control",
            default_size=(1000, 700),
            modal=False  # Non-modal to allow interaction with main window
        )

        self.batch_results = batch_results
        self.original_data_folder = original_data_folder
        self.selected_spectrum = None

        self.setup_ui()

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
                font-size: {FONT_SIZE_SECTION_LABEL}pt;
                font-weight: bold;
                color: {PRIMARY_TEXT};
            }}
        """)
        title_layout.addWidget(title_label)

        title_layout.addStretch()

        # Summary info
        if self.batch_results and hasattr(self.batch_results, 'get_summary'):
            summary = self.batch_results.get_summary()
            series_name = (getattr(self.batch_results, 'metadata', None) or {}).get('series_name')
            info_text = f"Series '{series_name}': " if series_name else "Series: "
            info_text += f"{summary.get('total_spectra', 0)} spectra, "
            info_text += f"{summary.get('successful', 0)} successful "
            info_text += f"({summary.get('success_rate', 0):.1f}%)"
            info_label = QLabel(info_text)
            info_label.setStyleSheet(f"""
                QLabel {{
                    font-size: {FONT_SIZE_BODY}pt;
                    color: {INFO_BLUE};
                }}
            """)
            title_layout.addWidget(info_label)

        layout.addLayout(title_layout)

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
        group = QGroupBox("Filter && Search")
        group.setStyleSheet(f"""
            QGroupBox {{
                font-size: {FONT_SIZE_BODY}pt;
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

        layout.addWidget(QLabel("Search:"))
        self.search_edit = QLineEdit()
        self.search_edit.setPlaceholderText("Filter by name...")
        self.search_edit.setMaximumWidth(200)
        self.search_edit.textChanged.connect(self.filter_spectra)
        layout.addWidget(self.search_edit)

        layout.addSpacing(SPACING_MD)

        # Spectrum score filter
        layout.addWidget(QLabel("Spectrum score:"))
        self.score_filter = QComboBox()
        self.score_filter.addItems(["All", "Excellent", "Good", "Fair", "Poor", "Failed"])
        self.score_filter.setMaximumWidth(120)
        self.score_filter.currentTextChanged.connect(self.filter_spectra)
        layout.addWidget(self.score_filter)

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
                font-size: {FONT_SIZE_BODY}pt;
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

        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(COLUMNS)

        # Every column as wide as its content, so that a spectrum's name is never cut
        self.tree.header().setSectionResizeMode(QHeaderView.ResizeToContents)

        self.tree.setRootIsDecorated(False)
        self.tree.setSortingEnabled(True)
        self.tree.sortByColumn(NAME, Qt.AscendingOrder)

        self.tree.itemDoubleClicked.connect(self.on_spectrum_double_click)
        self.tree.itemClicked.connect(self.on_spectrum_select)

        # Style
        self.tree.setStyleSheet(f"""
            QTreeWidget {{
                background-color: white;
                border: 1px solid #C7C7CC;
                border-radius: 6px;
            }}
            QTreeWidget::item {{
                padding: 4px;
            }}
            QTreeWidget::item:selected {{
                background-color: {PRIMARY_BUTTON_BG};
                color: white;
            }}
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

        open_viewer_btn = QPushButton("Open Spectrum Viewer")
        open_viewer_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        open_viewer_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {PRIMARY_BUTTON_BG};
                color: {PRIMARY_BUTTON_TEXT};
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px {SPACING_MD}px;
                font-size: {FONT_SIZE_BODY}pt;
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
                font-size: {FONT_SIZE_BODY}pt;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
        """)
        quick_analysis_btn.clicked.connect(self.show_quick_analysis)
        layout.addWidget(quick_analysis_btn)

        refresh_btn = QPushButton("Refresh List")
        refresh_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        refresh_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {SECONDARY_BUTTON_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px {SPACING_MD}px;
                font-size: {FONT_SIZE_BODY}pt;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
        """)
        refresh_btn.clicked.connect(self.refresh_list)
        layout.addWidget(refresh_btn)

        layout.addStretch()

        close_btn = QPushButton("Close")
        close_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        close_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {SECONDARY_BUTTON_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px {SPACING_MD}px;
                font-size: {FONT_SIZE_BODY}pt;
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

        for spectrum_name, result in self.batch_results.results.items():
            summary = spectrum_summary(result)
            score = self.spectrum_score(summary)
            band = self.score_band(score)

            # Status icon
            status_icon = "✓" if result.get('status') == 'success' else "✗"

            item = SpectrumRow([
                status_icon,
                spectrum_name,
                str(summary['detected']),
                str(summary['total']),
                f"{summary['detection_rate']:.1%}",
                f"{summary['mean_r_squared']:.4f}",
                f"{summary['mean_snr']:.1f}",
                f"{score:.0f} ({band})"
            ])

            item.setData(NAME, Qt.UserRole, spectrum_name)
            item.setData(DETECTED, Qt.UserRole, summary['detected'])
            item.setData(TOTAL, Qt.UserRole, summary['total'])
            item.setData(DETECTION_RATE, Qt.UserRole, summary['detection_rate'])
            item.setData(MEAN_R_SQUARED, Qt.UserRole, summary['mean_r_squared'])
            item.setData(MEAN_SNR, Qt.UserRole, summary['mean_snr'])
            item.setData(SCORE, Qt.UserRole, score)

            # Row tint by the band of the spectrum score
            if band in SCORE_BAND_COLORS:
                tint = QColor(SCORE_BAND_COLORS[band]).lighter(170)
                for col in range(self.tree.columnCount()):
                    item.setBackground(col, tint)

            self.tree.addTopLevelItem(item)

        logger.debug(f"Populated spectrum list with {self.tree.topLevelItemCount()} items")

    def spectrum_score(self, summary: Dict) -> float:
        """The score of a spectrum of the series, from what the series measured on it.

        Based on v0.9 calculate_quality_score() (spectrum_browser.py:301-325)

        Args:
            summary: spectrum_summary of the spectrum's result

        Returns:
            Spectrum score from 0-100
        """
        # Detection rate: 40% weight (scale 0-1 to 0-100)
        detection_score = summary['detection_rate'] * 100 * 0.4

        # R-squared: 40% weight (scale 0-1 to 0-100)
        r_squared_score = summary['mean_r_squared'] * 100 * 0.4

        # SNR: 20% weight (cap at 10 for 100% score)
        snr_score = min(summary['mean_snr'] / 10, 1.0) * 100 * 0.2

        total_score = detection_score + r_squared_score + snr_score
        return min(total_score, 100)

    def score_band(self, score: float) -> str:
        """The band a spectrum score falls in.

        Args:
            score: Spectrum score (0-100)

        Returns:
            Band name
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
        score_filter = self.score_filter.currentText()
        status_filter = self.status_filter.currentText()

        # Show/hide items based on filters
        for i in range(self.tree.topLevelItemCount()):
            item = self.tree.topLevelItem(i)
            spectrum_name = item.text(NAME).lower()
            band = self.score_band(item.data(SCORE, Qt.UserRole))
            status = "Success" if item.text(STATUS) == "✓" else "Failed"

            show = True

            if search_text and search_text not in spectrum_name:
                show = False

            # Spectrum score filter
            if score_filter != "All" and band != score_filter:
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
        spectrum_name = item.data(NAME, Qt.UserRole)
        self.selected_spectrum = spectrum_name
        logger.debug(f"Selected spectrum: {spectrum_name}")

    def on_spectrum_double_click(self, item: QTreeWidgetItem, column: int):
        """Handle double-click on spectrum - open viewer.

        Args:
            item: Double-clicked tree item
            column: Column that was clicked
        """
        self.selected_spectrum = item.data(NAME, Qt.UserRole)
        self.open_spectrum_viewer()

    def open_spectrum_viewer(self):
        """Open detailed spectrum viewer for selected spectrum."""
        if not self.selected_spectrum:
            QMessageBox.warning(
                self,
                "No Selection",
                "Please select a spectrum from the list first."
            )
            return

        result_data = self.batch_results.results.get(self.selected_spectrum)
        if not result_data:
            QMessageBox.warning(
                self,
                "Data Not Found",
                f"Could not find data for spectrum: {self.selected_spectrum}"
            )
            return

        data_folder = series_data_folder(self.batch_results, self.original_data_folder)
        viewer = SpectrumViewerDialog(
            parent=self,
            spectrum_name=self.selected_spectrum,
            spectrum_file_path=series_spectrum_file(data_folder, self.selected_spectrum, result_data),
            result_data=result_data
        )
        viewer.show()

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

        result_data = self.batch_results.results.get(self.selected_spectrum)
        if not result_data:
            return

        measured = spectrum_summary(result_data)

        summary = f"Quick Analysis: {self.selected_spectrum}\n"
        summary += "=" * 50 + "\n\n"

        summary += f"Status: {result_data.get('status', 'Unknown')}\n\n"

        summary += "Detection Statistics:\n"
        summary += f"  - Detected Peaks: {measured['detected']}\n"
        summary += f"  - Total Peaks: {measured['total']}\n"
        summary += f"  - Detection Rate: {measured['detection_rate']:.1%}\n\n"

        summary += "Fit Statistics:\n"
        summary += f"  - Mean R²: {measured['mean_r_squared']:.4f}\n"
        summary += f"  - Mean S/N: {measured['mean_snr']:.1f}\n\n"

        score = self.spectrum_score(measured)
        summary += f"Spectrum score: {score:.0f} ({self.score_band(score)})\n"

        QMessageBox.information(self, "Quick Analysis", summary)

    def refresh_list(self):
        """Refresh the spectrum list."""
        self.populate_spectrum_list()
        self.filter_spectra()
        logger.debug("Spectrum list refreshed")
