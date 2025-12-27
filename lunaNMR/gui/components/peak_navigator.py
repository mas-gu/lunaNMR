#!/usr/bin/env python3
# ABOUTME: Peak Navigator component for interactive peak list display and navigation
# ABOUTME: Qt/PySide6 port from CustomTkinter version with improved signal/slot architecture

"""
Peak Navigator Widget for lunaNMR

Provides an interactive table for browsing, selecting, and editing peaks.
Supports both reference peaks and detected peaks with quality indicators.

Features:
    - Peak list display with quality color coding
    - Selection with automatic spectrum centering
    - Previous/Next navigation through peaks
    - Peak editing capabilities (add, edit, delete, save)
    - Quality indicators (good fit, warning, failed)
    - Signal-based event handling for integration with main GUI

Author: Guillaume Mas
Date: 2025
"""

from typing import Optional, List, Dict, Any
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTableWidget, QTableWidgetItem,
    QLabel, QComboBox, QPushButton, QHeaderView, QAbstractItemView,
    QMessageBox, QDialog, QDialogButtonBox, QFormLayout, QLineEdit
)
from PySide6.QtCore import Signal, Qt, QSize
from PySide6.QtGui import QColor, QFont

from lunaNMR.gui.styles.design_system import (
    # Colors
    PRIMARY_TEXT, SECONDARY_TEXT, PANEL_BG_COLOR, FRAME_BG_COLOR,
    SUCCESS_GREEN, WARNING_ORANGE, ERROR_RED,
    PRIMARY_BUTTON_BG, PRIMARY_BUTTON_HOVER, PRIMARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BG, SECONDARY_BUTTON_HOVER, SECONDARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BORDER, DESTRUCTIVE_BUTTON_BG, DESTRUCTIVE_BUTTON_HOVER,
    DESTRUCTIVE_BUTTON_TEXT,
    # Spacing
    SPACING_XS, SPACING_SM, SPACING_MD,
    # Sizing
    BUTTON_WIDTH_ICON, BUTTON_WIDTH_STANDARD,
    # Typography
    FONT_SIZE_SECTION_LABEL, FONT_SIZE_BODY, FONT_SIZE_SMALL,
    FONT_WEIGHT_BOLD, FONT_WEIGHT_REGULAR
)


class PeakNavigator(QWidget):
    """Peak Navigator widget for interactive peak list display and selection.

    Displays peaks in a table with quality indicators and provides navigation
    controls. Supports both reference peaks and detected peaks.

    Signals:
        peak_selected(int): Emitted when a peak is selected (peak_id)
        peak_edited(int, dict): Emitted when a peak is edited (peak_id, new_values)
        navigation_requested(str): Emitted for prev/next navigation ('prev' or 'next')
        peak_analysis_requested(str, int): Emitted when analysis button clicked (peak_type, peak_id)
    """

    # Signals for communication with main GUI
    peak_selected = Signal(int)  # peak_id
    peak_edited = Signal(int, dict)  # peak_id, new_values
    navigation_requested = Signal(str)  # 'prev' or 'next'
    peak_analysis_requested = Signal(str, int)  # peak_type, peak_id

    def __init__(self, parent: Optional[QWidget] = None):
        """Initialize the Peak Navigator widget.

        Args:
            parent: Parent widget, or None for top-level widget
        """
        super().__init__(parent)

        # State variables
        self.selected_peak_index: Optional[int] = None
        self.selected_peak_type: str = "reference"
        self.reference_peaks: List[List[Any]] = []
        self.detected_peaks: List[List[Any]] = []
        self.navigation_enabled: bool = False
        self.current_navigation_index: Optional[int] = None

        # Setup UI
        self._setup_ui()

    def _setup_ui(self):
        """Create the peak navigator user interface."""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)
        layout.setSpacing(SPACING_SM)

        # Header section
        self._create_header(layout)

        # Peak table
        self._create_table(layout)

        # Status and legend
        self._create_status_section(layout)

        # Control buttons
        self._create_control_buttons(layout)

        # Edit buttons (for detected peaks)
        self._create_edit_buttons(layout)

    def _create_header(self, layout: QVBoxLayout):
        """Create header with title and peak type selector."""
        header_layout = QHBoxLayout()
        header_layout.setSpacing(SPACING_SM)

        # Title
        title = QLabel("🧭 Peak Navigator")
        title_font = QFont()
        title_font.setPointSize(FONT_SIZE_SECTION_LABEL)
        title_font.setWeight(QFont.Weight.Bold)
        title.setFont(title_font)
        header_layout.addWidget(title)

        header_layout.addStretch()

        # Peak type selector
        type_label = QLabel("Peak List:")
        type_label_font = QFont()
        type_label_font.setPointSize(FONT_SIZE_BODY)
        type_label.setFont(type_label_font)
        header_layout.addWidget(type_label)

        self.peak_type_combo = QComboBox()
        self.peak_type_combo.addItems(["Reference Peaks", "Detected Peaks"])
        self.peak_type_combo.setMinimumWidth(150)
        self.peak_type_combo.currentTextChanged.connect(self._on_peak_type_changed)
        header_layout.addWidget(self.peak_type_combo)

        layout.addLayout(header_layout)

    def _create_table(self, layout: QVBoxLayout):
        """Create the peak table widget."""
        self.table = QTableWidget()
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(["Assignment", "X (1H)", "Y", "R²"])

        # Configure table behavior
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setAlternatingRowColors(True)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)

        # Configure columns
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)  # Assignment
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)  # X
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)  # Y
        header.setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)  # R²

        # Connect selection signal
        self.table.itemSelectionChanged.connect(self._on_peak_selected)

        # Set minimum height
        self.table.setMinimumHeight(400)

        layout.addWidget(self.table)

    def _create_status_section(self, layout: QVBoxLayout):
        """Create status label and quality legend."""
        # Status label
        self.status_label = QLabel("No peaks loaded")
        status_font = QFont()
        status_font.setPointSize(FONT_SIZE_SMALL)
        self.status_label.setFont(status_font)
        self.status_label.setStyleSheet(f"color: {SECONDARY_TEXT};")
        layout.addWidget(self.status_label)

        # Quality legend
        legend = QLabel("Quality: ✅ R²≥0.85  🟠 R²<0.85  ❌ Failed")
        legend_font = QFont()
        legend_font.setPointSize(FONT_SIZE_SMALL)
        legend.setFont(legend_font)
        legend.setStyleSheet(f"color: {PRIMARY_TEXT};")
        layout.addWidget(legend)

    def _create_control_buttons(self, layout: QVBoxLayout):
        """Create navigation and analysis buttons."""
        button_layout = QHBoxLayout()
        button_layout.setSpacing(SPACING_XS)

        # Refresh button
        self.refresh_btn = self._create_button("🔄", BUTTON_WIDTH_ICON, is_secondary=True)
        self.refresh_btn.clicked.connect(self.refresh_peak_list)
        button_layout.addWidget(self.refresh_btn)

        button_layout.addSpacing(SPACING_SM)

        # Previous button
        self.prev_btn = self._create_button("◀", BUTTON_WIDTH_ICON, is_secondary=True)
        self.prev_btn.clicked.connect(self._navigate_previous)
        self.prev_btn.setEnabled(False)
        button_layout.addWidget(self.prev_btn)

        # Analyze button
        self.analyze_btn = self._create_button("🔬", BUTTON_WIDTH_ICON, is_primary=True)
        self.analyze_btn.clicked.connect(self._analyze_selected_peak)
        button_layout.addWidget(self.analyze_btn)

        # Next button
        self.next_btn = self._create_button("▶", BUTTON_WIDTH_ICON, is_secondary=True)
        self.next_btn.clicked.connect(self._navigate_next)
        self.next_btn.setEnabled(False)
        button_layout.addWidget(self.next_btn)

        button_layout.addStretch()

        layout.addLayout(button_layout)

    def _create_edit_buttons(self, layout: QVBoxLayout):
        """Create edit buttons for detected peaks."""
        edit_layout = QHBoxLayout()
        edit_layout.setSpacing(SPACING_XS)

        # Edit button
        self.edit_btn = self._create_button("✏️", BUTTON_WIDTH_ICON, is_secondary=True)
        self.edit_btn.clicked.connect(self._edit_selected_peak)
        edit_layout.addWidget(self.edit_btn)

        # Delete button
        self.delete_btn = self._create_button("🗑️", BUTTON_WIDTH_ICON, is_destructive=True)
        self.delete_btn.clicked.connect(self._delete_selected_peak)
        edit_layout.addWidget(self.delete_btn)

        # Add button
        self.add_btn = self._create_button("➕", BUTTON_WIDTH_ICON, is_primary=True)
        self.add_btn.clicked.connect(self._add_new_peak)
        edit_layout.addWidget(self.add_btn)

        # Save button
        self.save_btn = self._create_button("💾", BUTTON_WIDTH_ICON, is_primary=True)
        self.save_btn.clicked.connect(self._save_peak_changes)
        edit_layout.addWidget(self.save_btn)

        edit_layout.addStretch()

        layout.addLayout(edit_layout)

        # Update button states
        self._update_edit_button_states()

    def _create_button(
        self,
        text: str,
        width: int,
        is_primary: bool = False,
        is_secondary: bool = False,
        is_destructive: bool = False
    ) -> QPushButton:
        """Create a styled button.

        Args:
            text: Button text
            width: Button width in pixels
            is_primary: Whether this is a primary action button
            is_secondary: Whether this is a secondary action button
            is_destructive: Whether this is a destructive action button

        Returns:
            Configured QPushButton
        """
        btn = QPushButton(text)
        btn.setFixedWidth(width)
        btn.setMinimumHeight(28)

        # Apply styling based on button type
        if is_primary:
            btn.setStyleSheet(f"""
                QPushButton {{
                    background-color: {PRIMARY_BUTTON_BG};
                    color: {PRIMARY_BUTTON_TEXT};
                    border: none;
                    border-radius: 10px;
                    padding: 4px;
                    font-weight: bold;
                }}
                QPushButton:hover {{
                    background-color: {PRIMARY_BUTTON_HOVER};
                }}
                QPushButton:disabled {{
                    background-color: #CCCCCC;
                    color: #888888;
                }}
            """)
        elif is_destructive:
            btn.setStyleSheet(f"""
                QPushButton {{
                    background-color: {DESTRUCTIVE_BUTTON_BG};
                    color: {DESTRUCTIVE_BUTTON_TEXT};
                    border: none;
                    border-radius: 10px;
                    padding: 4px;
                    font-weight: bold;
                }}
                QPushButton:hover {{
                    background-color: {DESTRUCTIVE_BUTTON_HOVER};
                }}
                QPushButton:disabled {{
                    background-color: #CCCCCC;
                    color: #888888;
                }}
            """)
        else:  # Secondary
            btn.setStyleSheet(f"""
                QPushButton {{
                    background-color: {SECONDARY_BUTTON_BG};
                    color: {SECONDARY_BUTTON_TEXT};
                    border: 1px solid {SECONDARY_BUTTON_BORDER};
                    border-radius: 10px;
                    padding: 4px;
                }}
                QPushButton:hover {{
                    background-color: {SECONDARY_BUTTON_HOVER};
                }}
                QPushButton:disabled {{
                    background-color: #EEEEEE;
                    color: #AAAAAA;
                    border: 1px solid #DDDDDD;
                }}
            """)

        return btn

    # ============================================================================
    # Event Handlers
    # ============================================================================

    def _on_peak_type_changed(self, text: str):
        """Handle peak type dropdown selection change."""
        self.selected_peak_type = "reference" if "Reference" in text else "detected"
        self.refresh_peak_list()

    def _on_peak_selected(self):
        """Handle peak selection from table."""
        selected_rows = self.table.selectedIndexes()
        if not selected_rows:
            return

        row = selected_rows[0].row()

        # Get peak data
        current_peaks = self._get_current_peak_list()
        if row < 0 or row >= len(current_peaks):
            return

        peak = current_peaks[row]
        if len(peak) < 3:
            return

        # Update state
        self.selected_peak_index = row
        self.current_navigation_index = row

        # Extract coordinates
        assignment = peak[0]
        x_coord = float(peak[1])
        y_coord = float(peak[2])

        # Update status
        self.status_label.setText(f"Selected: {assignment} ({x_coord:.3f}, {y_coord:.1f})")

        # Emit signal for main GUI
        self.peak_selected.emit(row)

        # Update button states
        self._update_edit_button_states()
        self._update_navigation_button_states()

    def _navigate_previous(self):
        """Navigate to previous peak in list."""
        if not self.navigation_enabled:
            return

        current_peaks = self._get_current_peak_list()
        if not current_peaks or len(current_peaks) <= 1:
            return

        # Get current index
        current_index = self.current_navigation_index
        if current_index is None:
            current_index = self.selected_peak_index if self.selected_peak_index is not None else 0

        # Navigate to previous (wrap around)
        new_index = (current_index - 1) % len(current_peaks)
        self._navigate_to_peak_index(new_index)

    def _navigate_next(self):
        """Navigate to next peak in list."""
        if not self.navigation_enabled:
            return

        current_peaks = self._get_current_peak_list()
        if not current_peaks or len(current_peaks) <= 1:
            return

        # Get current index
        current_index = self.current_navigation_index
        if current_index is None:
            current_index = self.selected_peak_index if self.selected_peak_index is not None else 0

        # Navigate to next (wrap around)
        new_index = (current_index + 1) % len(current_peaks)
        self._navigate_to_peak_index(new_index)

    def _navigate_to_peak_index(self, peak_index: int):
        """Navigate to specific peak index and trigger analysis.

        Args:
            peak_index: Index of peak to navigate to
        """
        current_peaks = self._get_current_peak_list()
        if not current_peaks or peak_index < 0 or peak_index >= len(current_peaks):
            return

        # Update navigation state
        self.current_navigation_index = peak_index
        self.selected_peak_index = peak_index

        # Select the row in table
        self.table.selectRow(peak_index)
        self.table.scrollToItem(self.table.item(peak_index, 0))

        # Trigger analysis
        self._analyze_selected_peak()

    def _analyze_selected_peak(self):
        """Emit signal to request analysis of selected peak."""
        if self.selected_peak_index is None:
            QMessageBox.information(
                self,
                "No Selection",
                "Please select a peak from the table first."
            )
            return

        # Emit signal with peak type and index
        self.peak_analysis_requested.emit(self.selected_peak_type, self.selected_peak_index)

    def _edit_selected_peak(self):
        """Open dialog to edit selected peak."""
        if self.selected_peak_type != "detected":
            QMessageBox.information(
                self,
                "Cannot Edit",
                "You can only edit detected peaks, not reference peaks."
            )
            return

        if self.selected_peak_index is None:
            QMessageBox.information(
                self,
                "No Selection",
                "Please select a peak to edit."
            )
            return

        # Get current peak data
        peak = self.detected_peaks[self.selected_peak_index]
        assignment = str(peak[0])
        x_coord = float(peak[1])
        y_coord = float(peak[2])

        # Open edit dialog
        dialog = PeakEditDialog(self, assignment, x_coord, y_coord)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            new_assignment, new_x, new_y = dialog.get_values()

            # Update peak data
            self.detected_peaks[self.selected_peak_index][0] = new_assignment
            self.detected_peaks[self.selected_peak_index][1] = new_x
            self.detected_peaks[self.selected_peak_index][2] = new_y

            # Refresh display
            self.refresh_peak_list()

            # Emit signal
            new_values = {
                'assignment': new_assignment,
                'x_coord': new_x,
                'y_coord': new_y
            }
            self.peak_edited.emit(self.selected_peak_index, new_values)

    def _delete_selected_peak(self):
        """Delete selected peak from detected peaks list."""
        if self.selected_peak_type != "detected":
            QMessageBox.information(
                self,
                "Cannot Delete",
                "You can only delete detected peaks, not reference peaks."
            )
            return

        if self.selected_peak_index is None:
            QMessageBox.information(
                self,
                "No Selection",
                "Please select a peak to delete."
            )
            return

        # Get peak info for confirmation
        peak = self.detected_peaks[self.selected_peak_index]
        assignment = str(peak[0])

        # Confirm deletion
        reply = QMessageBox.question(
            self,
            "Confirm Deletion",
            f"Are you sure you want to delete peak '{assignment}'?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No
        )

        if reply == QMessageBox.StandardButton.Yes:
            # Remove from list
            del self.detected_peaks[self.selected_peak_index]
            self.selected_peak_index = None

            # Refresh display
            self.refresh_peak_list()

    def _add_new_peak(self):
        """Open dialog to add new peak."""
        if self.selected_peak_type != "detected":
            QMessageBox.information(
                self,
                "Cannot Add",
                "You can only add peaks to detected peaks list."
            )
            return

        # Open add dialog with default values
        new_id = len(self.detected_peaks) + 1
        dialog = PeakEditDialog(self, f"Det{new_id}", 8.0, 120.0)
        if dialog.exec() == QDialog.DialogCode.Accepted:
            new_assignment, new_x, new_y = dialog.get_values()

            # Add to detected peaks list
            self.detected_peaks.append([new_assignment, new_x, new_y, ""])

            # Refresh display
            self.refresh_peak_list()

    def _save_peak_changes(self):
        """Save detected peaks to file.

        Export detected peaks to text file in standard peak list format.
        Based on v0.9 gui_components.py:1564-1658 (save_peak_changes).
        """
        if self.selected_peak_type != "detected":
            QMessageBox.information(
                self,
                "Nothing to Save",
                "Only detected peaks can be saved."
            )
            return

        if not self.detected_peaks:
            QMessageBox.warning(
                self,
                "No Data",
                "No detected peaks to save."
            )
            return

        try:
            from datetime import datetime
            from PySide6.QtWidgets import QFileDialog

            # Generate default filename with timestamp
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            default_filename = f"detected_peaks_{timestamp}.txt"

            # Open file dialog for saving
            filename, _ = QFileDialog.getSaveFileName(
                self,
                "Export Detected Peak List",
                default_filename,
                "Text files (*.txt);;CSV files (*.csv);;All files (*.*)"
            )

            if not filename:
                return  # User cancelled

            # Prepare peak data for export
            peak_data = []
            for i, peak in enumerate(self.detected_peaks):
                assignment = peak[0]
                x_coord = peak[1]
                y_coord = peak[2]
                height = peak[3] if len(peak) > 3 else ""

                # Extract numeric height if available (v0.9 line 1604-1619)
                if isinstance(height, (int, float)):
                    height_value = height
                elif isinstance(height, str):
                    # Try to extract numeric value from formatted strings
                    if height and height not in ["", "Failed Fit", "No Signal", "Not Fitted"]:
                        try:
                            height_str = height.split()[0].split('(')[0]
                            height_value = float(height_str)
                        except (ValueError, IndexError):
                            height_value = ""
                    else:
                        height_value = ""
                else:
                    height_value = ""

                peak_data.append({
                    'Assignment': str(assignment),
                    'Position_X': float(x_coord),
                    'Position_Y': float(y_coord),
                    'Height': height_value
                })

            # Write to file in standard peak list format
            with open(filename, 'w', encoding='utf-8') as f:
                # Write header
                f.write("Assignment, Position_X, Position_Y, Height\n")

                # Write peak data
                for peak in peak_data:
                    assignment = peak['Assignment']
                    x = peak['Position_X']
                    y = peak['Position_Y']
                    h = peak['Height']

                    if h and h != "":
                        # Include height if available
                        f.write(f"{assignment}, {x:.6f}, {y:.6f}, {h:.6e}\n")
                    else:
                        # No height available
                        f.write(f"{assignment}, {x:.6f}, {y:.6f}\n")

            logger.info(f"Peak Navigator: Exported {len(peak_data)} peaks to {filename}")
            QMessageBox.information(
                self,
                "Export Successful",
                f"Exported {len(peak_data)} detected peaks to:\n{filename}\n\n"
                "This file can be loaded as a reference peak list in new projects."
            )

        except Exception as e:
            logger.error(f"Peak Navigator: Export error: {e}")
            import traceback
            traceback.print_exc()
            QMessageBox.critical(
                self,
                "Export Error",
                f"Failed to export peaks:\n{str(e)}"
            )

    # ============================================================================
    # Public API
    # ============================================================================

    def load_reference_peaks(self, peak_data: Any):
        """Load reference peaks data.

        Args:
            peak_data: Peak data in various formats (DataFrame, list, array)
        """
        self.reference_peaks = []

        if peak_data is None:
            return

        try:
            # Handle pandas DataFrame
            if hasattr(peak_data, 'iloc') and hasattr(peak_data, 'columns'):
                for i, row in peak_data.iterrows():
                    assignment = row.get('Assignment', f"Ref{i+1}")
                    x_coord = float(row.get('Position_X', 0))
                    y_coord = float(row.get('Position_Y', 0))
                    height = row.get('Height', row.get('Intensity', ''))
                    # R² at index 4 for reference peaks (used by refresh_peak_list)
                    r_squared = row.get('R_Squared', None)
                    self.reference_peaks.append([assignment, x_coord, y_coord, height, r_squared])

            # Handle numpy array or list
            # v0.9 format from main_window: [assignment, x_coord, y_coord, height, peak_id]
            elif hasattr(peak_data, 'shape') or isinstance(peak_data, (list, tuple)):
                for i, peak in enumerate(peak_data):
                    if len(peak) >= 3:
                        # Format: [assignment, x_coord, y_coord, height, peak_id]
                        assignment = str(peak[0]) if peak[0] else f"Ref{i+1}"
                        x_coord = float(peak[1])
                        y_coord = float(peak[2])
                        height = peak[3] if len(peak) > 3 else ""
                        self.reference_peaks.append([assignment, x_coord, y_coord, height])
                    elif len(peak) >= 2:
                        # Fallback for simple [x, y] format
                        assignment = f"Ref{i+1}"
                        x_coord = float(peak[0])
                        y_coord = float(peak[1])
                        height = ""
                        self.reference_peaks.append([assignment, x_coord, y_coord, height])

        except Exception as e:
            print(f"Peak Navigator: Error loading reference peaks: {e}")

        # Refresh if currently showing reference peaks
        if self.selected_peak_type == "reference":
            self.refresh_peak_list()

        self._update_navigation_button_states()

    def load_detected_peaks(self, fitted_peaks: Any):
        """Load detected peaks data.

        Args:
            fitted_peaks: Fitted peaks data in various formats (fit results list or DataFrame)
        """
        self.detected_peaks = []

        if fitted_peaks is None:
            return

        try:
            # Handle list of dictionaries (fit results format)
            if isinstance(fitted_peaks, list) and fitted_peaks and isinstance(fitted_peaks[0], dict):
                for i, peak in enumerate(fitted_peaks):
                    # Support multiple key formats for coordinates
                    x_coord = peak.get('center_x', peak.get('ppm_x', peak.get('Position_X', 0)))
                    y_coord = peak.get('center_y', peak.get('ppm_y', peak.get('Position_Y', 0)))
                    assignment = peak.get('assignment', peak.get('Assignment', f"Det{i+1}"))

                    if x_coord != 0 and y_coord != 0:
                        height = peak.get('height', peak.get('intensity', peak.get('amplitude', "")))
                        r_squared = peak.get('r_squared', None)
                        # Format: [assignment, x, y, height, peak_id, r_squared]
                        self.detected_peaks.append([assignment, float(x_coord), float(y_coord), height, i, r_squared])

            # Handle pandas DataFrame
            elif hasattr(fitted_peaks, 'iloc') and hasattr(fitted_peaks, 'columns'):
                for i, row in fitted_peaks.iterrows():
                    assignment = row.get('Assignment', f"Det{i+1}")
                    x_coord = float(row.get('Position_X', 0))
                    y_coord = float(row.get('Position_Y', 0))
                    if x_coord != 0 and y_coord != 0:
                        height = row.get('Height', row.get('Intensity', ''))
                        r_squared = row.get('R_Squared', None)
                        self.detected_peaks.append([assignment, x_coord, y_coord, height, i, r_squared])

        except Exception as e:
            print(f"Peak Navigator: Error loading detected peaks: {e}")

        # Auto-switch to detected peaks if loaded
        if len(self.detected_peaks) > 0:
            self.peak_type_combo.setCurrentText("Detected Peaks")
            self.selected_peak_type = "detected"
            self.refresh_peak_list()

        self._update_navigation_button_states()

    def refresh_peak_list(self):
        """Refresh the peak list display."""
        # Clear table
        self.table.setRowCount(0)

        # Get current peak list
        current_peaks = self._get_current_peak_list()
        peak_type_name = "Reference Peaks" if self.selected_peak_type == "reference" else "Detected Peaks"

        # Populate table
        for i, peak in enumerate(current_peaks):
            if len(peak) < 3:
                continue

            assignment = str(peak[0])
            x_coord = float(peak[1])
            y_coord = float(peak[2])
            height = peak[3] if len(peak) > 3 else ""

            # Extract R² value based on peak type
            # Reference peaks: [assignment, x, y, height, r_squared] - R² at index 4
            # Detected peaks: [assignment, x, y, height, peak_id, r_squared] - R² at index 5
            if self.selected_peak_type == "reference":
                r_squared = peak[4] if len(peak) > 4 else None
            else:
                r_squared = peak[5] if len(peak) > 5 else None

            # Format R² display
            r2_str = f"{r_squared:.3f}" if r_squared is not None else "—"

            # Add quality indicator to assignment (using R² thresholds)
            assignment_display = self._add_quality_indicator(assignment, r_squared, r2_str)

            # Insert row
            row = self.table.rowCount()
            self.table.insertRow(row)

            # Create items
            assignment_item = QTableWidgetItem(assignment_display)
            x_item = QTableWidgetItem(f"{x_coord:.2f}")
            y_item = QTableWidgetItem(f"{y_coord:.1f}")
            r2_item = QTableWidgetItem(r2_str)

            # Center align numeric columns
            x_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            y_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            r2_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)

            # Apply color coding based on R² quality thresholds
            row_color = self._get_quality_color(r_squared, r2_str)
            if row_color:
                for item in [assignment_item, x_item, y_item, r2_item]:
                    item.setBackground(row_color)

            # Add items to table
            self.table.setItem(row, 0, assignment_item)
            self.table.setItem(row, 1, x_item)
            self.table.setItem(row, 2, y_item)
            self.table.setItem(row, 3, r2_item)

        # Update status
        count = len(current_peaks)
        self.status_label.setText(f"{peak_type_name}: {count} peak{'s' if count != 1 else ''}")

        # Update button states
        self._update_edit_button_states()
        self._update_navigation_button_states()

    def update_heights_from_results(self, fitted_results: List[Dict[str, Any]]):
        """Update peak heights from fitting results.

        Args:
            fitted_results: List of fitting result dictionaries
        """
        if not fitted_results:
            return

        # Create mapping of assignments to results
        height_map = {}
        status_map = {}

        for result in fitted_results:
            assignment = str(result.get('assignment', ''))
            height = result.get('volume', result.get('intensity', ''))
            r_squared = result.get('r_squared', 0.0)
            fitted = result.get('fitted', True)

            if assignment:
                height_map[assignment] = height if height else "No Signal"
                status_map[assignment] = {
                    'fitted': fitted,
                    'r_squared': float(r_squared) if r_squared else 0.0,
                    'height_available': bool(height)
                }

        # Update detected peaks
        for i, peak in enumerate(self.detected_peaks):
            assignment = str(peak[0])
            if assignment in height_map:
                height_value = height_map[assignment]
                status = status_map[assignment]
                r_squared = status['r_squared']

                # Format height for display (volume with scientific notation)
                if not status['fitted']:
                    display_height = "Failed Fit"
                elif not status['height_available']:
                    display_height = "No Signal"
                elif isinstance(height_value, (int, float)):
                    display_height = f"{height_value:.2e}"
                else:
                    display_height = height_value

                # Update peak data: [assignment, x, y, height, peak_id?, r_squared]
                # Ensure we have at least 4 elements
                while len(self.detected_peaks[i]) < 4:
                    self.detected_peaks[i].append("")
                self.detected_peaks[i][3] = display_height

                # Store R² as element 5 (index 5, after peak_id at index 4)
                while len(self.detected_peaks[i]) < 6:
                    self.detected_peaks[i].append(None)
                self.detected_peaks[i][5] = r_squared

        # Update reference peaks similarly
        for i, peak in enumerate(self.reference_peaks):
            assignment = str(peak[0])
            if assignment in height_map:
                height_value = height_map[assignment]
                status = status_map[assignment]
                r_squared = status['r_squared']

                # Format height for display (volume with scientific notation)
                if not status['fitted']:
                    display_height = "Failed Fit"
                elif not status['height_available']:
                    display_height = "No Signal"
                elif isinstance(height_value, (int, float)):
                    display_height = f"{height_value:.2e}"
                else:
                    display_height = height_value

                # Update peak data: [assignment, x, y, height, r_squared]
                # Ensure we have at least 4 elements
                while len(self.reference_peaks[i]) < 4:
                    self.reference_peaks[i].append("")
                self.reference_peaks[i][3] = display_height

                # Store R² as element 4 for reference peaks (no peak_id)
                while len(self.reference_peaks[i]) < 5:
                    self.reference_peaks[i].append(None)
                self.reference_peaks[i][4] = r_squared

        # Refresh display
        self.refresh_peak_list()

    # ============================================================================
    # Helper Methods
    # ============================================================================

    def _get_current_peak_list(self) -> List[List[Any]]:
        """Get currently active peak list."""
        return self.reference_peaks if self.selected_peak_type == "reference" else self.detected_peaks

    def _add_quality_indicator(self, assignment: str, r_squared: Optional[float], r2_str: str) -> str:
        """Add quality indicator emoji to assignment based on R² thresholds.

        Unified thresholds (consistent with PeakNavigatorTable):
            - ❌ Red: R² == 0 or failed fit
            - 🟠 Orange: R² < 0.85 (low quality)
            - ✅ Green: R² >= 0.85 (good fit)

        Args:
            assignment: Peak assignment name
            r_squared: R² value from fitting (None if not fitted)
            r2_str: R² display string (for fallback quality detection)

        Returns:
            Assignment string with quality emoji prefix
        """
        # Use R² thresholds if available
        if r_squared is not None:
            if r_squared == 0:
                return f"❌ {assignment}"
            elif r_squared < 0.85:
                return f"🟠 {assignment}"
            else:  # r_squared >= 0.85
                return f"✅ {assignment}"

        # Fallback: no R² available, check if display string exists
        if r2_str and r2_str not in ["", "—"]:
            return f"✅ {assignment}"
        return assignment

    def _get_quality_color(self, r_squared: Optional[float], r2_str: str) -> Optional[QColor]:
        """Get background color based on R² quality thresholds.

        Unified thresholds (consistent with PeakNavigatorTable):
            - Red: R² == 0 or failed fit
            - Orange: R² < 0.85 (low quality)
            - Green: R² >= 0.85 (good fit)

        Args:
            r_squared: R² value from fitting (None if not fitted)
            r2_str: R² display string (for fallback quality detection)

        Returns:
            QColor for row background, or None for default
        """
        # Use R² thresholds if available
        if r_squared is not None:
            if r_squared == 0:
                return QColor(ERROR_RED).lighter(180)
            elif r_squared < 0.85:
                return QColor(WARNING_ORANGE).lighter(180)
            else:  # r_squared >= 0.85
                return QColor(SUCCESS_GREEN).lighter(180)

        # Fallback: no R² available, check if display string exists
        if r2_str and r2_str not in ["", "—"]:
            return QColor(SUCCESS_GREEN).lighter(180)
        return None

    def _update_edit_button_states(self):
        """Update edit button states based on current selection and peak type."""
        is_detected = (self.selected_peak_type == "detected")
        has_selection = (self.selected_peak_index is not None)

        self.edit_btn.setEnabled(is_detected and has_selection)
        self.delete_btn.setEnabled(is_detected and has_selection)
        self.add_btn.setEnabled(is_detected)
        self.save_btn.setEnabled(is_detected)

    def _update_navigation_button_states(self):
        """Update navigation button states based on available peaks."""
        current_peaks = self._get_current_peak_list()
        has_multiple = len(current_peaks) > 1

        self.navigation_enabled = has_multiple
        self.prev_btn.setEnabled(has_multiple)
        self.next_btn.setEnabled(has_multiple)


class PeakEditDialog(QDialog):
    """Dialog for editing peak properties."""

    def __init__(
        self,
        parent: Optional[QWidget],
        assignment: str = "",
        x_coord: float = 0.0,
        y_coord: float = 0.0
    ):
        """Initialize peak edit dialog.

        Args:
            parent: Parent widget
            assignment: Current peak assignment
            x_coord: Current X coordinate
            y_coord: Current Y coordinate
        """
        super().__init__(parent)

        self.setWindowTitle("Edit Peak")
        self.setModal(True)
        self.setMinimumWidth(300)

        # Create form
        layout = QVBoxLayout(self)
        form = QFormLayout()

        # Assignment field
        self.assignment_edit = QLineEdit(assignment)
        form.addRow("Assignment:", self.assignment_edit)

        # X coordinate field
        self.x_edit = QLineEdit(f"{x_coord:.6f}")
        form.addRow("X (1H ppm):", self.x_edit)

        # Y coordinate field
        self.y_edit = QLineEdit(f"{y_coord:.6f}")
        form.addRow("Y (15N/13C ppm):", self.y_edit)

        layout.addLayout(form)

        # Buttons
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def get_values(self) -> tuple[str, float, float]:
        """Get edited values.

        Returns:
            Tuple of (assignment, x_coord, y_coord)
        """
        assignment = self.assignment_edit.text().strip()
        if not assignment:
            assignment = "Peak"

        try:
            x_coord = float(self.x_edit.text())
            y_coord = float(self.y_edit.text())
        except ValueError:
            x_coord = 0.0
            y_coord = 0.0

        return assignment, x_coord, y_coord
