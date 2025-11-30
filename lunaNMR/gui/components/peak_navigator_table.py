# ABOUTME: Read-only table-based peak navigator for multi-spectrum viewer
# ABOUTME: Simplified version of PeakNavigator with navigation but no edit features

"""
Peak Navigator Table Widget

A simplified, read-only table widget for displaying and navigating peaks.
Used in multi-spectrum overlay viewer for consistent UI with main GUI.

Features:
    - QTableWidget with columns: Assignment, X (1H), Y, R²
    - Quality color coding (green/orange/red rows)
    - Quality emoji indicators
    - Prev/Next navigation buttons
    - peak_selected signal for integration
"""

from typing import Optional, List, Dict, Any
from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTableWidget, QTableWidgetItem,
    QLabel, QPushButton, QHeaderView, QAbstractItemView
)
from PySide6.QtCore import Signal, Qt
from PySide6.QtGui import QColor, QFont

from lunaNMR.gui.styles.design_system import (
    PRIMARY_TEXT, SECONDARY_TEXT, FRAME_BG_COLOR,
    SUCCESS_GREEN, WARNING_ORANGE, ERROR_RED,
    PRIMARY_BUTTON_BG, PRIMARY_BUTTON_HOVER, PRIMARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BG, SECONDARY_BUTTON_HOVER, SECONDARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BORDER,
    SPACING_XS, SPACING_SM, SPACING_MD,
    BUTTON_CORNER_RADIUS,
    FONT_SIZE_BODY, FONT_SIZE_SMALL
)


class PeakNavigatorTable(QWidget):
    """Read-only table-based peak navigator for multi-spectrum viewer.

    Displays peaks in a table with quality indicators and provides navigation.
    This is a simplified version of PeakNavigator without edit capabilities.

    Signals:
        peak_selected(int): Emitted when a peak is selected (peak index)
    """

    peak_selected = Signal(int)

    def __init__(self, parent: Optional[QWidget] = None):
        """Initialize the Peak Navigator Table widget."""
        super().__init__(parent)

        self.peaks: List[Dict] = []
        self.selected_index: Optional[int] = None

        self._setup_ui()

    def _setup_ui(self):
        """Create the widget UI."""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(SPACING_SM)

        # Peak table
        self._create_table(layout)

        # Status label
        self.status_label = QLabel("No peaks loaded")
        self.status_label.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}px; color: {SECONDARY_TEXT};")
        layout.addWidget(self.status_label)

        # Quality legend
        legend = QLabel("Quality: ✅ R²≥0.85  🟠 R²<0.85  ❌ Failed")
        legend.setStyleSheet(f"font-size: {FONT_SIZE_SMALL}px; color: {PRIMARY_TEXT};")
        layout.addWidget(legend)

        # Navigation buttons
        self._create_nav_buttons(layout)

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
        self.table.itemSelectionChanged.connect(self._on_selection_changed)

        layout.addWidget(self.table)

    def _create_nav_buttons(self, layout: QVBoxLayout):
        """Create compact navigation buttons with symbols only."""
        btn_layout = QHBoxLayout()
        btn_layout.setSpacing(SPACING_XS)

        # Previous button (symbol only)
        self.prev_btn = QPushButton("◀")
        self.prev_btn.setToolTip("Previous peak")
        self.prev_btn.setStyleSheet(self._get_compact_button_style())
        self.prev_btn.clicked.connect(self._navigate_prev)
        self.prev_btn.setEnabled(False)
        btn_layout.addWidget(self.prev_btn)

        # Next button (symbol only)
        self.next_btn = QPushButton("▶")
        self.next_btn.setToolTip("Next peak")
        self.next_btn.setStyleSheet(self._get_compact_button_style())
        self.next_btn.clicked.connect(self._navigate_next)
        self.next_btn.setEnabled(False)
        btn_layout.addWidget(self.next_btn)

        # Analyze button (symbol only)
        self.analyze_btn = QPushButton("📊")
        self.analyze_btn.setToolTip("Analyze peak")
        self.analyze_btn.setStyleSheet(self._get_compact_primary_button_style())
        self.analyze_btn.clicked.connect(self._on_analyze_clicked)
        btn_layout.addWidget(self.analyze_btn)

        layout.addLayout(btn_layout)

    def _get_primary_button_style(self) -> str:
        """Get primary button style."""
        return f"""
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
            QPushButton:disabled {{
                background-color: #cccccc;
                color: #888888;
            }}
        """

    def _get_secondary_button_style(self) -> str:
        """Get secondary button style."""
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
            QPushButton:disabled {{
                background-color: #f0f0f0;
                color: #aaaaaa;
                border-color: #cccccc;
            }}
        """

    def _get_compact_button_style(self) -> str:
        """Get compact button style for symbol-only buttons."""
        return f"""
            QPushButton {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {SECONDARY_BUTTON_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_XS}px {SPACING_SM}px;
                font-size: 16px;
                min-width: 32px;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
            QPushButton:disabled {{
                background-color: #f0f0f0;
                color: #aaaaaa;
                border-color: #cccccc;
            }}
        """

    def _get_compact_primary_button_style(self) -> str:
        """Get compact primary button style for symbol-only buttons."""
        return f"""
            QPushButton {{
                background-color: {PRIMARY_BUTTON_BG};
                color: {PRIMARY_BUTTON_TEXT};
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_XS}px {SPACING_SM}px;
                font-size: 16px;
                min-width: 32px;
            }}
            QPushButton:hover {{
                background-color: {PRIMARY_BUTTON_HOVER};
            }}
            QPushButton:disabled {{
                background-color: #cccccc;
                color: #888888;
            }}
        """

    def load_peaks(self, fitted_peaks: List[Dict]):
        """Load peaks from fitted_peaks list.

        Args:
            fitted_peaks: List of peak dictionaries with keys:
                - assignment: Peak assignment name
                - center_x/peak_x/ppm_x/pos_f2: X coordinate (1H)
                - center_y/peak_y/ppm_y/pos_f1: Y coordinate (15N)
                - r_squared/r2: R² value
                - fitting_quality/quality: Quality rating
        """
        self.peaks = fitted_peaks or []
        self.selected_index = None

        # Clear and populate table
        self.table.setRowCount(0)
        self.table.setRowCount(len(self.peaks))

        for row, peak in enumerate(self.peaks):
            self._populate_row(row, peak)

        # Update status
        self.status_label.setText(f"Detected Peaks: {len(self.peaks)} peaks")

        # Update button states
        has_multiple = len(self.peaks) > 1
        self.prev_btn.setEnabled(has_multiple)
        self.next_btn.setEnabled(has_multiple)

    def _populate_row(self, row: int, peak: Dict):
        """Populate a single table row with peak data."""
        # Extract values with fallbacks
        assignment = peak.get('assignment', peak.get('Assignment', f'Peak_{row+1}'))
        x_ppm = (peak.get('center_x') or peak.get('peak_x') or peak.get('ppm_x') or
                 peak.get('pos_f2') or peak.get('Position_X') or 0)
        y_ppm = (peak.get('center_y') or peak.get('peak_y') or peak.get('ppm_y') or
                 peak.get('pos_f1') or peak.get('Position_Y') or 0)
        r_squared = peak.get('r_squared', peak.get('r2', peak.get('avg_r_squared', 0)))
        quality = peak.get('fitting_quality', peak.get('quality', peak.get('Quality', '')))

        # Determine quality indicator and color
        quality_indicator, row_color = self._get_quality_display(r_squared, quality)

        # Create table items
        assignment_item = QTableWidgetItem(f"{quality_indicator} {assignment}")
        x_item = QTableWidgetItem(f"{float(x_ppm):.3f}")
        y_item = QTableWidgetItem(f"{float(y_ppm):.2f}")
        r2_item = QTableWidgetItem(f"{float(r_squared):.3f}" if r_squared else "—")

        # Center align numeric columns
        x_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        y_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        r2_item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)

        # Apply background color
        if row_color:
            for item in [assignment_item, x_item, y_item, r2_item]:
                item.setBackground(row_color)

        # Set items in table
        self.table.setItem(row, 0, assignment_item)
        self.table.setItem(row, 1, x_item)
        self.table.setItem(row, 2, y_item)
        self.table.setItem(row, 3, r2_item)

    def _get_quality_display(self, r_squared: float, quality: str) -> tuple:
        """Get quality indicator emoji and row color.

        Unified thresholds (consistent with PeakNavigator):
            - ❌ Red: R² == 0 or failed fit
            - 🟠 Orange: R² < 0.85 (low quality)
            - ✅ Green: R² >= 0.85 (good fit)

        Returns:
            Tuple of (emoji, QColor or None)
        """
        # Check quality string first
        quality_lower = str(quality).lower() if quality else ''

        # Failed fit
        if 'failed' in quality_lower or r_squared == 0:
            return "❌", QColor(ERROR_RED).lighter(180)

        # Good fit (R² >= 0.85 or quality string indicates good/excellent)
        if 'excellent' in quality_lower or 'good' in quality_lower or (r_squared and r_squared >= 0.85):
            return "✅", QColor(SUCCESS_GREEN).lighter(180)

        # Low quality (R² < 0.85 or quality string indicates poor/fair)
        if 'poor' in quality_lower or 'fair' in quality_lower or (r_squared and r_squared < 0.85):
            return "🟠", QColor(WARNING_ORANGE).lighter(180)

        # Default for unknown/no data
        if r_squared and r_squared > 0:
            return "✅", QColor(SUCCESS_GREEN).lighter(180)

        return "", None

    def _on_selection_changed(self):
        """Handle table selection change."""
        selected_rows = self.table.selectionModel().selectedRows()
        if selected_rows:
            row = selected_rows[0].row()
            self.selected_index = row
            self.peak_selected.emit(row)

    def _navigate_prev(self):
        """Navigate to previous peak (with wrap-around)."""
        if not self.peaks:
            return

        if self.selected_index is None:
            new_index = len(self.peaks) - 1
        elif self.selected_index <= 0:
            new_index = len(self.peaks) - 1
        else:
            new_index = self.selected_index - 1

        self.select_peak(new_index)

    def _navigate_next(self):
        """Navigate to next peak (with wrap-around)."""
        if not self.peaks:
            return

        if self.selected_index is None:
            new_index = 0
        elif self.selected_index >= len(self.peaks) - 1:
            new_index = 0
        else:
            new_index = self.selected_index + 1

        self.select_peak(new_index)

    def _on_analyze_clicked(self):
        """Handle analyze button click."""
        if self.selected_index is not None:
            self.peak_selected.emit(self.selected_index)

    def select_peak(self, index: int):
        """Programmatically select a peak by index.

        Args:
            index: Peak index to select
        """
        if 0 <= index < len(self.peaks):
            self.table.selectRow(index)
            self.table.scrollToItem(self.table.item(index, 0))

    def get_selected_index(self) -> Optional[int]:
        """Return currently selected peak index."""
        return self.selected_index

    def clear(self):
        """Clear all peaks from the table."""
        self.peaks = []
        self.selected_index = None
        self.table.setRowCount(0)
        self.status_label.setText("No peaks loaded")
        self.prev_btn.setEnabled(False)
        self.next_btn.setEnabled(False)
