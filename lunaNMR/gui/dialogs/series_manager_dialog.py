# ABOUTME: Qt dialog for managing saved series integrations (view, rename, delete)
# ABOUTME: Provides a centralized interface for series management operations

import logging
from typing import Optional, Dict, Any

from PySide6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QListWidget, QListWidgetItem, QMessageBox, QInputDialog
)
from PySide6.QtCore import Qt, Signal

from lunaNMR.gui.base.base_dialog import BaseDialog
from lunaNMR.gui.styles.design_system import (
    SPACING_SM, SPACING_MD, SPACING_LG,
    FONT_SIZE_BODY, FONT_SIZE_SECTION_LABEL,
    PRIMARY_BUTTON_BG, PRIMARY_BUTTON_HOVER, PRIMARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BG, SECONDARY_BUTTON_HOVER, SECONDARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BORDER, BUTTON_CORNER_RADIUS, BUTTON_HEIGHT_DIALOG,
    DESTRUCTIVE_BUTTON_BG, DESTRUCTIVE_BUTTON_HOVER, DESTRUCTIVE_BUTTON_TEXT
)

logger = logging.getLogger(__name__)


class SeriesManagerDialog(BaseDialog):
    """Dialog for managing saved series integrations."""

    series_updated = Signal()  # Emitted when series are modified

    def __init__(self, parent=None, main_window=None):
        super().__init__(
            parent=parent,
            title="Series Manager",
            default_size=(400, 350),
            modal=True
        )

        self.main_window = main_window
        self.setup_ui()
        self._populate_list()

    def setup_ui(self):
        """Set up the dialog UI."""
        layout = QVBoxLayout()
        layout.setSpacing(SPACING_MD)
        layout.setContentsMargins(SPACING_LG, SPACING_LG, SPACING_LG, SPACING_LG)

        # Header
        header = QLabel("Saved Series Integrations")
        header.setStyleSheet(f"font-size: {FONT_SIZE_SECTION_LABEL}px; font-weight: bold;")
        layout.addWidget(header)

        # Series list
        self.series_list = QListWidget()
        self.series_list.setStyleSheet(f"""
            QListWidget {{
                font-size: {FONT_SIZE_BODY}px;
                border: 1px solid #ccc;
                border-radius: 4px;
            }}
            QListWidget::item {{
                padding: 8px;
            }}
            QListWidget::item:selected {{
                background-color: {PRIMARY_BUTTON_BG};
                color: white;
            }}
        """)
        self.series_list.itemSelectionChanged.connect(self._on_selection_changed)
        layout.addWidget(self.series_list)

        # Button row
        button_layout = QHBoxLayout()
        button_layout.setSpacing(SPACING_SM)

        # Rename button
        self.rename_btn = QPushButton("Rename")
        self.rename_btn.setEnabled(False)
        self.rename_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {SECONDARY_BUTTON_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: 8px 16px;
                min-height: {BUTTON_HEIGHT_DIALOG}px;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
            QPushButton:disabled {{
                opacity: 0.5;
            }}
        """)
        self.rename_btn.clicked.connect(self._rename_series)
        button_layout.addWidget(self.rename_btn)

        # Delete button
        self.delete_btn = QPushButton("Delete")
        self.delete_btn.setEnabled(False)
        self.delete_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {DESTRUCTIVE_BUTTON_BG};
                color: {DESTRUCTIVE_BUTTON_TEXT};
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: 8px 16px;
                min-height: {BUTTON_HEIGHT_DIALOG}px;
            }}
            QPushButton:hover {{
                background-color: {DESTRUCTIVE_BUTTON_HOVER};
            }}
            QPushButton:disabled {{
                opacity: 0.5;
            }}
        """)
        self.delete_btn.clicked.connect(self._delete_series)
        button_layout.addWidget(self.delete_btn)

        button_layout.addStretch()

        # Close button
        close_btn = QPushButton("Close")
        close_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {PRIMARY_BUTTON_BG};
                color: {PRIMARY_BUTTON_TEXT};
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: 8px 16px;
                min-height: {BUTTON_HEIGHT_DIALOG}px;
            }}
            QPushButton:hover {{
                background-color: {PRIMARY_BUTTON_HOVER};
            }}
        """)
        close_btn.clicked.connect(self.accept)
        button_layout.addWidget(close_btn)

        layout.addLayout(button_layout)

        self.setLayout(layout)

    def _populate_list(self):
        """Populate the series list from main_window.saved_series."""
        self.series_list.clear()

        if not self.main_window:
            return

        saved_series = getattr(self.main_window, 'saved_series', {}) or {}

        if not saved_series:
            item = QListWidgetItem("No series integrations saved")
            item.setFlags(item.flags() & ~Qt.ItemIsSelectable)
            self.series_list.addItem(item)
            return

        for name, batch in saved_series.items():
            count = len(batch.results) if hasattr(batch, 'results') else 0
            item = QListWidgetItem(f"{name}  ({count} spectra)")
            item.setData(Qt.UserRole, name)  # Store actual name
            self.series_list.addItem(item)

    def _on_selection_changed(self):
        """Handle selection change in the list."""
        has_selection = len(self.series_list.selectedItems()) > 0
        selected_item = self.series_list.currentItem()

        # Check if it's a real selectable item
        if selected_item and selected_item.data(Qt.UserRole) is None:
            has_selection = False

        self.rename_btn.setEnabled(has_selection)
        self.delete_btn.setEnabled(has_selection)

    def _get_selected_name(self) -> Optional[str]:
        """Get the name of the currently selected series."""
        selected = self.series_list.currentItem()
        if selected:
            return selected.data(Qt.UserRole)
        return None

    def _rename_series(self):
        """Rename the selected series."""
        old_name = self._get_selected_name()
        if not old_name:
            return

        saved_series = getattr(self.main_window, 'saved_series', {})

        new_name, ok = QInputDialog.getText(
            self,
            "Rename Series",
            f"Enter new name for '{old_name}':",
            text=old_name
        )

        if not ok or not new_name.strip():
            return

        new_name = new_name.strip()

        # Validate name
        invalid_chars = '<>:"/\\|?*'
        for char in invalid_chars:
            new_name = new_name.replace(char, '_')

        # Check for duplicate
        if new_name != old_name and new_name in saved_series:
            QMessageBox.warning(
                self,
                "Name Exists",
                f"A series named '{new_name}' already exists."
            )
            return

        # Perform rename
        batch = saved_series.pop(old_name)
        if hasattr(batch, 'metadata'):
            batch.metadata['series_name'] = new_name
        saved_series[new_name] = batch

        self._populate_list()
        self.series_updated.emit()
        logger.info(f"Renamed series '{old_name}' to '{new_name}'")

    def _delete_series(self):
        """Delete the selected series."""
        name = self._get_selected_name()
        if not name:
            return

        result = QMessageBox.question(
            self,
            "Delete Series",
            f"Are you sure you want to delete the series '{name}'?\n\n"
            "This action cannot be undone.",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )

        if result != QMessageBox.Yes:
            return

        saved_series = getattr(self.main_window, 'saved_series', {})
        if name in saved_series:
            del saved_series[name]

            # If deleted series was current batch_results, clear it
            if hasattr(self.main_window, 'batch_results'):
                current = self.main_window.batch_results
                if current and hasattr(current, 'metadata'):
                    if current.metadata.get('series_name') == name:
                        self.main_window.batch_results = None

            self._populate_list()
            self.series_updated.emit()
            logger.info(f"Deleted series '{name}'")
