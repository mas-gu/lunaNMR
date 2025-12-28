# ABOUTME: Dialog for resolving missing file references when loading a project
# ABOUTME: Allows browsing to locate files or skipping missing files

from pathlib import Path
from typing import Dict, Set, List
import logging

from PySide6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QFileDialog, QScrollArea, QWidget,
    QFrame
)

from lunaNMR.gui.base.base_dialog import BaseDialog
from lunaNMR.gui.styles.design_system import (
    SPACING_SM, SPACING_MD, SPACING_LG,
    FONT_SIZE_BODY, FONT_SIZE_SECTION_LABEL, FONT_SIZE_SMALL,
    PRIMARY_TEXT, SECONDARY_TEXT, DISABLED_TEXT,
    PRIMARY_BUTTON_BG, PRIMARY_BUTTON_HOVER, PRIMARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BG, SECONDARY_BUTTON_HOVER, SECONDARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BORDER, BUTTON_CORNER_RADIUS, BUTTON_HEIGHT_DIALOG,
    FRAME_BG_COLOR, WARNING_ORANGE
)

logger = logging.getLogger(__name__)


# File type labels for display
FILE_TYPE_LABELS = {
    "nmr_spectrum": "NMR Spectrum",
    "peak_list": "Peak List",
    "series_spectrum": "Series Spectrum",
    "dynamixs_result": "DynamiXs Result",
}


class MissingFilesDialog(BaseDialog):
    """Dialog for resolving missing file references when loading a project.

    Displays a list of missing files and allows the user to:
    - Browse to locate each file
    - Skip individual files
    - Skip all missing files
    - Proceed with loading (using remapped or skipped files)

    Example:
        ```python
        missing = [
            {"path": "/old/path/spectrum.ft", "type": "nmr_spectrum"},
            {"path": "/old/path/peaks.txt", "type": "peak_list"},
        ]
        dialog = MissingFilesDialog(parent, missing_files=missing)
        if dialog.exec():
            remapped = dialog.get_remapped_paths()
            skipped = dialog.get_skipped_files()
        ```
    """

    def __init__(self, parent=None, missing_files: List[Dict] = None):
        """Initialize the missing files dialog.

        Args:
            parent: Parent widget
            missing_files: List of dicts with 'path' and 'type' keys
        """
        self.missing_files = missing_files or []
        self.remapped_paths: Dict[str, str] = {}
        self.skipped_files: Set[str] = set()
        self._file_rows: Dict[str, Dict] = {}  # Track UI elements per file

        super().__init__(
            parent=parent,
            title="Missing Files - lunaNMR",
            default_size=(700, 400),
            modal=True
        )

        self.setup_ui()

        if parent:
            self.center_on_parent()
        else:
            self.center_on_screen()

        logger.debug(f"MissingFilesDialog initialized with {len(self.missing_files)} missing files")

    @staticmethod
    def get_file_type_label(file_type: str) -> str:
        """Get human-readable label for file type.

        Args:
            file_type: Internal file type string

        Returns:
            Human-readable label
        """
        if file_type in FILE_TYPE_LABELS:
            return FILE_TYPE_LABELS[file_type]
        # Convert snake_case to Title Case for unknown types
        return file_type.replace("_", " ").title()

    def setup_ui(self):
        """Setup the dialog user interface."""
        layout = QVBoxLayout()
        layout.setSpacing(SPACING_MD)
        layout.setContentsMargins(SPACING_LG, SPACING_LG, SPACING_LG, SPACING_LG)

        # Warning header
        header_layout = QHBoxLayout()
        warning_icon = QLabel("\u26A0")  # Warning sign unicode
        warning_icon.setStyleSheet(f"""
            QLabel {{
                font-size: 24px;
                color: {WARNING_ORANGE};
                padding-right: {SPACING_SM}px;
            }}
        """)
        header_layout.addWidget(warning_icon)

        header_text = QLabel("Some files referenced by this project could not be found.")
        header_text.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_SECTION_LABEL}px;
                font-weight: bold;
                color: {PRIMARY_TEXT};
            }}
        """)
        header_layout.addWidget(header_text)
        header_layout.addStretch()
        layout.addLayout(header_layout)

        # Instructions
        instructions = QLabel(
            "Use the Browse button to locate each file, or Skip to proceed without it."
        )
        instructions.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_BODY}px;
                color: {SECONDARY_TEXT};
                padding-bottom: {SPACING_SM}px;
            }}
        """)
        layout.addWidget(instructions)

        # Scrollable area for file list
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.Shape.NoFrame)

        scroll_content = QWidget()
        scroll_layout = QVBoxLayout(scroll_content)
        scroll_layout.setSpacing(SPACING_SM)
        scroll_layout.setContentsMargins(0, 0, SPACING_SM, 0)

        # Create a row for each missing file
        for file_info in self.missing_files:
            file_path = file_info.get("path", "")
            file_type = file_info.get("type", "unknown")
            row = self._create_file_row(file_path, file_type)
            scroll_layout.addWidget(row)

        scroll_layout.addStretch()
        scroll_area.setWidget(scroll_content)
        layout.addWidget(scroll_area, stretch=1)

        # Bottom buttons
        button_layout = QHBoxLayout()
        button_layout.setSpacing(SPACING_MD)

        # Skip All button
        skip_all_btn = QPushButton("Skip All")
        skip_all_btn.setStyleSheet(f"""
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
        """)
        skip_all_btn.clicked.connect(self._on_skip_all_clicked)
        button_layout.addWidget(skip_all_btn)

        button_layout.addStretch()

        # Cancel button
        cancel_btn = QPushButton("Cancel")
        cancel_btn.setStyleSheet(f"""
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
        """)
        cancel_btn.clicked.connect(self.reject)
        button_layout.addWidget(cancel_btn)

        # Load Anyway button
        load_btn = QPushButton("Load Anyway")
        load_btn.setStyleSheet(f"""
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
        load_btn.clicked.connect(self.accept)
        button_layout.addWidget(load_btn)

        layout.addLayout(button_layout)
        self.setLayout(layout)

    def _create_file_row(self, file_path: str, file_type: str) -> QFrame:
        """Create a row widget for a missing file.

        Args:
            file_path: Original path of the missing file
            file_type: Type of file (for label)

        Returns:
            QFrame containing the row widgets
        """
        row = QFrame()
        row.setStyleSheet(f"""
            QFrame {{
                background-color: {FRAME_BG_COLOR};
                border-radius: 4px;
                padding: {SPACING_SM}px;
            }}
        """)

        row_layout = QHBoxLayout(row)
        row_layout.setContentsMargins(SPACING_SM, SPACING_SM, SPACING_SM, SPACING_SM)
        row_layout.setSpacing(SPACING_MD)

        # File info section
        info_layout = QVBoxLayout()
        info_layout.setSpacing(2)

        # Type label
        type_label = QLabel(self.get_file_type_label(file_type))
        type_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_SMALL}px;
                color: {DISABLED_TEXT};
                font-weight: bold;
            }}
        """)
        info_layout.addWidget(type_label)

        # Path label (shortened to filename)
        path_obj = Path(file_path)
        path_label = QLabel(path_obj.name)
        path_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_BODY}px;
                color: {PRIMARY_TEXT};
            }}
        """)
        path_label.setToolTip(str(file_path))
        info_layout.addWidget(path_label)

        # Status label (shows remapped path or "Skipped")
        status_label = QLabel("")
        status_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_SMALL}px;
                color: {SECONDARY_TEXT};
                font-style: italic;
            }}
        """)
        info_layout.addWidget(status_label)

        row_layout.addLayout(info_layout, stretch=1)

        # Buttons
        browse_btn = QPushButton("Browse...")
        browse_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: {SECONDARY_BUTTON_BG};
                color: {SECONDARY_BUTTON_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: 4px 12px;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
        """)
        browse_btn.clicked.connect(lambda: self.browse_for_file(file_path))
        row_layout.addWidget(browse_btn)

        skip_btn = QPushButton("Skip")
        skip_btn.setStyleSheet(f"""
            QPushButton {{
                background-color: transparent;
                color: {SECONDARY_TEXT};
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: 4px 12px;
            }}
            QPushButton:hover {{
                background-color: {SECONDARY_BUTTON_HOVER};
            }}
        """)
        skip_btn.clicked.connect(lambda: self._on_skip_clicked(file_path))
        row_layout.addWidget(skip_btn)

        # Store references for updating later
        self._file_rows[file_path] = {
            "row": row,
            "status_label": status_label,
            "browse_btn": browse_btn,
            "skip_btn": skip_btn,
        }

        return row

    def browse_for_file(self, original_path: str):
        """Open file dialog to locate a missing file.

        Args:
            original_path: Original path of the missing file
        """
        path_obj = Path(original_path)
        file_filter = self._get_file_filter(path_obj.suffix)

        new_path, _ = QFileDialog.getOpenFileName(
            self,
            f"Locate {path_obj.name}",
            str(path_obj.parent) if path_obj.parent.exists() else "",
            file_filter
        )

        if new_path:
            self.remap_path(original_path, new_path)
            self._update_file_row_status(original_path, f"Found: {Path(new_path).name}")

    def _get_file_filter(self, suffix: str) -> str:
        """Get file dialog filter based on file suffix.

        Args:
            suffix: File extension including dot

        Returns:
            File filter string for QFileDialog
        """
        suffix_lower = suffix.lower()
        if suffix_lower in [".ft", ".ft2", ".ft3"]:
            return "NMR Files (*.ft *.ft2 *.ft3);;All Files (*)"
        elif suffix_lower in [".txt", ".csv"]:
            return "Text Files (*.txt *.csv);;All Files (*)"
        elif suffix_lower == ".json":
            return "JSON Files (*.json);;All Files (*)"
        else:
            return "All Files (*)"

    def remap_path(self, original_path: str, new_path: str):
        """Record a path remapping.

        Args:
            original_path: Original missing path
            new_path: New path to use instead
        """
        self.remapped_paths[original_path] = new_path
        # Remove from skipped if it was there
        self.skipped_files.discard(original_path)
        logger.debug(f"Remapped path: {original_path} -> {new_path}")

    def skip_file(self, file_path: str):
        """Mark a file as skipped.

        Args:
            file_path: Path of file to skip
        """
        self.skipped_files.add(file_path)
        # Remove from remapped if it was there
        self.remapped_paths.pop(file_path, None)
        logger.debug(f"Skipped file: {file_path}")

    def skip_all(self):
        """Skip all missing files."""
        for file_info in self.missing_files:
            file_path = file_info.get("path", "")
            if file_path and file_path not in self.remapped_paths:
                self.skipped_files.add(file_path)
        logger.debug(f"Skipped all remaining files ({len(self.skipped_files)} total)")

    def _on_skip_clicked(self, file_path: str):
        """Handle skip button click for a file.

        Args:
            file_path: Path of file to skip
        """
        self.skip_file(file_path)
        self._update_file_row_status(file_path, "Skipped")

    def _on_skip_all_clicked(self):
        """Handle skip all button click."""
        self.skip_all()
        # Update all row statuses
        for file_info in self.missing_files:
            file_path = file_info.get("path", "")
            if file_path in self.skipped_files:
                self._update_file_row_status(file_path, "Skipped")

    def _update_file_row_status(self, file_path: str, status_text: str):
        """Update the status label for a file row.

        Args:
            file_path: Path of the file
            status_text: New status text to display
        """
        if file_path in self._file_rows:
            row_info = self._file_rows[file_path]
            row_info["status_label"].setText(status_text)

    def get_remapped_paths(self) -> Dict[str, str]:
        """Get the mapping of original paths to new paths.

        Returns:
            Dict mapping original paths to new paths
        """
        return self.remapped_paths.copy()

    def get_skipped_files(self) -> Set[str]:
        """Get the set of skipped file paths.

        Returns:
            Set of skipped file paths
        """
        return self.skipped_files.copy()
