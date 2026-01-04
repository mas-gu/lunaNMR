"""ABOUTME: Data loading dialog for selecting NMR spectrum and peak list files
ABOUTME: Provides workflow-aware file selection with two modes: peak_list and sn_threshold
"""

import logging
from pathlib import Path

from PySide6.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QFileDialog, QLineEdit, QGroupBox, QMessageBox
)

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


class DataLoadingDialog(BaseDialog):
    """Dialog for loading NMR spectrum and peak list files.

    This dialog provides a workflow-aware interface for selecting data files:

    Workflow Modes:
        - peak_list: Both NMR spectrum and peak list required
        - sn_threshold: Only NMR spectrum required, peak list optional

    Features:
        - File browser with filters for NMR file formats
        - Optional peak list selection
        - Validation before accepting
        - Returns file paths on success

    Example:
        ```python
        dialog = DataLoadingDialog(parent, workflow_mode='peak_list')
        if dialog.exec():
            result = dialog.get_result()
            nmr_file = result['nmr_file']
            peak_file = result['peak_file']
        ```

    Based on original tkinter DataLoadingDialog (main_gui.py:163-432)
    """

    def __init__(self, parent=None, workflow_mode='peak_list',
                 current_nmr_folder=None, current_peak_folder=None):
        """Initialize the data loading dialog.

        Args:
            parent: Parent widget
            workflow_mode: 'peak_list' or 'sn_threshold'
            current_nmr_folder: Initial folder for NMR file browser
            current_peak_folder: Initial folder for peak list browser
        """
        super().__init__(
            parent=parent,
            title="Load Data - lunaNMR v1.0",
            default_size=(700, 400),
            modal=True
        )

        self.workflow_mode = workflow_mode
        self.current_nmr_folder = current_nmr_folder
        self.current_peak_folder = current_peak_folder

        # Store selections
        self.nmr_file = None
        self.peak_file = None

        # Build UI
        self.setup_ui()

        # Center on parent
        if parent:
            self.center_on_parent()
        else:
            self.center_on_screen()

        logger.debug(f"DataLoadingDialog initialized (workflow_mode={workflow_mode})")

    def setup_ui(self):
        """Setup the dialog user interface."""
        layout = QVBoxLayout()
        layout.setSpacing(SPACING_MD)
        layout.setContentsMargins(SPACING_LG, SPACING_LG, SPACING_LG, SPACING_LG)

        # Title label (workflow-aware)
        if self.workflow_mode == 'sn_threshold':
            title_text = "Select NMR Spectrum (Peak List Optional)"
        else:
            title_text = "Select NMR Spectrum and Peak List"

        title_label = QLabel(title_text)
        title_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_SECTION_LABEL}px;
                font-weight: bold;
                color: {PRIMARY_TEXT};
                padding-bottom: {SPACING_SM}px;
            }}
        """)
        layout.addWidget(title_label)

        # NMR File section
        nmr_group = self.create_nmr_section()
        layout.addWidget(nmr_group)

        # Peak File section (workflow-aware)
        peak_group = self.create_peak_section()
        layout.addWidget(peak_group)

        # Status label
        self.status_label = QLabel("")
        self.status_label.setStyleSheet(f"""
            QLabel {{
                font-size: {FONT_SIZE_SMALL}px;
                color: {SECONDARY_TEXT};
                padding: {SPACING_SM}px 0;
            }}
        """)
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)

        # Add stretch to push buttons to bottom
        layout.addStretch()

        # Button row
        button_layout = self.create_button_row()
        layout.addLayout(button_layout)

        self.setLayout(layout)

    def create_nmr_section(self) -> QGroupBox:
        """Create the NMR spectrum file selection section.

        Returns:
            QGroupBox containing NMR file browser controls
        """
        group = QGroupBox("NMR Spectrum File")
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
        layout.setSpacing(SPACING_SM)
        layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)

        # Path display (read-only)
        self.nmr_path_edit = QLineEdit()
        self.nmr_path_edit.setReadOnly(True)
        self.nmr_path_edit.setPlaceholderText("No file selected")
        self.nmr_path_edit.setStyleSheet(f"""
            QLineEdit {{
                font-size: {FONT_SIZE_BODY}px;
                color: {PRIMARY_TEXT};
                background-color: white;
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px;
            }}
        """)
        layout.addWidget(self.nmr_path_edit, stretch=1)

        # Browse button
        nmr_browse_btn = QPushButton("Browse...")
        nmr_browse_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        nmr_browse_btn.setMinimumWidth(100)
        nmr_browse_btn.setStyleSheet(f"""
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
            QPushButton:pressed {{
                background-color: {SECONDARY_BUTTON_BORDER};
            }}
        """)
        nmr_browse_btn.clicked.connect(self.browse_nmr)
        layout.addWidget(nmr_browse_btn)

        group.setLayout(layout)
        return group

    def create_peak_section(self) -> QGroupBox:
        """Create the peak list file selection section.

        Returns:
            QGroupBox containing peak list file browser controls
        """
        # Workflow-aware title
        if self.workflow_mode == 'sn_threshold':
            title = "Peak List File (Optional)"
        else:
            title = "Peak List File"

        group = QGroupBox(title)
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
        layout.setSpacing(SPACING_SM)
        layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)

        # Path display (read-only)
        self.peak_path_edit = QLineEdit()
        self.peak_path_edit.setReadOnly(True)
        self.peak_path_edit.setPlaceholderText("No file selected (optional)" if self.workflow_mode == 'sn_threshold' else "No file selected")
        self.peak_path_edit.setStyleSheet(f"""
            QLineEdit {{
                font-size: {FONT_SIZE_BODY}px;
                color: {PRIMARY_TEXT};
                background-color: white;
                border: 1px solid {SECONDARY_BUTTON_BORDER};
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px;
            }}
        """)
        layout.addWidget(self.peak_path_edit, stretch=1)

        # Browse button
        peak_browse_btn = QPushButton("Browse...")
        peak_browse_btn.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        peak_browse_btn.setMinimumWidth(100)
        peak_browse_btn.setStyleSheet(f"""
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
            QPushButton:pressed {{
                background-color: {SECONDARY_BUTTON_BORDER};
            }}
        """)
        peak_browse_btn.clicked.connect(self.browse_peak)
        layout.addWidget(peak_browse_btn)

        group.setLayout(layout)
        return group

    def create_button_row(self) -> QHBoxLayout:
        """Create the button row with Load and Cancel buttons.

        Returns:
            QHBoxLayout containing centered buttons
        """
        layout = QHBoxLayout()
        layout.setSpacing(SPACING_SM)

        # Add stretch to center buttons
        layout.addStretch()

        # Load button (primary - enabled based on validation)
        self.load_button = QPushButton("✓ Load Selected Data")
        self.load_button.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        self.load_button.setMinimumWidth(180)
        self.load_button.setEnabled(False)  # Initially disabled
        self.load_button.setStyleSheet(f"""
            QPushButton {{
                background-color: {PRIMARY_BUTTON_BG};
                color: {PRIMARY_BUTTON_TEXT};
                border: none;
                border-radius: {BUTTON_CORNER_RADIUS}px;
                padding: {SPACING_SM}px {SPACING_MD}px;
                font-size: {FONT_SIZE_BODY}px;
                font-weight: bold;
            }}
            QPushButton:hover:enabled {{
                background-color: {PRIMARY_BUTTON_HOVER};
            }}
            QPushButton:pressed:enabled {{
                background-color: #3A7CC3;
            }}
            QPushButton:disabled {{
                background-color: #C7C7CC;
                color: #8E8E93;
            }}
        """)
        self.load_button.clicked.connect(self.on_load)
        layout.addWidget(self.load_button)

        # Cancel button (secondary)
        cancel_button = QPushButton("Cancel")
        cancel_button.setMinimumHeight(BUTTON_HEIGHT_DIALOG)
        cancel_button.setMinimumWidth(120)
        cancel_button.setStyleSheet(f"""
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
            QPushButton:pressed {{
                background-color: {SECONDARY_BUTTON_BORDER};
            }}
        """)
        cancel_button.clicked.connect(self.reject)
        layout.addWidget(cancel_button)

        # Add stretch to center buttons
        layout.addStretch()

        return layout

    def browse_nmr(self):
        """Open file dialog to select NMR spectrum file."""
        initial_dir = str(self.current_nmr_folder) if self.current_nmr_folder else ""

        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select NMR Spectrum",
            initial_dir,
            "NMR Files (*.ft2 *.ft *.fid *.pipe *.ucsf 2rr 2ri 2ir 2ii 1r 1i);;Bruker Processed (2rr 2ri 2ir 2ii);;NMRPipe (*.ft2 *.ft *.pipe);;All Files (*)"
        )

        if file_path:
            self.nmr_file = file_path
            self.nmr_path_edit.setText(file_path)

            # Update current folder
            self.current_nmr_folder = str(Path(file_path).parent)

            logger.debug(f"NMR file selected: {file_path}")
            self.update_status()

    def browse_peak(self):
        """Open file dialog to select peak list file."""
        initial_dir = str(self.current_peak_folder) if self.current_peak_folder else ""

        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Peak List",
            initial_dir,
            "Peak Lists (*.csv *.txt *.peaks);;All Files (*)"
        )

        if file_path:
            self.peak_file = file_path
            self.peak_path_edit.setText(file_path)

            # Update current folder
            self.current_peak_folder = str(Path(file_path).parent)

            logger.debug(f"Peak list file selected: {file_path}")
            self.update_status()

    def update_status(self):
        """Update status label and enable/disable load button based on selections.

        Validation logic (workflow-aware):
            - peak_list mode: Requires both NMR file and peak list
            - sn_threshold mode: Requires only NMR file (peak list optional)
        """
        status_parts = []

        # Build status message
        if self.nmr_file:
            nmr_name = Path(self.nmr_file).name
            status_parts.append(f"📊 NMR: {nmr_name}")

        if self.peak_file:
            peak_name = Path(self.peak_file).name
            status_parts.append(f"📌 Peaks: {peak_name}")

        # Update status label
        if status_parts:
            self.status_label.setText(" | ".join(status_parts))
            self.status_label.setStyleSheet(f"""
                QLabel {{
                    font-size: {FONT_SIZE_SMALL}px;
                    color: #007AFF;
                    padding: {SPACING_SM}px 0;
                }}
            """)
        else:
            self.status_label.setText("")

        # Validate and enable/disable load button (workflow-aware)
        if self.workflow_mode == 'sn_threshold':
            # S/N Threshold mode: Only NMR file required
            if self.nmr_file:
                self.load_button.setEnabled(True)
            else:
                self.load_button.setEnabled(False)
        else:
            # Peak List mode: Both files required
            if self.nmr_file and self.peak_file:
                self.load_button.setEnabled(True)
            else:
                self.load_button.setEnabled(False)

    def on_load(self):
        """Handle Load button click - validate and accept dialog."""
        # Final validation before accepting
        if not self.nmr_file:
            QMessageBox.critical(
                self,
                "No File Selected",
                "Please select an NMR spectrum file."
            )
            return

        # Workflow-aware validation
        if self.workflow_mode == 'peak_list' and not self.peak_file:
            QMessageBox.critical(
                self,
                "No Peak List Selected",
                "Please select a peak list file.\n\n"
                "Tip: Switch to 'S/N Threshold' workflow mode if you want "
                "to detect peaks automatically without a reference list."
            )
            return

        # Verify files exist
        if not Path(self.nmr_file).exists():
            QMessageBox.critical(
                self,
                "File Not Found",
                f"NMR spectrum file not found:\n{self.nmr_file}"
            )
            return

        if self.peak_file and not Path(self.peak_file).exists():
            QMessageBox.critical(
                self,
                "File Not Found",
                f"Peak list file not found:\n{self.peak_file}"
            )
            return

        # All validations passed - accept dialog
        logger.info(f"Data loading accepted: NMR={self.nmr_file}, Peaks={self.peak_file}")
        self.accept()

    def get_result(self) -> dict:
        """Get the selected file paths.

        Returns:
            Dictionary with keys:
                - nmr_file: Path to NMR spectrum file
                - peak_file: Path to peak list file (or None)
                - nmr_folder: Current NMR folder path
                - peak_folder: Current peak folder path
        """
        return {
            'nmr_file': self.nmr_file,
            'peak_file': self.peak_file,
            'nmr_folder': self.current_nmr_folder,
            'peak_folder': self.current_peak_folder
        }
