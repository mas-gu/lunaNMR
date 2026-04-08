#!/usr/bin/env python3
"""
GUI Components and Style Utilities for DynamiXs v2.0
Based on LunaNMR v1.0 Qt Style Guide

This module provides reusable Qt components and styling utilities
for the DynamiXs NMR analysis GUI, ensuring consistency with the
LunaNMR software suite design language.
"""

from pathlib import Path
from typing import Optional, Callable, List

from PySide6.QtWidgets import (
    QWidget, QFrame, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QPushButton, QLineEdit, QTextEdit, QPlainTextEdit,
    QGroupBox, QScrollArea, QCheckBox, QRadioButton, QButtonGroup,
    QComboBox, QSpinBox, QDoubleSpinBox, QSlider, QProgressBar,
    QFileDialog, QMessageBox, QApplication, QSizePolicy
)
from PySide6.QtCore import Qt, Signal, QSize
from PySide6.QtGui import QFont, QFontDatabase

# Import design constants
from constants import (
    BG_COLOR, PANEL_BG_COLOR, FRAME_BG_COLOR,
    PRIMARY_TEXT, SECONDARY_TEXT, DISABLED_TEXT,
    PRIMARY_BUTTON_BG, PRIMARY_BUTTON_HOVER, PRIMARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BG, SECONDARY_BUTTON_HOVER, SECONDARY_BUTTON_TEXT,
    DESTRUCTIVE_BUTTON_BG, DESTRUCTIVE_BUTTON_HOVER, DESTRUCTIVE_BUTTON_TEXT,
    SUCCESS_GREEN, WARNING_ORANGE, ERROR_RED, INFO_BLUE,
    BORDER_COLOR, SEPARATOR_COLOR,
    BUTTON_CORNER_RADIUS, FRAME_CORNER_RADIUS, CARD_CORNER_RADIUS,
    SPACING_XS, SPACING_SM, SPACING_MD, SPACING_LG, SPACING_XL,
    FONT_SIZE_LARGE_HEADER, FONT_SIZE_MEDIUM_HEADER, FONT_SIZE_SECTION_LABEL,
    FONT_SIZE_BODY, FONT_SIZE_SMALL, FONT_SIZE_TINY, FONT_SIZE_MONOSPACE,
    BUTTON_HEIGHT_STANDARD, BUTTON_WIDTH_STANDARD, BUTTON_WIDTH_LARGE,
    get_system_font
)


# =============================================================================
# STYLESHEET LOADING
# =============================================================================

def load_stylesheet(app: QApplication) -> bool:
    """
    Load the LunaNMR Qt stylesheet.

    Args:
        app: The QApplication instance

    Returns:
        True if stylesheet was loaded successfully, False otherwise
    """
    qss_path = Path(__file__).parent / "styles" / "main.qss"

    if qss_path.exists():
        with open(qss_path, 'r') as f:
            app.setStyleSheet(f.read())
        return True
    else:
        print(f"Warning: Stylesheet not found at {qss_path}")
        return False


# =============================================================================
# FONT UTILITIES
# =============================================================================

def get_font(size: int = FONT_SIZE_BODY, bold: bool = False) -> QFont:
    """
    Create a QFont with the system font family.

    Args:
        size: Font size in points
        bold: Whether to use bold weight

    Returns:
        QFont instance
    """
    font = QFont(get_system_font(), size)
    if bold:
        font.setBold(True)
    return font


def get_monospace_font(size: int = FONT_SIZE_MONOSPACE) -> QFont:
    """
    Create a monospace QFont.

    Args:
        size: Font size in points

    Returns:
        QFont instance with monospace family
    """
    return QFont("Courier New", size)


# =============================================================================
# CUSTOM COMPONENTS
# =============================================================================

class CollapsibleGroupBox(QFrame):
    """
    Custom collapsible group box with LunaNMR styling.

    This component creates a frame with a header label and optional
    collapse/expand functionality.
    """

    collapsed_changed = Signal(bool)

    def __init__(
        self,
        title: str = "",
        collapsible: bool = False,
        initially_collapsed: bool = False,
        parent: Optional[QWidget] = None
    ):
        """
        Initialize the collapsible group box.

        Args:
            title: Header text for the group box
            collapsible: Whether the group can be collapsed
            initially_collapsed: Start in collapsed state
            parent: Parent widget
        """
        super().__init__(parent)

        self._collapsed = initially_collapsed
        self._collapsible = collapsible

        # Set up styling
        self.setProperty("class", "panel")

        # Main layout
        self._layout = QVBoxLayout(self)
        self._layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)
        self._layout.setSpacing(SPACING_SM)

        # Header
        if title or collapsible:
            self._header = QFrame()
            self._header.setProperty("class", "")  # Transparent
            header_layout = QHBoxLayout(self._header)
            header_layout.setContentsMargins(0, 0, 0, 0)
            header_layout.setSpacing(SPACING_SM)

            # Collapse button (if collapsible)
            if collapsible:
                self._collapse_btn = QPushButton("▼" if not initially_collapsed else "▶")
                self._collapse_btn.setProperty("class", "icon")
                self._collapse_btn.setFixedSize(24, 24)
                self._collapse_btn.clicked.connect(self._toggle_collapse)
                header_layout.addWidget(self._collapse_btn)

            # Title label
            self._title_label = QLabel(title)
            self._title_label.setProperty("class", "section-label")
            self._title_label.setFont(get_font(FONT_SIZE_SECTION_LABEL, bold=True))
            header_layout.addWidget(self._title_label)
            header_layout.addStretch()

            self._layout.addWidget(self._header)

        # Content container
        self._content = QFrame()
        self._content.setProperty("class", "")  # Transparent
        self._content_layout = QVBoxLayout(self._content)
        self._content_layout.setContentsMargins(0, 0, 0, 0)
        self._content_layout.setSpacing(SPACING_SM)
        self._layout.addWidget(self._content)

        # Apply initial collapsed state
        if initially_collapsed:
            self._content.setVisible(False)

    @property
    def content(self) -> QFrame:
        """Get the content container widget."""
        return self._content

    @property
    def content_layout(self) -> QVBoxLayout:
        """Get the content layout for adding widgets."""
        return self._content_layout

    @property
    def header(self) -> Optional[QFrame]:
        """Get the header widget if it exists."""
        return getattr(self, '_header', None)

    def add_widget(self, widget: QWidget):
        """Add a widget to the content area."""
        self._content_layout.addWidget(widget)

    def add_layout(self, layout):
        """Add a layout to the content area."""
        self._content_layout.addLayout(layout)

    def add_header_widget(self, widget: QWidget):
        """Add a widget to the header area (e.g., MHz input)."""
        if hasattr(self, '_header'):
            self._header.layout().insertWidget(
                self._header.layout().count() - 1,  # Before stretch
                widget
            )

    def set_collapsed(self, collapsed: bool):
        """Set the collapsed state."""
        if self._collapsible and self._collapsed != collapsed:
            self._collapsed = collapsed
            self._content.setVisible(not collapsed)
            if hasattr(self, '_collapse_btn'):
                self._collapse_btn.setText("▶" if collapsed else "▼")
            self.collapsed_changed.emit(collapsed)

    def is_collapsed(self) -> bool:
        """Check if the group is collapsed."""
        return self._collapsed

    def _toggle_collapse(self):
        """Toggle the collapsed state."""
        self.set_collapsed(not self._collapsed)


class StyledGroupBox(QGroupBox):
    """
    Standard QGroupBox with LunaNMR styling applied.
    """

    def __init__(self, title: str = "", parent: Optional[QWidget] = None):
        super().__init__(title, parent)
        self.setFont(get_font(FONT_SIZE_SECTION_LABEL, bold=True))


# =============================================================================
# BUTTON FACTORIES
# =============================================================================

def create_primary_button(
    text: str,
    clicked: Optional[Callable] = None,
    width: Optional[int] = None,
    enabled: bool = True,
    tooltip: str = ""
) -> QPushButton:
    """
    Create a primary action button with LunaNMR styling.

    Args:
        text: Button text
        clicked: Click handler callback
        width: Fixed width (optional)
        enabled: Initial enabled state
        tooltip: Tooltip text

    Returns:
        Styled QPushButton
    """
    btn = QPushButton(text)
    btn.setProperty("class", "primary")
    btn.setFont(get_font(FONT_SIZE_LARGE_HEADER, bold=True))
    btn.setEnabled(enabled)

    if width:
        btn.setFixedWidth(width)

    if tooltip:
        btn.setToolTip(tooltip)

    if clicked:
        btn.clicked.connect(clicked)

    # Force style refresh
    btn.style().unpolish(btn)
    btn.style().polish(btn)

    return btn


def create_secondary_button(
    text: str,
    clicked: Optional[Callable] = None,
    width: Optional[int] = None,
    enabled: bool = True,
    tooltip: str = ""
) -> QPushButton:
    """
    Create a secondary action button with LunaNMR styling.

    Args:
        text: Button text
        clicked: Click handler callback
        width: Fixed width (optional)
        enabled: Initial enabled state
        tooltip: Tooltip text

    Returns:
        Styled QPushButton
    """
    btn = QPushButton(text)
    btn.setProperty("class", "secondary")
    btn.setFont(get_font(FONT_SIZE_BODY))
    btn.setEnabled(enabled)

    if width:
        btn.setFixedWidth(width)

    if tooltip:
        btn.setToolTip(tooltip)

    if clicked:
        btn.clicked.connect(clicked)

    # Force style refresh
    btn.style().unpolish(btn)
    btn.style().polish(btn)

    return btn


def create_destructive_button(
    text: str,
    clicked: Optional[Callable] = None,
    width: Optional[int] = None,
    enabled: bool = True,
    tooltip: str = ""
) -> QPushButton:
    """
    Create a destructive action button with LunaNMR styling.

    Args:
        text: Button text
        clicked: Click handler callback
        width: Fixed width (optional)
        enabled: Initial enabled state
        tooltip: Tooltip text

    Returns:
        Styled QPushButton
    """
    btn = QPushButton(text)
    btn.setProperty("class", "destructive")
    btn.setFont(get_font(FONT_SIZE_BODY))
    btn.setEnabled(enabled)

    if width:
        btn.setFixedWidth(width)

    if tooltip:
        btn.setToolTip(tooltip)

    if clicked:
        btn.clicked.connect(clicked)

    # Force style refresh
    btn.style().unpolish(btn)
    btn.style().polish(btn)

    return btn


def create_success_button(
    text: str,
    clicked: Optional[Callable] = None,
    width: Optional[int] = None,
    enabled: bool = True,
    tooltip: str = ""
) -> QPushButton:
    """
    Create a success action button with LunaNMR styling.

    Args:
        text: Button text
        clicked: Click handler callback
        width: Fixed width (optional)
        enabled: Initial enabled state
        tooltip: Tooltip text

    Returns:
        Styled QPushButton
    """
    btn = QPushButton(text)
    btn.setProperty("class", "success")
    btn.setFont(get_font(FONT_SIZE_LARGE_HEADER, bold=True))
    btn.setEnabled(enabled)

    if width:
        btn.setFixedWidth(width)

    if tooltip:
        btn.setToolTip(tooltip)

    if clicked:
        btn.clicked.connect(clicked)

    # Force style refresh
    btn.style().unpolish(btn)
    btn.style().polish(btn)

    return btn


def create_icon_button(
    icon_text: str,
    clicked: Optional[Callable] = None,
    size: int = 36,
    tooltip: str = ""
) -> QPushButton:
    """
    Create an icon-only button.

    Args:
        icon_text: Icon character or emoji
        clicked: Click handler callback
        size: Button size (square)
        tooltip: Tooltip text

    Returns:
        Styled QPushButton
    """
    btn = QPushButton(icon_text)
    btn.setProperty("class", "icon")
    btn.setFixedSize(size, size)

    if tooltip:
        btn.setToolTip(tooltip)

    if clicked:
        btn.clicked.connect(clicked)

    # Force style refresh
    btn.style().unpolish(btn)
    btn.style().polish(btn)

    return btn


# =============================================================================
# LABEL FACTORIES
# =============================================================================

def create_label(
    text: str,
    style_class: str = "",
    font_size: int = FONT_SIZE_BODY,
    bold: bool = False,
    color: str = ""
) -> QLabel:
    """
    Create a label with LunaNMR styling.

    Args:
        text: Label text
        style_class: CSS class name (header-large, header-medium, section-label, secondary)
        font_size: Font size in points
        bold: Whether to use bold font
        color: Custom text color (hex)

    Returns:
        Styled QLabel
    """
    label = QLabel(text)
    label.setFont(get_font(font_size, bold))

    if style_class:
        label.setProperty("class", style_class)
        # Force style refresh
        label.style().unpolish(label)
        label.style().polish(label)

    if color:
        label.setStyleSheet(f"color: {color};")

    return label


def create_header_label(text: str, level: int = 1) -> QLabel:
    """
    Create a header label.

    Args:
        text: Header text
        level: Header level (1=large, 2=medium, 3=section)

    Returns:
        Styled QLabel
    """
    if level == 1:
        return create_label(text, "header-large", FONT_SIZE_LARGE_HEADER, bold=True)
    elif level == 2:
        return create_label(text, "header-medium", FONT_SIZE_MEDIUM_HEADER, bold=True)
    else:
        return create_label(text, "section-label", FONT_SIZE_SECTION_LABEL, bold=True)


def create_secondary_label(text: str) -> QLabel:
    """
    Create a secondary/help text label.

    Args:
        text: Label text

    Returns:
        Styled QLabel with secondary color
    """
    return create_label(text, "secondary", FONT_SIZE_SMALL)


# =============================================================================
# INPUT FACTORIES
# =============================================================================

def create_line_edit(
    text: str = "",
    placeholder: str = "",
    width: Optional[int] = None,
    read_only: bool = False,
    enabled: bool = True
) -> QLineEdit:
    """
    Create a line edit with LunaNMR styling.

    Args:
        text: Initial text
        placeholder: Placeholder text
        width: Fixed width (optional)
        read_only: Whether the field is read-only
        enabled: Initial enabled state

    Returns:
        Styled QLineEdit
    """
    edit = QLineEdit(text)
    edit.setFont(get_font(FONT_SIZE_BODY))
    edit.setEnabled(enabled)
    edit.setReadOnly(read_only)

    if placeholder:
        edit.setPlaceholderText(placeholder)

    if width:
        edit.setFixedWidth(width)

    return edit


def create_text_edit(
    text: str = "",
    placeholder: str = "",
    read_only: bool = False,
    monospace: bool = False
) -> QTextEdit:
    """
    Create a multi-line text edit with LunaNMR styling.

    Args:
        text: Initial text
        placeholder: Placeholder text
        read_only: Whether the field is read-only
        monospace: Use monospace font

    Returns:
        Styled QTextEdit
    """
    edit = QTextEdit()
    edit.setText(text)
    edit.setReadOnly(read_only)

    if monospace:
        edit.setFont(get_monospace_font())
    else:
        edit.setFont(get_font(FONT_SIZE_BODY))

    if placeholder:
        edit.setPlaceholderText(placeholder)

    return edit


def create_plain_text_edit(
    text: str = "",
    read_only: bool = False
) -> QPlainTextEdit:
    """
    Create a plain text edit (for logs/output) with LunaNMR styling.

    Args:
        text: Initial text
        read_only: Whether the field is read-only

    Returns:
        Styled QPlainTextEdit with monospace font
    """
    edit = QPlainTextEdit()
    edit.setPlainText(text)
    edit.setReadOnly(read_only)
    edit.setFont(get_monospace_font())

    return edit


def create_spin_box(
    value: int = 0,
    minimum: int = 0,
    maximum: int = 100,
    suffix: str = "",
    width: Optional[int] = None
) -> QSpinBox:
    """
    Create a spin box with LunaNMR styling.

    Args:
        value: Initial value
        minimum: Minimum value
        maximum: Maximum value
        suffix: Suffix text (e.g., " MHz")
        width: Fixed width (optional)

    Returns:
        Styled QSpinBox
    """
    spin = QSpinBox()
    spin.setFont(get_font(FONT_SIZE_BODY))
    spin.setRange(minimum, maximum)
    spin.setValue(value)

    if suffix:
        spin.setSuffix(suffix)

    if width:
        spin.setFixedWidth(width)

    return spin


def create_double_spin_box(
    value: float = 0.0,
    minimum: float = 0.0,
    maximum: float = 100.0,
    decimals: int = 2,
    suffix: str = "",
    width: Optional[int] = None
) -> QDoubleSpinBox:
    """
    Create a double spin box with LunaNMR styling.

    Args:
        value: Initial value
        minimum: Minimum value
        maximum: Maximum value
        decimals: Number of decimal places
        suffix: Suffix text
        width: Fixed width (optional)

    Returns:
        Styled QDoubleSpinBox
    """
    spin = QDoubleSpinBox()
    spin.setFont(get_font(FONT_SIZE_BODY))
    spin.setRange(minimum, maximum)
    spin.setDecimals(decimals)
    spin.setValue(value)

    if suffix:
        spin.setSuffix(suffix)

    if width:
        spin.setFixedWidth(width)

    return spin


def create_combo_box(
    items: List[str] = None,
    current_index: int = 0,
    width: Optional[int] = None
) -> QComboBox:
    """
    Create a combo box with LunaNMR styling.

    Args:
        items: List of items to add
        current_index: Initially selected index
        width: Fixed width (optional)

    Returns:
        Styled QComboBox
    """
    combo = QComboBox()
    combo.setFont(get_font(FONT_SIZE_BODY))

    if items:
        combo.addItems(items)
        combo.setCurrentIndex(current_index)

    if width:
        combo.setFixedWidth(width)

    return combo


def create_check_box(
    text: str = "",
    checked: bool = False,
    enabled: bool = True
) -> QCheckBox:
    """
    Create a check box with LunaNMR styling.

    Args:
        text: Checkbox label text
        checked: Initial checked state
        enabled: Initial enabled state

    Returns:
        Styled QCheckBox
    """
    checkbox = QCheckBox(text)
    checkbox.setFont(get_font(FONT_SIZE_BODY))
    checkbox.setChecked(checked)
    checkbox.setEnabled(enabled)

    return checkbox


def create_radio_button(
    text: str = "",
    checked: bool = False,
    enabled: bool = True
) -> QRadioButton:
    """
    Create a radio button with LunaNMR styling.

    Args:
        text: Radio button label text
        checked: Initial checked state
        enabled: Initial enabled state

    Returns:
        Styled QRadioButton
    """
    radio = QRadioButton(text)
    radio.setFont(get_font(FONT_SIZE_BODY))
    radio.setChecked(checked)
    radio.setEnabled(enabled)

    return radio


# =============================================================================
# LAYOUT UTILITIES
# =============================================================================

def create_h_layout(
    spacing: int = SPACING_SM,
    margins: tuple = (0, 0, 0, 0)
) -> QHBoxLayout:
    """
    Create a horizontal layout with standard spacing.

    Args:
        spacing: Space between widgets
        margins: Layout margins (left, top, right, bottom)

    Returns:
        QHBoxLayout
    """
    layout = QHBoxLayout()
    layout.setSpacing(spacing)
    layout.setContentsMargins(*margins)
    return layout


def create_v_layout(
    spacing: int = SPACING_SM,
    margins: tuple = (0, 0, 0, 0)
) -> QVBoxLayout:
    """
    Create a vertical layout with standard spacing.

    Args:
        spacing: Space between widgets
        margins: Layout margins (left, top, right, bottom)

    Returns:
        QVBoxLayout
    """
    layout = QVBoxLayout()
    layout.setSpacing(spacing)
    layout.setContentsMargins(*margins)
    return layout


def create_grid_layout(
    h_spacing: int = SPACING_SM,
    v_spacing: int = SPACING_SM,
    margins: tuple = (0, 0, 0, 0)
) -> QGridLayout:
    """
    Create a grid layout with standard spacing.

    Args:
        h_spacing: Horizontal space between widgets
        v_spacing: Vertical space between widgets
        margins: Layout margins (left, top, right, bottom)

    Returns:
        QGridLayout
    """
    layout = QGridLayout()
    layout.setHorizontalSpacing(h_spacing)
    layout.setVerticalSpacing(v_spacing)
    layout.setContentsMargins(*margins)
    return layout


def create_scroll_area(widget: QWidget) -> QScrollArea:
    """
    Create a scroll area containing the given widget.

    Args:
        widget: Widget to make scrollable

    Returns:
        QScrollArea containing the widget
    """
    scroll = QScrollArea()
    scroll.setWidget(widget)
    scroll.setWidgetResizable(True)
    scroll.setFrameShape(QFrame.NoFrame)
    return scroll


# =============================================================================
# DIALOG UTILITIES
# =============================================================================

def show_info(parent: QWidget, title: str, message: str):
    """Show an information message box."""
    QMessageBox.information(parent, title, message)


def show_warning(parent: QWidget, title: str, message: str):
    """Show a warning message box."""
    QMessageBox.warning(parent, title, message)


def show_error(parent: QWidget, title: str, message: str):
    """Show an error message box."""
    QMessageBox.critical(parent, title, message)


def ask_yes_no(parent: QWidget, title: str, message: str) -> bool:
    """
    Show a yes/no question dialog.

    Returns:
        True if user clicked Yes, False otherwise
    """
    result = QMessageBox.question(
        parent, title, message,
        QMessageBox.Yes | QMessageBox.No
    )
    return result == QMessageBox.Yes


def open_file_dialog(
    parent: QWidget,
    title: str = "Open File",
    directory: str = "",
    filter: str = "All Files (*)"
) -> str:
    """
    Open a file selection dialog.

    Returns:
        Selected file path or empty string if cancelled
    """
    file_path, _ = QFileDialog.getOpenFileName(parent, title, directory, filter)
    return file_path


def open_files_dialog(
    parent: QWidget,
    title: str = "Open Files",
    directory: str = "",
    filter: str = "All Files (*)"
) -> List[str]:
    """
    Open a multiple file selection dialog.

    Returns:
        List of selected file paths or empty list if cancelled
    """
    file_paths, _ = QFileDialog.getOpenFileNames(parent, title, directory, filter)
    return file_paths


def save_file_dialog(
    parent: QWidget,
    title: str = "Save File",
    directory: str = "",
    filter: str = "All Files (*)"
) -> str:
    """
    Open a save file dialog.

    Returns:
        Selected file path or empty string if cancelled
    """
    file_path, _ = QFileDialog.getSaveFileName(parent, title, directory, filter)
    return file_path


def open_directory_dialog(
    parent: QWidget,
    title: str = "Select Directory",
    directory: str = ""
) -> str:
    """
    Open a directory selection dialog.

    Returns:
        Selected directory path or empty string if cancelled
    """
    return QFileDialog.getExistingDirectory(parent, title, directory)


# =============================================================================
# SPECTRA INPUT WIDGET (LunaNMR Integration)
# =============================================================================

class SpectraInputWidget(QFrame):
    """
    Widget for NMR spectra folder input with auto-detection.

    Supports master folder auto-detect (scans for T1/, T2/, hetNOE_sat/,
    hetNOE_unsat/ subfolders) and individual folder browse overrides.

    Modes:
    - mode='full': Full dual-field support with all experiments (T1, T2, hetNOE) - for Model-Free
    - mode='t1t2': Multi-field support for T1/T2 only (no hetNOE) - for T1/T2 Fitting
    - mode='single': Simple single-folder input - legacy mode

    Legacy parameter single_experiment_mode=True maps to mode='single'.
    """

    folders_changed = Signal()
    peak_list_changed = Signal(str)

    # Expected subfolder names for auto-detection
    EXPERIMENT_FOLDERS = {
        'T1': ['T1', 't1', 'T1_relaxation', 'T1_data'],
        'T2': ['T2', 't2', 'T2_relaxation', 'T2_data'],
        'hetNOE_sat': ['hetNOE_sat', 'NOE_sat', 'noe_sat', 'saturated'],
        'hetNOE_unsat': ['hetNOE_unsat', 'NOE_unsat', 'noe_unsat', 'unsaturated', 'reference'],
    }

    def __init__(self, parent: Optional[QWidget] = None, single_experiment_mode: bool = False,
                 mode: str = None):
        super().__init__(parent)

        # Determine mode
        if mode is not None:
            self._mode = mode
        elif single_experiment_mode:
            self._mode = 'single'
        else:
            self._mode = 'full'

        # Legacy compatibility
        self._single_experiment_mode = (self._mode == 'single')

        self._peak_list_path = ""

        # Field 1 data
        self._field1_folder = ""
        self._field1_experiments = {}  # {experiment_type: folder_path}
        self._field1_delays = {}       # {experiment_type: [delays]}
        self._field1_freq = 600.0

        # Field 2 data (for dual-field)
        self._field2_enabled = False
        self._field2_folder = ""
        self._field2_experiments = {}
        self._field2_delays = {}
        self._field2_freq = 700.0

        self._setup_ui()

    def _setup_ui(self):
        """Set up the widget UI."""
        self.setProperty("class", "panel")

        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(SPACING_MD, SPACING_MD, SPACING_MD, SPACING_MD)
        main_layout.setSpacing(SPACING_MD)

        # --- Peak List Section (shared) ---
        peak_group = CollapsibleGroupBox("Peak List", collapsible=False)

        peak_row = QHBoxLayout()
        peak_row.setSpacing(SPACING_SM)

        self._peak_list_edit = create_line_edit(
            placeholder="Select peak list file (CSV, TXT, or SPARKY peaks)",
            read_only=True
        )
        peak_row.addWidget(self._peak_list_edit, stretch=1)

        self._peak_browse_btn = create_secondary_button("Browse...", self._browse_peak_list)
        peak_row.addWidget(self._peak_browse_btn)

        peak_group.add_layout(peak_row)
        main_layout.addWidget(peak_group)

        if self._mode == 'single':
            # Simple single-folder mode (legacy)
            self._setup_single_mode(main_layout)
        elif self._mode == 't1t2':
            # Multi-field T1/T2 mode (no hetNOE)
            self._setup_t1t2_multi_field_mode(main_layout)
        else:
            # Full dual-field mode for Model-Free analysis
            self._setup_dual_field_mode(main_layout)

        main_layout.addStretch()

    def _setup_single_mode(self, main_layout):
        """Setup simple single-folder mode for T1/T2 fitting."""
        folder_group = CollapsibleGroupBox("Spectra Folder", collapsible=False)

        folder_row = QHBoxLayout()
        folder_row.setSpacing(SPACING_SM)

        self._field1_edit = create_line_edit(
            placeholder="Select folder containing T1 or T2 spectra (with delay in filenames)",
            read_only=True
        )
        folder_row.addWidget(self._field1_edit, stretch=1)

        self._field1_browse_btn = create_secondary_button("Browse...", lambda: self._browse_field_folder(1))
        folder_row.addWidget(self._field1_browse_btn)

        folder_group.add_layout(folder_row)
        main_layout.addWidget(folder_group)

        # Detected spectra display
        self._field1_exp_widgets = {}
        experiments_group = CollapsibleGroupBox("Detected Spectra", collapsible=True)
        self._setup_experiments_grid(experiments_group, self._field1_exp_widgets, ['T1', 'T2'], field=1)
        main_layout.addWidget(experiments_group)

    def _setup_t1t2_multi_field_mode(self, main_layout):
        """Setup multi-field mode for T1/T2 analysis (supports multiple fields, T1/T2 only)."""
        experiment_types = ['T1', 'T2']  # No hetNOE in T1/T2 mode

        # --- Field 1 Section ---
        field1_group = CollapsibleGroupBox("Field 1 Data (Required)", collapsible=False)

        # Frequency input
        freq1_row = QHBoxLayout()
        freq1_row.setSpacing(SPACING_SM)
        freq1_row.addWidget(create_label("Frequency:"))
        self._field1_freq_spin = create_double_spin_box(600.0, 100.0, 1200.0, 1, " MHz", width=100)
        self._field1_freq_spin.valueChanged.connect(lambda v: setattr(self, '_field1_freq', v))
        freq1_row.addWidget(self._field1_freq_spin)
        freq1_row.addStretch()
        field1_group.add_layout(freq1_row)

        # Master folder row
        folder1_row = QHBoxLayout()
        folder1_row.setSpacing(SPACING_SM)
        folder1_row.addWidget(create_label("Spectra Folder:"))
        self._field1_edit = create_line_edit(
            placeholder="Select folder containing T1 or T2 spectra (with delay in filenames)",
            read_only=True
        )
        folder1_row.addWidget(self._field1_edit, stretch=1)
        self._field1_browse_btn = create_secondary_button("Browse...", lambda: self._browse_field_folder(1))
        folder1_row.addWidget(self._field1_browse_btn)
        field1_group.add_layout(folder1_row)

        # Experiments grid for Field 1
        self._field1_exp_widgets = {}
        self._setup_experiments_grid(field1_group, self._field1_exp_widgets, experiment_types, field=1)

        main_layout.addWidget(field1_group)

        # --- Field 2 Section (Optional) ---
        field2_group = CollapsibleGroupBox("Field 2 Data (Optional)", collapsible=True, initially_collapsed=True)

        # Enable toggle and frequency
        enable_row = QHBoxLayout()
        enable_row.setSpacing(SPACING_SM)

        self._field2_enable_check = QCheckBox("Enable Field 2")
        self._field2_enable_check.setStyleSheet(f"color: {PRIMARY_TEXT};")
        self._field2_enable_check.toggled.connect(self._on_field2_toggle)
        enable_row.addWidget(self._field2_enable_check)

        enable_row.addWidget(create_label("Frequency:"))
        self._field2_freq_spin = create_double_spin_box(700.0, 100.0, 1200.0, 1, " MHz", width=100)
        self._field2_freq_spin.valueChanged.connect(lambda v: setattr(self, '_field2_freq', v))
        self._field2_freq_spin.setEnabled(False)
        enable_row.addWidget(self._field2_freq_spin)
        enable_row.addStretch()
        field2_group.add_layout(enable_row)

        # Master folder row for Field 2
        folder2_row = QHBoxLayout()
        folder2_row.setSpacing(SPACING_SM)
        folder2_row.addWidget(create_label("Spectra Folder:"))
        self._field2_edit = create_line_edit(
            placeholder="Select folder containing T1 or T2 spectra (with delay in filenames)",
            read_only=True
        )
        self._field2_edit.setEnabled(False)
        folder2_row.addWidget(self._field2_edit, stretch=1)
        self._field2_browse_btn = create_secondary_button("Browse...", lambda: self._browse_field_folder(2))
        self._field2_browse_btn.setEnabled(False)
        folder2_row.addWidget(self._field2_browse_btn)
        field2_group.add_layout(folder2_row)

        # Experiments grid for Field 2
        self._field2_exp_widgets = {}
        self._setup_experiments_grid(field2_group, self._field2_exp_widgets, experiment_types, field=2, enabled=False)

        main_layout.addWidget(field2_group)

        # --- Additional Fields Section (Optional) ---
        # Note: For more than 2 fields, users can add more via a future enhancement
        more_fields_label = create_secondary_label(
            "Note: For analysis with more than 2 fields, process each field separately and combine results."
        )
        main_layout.addWidget(more_fields_label)

    def _setup_dual_field_mode(self, main_layout):
        """Setup full dual-field mode for Model-Free analysis."""
        experiment_types = ['T1', 'T2', 'hetNOE_sat', 'hetNOE_unsat']

        # --- Field 1 Section ---
        field1_group = CollapsibleGroupBox("Field 1 Data (Required)", collapsible=False)

        # Frequency input in header
        freq1_row = QHBoxLayout()
        freq1_row.setSpacing(SPACING_SM)
        freq1_row.addWidget(create_label("Frequency:"))
        self._field1_freq_spin = create_double_spin_box(600.0, 100.0, 1200.0, 1, " MHz", width=100)
        self._field1_freq_spin.valueChanged.connect(lambda v: setattr(self, '_field1_freq', v))
        freq1_row.addWidget(self._field1_freq_spin)
        freq1_row.addStretch()
        field1_group.add_layout(freq1_row)

        # Master folder row
        folder1_row = QHBoxLayout()
        folder1_row.setSpacing(SPACING_SM)
        folder1_row.addWidget(create_label("Master Folder:"))
        self._field1_edit = create_line_edit(
            placeholder="Select folder containing T1/, T2/, hetNOE_sat/, hetNOE_unsat/",
            read_only=True
        )
        folder1_row.addWidget(self._field1_edit, stretch=1)
        self._field1_browse_btn = create_secondary_button("Browse...", lambda: self._browse_field_folder(1))
        folder1_row.addWidget(self._field1_browse_btn)
        field1_group.add_layout(folder1_row)

        # Experiments grid for Field 1
        self._field1_exp_widgets = {}
        self._setup_experiments_grid(field1_group, self._field1_exp_widgets, experiment_types, field=1)

        main_layout.addWidget(field1_group)

        # --- Field 2 Section (Optional) ---
        field2_group = CollapsibleGroupBox("Field 2 Data (Optional - for dual-field)", collapsible=True, initially_collapsed=True)

        # Enable toggle and frequency in header
        enable_row = QHBoxLayout()
        enable_row.setSpacing(SPACING_SM)

        self._field2_enable_check = QCheckBox("Enable Field 2")
        self._field2_enable_check.setStyleSheet(f"color: {PRIMARY_TEXT};")
        self._field2_enable_check.toggled.connect(self._on_field2_toggle)
        enable_row.addWidget(self._field2_enable_check)

        enable_row.addWidget(create_label("Frequency:"))
        self._field2_freq_spin = create_double_spin_box(700.0, 100.0, 1200.0, 1, " MHz", width=100)
        self._field2_freq_spin.valueChanged.connect(lambda v: setattr(self, '_field2_freq', v))
        self._field2_freq_spin.setEnabled(False)
        enable_row.addWidget(self._field2_freq_spin)
        enable_row.addStretch()
        field2_group.add_layout(enable_row)

        # Master folder row for Field 2
        folder2_row = QHBoxLayout()
        folder2_row.setSpacing(SPACING_SM)
        folder2_row.addWidget(create_label("Master Folder:"))
        self._field2_edit = create_line_edit(
            placeholder="Select folder containing T1/, T2/, hetNOE_sat/, hetNOE_unsat/",
            read_only=True
        )
        self._field2_edit.setEnabled(False)
        folder2_row.addWidget(self._field2_edit, stretch=1)
        self._field2_browse_btn = create_secondary_button("Browse...", lambda: self._browse_field_folder(2))
        self._field2_browse_btn.setEnabled(False)
        folder2_row.addWidget(self._field2_browse_btn)
        field2_group.add_layout(folder2_row)

        # Experiments grid for Field 2
        self._field2_exp_widgets = {}
        self._setup_experiments_grid(field2_group, self._field2_exp_widgets, experiment_types, field=2, enabled=False)

        main_layout.addWidget(field2_group)

    def _setup_experiments_grid(self, parent_group, widgets_dict, experiment_types, field, enabled=True):
        """Setup experiments detection grid for a field."""
        grid = QGridLayout()
        grid.setHorizontalSpacing(SPACING_MD)
        grid.setVerticalSpacing(SPACING_SM)

        # Headers
        grid.addWidget(create_label("Experiment", bold=True), 0, 0)
        grid.addWidget(create_label("Status", bold=True), 0, 1)
        grid.addWidget(create_label("Details", bold=True), 0, 2)
        if not self._single_experiment_mode:
            grid.addWidget(create_label(""), 0, 3)

        row = 1
        for exp_type in experiment_types:
            label = create_label(exp_type)
            grid.addWidget(label, row, 0)

            status_label = create_label("Not found", color=ERROR_RED if enabled else SECONDARY_TEXT)
            grid.addWidget(status_label, row, 1)

            details_label = create_secondary_label("")
            grid.addWidget(details_label, row, 2)

            if not self._single_experiment_mode:
                browse_btn = create_secondary_button("...", width=40)
                browse_btn.setToolTip(f"Browse for {exp_type} folder")
                browse_btn.clicked.connect(
                    lambda checked, et=exp_type, f=field: self._browse_experiment_folder(et, f)
                )
                browse_btn.setEnabled(enabled)
                grid.addWidget(browse_btn, row, 3)
            else:
                browse_btn = None

            widgets_dict[exp_type] = {
                'status': status_label,
                'details': details_label,
                'browse': browse_btn
            }
            row += 1

        parent_group.add_layout(grid)

    def _on_field2_toggle(self, enabled: bool):
        """Handle Field 2 enable/disable toggle."""
        self._field2_enabled = enabled
        self._field2_freq_spin.setEnabled(enabled)
        self._field2_edit.setEnabled(enabled)
        self._field2_browse_btn.setEnabled(enabled)

        for exp_type, widgets in self._field2_exp_widgets.items():
            if widgets['browse']:
                widgets['browse'].setEnabled(enabled)
            if not enabled:
                widgets['status'].setText("Disabled")
                widgets['status'].setStyleSheet(f"color: {SECONDARY_TEXT};")
                widgets['details'].setText("")
            else:
                widgets['status'].setText("Not found")
                widgets['status'].setStyleSheet(f"color: {ERROR_RED};")

        self.folders_changed.emit()

    def _browse_peak_list(self):
        """Browse for peak list file."""
        file_path = open_file_dialog(
            self,
            "Select Peak List",
            "",
            "Peak Lists (*.csv *.txt *.peaks);;CSV Files (*.csv);;Text Files (*.txt);;SPARKY Peaks (*.peaks);;All Files (*)"
        )
        if file_path:
            self._peak_list_path = file_path
            self._peak_list_edit.setText(file_path)
            self.peak_list_changed.emit(file_path)

    def _browse_field_folder(self, field: int):
        """Browse for field master folder and auto-detect experiments."""
        current_folder = self._field1_folder if field == 1 else self._field2_folder
        folder = open_directory_dialog(
            self,
            f"Select Field {field} Master Folder",
            current_folder
        )
        if folder:
            self._set_field_folder(field, folder)

    def _browse_experiment_folder(self, experiment_type: str, field: int):
        """Browse for individual experiment folder."""
        current_folder = self._field1_folder if field == 1 else self._field2_folder
        folder = open_directory_dialog(
            self,
            f"Select {experiment_type} Folder (Field {field})",
            current_folder
        )
        if folder:
            experiments = self._field1_experiments if field == 1 else self._field2_experiments
            experiments[experiment_type] = folder
            widgets = self._field1_exp_widgets if field == 1 else self._field2_exp_widgets
            delays = self._field1_delays if field == 1 else self._field2_delays
            self._update_experiment_display(experiment_type, folder, widgets, delays)
            self.folders_changed.emit()

    def _set_field_folder(self, field: int, folder_path: str):
        """Set master folder for a field and auto-detect experiments."""
        import os

        if field == 1:
            self._field1_folder = folder_path
            self._field1_edit.setText(folder_path)
            self._field1_experiments = {}
            self._field1_delays = {}
            experiments = self._field1_experiments
            delays = self._field1_delays
            widgets = self._field1_exp_widgets
        else:
            self._field2_folder = folder_path
            self._field2_edit.setText(folder_path)
            self._field2_experiments = {}
            self._field2_delays = {}
            experiments = self._field2_experiments
            delays = self._field2_delays
            widgets = self._field2_exp_widgets

        if not os.path.isdir(folder_path):
            return

        if self._single_experiment_mode:
            # In single experiment mode, scan the folder directly for delays
            for exp_type in ['T1', 'T2']:
                experiments[exp_type] = folder_path
                self._update_experiment_display(exp_type, folder_path, widgets, delays)
        else:
            # Full mode - scan for experiment subfolders
            for exp_type, possible_names in self.EXPERIMENT_FOLDERS.items():
                for name in possible_names:
                    potential_path = os.path.join(folder_path, name)
                    if os.path.isdir(potential_path):
                        experiments[exp_type] = potential_path
                        self._update_experiment_display(exp_type, potential_path, widgets, delays)
                        break
                else:
                    # Not found - clear display
                    self._clear_experiment_display(exp_type, widgets, delays)

        self.folders_changed.emit()

    def _update_experiment_display(self, experiment_type: str, folder_path: str,
                                    widgets_dict: dict, delays_dict: dict):
        """Update the display for an experiment type."""
        import os
        import sys

        widgets = widgets_dict[experiment_type]

        # Import delay extractor (add lunaNMR to path if needed)
        try:
            lunaNMR_root = Path(__file__).parent.parent.parent
            if str(lunaNMR_root) not in sys.path:
                sys.path.insert(0, str(lunaNMR_root))
            from lunaNMR.utils.delay_extractor import DelayExtractor

            extractor = DelayExtractor()
            result = extractor.scan_folder(folder_path)

            if result:
                delay_values = [d for _, d in result]
                delays_dict[experiment_type] = delay_values

                widgets['status'].setText("Found")
                widgets['status'].setStyleSheet(f"color: {SUCCESS_GREEN};")

                # Format details
                if experiment_type in ['hetNOE_sat', 'hetNOE_unsat']:
                    details = f"{len(result)} spectrum"
                else:
                    min_d = min(delay_values)
                    max_d = max(delay_values)
                    details = f"{len(result)} spectra, {min_d:.0f}-{max_d:.0f}ms"
                widgets['details'].setText(details)
            else:
                # Folder exists but no NMR files with delays
                # For hetNOE, this is OK (no delay expected)
                if experiment_type in ['hetNOE_sat', 'hetNOE_unsat']:
                    # Count any NMR files
                    nmr_files = [f for f in os.listdir(folder_path)
                                 if os.path.splitext(f)[1].lower() in
                                 {'.ft', '.ucsf', '.pipe', '.2rr', '.2ii'}]
                    if nmr_files:
                        widgets['status'].setText("Found")
                        widgets['status'].setStyleSheet(f"color: {SUCCESS_GREEN};")
                        widgets['details'].setText(f"{len(nmr_files)} spectrum")
                        delays_dict[experiment_type] = []
                    else:
                        self._clear_experiment_display(experiment_type, widgets_dict, delays_dict)
                else:
                    self._clear_experiment_display(experiment_type, widgets_dict, delays_dict)

        except ImportError as e:
            widgets['status'].setText("Error")
            widgets['status'].setStyleSheet(f"color: {WARNING_ORANGE};")
            widgets['details'].setText(str(e))

    def _clear_experiment_display(self, experiment_type: str, widgets_dict: dict, delays_dict: dict):
        """Clear the display for an experiment type."""
        widgets = widgets_dict[experiment_type]
        widgets['status'].setText("Not found")
        widgets['status'].setStyleSheet(f"color: {ERROR_RED};")
        widgets['details'].setText("")

        if experiment_type in delays_dict:
            del delays_dict[experiment_type]

    # -------------------------------------------------------------------------
    # Public API
    # -------------------------------------------------------------------------

    def get_peak_list_path(self) -> str:
        """Get the selected peak list path."""
        return self._peak_list_path

    def get_field1_experiments(self) -> dict:
        """Get Field 1 experiment folders."""
        return dict(self._field1_experiments)

    def get_field2_experiments(self) -> dict:
        """Get Field 2 experiment folders (empty if not enabled)."""
        if not self._field2_enabled:
            return {}
        return dict(self._field2_experiments)

    def get_field1_freq(self) -> float:
        """Get Field 1 frequency in MHz."""
        return self._field1_freq

    def get_field2_freq(self) -> float:
        """Get Field 2 frequency in MHz."""
        return self._field2_freq

    def is_dual_field_enabled(self) -> bool:
        """Check if dual-field mode is enabled."""
        return self._field2_enabled

    def get_field1_delays(self) -> dict:
        """Get Field 1 delay summaries."""
        return dict(self._field1_delays)

    def get_field2_delays(self) -> dict:
        """Get Field 2 delay summaries."""
        return dict(self._field2_delays)

    def is_ready_for_t1t2(self) -> bool:
        """Check if enough data is available for T1/T2 analysis."""
        has_t1 = 'T1' in self._field1_experiments
        has_t2 = 'T2' in self._field1_experiments
        has_peak_list = bool(self._peak_list_path)
        return (has_t1 or has_t2) and has_peak_list

    def is_ready_for_model_free(self) -> bool:
        """Check if all data is available for full Model-Free analysis (Field 1)."""
        has_t1 = 'T1' in self._field1_experiments
        has_t2 = 'T2' in self._field1_experiments
        has_noe_sat = 'hetNOE_sat' in self._field1_experiments
        has_noe_unsat = 'hetNOE_unsat' in self._field1_experiments
        has_peak_list = bool(self._peak_list_path)
        return has_t1 and has_t2 and has_noe_sat and has_noe_unsat and has_peak_list

    def is_ready_for_dual_field_model_free(self) -> bool:
        """Check if all data is available for dual-field Model-Free analysis."""
        if not self._field2_enabled:
            return False
        # Check Field 1
        if not self.is_ready_for_model_free():
            return False
        # Check Field 2
        has_t1 = 'T1' in self._field2_experiments
        has_t2 = 'T2' in self._field2_experiments
        has_noe_sat = 'hetNOE_sat' in self._field2_experiments
        has_noe_unsat = 'hetNOE_unsat' in self._field2_experiments
        return has_t1 and has_t2 and has_noe_sat and has_noe_unsat

    def get_validation_message(self) -> str:
        """Get a message describing what's missing."""
        missing = []
        if not self._peak_list_path:
            missing.append("Peak list")

        if self._single_experiment_mode:
            # Only need the spectra folder in single mode
            if 'T1' not in self._field1_experiments and 'T2' not in self._field1_experiments:
                missing.append("Spectra folder")
        else:
            # Field 1 requirements
            if 'T1' not in self._field1_experiments:
                missing.append("Field 1: T1 folder")
            if 'T2' not in self._field1_experiments:
                missing.append("Field 1: T2 folder")
            if 'hetNOE_sat' not in self._field1_experiments:
                missing.append("Field 1: hetNOE saturated folder")
            if 'hetNOE_unsat' not in self._field1_experiments:
                missing.append("Field 1: hetNOE unsaturated folder")

            # Field 2 requirements (only if enabled)
            if self._field2_enabled:
                if 'T1' not in self._field2_experiments:
                    missing.append("Field 2: T1 folder")
                if 'T2' not in self._field2_experiments:
                    missing.append("Field 2: T2 folder")
                if 'hetNOE_sat' not in self._field2_experiments:
                    missing.append("Field 2: hetNOE saturated folder")
                if 'hetNOE_unsat' not in self._field2_experiments:
                    missing.append("Field 2: hetNOE unsaturated folder")

        if not missing:
            return "Ready for analysis"
        return f"Missing: {', '.join(missing)}"

    def get_spectra_folder(self) -> str:
        """
        Get the spectra folder path (for single_experiment_mode).

        In single_experiment_mode, returns the master folder directly since
        it contains the spectra. Returns empty string if not set.
        """
        return self._field1_folder

    # Legacy compatibility methods
    def get_experiment_folders(self) -> dict:
        """Get Field 1 experiment folders (legacy compatibility)."""
        return self.get_field1_experiments()

    def get_delay_summary(self) -> dict:
        """Get Field 1 delay summaries (legacy compatibility)."""
        return self.get_field1_delays()
