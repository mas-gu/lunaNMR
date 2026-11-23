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
    QLabel, QPushButton, QLineEdit, QPlainTextEdit,
    QGroupBox, QCheckBox, QRadioButton,
    QComboBox, QSpinBox, QDoubleSpinBox,
    QFileDialog, QMessageBox, QApplication, QInputDialog
)
from PySide6.QtCore import Signal
from PySide6.QtGui import QFont

from constants import (
    SPACING_SM, SPACING_MD,
    FONT_SIZE_LARGE_HEADER, FONT_SIZE_MEDIUM_HEADER, FONT_SIZE_SECTION_LABEL,
    FONT_SIZE_BODY, FONT_SIZE_SMALL, FONT_SIZE_MONOSPACE,
    get_system_font
)


# STYLESHEET LOADING

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


# FONT UTILITIES

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


# CUSTOM COMPONENTS

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


# BUTTON FACTORIES

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


# LABEL FACTORIES

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


# INPUT FACTORIES

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


# LAYOUT UTILITIES

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


# DIALOG UTILITIES

def show_info(parent: QWidget, title: str, message: str):
    """Show an information message box."""
    QMessageBox.information(parent, title, message)


def show_warning(parent: QWidget, title: str, message: str):
    """Show a warning message box."""
    QMessageBox.warning(parent, title, message)


def show_error(parent: QWidget, title: str, message: str):
    """Show an error message box."""
    QMessageBox.critical(parent, title, message)


# NAMING A SAVED ANALYSIS

def unique_name(base: str, existing) -> str:
    """`base`, or the first of base_2, base_3, ... not already in `existing`."""
    base = base or "analysis"
    if base not in existing:
        return base
    i = 2
    while f"{base}_{i}" in existing:
        i += 1
    return f"{base}_{i}"


def prompt_name_collision(parent: QWidget, base: str, existing, kind: str) -> Optional[str]:
    """A saved `kind` (e.g. "Kd analysis") named `base` already exists: ask whether to
    replace it or save under a different name. Returns the name to store under (the base
    to replace, or a new name), or None to cancel."""
    box = QMessageBox(parent)
    box.setWindowTitle("Analysis already exists")
    box.setText(f"A saved {kind} named '{base}' already exists in this project.")
    box.setInformativeText("Replace it, or save under a different name?")
    replace_btn = box.addButton("Replace", QMessageBox.AcceptRole)
    rename_btn = box.addButton("Save as…", QMessageBox.ActionRole)
    box.addButton(QMessageBox.Cancel)
    box.exec()
    clicked = box.clickedButton()
    if clicked is replace_btn:
        return base
    if clicked is rename_btn:
        return prompt_new_name(parent, base, existing, kind)
    return None


def prompt_new_name(parent: QWidget, base: str, existing, kind: str) -> Optional[str]:
    """Ask for a name for a `kind`, pre-filled with a collision-free suggestion.
    Re-prompts on a blank name and confirms before overwriting an existing one.
    Returns the chosen name, or None if the user cancels."""
    suggestion = unique_name(base, existing)
    while True:
        name, ok = QInputDialog.getText(
            parent, "Save analysis as", f"Name for this {kind}:", text=suggestion)
        if not ok:
            return None
        name = name.strip()
        if not name:
            continue
        if name in existing and QMessageBox.question(
                parent, "Name in use", f"'{name}' also exists. Replace it?") != QMessageBox.Yes:
            continue
        return name


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
