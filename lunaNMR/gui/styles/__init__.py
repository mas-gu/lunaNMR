#!/usr/bin/env python3
# ABOUTME: Styles module initialization - provides easy access to design system and stylesheet loading
# ABOUTME: Centralized import point for LunaNMR Qt design system constants and QSS loading utilities

"""
LunaNMR Styles Module

This module provides centralized access to the design system and stylesheet loading
for the Qt/PySide6 version of LunaNMR.

Usage:
    # Import design constants
    from lunaNMR.gui.styles import (
        BG_COLOR, PRIMARY_BUTTON_BG, SPACING_MD, load_stylesheet
    )

    # Load and apply stylesheet
    app = QApplication(sys.argv)
    app.setStyleSheet(load_stylesheet())

Author: Guillaume Mas
Date: 2025
"""

import os
from pathlib import Path
from typing import Optional

# Import all design system constants for easy access
from .design_system import (
    # Background Colors
    BG_COLOR,
    PANEL_BG_COLOR,
    FRAME_BG_COLOR,

    # Text Colors
    PRIMARY_TEXT,
    SECONDARY_TEXT,
    DISABLED_TEXT,
    LABEL_TEXT,

    # Primary Button Colors
    PRIMARY_BUTTON_BG,
    PRIMARY_BUTTON_HOVER,
    PRIMARY_BUTTON_TEXT,

    # Secondary Button Colors
    SECONDARY_BUTTON_BG,
    SECONDARY_BUTTON_HOVER,
    SECONDARY_BUTTON_TEXT,
    SECONDARY_BUTTON_BORDER,

    # Destructive Button Colors
    DESTRUCTIVE_BUTTON_BG,
    DESTRUCTIVE_BUTTON_HOVER,
    DESTRUCTIVE_BUTTON_TEXT,

    # Accent Colors
    SUCCESS_GREEN,
    WARNING_ORANGE,
    ERROR_RED,
    INFO_BLUE,

    # Border and Separator Colors
    BORDER_COLOR,
    SEPARATOR_COLOR,

    # Corner Radius Values
    BUTTON_CORNER_RADIUS,
    FRAME_CORNER_RADIUS,
    DIALOG_CORNER_RADIUS,
    CARD_CORNER_RADIUS,

    # Spacing System
    SPACING_XS,
    SPACING_SM,
    SPACING_MD,
    SPACING_LG,
    SPACING_XL,

    # Typography
    FONT_SIZE_LARGE_HEADER,
    FONT_SIZE_MEDIUM_HEADER,
    FONT_SIZE_SECTION_LABEL,
    FONT_SIZE_BODY,
    FONT_SIZE_SMALL,
    FONT_SIZE_TINY,
    FONT_SIZE_MONOSPACE,
    FONT_WEIGHT_REGULAR,
    FONT_WEIGHT_BOLD,

    # Button Sizing
    BUTTON_WIDTH_ICON,
    BUTTON_WIDTH_COMPACT,
    BUTTON_WIDTH_SMALL,
    BUTTON_WIDTH_STANDARD,
    BUTTON_WIDTH_MEDIUM,
    BUTTON_WIDTH_LARGE,
    BUTTON_WIDTH_DIALOG,
    BUTTON_HEIGHT_STANDARD,
    BUTTON_HEIGHT_DIALOG,

    # Utility Functions
    hex_to_rgb,
    rgb_to_hex,
    color_with_alpha,
    validate_design_system,
)


def get_styles_dir() -> Path:
    """
    Get the absolute path to the styles directory.

    Returns:
        Path object pointing to the styles directory
    """
    return Path(__file__).parent


def load_stylesheet(stylesheet_name: str = "main.qss") -> str:
    """
    Load a Qt stylesheet (QSS) file from the styles directory.

    Args:
        stylesheet_name: Name of the stylesheet file (default: "main.qss")

    Returns:
        String containing the QSS stylesheet content

    Raises:
        FileNotFoundError: If the stylesheet file doesn't exist
        IOError: If there's an error reading the file

    Example:
        >>> from lunaNMR.gui.styles import load_stylesheet
        >>> app = QApplication(sys.argv)
        >>> app.setStyleSheet(load_stylesheet())
    """
    stylesheet_path = get_styles_dir() / stylesheet_name

    if not stylesheet_path.exists():
        raise FileNotFoundError(
            f"Stylesheet '{stylesheet_name}' not found at {stylesheet_path}"
        )

    try:
        with open(stylesheet_path, 'r', encoding='utf-8') as f:
            stylesheet = f.read()
        return stylesheet
    except Exception as e:
        raise IOError(f"Error reading stylesheet '{stylesheet_name}': {e}")


def get_stylesheet_path(stylesheet_name: str = "main.qss") -> str:
    """
    Get the absolute path to a stylesheet file as a string.

    This is useful for Qt applications that need a file path rather than
    the stylesheet content directly.

    Args:
        stylesheet_name: Name of the stylesheet file (default: "main.qss")

    Returns:
        Absolute path to the stylesheet file as a string

    Example:
        >>> from lunaNMR.gui.styles import get_stylesheet_path
        >>> with open(get_stylesheet_path()) as f:
        ...     stylesheet = f.read()
    """
    return str(get_styles_dir() / stylesheet_name)


def apply_stylesheet_to_widget(widget, stylesheet_name: str = "main.qss"):
    """
    Load and apply a stylesheet directly to a Qt widget.

    Args:
        widget: Qt widget (QWidget, QMainWindow, QDialog, etc.)
        stylesheet_name: Name of the stylesheet file (default: "main.qss")

    Example:
        >>> from lunaNMR.gui.styles import apply_stylesheet_to_widget
        >>> window = QMainWindow()
        >>> apply_stylesheet_to_widget(window)
    """
    stylesheet = load_stylesheet(stylesheet_name)
    widget.setStyleSheet(stylesheet)


def validate_styles():
    """
    Validate that all style resources are available and valid.

    This function checks:
    - Design system constants are valid
    - Stylesheet files exist
    - Color values are properly formatted

    Returns:
        bool: True if all validations pass

    Raises:
        ValueError: If any validation fails
    """
    # Validate design system constants
    validate_design_system()

    # Check that main stylesheet exists
    main_stylesheet = get_styles_dir() / "main.qss"
    if not main_stylesheet.exists():
        raise FileNotFoundError(f"Main stylesheet not found at {main_stylesheet}")

    # Try to load the stylesheet
    try:
        load_stylesheet()
    except Exception as e:
        raise ValueError(f"Failed to load main stylesheet: {e}")

    return True


# Module-level constants for convenient access
__all__ = [
    # Background Colors
    'BG_COLOR',
    'PANEL_BG_COLOR',
    'FRAME_BG_COLOR',

    # Text Colors
    'PRIMARY_TEXT',
    'SECONDARY_TEXT',
    'DISABLED_TEXT',
    'LABEL_TEXT',

    # Button Colors
    'PRIMARY_BUTTON_BG',
    'PRIMARY_BUTTON_HOVER',
    'PRIMARY_BUTTON_TEXT',
    'SECONDARY_BUTTON_BG',
    'SECONDARY_BUTTON_HOVER',
    'SECONDARY_BUTTON_TEXT',
    'SECONDARY_BUTTON_BORDER',
    'DESTRUCTIVE_BUTTON_BG',
    'DESTRUCTIVE_BUTTON_HOVER',
    'DESTRUCTIVE_BUTTON_TEXT',

    # Accent Colors
    'SUCCESS_GREEN',
    'WARNING_ORANGE',
    'ERROR_RED',
    'INFO_BLUE',

    # Borders
    'BORDER_COLOR',
    'SEPARATOR_COLOR',

    # Corner Radius
    'BUTTON_CORNER_RADIUS',
    'FRAME_CORNER_RADIUS',
    'DIALOG_CORNER_RADIUS',
    'CARD_CORNER_RADIUS',

    # Spacing
    'SPACING_XS',
    'SPACING_SM',
    'SPACING_MD',
    'SPACING_LG',
    'SPACING_XL',

    # Typography
    'FONT_SIZE_LARGE_HEADER',
    'FONT_SIZE_MEDIUM_HEADER',
    'FONT_SIZE_SECTION_LABEL',
    'FONT_SIZE_BODY',
    'FONT_SIZE_SMALL',
    'FONT_SIZE_TINY',
    'FONT_SIZE_MONOSPACE',
    'FONT_WEIGHT_REGULAR',
    'FONT_WEIGHT_BOLD',

    # Button Sizing
    'BUTTON_WIDTH_ICON',
    'BUTTON_WIDTH_COMPACT',
    'BUTTON_WIDTH_SMALL',
    'BUTTON_WIDTH_STANDARD',
    'BUTTON_WIDTH_MEDIUM',
    'BUTTON_WIDTH_LARGE',
    'BUTTON_WIDTH_DIALOG',
    'BUTTON_HEIGHT_STANDARD',
    'BUTTON_HEIGHT_DIALOG',

    # Functions
    'load_stylesheet',
    'get_stylesheet_path',
    'get_styles_dir',
    'apply_stylesheet_to_widget',
    'validate_styles',
    'hex_to_rgb',
    'rgb_to_hex',
    'color_with_alpha',
    'validate_design_system',
]


if __name__ == "__main__":
    # Run validation when module is executed directly
    try:
        validate_styles()
        print("✓ All styles validated successfully")
        print(f"✓ Styles directory: {get_styles_dir()}")
        print(f"✓ Main stylesheet: {get_stylesheet_path()}")
        print(f"✓ Design constants exported: {len([k for k in __all__ if not k.startswith('_')])}")
    except Exception as e:
        print(f"✗ Style validation failed: {e}")
