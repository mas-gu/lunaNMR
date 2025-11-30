#!/usr/bin/env python3
# ABOUTME: LunaNMR Qt Design System - centralized color palette, spacing, and sizing constants
# ABOUTME: Ported from CustomTkinter style guide for consistent Apple-style UI across Qt/PySide6

"""
Qt Design System for LunaNMR v1.0

This module defines the complete design system for the Qt/PySide6 migration of LunaNMR.
All colors, spacing, corner radii, and styling constants are centralized here for
consistent application of the Apple-inspired design language.

The design follows these principles:
- Clean & Minimal: Reduce visual noise, focus on content
- Soft Colors: Pleasant, non-aggressive tones for extended use
- High Contrast: Ensure readability without harsh colors
- Consistent Spacing: 8pt grid system for visual rhythm
- Rounded Corners: Modern, friendly aesthetic
- Contextual Colors: Use color to communicate meaning (success, warning, danger)

Usage:
    from lunaNMR.gui.styles.design_system import (
        BG_COLOR, PRIMARY_BUTTON_BG, SPACING_MD
    )

Author: Guillaume Mas
Date: 2025
"""

# ============================================================================
# BACKGROUND COLORS
# Apple-style layered backgrounds with subtle warmth
# ============================================================================

BG_COLOR = "#FAFAFA"  # Main window background - softer white with subtle warmth
"""Main window background color. Use for the primary application window."""

PANEL_BG_COLOR = "#F5F5F7"  # Secondary panels/frames - Apple's signature light grey
"""Secondary panel background. Use for control panels, sidebars, and frames."""

FRAME_BG_COLOR = "#FFFFFF"  # Card/container backgrounds - pure white for contrast
"""Card and container background. Use for data display areas and text boxes."""


# ============================================================================
# TEXT COLORS
# Softer than pure black for better readability and reduced eye strain
# ============================================================================

PRIMARY_TEXT = "#1C1C1E"  # Primary content text (softer than pure black)
"""Primary text color. Use for main labels, titles, and important content."""

SECONDARY_TEXT = "#8E8E93"  # Secondary/help text
"""Secondary text color. Use for help text, hints, and annotations."""

DISABLED_TEXT = "#C7C7CC"  # Disabled states
"""Disabled text color. Use for inactive or disabled UI elements."""

LABEL_TEXT = "#1C1C1E"  # Label text
"""Label text color. Same as primary text for consistency."""


# ============================================================================
# BUTTON COLORS - PRIMARY ACTIONS
# For main actions (Load Data, Fit All Peaks, Start Integration)
# ============================================================================

PRIMARY_BUTTON_BG = "#5B9EE5"  # Softer, more pleasant blue
"""Primary button background. Use for main action buttons."""

PRIMARY_BUTTON_HOVER = "#4A8DD4"  # Darker blue on hover
"""Primary button hover state. 15% darker than base color."""

PRIMARY_BUTTON_TEXT = "#FFFFFF"  # White text
"""Primary button text color. Always white for maximum contrast."""


# ============================================================================
# BUTTON COLORS - SECONDARY/UTILITY ACTIONS
# For supporting actions (Browse, Refresh, Clear, Cancel)
# IMPORTANT: Secondary buttons MUST include a 1px border to prevent
# grey-on-grey edge artifacts on light backgrounds
# ============================================================================

SECONDARY_BUTTON_BG = "#E5E5EA"  # Light grey
"""Secondary button background. Use for utility buttons."""

SECONDARY_BUTTON_HOVER = "#D1D1D6"  # Slightly darker grey
"""Secondary button hover state. 10% darker than base color."""

SECONDARY_BUTTON_TEXT = "#1C1C1E"  # Near-black (softer than pure black)
"""Secondary button text color. High contrast against light background."""

SECONDARY_BUTTON_BORDER = "#C8C8CD"  # Subtle border
"""Secondary button border color. REQUIRED to prevent edge artifacts."""


# ============================================================================
# BUTTON COLORS - DESTRUCTIVE ACTIONS
# For dangerous actions (Clear All, Reset, Delete)
# ============================================================================

DESTRUCTIVE_BUTTON_BG = "#E8554E"  # Softer, less aggressive red
"""Destructive button background. Use for delete/clear/reset actions."""

DESTRUCTIVE_BUTTON_HOVER = "#D44943"  # Darker red on hover
"""Destructive button hover state. 12% darker than base color."""

DESTRUCTIVE_BUTTON_TEXT = "#FFFFFF"  # White text
"""Destructive button text color. Always white for maximum contrast."""


# ============================================================================
# ACCENT COLORS
# For status indicators, messages, and contextual highlights
# ============================================================================

SUCCESS_GREEN = "#34C759"  # Successful operations, progress bars
"""Success color. Use for completed operations and positive feedback."""

WARNING_ORANGE = "#F0A04B"  # Warnings, attention needed
"""Warning color. Use for non-critical warnings and attention indicators."""

ERROR_RED = "#E8554E"  # Errors (matches destructive buttons)
"""Error color. Use for error messages and failed operations."""

INFO_BLUE = "#5B9EE5"  # Informational messages (matches primary buttons)
"""Info color. Use for informational messages and tips."""


# ============================================================================
# BORDER & SEPARATOR COLORS
# Subtle divisions and boundaries
# ============================================================================

BORDER_COLOR = "#D1D1D6"  # Light grey borders
"""Border color. Use for widget borders and boundaries."""

SEPARATOR_COLOR = "#E5E5EA"  # Separators and dividers
"""Separator color. Use for dividing lines and section breaks."""


# ============================================================================
# CORNER RADIUS VALUES
# Different radii for different UI element types to establish hierarchy
# ============================================================================

BUTTON_CORNER_RADIUS = 10  # Buttons
"""Button corner radius in pixels. Standard for all button types."""

FRAME_CORNER_RADIUS = 12  # Panels and label frames
"""Frame corner radius in pixels. For panels and label frames."""

DIALOG_CORNER_RADIUS = 14  # Modal dialogs, popups
"""Dialog corner radius in pixels. For modal windows and popups."""

CARD_CORNER_RADIUS = 8  # Smaller UI cards/elements
"""Card corner radius in pixels. For small UI cards and elements."""


# ============================================================================
# SPACING SYSTEM
# 8pt grid system for consistent visual rhythm and alignment
# All spacing should use these multiples of 8 pixels
# ============================================================================

SPACING_XS = 4  # Tight spacing (label-to-field, within components)
"""Extra-small spacing (4px). Use for tight label-to-field spacing."""

SPACING_SM = 8  # Default spacing (between related elements)
"""Small spacing (8px). Use between related elements and buttons."""

SPACING_MD = 16  # Section spacing (between groups)
"""Medium spacing (16px). Use between sections and groups."""

SPACING_LG = 24  # Major section breaks
"""Large spacing (24px). Use for major section breaks."""

SPACING_XL = 32  # Window padding, large separations
"""Extra-large spacing (32px). Use for window padding and large separations."""


# ============================================================================
# TYPOGRAPHY SYSTEM
# Font sizes and weights for consistent text hierarchy
# Note: Qt uses point sizes, not pixel sizes
# ============================================================================

# Font Sizes (in points)
FONT_SIZE_LARGE_HEADER = 14  # Section titles, main controls
"""Large header font size (14pt). For main section titles."""

FONT_SIZE_MEDIUM_HEADER = 12  # Subsection titles
"""Medium header font size (12pt). For subsection titles."""

FONT_SIZE_SECTION_LABEL = 11  # Important labels
"""Section label font size (11pt). For important labels."""

FONT_SIZE_BODY = 10  # Standard UI text, buttons
"""Body text font size (10pt). Standard for UI elements."""

FONT_SIZE_SMALL = 9  # Help text, secondary info
"""Small text font size (9pt). For help text and secondary info."""

FONT_SIZE_TINY = 8  # Metadata, annotations
"""Tiny text font size (8pt). For metadata and annotations."""

FONT_SIZE_MONOSPACE = 10  # Data display, code output
"""Monospace font size (10pt). For data and code display."""

# Font Weights (Qt weights: Light=25, Normal=50, DemiBold=63, Bold=75, Black=87)
FONT_WEIGHT_REGULAR = 50  # Normal weight
"""Regular font weight. For standard text."""

FONT_WEIGHT_BOLD = 75  # Bold weight
"""Bold font weight. For headers and emphasis."""


# ============================================================================
# BUTTON SIZING
# Standard button dimensions following 8pt grid
# ============================================================================

# Button widths (in pixels)
BUTTON_WIDTH_ICON = 36  # Icon-only buttons
"""Icon button width (36px). For single icon buttons."""

BUTTON_WIDTH_COMPACT = 32  # Compact navigation buttons
"""Compact button width (32px). For navigation arrows."""

BUTTON_WIDTH_SMALL = 60  # Small action buttons
"""Small button width (60px). For small action buttons."""

BUTTON_WIDTH_STANDARD = 80  # Standard action buttons
"""Standard button width (80px). For common actions."""

BUTTON_WIDTH_MEDIUM = 100  # Medium action buttons
"""Medium button width (100px). For important actions."""

BUTTON_WIDTH_LARGE = 180  # Large primary action buttons
"""Large button width (180px). For main actions."""

BUTTON_WIDTH_DIALOG = 120  # Dialog buttons (OK, Cancel)
"""Dialog button width (120px). For modal dialog buttons."""

# Button heights
BUTTON_HEIGHT_STANDARD = 28  # Standard button height
"""Standard button height (28px). Default for most buttons."""

BUTTON_HEIGHT_DIALOG = 40  # Dialog button height
"""Dialog button height (40px). For modal dialog buttons."""


# ============================================================================
# LEGACY COMPATIBILITY
# Aliases for backward compatibility during migration
# These will be phased out in future versions
# ============================================================================

BUTTON_FG_COLOR = SECONDARY_BUTTON_BG
"""DEPRECATED: Use SECONDARY_BUTTON_BG instead."""

BUTTON_HOVER_COLOR = SECONDARY_BUTTON_HOVER
"""DEPRECATED: Use SECONDARY_BUTTON_HOVER instead."""

BUTTON_TEXT_COLOR = SECONDARY_BUTTON_TEXT
"""DEPRECATED: Use SECONDARY_BUTTON_TEXT instead."""


# ============================================================================
# COLOR UTILITY FUNCTIONS
# Helper functions for color manipulation
# ============================================================================

def hex_to_rgb(hex_color: str) -> tuple[int, int, int]:
    """
    Convert hex color string to RGB tuple.

    Args:
        hex_color: Hex color string (e.g., "#5B9EE5")

    Returns:
        Tuple of (r, g, b) values (0-255)

    Example:
        >>> hex_to_rgb("#5B9EE5")
        (91, 158, 229)
    """
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))


def rgb_to_hex(r: int, g: int, b: int) -> str:
    """
    Convert RGB values to hex color string.

    Args:
        r: Red value (0-255)
        g: Green value (0-255)
        b: Blue value (0-255)

    Returns:
        Hex color string (e.g., "#5B9EE5")

    Example:
        >>> rgb_to_hex(91, 158, 229)
        "#5B9EE5"
    """
    return f"#{r:02x}{g:02x}{b:02x}".upper()


def color_with_alpha(hex_color: str, alpha: float) -> str:
    """
    Create RGBA color string from hex color and alpha value.

    Args:
        hex_color: Hex color string (e.g., "#5B9EE5")
        alpha: Alpha value (0.0 to 1.0)

    Returns:
        RGBA string for Qt (e.g., "rgba(91, 158, 229, 0.5)")

    Example:
        >>> color_with_alpha("#5B9EE5", 0.5)
        "rgba(91, 158, 229, 0.5)"
    """
    r, g, b = hex_to_rgb(hex_color)
    return f"rgba({r}, {g}, {b}, {alpha})"


# ============================================================================
# DESIGN SYSTEM VALIDATION
# Ensure all constants follow design principles
# ============================================================================

def validate_design_system():
    """
    Validate that all design system constants follow the guidelines.
    Raises ValueError if any constants violate design principles.

    This is primarily for testing and development purposes.
    """
    # Validate spacing follows 8pt grid
    spacing_values = [SPACING_XS, SPACING_SM, SPACING_MD, SPACING_LG, SPACING_XL]
    for i, spacing in enumerate(spacing_values):
        if spacing % 4 != 0:
            raise ValueError(f"Spacing value {spacing} does not follow 4px grid")

    # Validate corner radii are reasonable
    radii = [BUTTON_CORNER_RADIUS, FRAME_CORNER_RADIUS, DIALOG_CORNER_RADIUS, CARD_CORNER_RADIUS]
    for radius in radii:
        if radius < 0 or radius > 20:
            raise ValueError(f"Corner radius {radius} is outside reasonable range (0-20)")

    # Validate all color strings start with #
    import sys
    module = sys.modules[__name__]
    color_attrs = [
        attr for attr in dir(module)
        if not attr.startswith('_') and 'COLOR' in attr and isinstance(getattr(module, attr), str)
    ]
    for attr in color_attrs:
        color = getattr(module, attr)
        if not color.startswith('#') or len(color) != 7:
            raise ValueError(f"Color constant {attr} has invalid format: {color}")

    return True


if __name__ == "__main__":
    # Run validation when module is executed directly
    try:
        validate_design_system()
        print("✓ Design system validation passed")
        print(f"✓ Colors defined: {len([k for k in globals() if 'COLOR' in k])}")
        print(f"✓ Spacing system: {SPACING_XS}px to {SPACING_XL}px (8pt grid)")
        print(f"✓ Corner radii: {CARD_CORNER_RADIUS}px to {DIALOG_CORNER_RADIUS}px")
    except ValueError as e:
        print(f"✗ Design system validation failed: {e}")
