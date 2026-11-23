"""
DynamiXs v2.0 - Design Constants
Based on LunaNMR v1.0 Qt Style Guide

This module contains all color, spacing, typography, and sizing constants
used throughout the application for consistent styling.
"""

# Background Colors
BG_COLOR = "#FAFAFA"            # Main window background (soft white)
PANEL_BG_COLOR = "#F5F5F7"      # Secondary panels, frames, sidebars (Apple grey)
FRAME_BG_COLOR = "#FFFFFF"      # Elevated cards, input fields (pure white)

# Text Colors
PRIMARY_TEXT = "#1C1C1E"        # Main content, labels (near black)
SECONDARY_TEXT = "#8E8E93"      # Help text, descriptions (grey)

# Primary Button Colors
PRIMARY_BUTTON_BG = "#5B9EE5"           # Main action buttons (pleasant blue)
PRIMARY_BUTTON_HOVER = "#4A8DD4"        # Primary button hover state

# Status Colors
SUCCESS_GREEN = "#34C759"       # Success states, progress
WARNING_ORANGE = "#F0A04B"      # Warnings, attention

# Spacing System (8-point grid)
SPACING_XS = 4                  # Tight spacing (label-to-field)
SPACING_SM = 8                  # Default spacing (related elements)
SPACING_MD = 16                 # Section spacing (between groups)
SPACING_XL = 32                 # Window padding, large separations

# Font Sizes (in points)
FONT_SIZE_LARGE_HEADER = 14     # Section titles, main action buttons
FONT_SIZE_MEDIUM_HEADER = 12    # Subsection titles
FONT_SIZE_SECTION_LABEL = 11    # Important labels
FONT_SIZE_BODY = 10             # Standard UI text, buttons
FONT_SIZE_SMALL = 9             # Help text, secondary info
FONT_SIZE_MONOSPACE = 10        # Data display, code output

# Window Sizes
WINDOW_MIN_WIDTH = 900          # Minimum main window width
WINDOW_MIN_HEIGHT = 700         # Minimum main window height

# Platform-specific settings
import platform

def get_system_font() -> str:
    """Get the appropriate system font for the current platform."""
    system = platform.system()
    if system == "Darwin":  # macOS
        return "Helvetica Neue"
    elif system == "Windows":
        return "Segoe UI"
    else:  # Linux
        return "Ubuntu"
