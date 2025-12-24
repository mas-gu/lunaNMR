"""
DynamiXs v2.0 - Design Constants
Based on LunaNMR v1.0 Qt Style Guide

This module contains all color, spacing, typography, and sizing constants
used throughout the application for consistent styling.
"""

# =============================================================================
# Background Colors
# =============================================================================
BG_COLOR = "#FAFAFA"            # Main window background (soft white)
PANEL_BG_COLOR = "#F5F5F7"      # Secondary panels, frames, sidebars (Apple grey)
FRAME_BG_COLOR = "#FFFFFF"      # Elevated cards, input fields (pure white)

# =============================================================================
# Text Colors
# =============================================================================
PRIMARY_TEXT = "#1C1C1E"        # Main content, labels (near black)
SECONDARY_TEXT = "#8E8E93"      # Help text, descriptions (grey)
DISABLED_TEXT = "#C7C7CC"       # Disabled states (light grey)

# =============================================================================
# Primary Button Colors
# =============================================================================
PRIMARY_BUTTON_BG = "#5B9EE5"           # Main action buttons (pleasant blue)
PRIMARY_BUTTON_HOVER = "#4A8DD4"        # Primary button hover state
PRIMARY_BUTTON_PRESSED = "#3A7DC4"      # Primary button pressed state
PRIMARY_BUTTON_TEXT = "#FFFFFF"         # Primary button text (white)

# =============================================================================
# Secondary Button Colors
# =============================================================================
SECONDARY_BUTTON_BG = "#E5E5EA"         # Secondary actions (light grey)
SECONDARY_BUTTON_HOVER = "#D1D1D6"      # Secondary hover state
SECONDARY_BUTTON_PRESSED = "#C8C8CD"    # Secondary pressed state
SECONDARY_BUTTON_TEXT = "#1C1C1E"       # Secondary button text (dark)
SECONDARY_BUTTON_BORDER = "#C8C8CD"     # Required 1px border

# =============================================================================
# Destructive Button Colors
# =============================================================================
DESTRUCTIVE_BUTTON_BG = "#E8554E"       # Errors, destructive actions (soft red)
DESTRUCTIVE_BUTTON_HOVER = "#D44943"    # Destructive hover state
DESTRUCTIVE_BUTTON_PRESSED = "#C43D37"  # Destructive pressed state
DESTRUCTIVE_BUTTON_TEXT = "#FFFFFF"     # Destructive button text (white)

# =============================================================================
# Success Button Colors
# =============================================================================
SUCCESS_BUTTON_BG = "#34C759"           # Success states (green)
SUCCESS_BUTTON_HOVER = "#2DB84E"        # Success hover state
SUCCESS_BUTTON_PRESSED = "#28A745"      # Success pressed state
SUCCESS_BUTTON_TEXT = "#FFFFFF"         # Success button text (white)

# =============================================================================
# Status Colors
# =============================================================================
SUCCESS_GREEN = "#34C759"       # Success states, progress
WARNING_ORANGE = "#F0A04B"      # Warnings, attention
ERROR_RED = "#E8554E"           # Errors, destructive actions
INFO_BLUE = "#5B9EE5"           # Informational messages

# =============================================================================
# Border & Separator Colors
# =============================================================================
BORDER_COLOR = "#D1D1D6"        # Input borders, frame borders
SEPARATOR_COLOR = "#E5E5EA"     # Dividers, grid lines

# =============================================================================
# Corner Radii (in pixels)
# =============================================================================
BUTTON_CORNER_RADIUS = 10       # All button types
FRAME_CORNER_RADIUS = 12        # Panels, group boxes
DIALOG_CORNER_RADIUS = 14       # Modal dialogs, popups
CARD_CORNER_RADIUS = 8          # Smaller UI elements, inputs

# =============================================================================
# Spacing System (8-point grid)
# =============================================================================
SPACING_XS = 4                  # Tight spacing (label-to-field)
SPACING_SM = 8                  # Default spacing (related elements)
SPACING_MD = 16                 # Section spacing (between groups)
SPACING_LG = 24                 # Major section breaks
SPACING_XL = 32                 # Window padding, large separations

# =============================================================================
# Font Sizes (in points)
# =============================================================================
FONT_SIZE_LARGE_HEADER = 14     # Section titles, main action buttons
FONT_SIZE_MEDIUM_HEADER = 12    # Subsection titles
FONT_SIZE_SECTION_LABEL = 11    # Important labels
FONT_SIZE_BODY = 10             # Standard UI text, buttons
FONT_SIZE_SMALL = 9             # Help text, secondary info
FONT_SIZE_TINY = 8              # Metadata, annotations
FONT_SIZE_MONOSPACE = 10        # Data display, code output

# =============================================================================
# Font Families
# =============================================================================
FONT_FAMILY = '"Segoe UI", "Helvetica Neue", "Arial", sans-serif'
MONOSPACE_FONT = "Courier New"
MONOSPACE_SIZE = 9

# =============================================================================
# Font Weights
# =============================================================================
FONT_WEIGHT_REGULAR = 50        # QFont.Weight.Normal
FONT_WEIGHT_BOLD = 75           # QFont.Weight.Bold

# =============================================================================
# Button Sizing
# =============================================================================
BUTTON_HEIGHT_STANDARD = 28     # Standard button height
BUTTON_HEIGHT_DIALOG = 40       # Dialog OK/Cancel button height
BUTTON_WIDTH_ICON = 36          # Icon-only button size
BUTTON_WIDTH_COMPACT = 32       # Navigation arrows
BUTTON_WIDTH_SMALL = 60         # Small actions
BUTTON_WIDTH_STANDARD = 80      # Standard actions
BUTTON_WIDTH_MEDIUM = 100       # Medium actions
BUTTON_WIDTH_LARGE = 180        # Primary actions
BUTTON_WIDTH_DIALOG = 120       # Dialog buttons

# =============================================================================
# Input Sizing
# =============================================================================
INPUT_HEIGHT_STANDARD = 20      # Standard input minimum height
INPUT_PADDING_H = 8             # Horizontal padding
INPUT_PADDING_V = 6             # Vertical padding

# =============================================================================
# Window Sizes
# =============================================================================
WINDOW_MIN_WIDTH = 900          # Minimum main window width
WINDOW_MIN_HEIGHT = 700         # Minimum main window height
DIALOG_MIN_WIDTH = 400          # Minimum dialog width

# =============================================================================
# Matplotlib Integration Colors
# =============================================================================
MPL_FIGURE_BG = PANEL_BG_COLOR  # Figure background (matches Qt panels)
MPL_AXES_BG = FRAME_BG_COLOR    # Axes background (white)
MPL_TEXT_COLOR = PRIMARY_TEXT   # Text color
MPL_GRID_COLOR = SEPARATOR_COLOR  # Grid lines
MPL_EDGE_COLOR = BORDER_COLOR   # Axes edge color

# Quality-based color coding for plots
def get_quality_color(r_squared: float) -> str:
    """
    Get color for quality indicator based on R² value.

    Args:
        r_squared: R² fit quality value (0-1)

    Returns:
        Hex color string
    """
    if r_squared >= 0.95:
        return SUCCESS_GREEN    # Excellent
    elif r_squared >= 0.85:
        return INFO_BLUE        # Good
    elif r_squared >= 0.75:
        return WARNING_ORANGE   # Fair
    else:
        return ERROR_RED        # Poor

# =============================================================================
# Platform-specific settings
# =============================================================================
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

def get_emoji_font() -> str:
    """Get the appropriate emoji font for the current platform."""
    system = platform.system()
    if system == "Darwin":
        return "Apple Color Emoji"
    elif system == "Windows":
        return "Segoe UI Emoji"
    else:
        return "Noto Color Emoji"

# =============================================================================
# Emoji Fallbacks (for systems without emoji support)
# =============================================================================
EMOJI_FALLBACKS = {
    "✅": "[OK]",
    "❌": "[FAIL]",
    "📍": "[REF]",
    "📄": "[DOC]",
    "📁": "[DIR]",
    "🔍": "[SEARCH]",
    "📊": "[CHART]",
    "📈": "[UP]",
    "📉": "[DOWN]",
    "💾": "[SAVE]",
    "🔧": "[TOOL]",
    "⚙️": "[SETTINGS]",
    "⚠️": "[WARN]",
    "🎛️": "[CTRL]",
    "⏱️": "[TIME]",
    "🔬": "[LAB]",
    "✔️": "[CHECK]",
    "⏹️": "[STOP]",
}

# =============================================================================
# Hex to RGB/Matplotlib conversion utilities
# =============================================================================
def hex_to_rgb(hex_color: str) -> tuple:
    """
    Convert hex color to RGB tuple (0-255 range).

    Args:
        hex_color: Hex color string (e.g., "#FFFFFF")

    Returns:
        Tuple of (R, G, B) values (0-255)
    """
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

def hex_to_mpl(hex_color: str) -> tuple:
    """
    Convert hex color to matplotlib-compatible tuple (0-1 range).

    Args:
        hex_color: Hex color string (e.g., "#FFFFFF")

    Returns:
        Tuple of (R, G, B) values (0.0-1.0)
    """
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) / 255 for i in (0, 2, 4))
