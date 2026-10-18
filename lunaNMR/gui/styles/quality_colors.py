# ABOUTME: The colour and symbol a fitted peak is shown with, from its fit's quality label.
# ABOUTME: Shared by the spectrum markers and both peak navigators so the three always agree.
"""
The label is the fit's own verdict (Excellent / Good / Fair / Poor / Failed). It is
read, never recomputed from R-squared: the fitters assign it with their own
thresholds and special cases (a near-noise peak is Poor at any R-squared).
"""

from typing import Mapping, Optional

from PySide6.QtGui import QColor

from lunaNMR.gui.styles.design_system import ERROR_RED, SUCCESS_GREEN, WARNING_ORANGE

QUALITY_COLORS = {
    'excellent': 'lime',
    'good': SUCCESS_GREEN,
    'fair': WARNING_ORANGE,
    'poor': ERROR_RED,
}

QUALITY_SYMBOLS = {
    'excellent': "✅",
    'good': "✅",
    'fair': "🟠",
    'poor': "🔴",
}

# A label outside the four above: Failed, Unknown, or a peak detected but not fitted.
UNRATED_COLOR = ERROR_RED
UNRATED_SYMBOL = "❌"


def peak_quality_label(peak: Mapping) -> Optional[str]:
    """The quality label a peak record carries, under any of the keys the pipeline uses."""
    for key in ('quality', 'Quality', 'fitting_quality'):
        label = peak.get(key)
        if isinstance(label, str) and label:
            return label
    return None


def quality_color(label: Optional[str]) -> str:
    """Colour of a peak with this quality label; a missing label is unrated."""
    return QUALITY_COLORS.get(str(label).lower(), UNRATED_COLOR)


def quality_symbol(label: Optional[str]) -> str:
    """Symbol prefixed to a peak's name in a navigator, matching quality_color."""
    return QUALITY_SYMBOLS.get(str(label).lower(), UNRATED_SYMBOL)


def navigator_row_display(label: Optional[str]) -> tuple:
    """(symbol, row background) of a navigator row whose peak carries this label.

    A peak with no label has no fit to rate, and its row is left plain.
    """
    if label is None:
        return "", None
    return quality_symbol(label), QColor(quality_color(label)).lighter(180)
