# ABOUTME: GUI module for LunaNMR Qt/PySide6 interface
# ABOUTME: Exports main window and dialog components for application launch

"""
GUI Components for LunaNMR v1.0 (Qt/PySide6)

This module contains all graphical user interface components including
the main window, dialogs, and visualization tools.

Note: Legacy Tkinter files have been moved to gui/_deprecated/

Author: Guillaume Mas
Date: 2025
"""

# Main Qt window
from .main_window import LunaNMRMainWindow, main

# Qt Dialogs
from .dialogs import (
    DataLoadingDialog,
    ConfigurationDialog,
    SpectrumBrowserDialog,
    MultiSpectrumViewerDialog,
    SpectrumViewerDialog,
    SeriesIntegrationDialog
)

# Qt Components
from .components import (
    MatplotlibWidget,
    MatplotlibMultiAxesWidget,
    PeakNavigator
)

__all__ = [
    # Main window
    'LunaNMRMainWindow',
    'main',
    # Dialogs
    'DataLoadingDialog',
    'ConfigurationDialog',
    'SpectrumBrowserDialog',
    'MultiSpectrumViewerDialog',
    'SpectrumViewerDialog',
    'SeriesIntegrationDialog',
    # Components
    'MatplotlibWidget',
    'MatplotlibMultiAxesWidget',
    'PeakNavigator',
]
