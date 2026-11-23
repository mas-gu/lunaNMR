"""ABOUTME: Qt reusable GUI components for lunaNMR
ABOUTME: Modular widgets and dialogs following the design system
"""

from lunaNMR.gui.components.progress_dialog import ProgressDialog
from lunaNMR.gui.components.dialogs import show_error
from lunaNMR.gui.components.peak_navigator import PeakNavigator, PeakEditDialog
from lunaNMR.gui.components.peak_navigator_table import PeakNavigatorTable
from lunaNMR.gui.components.matplotlib_widget import (
    MatplotlibWidget,
    MatplotlibMultiAxesWidget
)
from lunaNMR.gui.components.spectrum_plotter import SpectrumPlotter
from lunaNMR.gui.components.nmr_navigation_handler import NMRNavigationHandler

__all__ = [
    "ProgressDialog",
    "show_error",
    "PeakNavigator",
    "PeakEditDialog",
    "PeakNavigatorTable",
    "MatplotlibWidget",
    "MatplotlibMultiAxesWidget",
    "SpectrumPlotter",
    "NMRNavigationHandler"
]
