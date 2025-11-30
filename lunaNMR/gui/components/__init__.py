"""ABOUTME: Qt reusable GUI components for lunaNMR
ABOUTME: Modular widgets and dialogs following the design system
"""

from lunaNMR.gui.components.progress_dialog import ProgressDialog
from lunaNMR.gui.components.dialogs import (
    show_question,
    show_confirmation,
    show_error,
    show_warning,
    show_info,
    ConfirmationDialog,
    InputDialog,
    ask_yes_no,
    ask_ok_cancel
)
from lunaNMR.gui.components.peak_navigator import PeakNavigator, PeakEditDialog
from lunaNMR.gui.components.peak_navigator_table import PeakNavigatorTable
from lunaNMR.gui.components.matplotlib_widget import (
    MatplotlibWidget,
    MatplotlibMultiAxesWidget
)
from lunaNMR.gui.components.series_plotter import SeriesPlotter
from lunaNMR.gui.components.spectrum_plotter import SpectrumPlotter

__all__ = [
    "ProgressDialog",
    "show_question",
    "show_confirmation",
    "show_error",
    "show_warning",
    "show_info",
    "ConfirmationDialog",
    "InputDialog",
    "ask_yes_no",
    "ask_ok_cancel",
    "PeakNavigator",
    "PeakEditDialog",
    "PeakNavigatorTable",
    "MatplotlibWidget",
    "MatplotlibMultiAxesWidget",
    "SeriesPlotter",
    "SpectrumPlotter"
]
