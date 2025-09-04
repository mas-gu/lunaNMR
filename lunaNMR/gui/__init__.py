"""
GUI Components for LunaNMR

This module contains all graphical user interface components including
the main GUI, spectrum browser, and visualization tools.

Author: Guillaume Mas
Date: 2025
"""

try:
    from .main_gui import NMRPeaksSeriesGUI
except ImportError:
    pass

try:
    from .spectrum_browser import SpectrumBrowser
except ImportError:
    pass

try:
    from .visualization import VoigtPlotter, create_residual_plot, create_comparison_plot
except ImportError:
    pass

try:
    from .gui_components import *
except ImportError:
    pass

__all__ = [
    'NMRPeaksSeriesGUI',
    'SpectrumBrowser',
    'VoigtPlotter',
    'create_residual_plot',
    'create_comparison_plot'
]
