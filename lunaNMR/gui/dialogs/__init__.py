"""ABOUTME: Dialog modules for lunaNMR Qt interface providing specialized popup windows
ABOUTME: Includes data loading, configuration, spectrum browsers, and series integration dialogs
"""

from .data_loading_dialog import DataLoadingDialog
from .configuration_dialog import ConfigurationDialog
from .spectrum_browser_dialog import SpectrumBrowserDialog
from .multi_spectrum_viewer_dialog import MultiSpectrumViewerDialog
from .spectrum_viewer_dialog import SpectrumViewerDialog
from .series_integration_dialog import SeriesIntegrationDialog
from .results_browser_dialog import ResultsBrowserDialog

__all__ = [
    'DataLoadingDialog',
    'ConfigurationDialog',
    'SpectrumBrowserDialog',
    'MultiSpectrumViewerDialog',
    'SpectrumViewerDialog',
    'SeriesIntegrationDialog',
    'ResultsBrowserDialog'
]
