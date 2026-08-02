"""ABOUTME: Dialog modules for lunaNMR Qt interface providing specialized popup windows
ABOUTME: Includes data loading, configuration, spectrum browsers, and series integration dialogs
"""

from .data_loading_dialog import DataLoadingDialog
from .configuration_dialog import ConfigurationDialog
from .spectrum_browser_dialog import SpectrumBrowserDialog
from .multi_spectrum_viewer_dialog import MultiSpectrumViewerDialog
from .spectrum_viewer_dialog import SpectrumViewerDialog
from .series_integration_dialog import SeriesIntegrationDialog
from .series_manager_dialog import SeriesManagerDialog
from .results_browser_dialog import ResultsBrowserDialog
from .missing_files_dialog import MissingFilesDialog
from .dynamixs_dialog import DynamiXsDialog
from .kd_titration_dialog import KdTitrationDialog
from .series_qc_dialog import SeriesQCDialog
from .spectral_inspector import SpectralInspector
from .integration_1d_dialog import Integration1DDialog

__all__ = [
    'DataLoadingDialog',
    'ConfigurationDialog',
    'SpectrumBrowserDialog',
    'MultiSpectrumViewerDialog',
    'SpectrumViewerDialog',
    'SeriesIntegrationDialog',
    'SeriesManagerDialog',
    'ResultsBrowserDialog',
    'MissingFilesDialog',
    'DynamiXsDialog',
    'KdTitrationDialog',
    'SeriesQCDialog',
    'SpectralInspector',
    'Integration1DDialog',
]
