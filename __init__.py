"""
NMR Peak Series Analysis Package

This package provides a comprehensive NMR peak detection, integration, and series analysis toolkit
with both full peak detection and in-place fitting capabilities.

Modules:
- core_integrator: Core VoigtIntegrator class with advanced fitting capabilities
- gui_components: Reusable GUI components (frames, dialogs, widgets)
- file_manager: File handling and data loading utilities
- data_processor: Data processing and analysis functions
- series_processor: Series integration workflow management
- visualization: Plotting and visualization utilities
- config_manager: Configuration management and settings
- main_gui: Main GUI application
"""

__version__ = "1.0.0"
__author__ = "Guillaume Mas"

# Import key classes for easy access
from .core.core_integrator import VoigtIntegrator, EnhancedVoigtIntegrator
from .gui.gui_components import (
    ScrollableFrame,
    EnhancedFileListFrame,
    AdvancedProgressDialog,
    StatisticsPanel
)
from .processors.series_processor import SeriesProcessor
from .utils.config_manager import ConfigurationManager

__all__ = [
    'VoigtIntegrator',
    'EnhancedVoigtIntegrator',
    'ScrollableFrame',
    'EnhancedFileListFrame',
    'AdvancedProgressDialog',
    'StatisticsPanel',
    'SeriesProcessor',
    'ConfigurationManager'



]
