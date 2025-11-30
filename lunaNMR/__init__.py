# ABOUTME: Main package init for LunaNMR v1.0 Qt/PySide6 NMR analysis suite
# ABOUTME: Exports main window, core processing, and utility classes

"""
LunaNMR v1.0: Advanced NMR Peak Analysis and Integration

A comprehensive toolkit for NMR peak detection, fitting, and integration
with advanced Voigt profile analysis and multi-peak deconvolution.

This version uses Qt/PySide6 for the GUI framework.
"""

__version__ = "1.0.0"
__author__ = "Guillaume Mas"
__description__ = "Advanced NMR Peak Analysis and Integration"

# Main Qt GUI access
try:
    from .gui.main_window import LunaNMRMainWindow, main
except ImportError:
    pass

# Core functionality
try:
    from .core.core_integrator import EnhancedVoigtIntegrator
    from .core.enhanced_voigt_fitter import EnhancedVoigtFitter
    from .core.enhanced_peak_picker import EnhancedPeakPicker
except ImportError:
    pass

# Processors
try:
    from .processors.multi_spectrum_processor import MultiSpectrumProcessor
    from .processors.single_spectrum_processor import SingleSpectrumProcessor
    from .processors.parallel_fitting import ParallelVoigtProcessor
except ImportError:
    pass

# Utilities
try:
    from .utils.config_manager import ConfigurationManager
    from .utils.file_manager import NMRFileManager
    from .utils.parameter_manager import NMRParameterManager
except ImportError:
    pass

# DynamiXs Integration (optional submodule)
try:
    import sys
    import os
    # Add modules directory to path
    modules_path = os.path.join(os.path.dirname(__file__), '..', 'modules')
    if modules_path not in sys.path:
        sys.path.append(modules_path)

    from dynamiXs import DynamiXsGUI, run_dynamixs
except ImportError:
    pass

__all__ = [
    # Main window
    'LunaNMRMainWindow',
    'main',
    # Core
    'EnhancedVoigtIntegrator',
    'EnhancedVoigtFitter',
    'EnhancedPeakPicker',
    # Processors
    'MultiSpectrumProcessor',
    'SingleSpectrumProcessor',
    'ParallelVoigtProcessor',
    # Utilities
    'ConfigurationManager',
    'NMRFileManager',
    'NMRParameterManager',
]
