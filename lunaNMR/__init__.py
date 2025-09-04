"""
LunaNMR: Advanced NMR Peak Analysis and Integration

A comprehensive toolkit for NMR peak detection, fitting, and integration
with advanced Voigt profile analysis and multi-peak deconvolution.
"""

__version__ = "0.9.0"
__author__ = "Guillaume"
__description__ = "Advanced NMR Peak Analysis and Integration"

# Main GUI access
try:
    from .gui.main_gui import NMRPeaksSeriesGUI
except ImportError:
    pass

# Core functionality
try:
    from .core.core_integrator import CoreIntegrator
    from .core.enhanced_voigt_fitter import EnhancedVoigtFitter
    from .core.enhanced_peak_picker import EnhancedPeakPicker
    from .core.integrated_detection_fitter import IntegratedDetectionFitter
except ImportError:
    pass

# Processors
try:
    from .processors.series_processor import SeriesProcessor
    from .processors.multi_spectrum_processor import MultiSpectrumProcessor
    from .processors.single_spectrum_processor import SingleSpectrumProcessor
    from .processors.parallel_fitting import ParallelFitting
except ImportError:
    pass

# Utilities
try:
    from .utils.config_manager import ConfigManager
    from .utils.file_manager import FileManager
    from .utils.parameter_manager import ParameterManager
    from .utils.global_optimization_manager import GlobalOptimizationManager
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
    __all__.extend(['DynamiXsGUI', 'run_dynamixs'])
except ImportError:
    pass

# Backwards compatibility
try:
    from .gui.main_gui import NMRPeaksSeriesGUI as lunaNMRv0o9_gui
except ImportError:
    pass

__all__ = [
    'NMRPeaksSeriesGUI',
    'CoreIntegrator', 
    'EnhancedVoigtFitter',
    'EnhancedPeakPicker',
    'IntegratedDetectionFitter',
    'SeriesProcessor',
    'MultiSpectrumProcessor', 
    'SingleSpectrumProcessor',
    'ParallelFitting',
    'ConfigManager',
    'FileManager',
    'ParameterManager',
    'GlobalOptimizationManager'
]