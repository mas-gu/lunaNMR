# ABOUTME: Main package init for LunaNMR v1.0 Qt/PySide6 NMR analysis suite
# ABOUTME: Lazily exports main window, core processing, and utility classes on first access

"""
LunaNMR v1.0: Advanced NMR Peak Analysis and Integration

A comprehensive toolkit for NMR peak detection, fitting, and integration
with advanced Voigt profile analysis and multi-peak deconvolution.

This version uses Qt/PySide6 for the GUI framework.

The public classes are imported lazily (PEP 562 ``__getattr__``): ``import lunaNMR``
stays cheap and headless — the GUI, core, and nmrglue-backed integrators are only
loaded when their name is first accessed. This lets ``python -m lunaNMR`` run without
a display or the full scientific stack for subcommands that don't need them.
"""

import importlib
import os
import sys

__version__ = "1.0.0"
__author__ = "Guillaume Mas"
__description__ = "Advanced NMR Peak Analysis and Integration"

# Make the optional dynamiXs submodule importable (path side-effect only, no heavy import).
_modules_path = os.path.join(os.path.dirname(__file__), '..', 'modules')
if _modules_path not in sys.path:
    sys.path.append(_modules_path)

# Public name -> owning module. Leading '.' means a lunaNMR subpackage; otherwise a
# top-level module reachable via the dynamiXs path appended above.
_LAZY_EXPORTS = {
    'LunaNMRMainWindow': '.gui.main_window',
    'main': '.gui.main_window',
    'EnhancedVoigtIntegrator': '.core.core_integrator',
    'EnhancedVoigtFitter': '.core.enhanced_voigt_fitter',
    'EnhancedPeakPicker': '.core.enhanced_peak_picker',
    'MultiSpectrumProcessor': '.processors.multi_spectrum_processor',
    'SingleSpectrumProcessor': '.processors.single_spectrum_processor',
    'ParallelVoigtProcessor': '.processors.parallel_fitting',
    'ConfigurationManager': '.utils.config_manager',
    'NMRFileManager': '.utils.file_manager',
    'NMRParameterManager': '.utils.parameter_manager',
    'DynamiXsGUI': 'dynamiXs',
    'run_dynamixs': 'dynamiXs',
}


def __getattr__(name):
    module_path = _LAZY_EXPORTS.get(name)
    if module_path is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    package = __name__ if module_path.startswith('.') else None
    module = importlib.import_module(module_path, package)
    value = getattr(module, name)
    globals()[name] = value  # cache so subsequent access skips __getattr__
    return value


def __dir__():
    return sorted(list(globals().keys()) + list(_LAZY_EXPORTS.keys()))


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
