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

__version__ = "1.1.0"
__author__ = "Guillaume Mas"
__description__ = "Advanced NMR Peak Analysis and Integration"

# Public name -> owning lunaNMR subpackage module (leading '.').
_LAZY_EXPORTS = {
    'LunaNMRMainWindow': '.gui.main_window',
    'main': '.gui.main_window',
    'EnhancedVoigtIntegrator': '.core.core_integrator',
    'EnhancedVoigtFitter': '.core.enhanced_voigt_fitter',
    'MultiSpectrumProcessor': '.processors.multi_spectrum_processor',
    'SingleSpectrumProcessor': '.processors.single_spectrum_processor',
    'NMRFileManager': '.utils.file_manager',
    'NMRParameterManager': '.utils.parameter_manager',
}


def __getattr__(name):
    module_path = _LAZY_EXPORTS.get(name)
    if module_path is not None:
        try:
            module = importlib.import_module(module_path, __name__)
            value = getattr(module, name)
        except ImportError as exc:
            # An optional dependency (e.g. PySide6) is absent: surface the name as
            # simply unavailable so hasattr() works, as the old eager try/except did.
            raise AttributeError(f"{name!r} is unavailable: {exc}") from exc
        globals()[name] = value  # cache so subsequent access skips __getattr__
        return value
    # Otherwise treat `name` as a submodule so `lunaNMR.core`, `lunaNMR.utils`, etc.
    # resolve after a bare `import lunaNMR` (the old eager imports bound these).
    try:
        return importlib.import_module(f'.{name}', __name__)
    except ImportError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(list(globals().keys()) + list(_LAZY_EXPORTS.keys()))


__all__ = [
    # Main window
    'LunaNMRMainWindow',
    'main',
    # Core
    'EnhancedVoigtIntegrator',
    'EnhancedVoigtFitter',
    # Processors
    'MultiSpectrumProcessor',
    'SingleSpectrumProcessor',
    # Utilities
    'NMRFileManager',
    'NMRParameterManager',
]
