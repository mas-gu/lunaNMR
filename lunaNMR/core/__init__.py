# ABOUTME: Core processing engines for LunaNMR: integrator, Voigt fitter, peak picker.
# ABOUTME: Exported lazily (PEP 562) so importing a light submodule never loads the integrator.

"""
Core Processing Engines for LunaNMR

This module contains the core processing engines including
integrators, fitters, and peak pickers.

The engines are imported on first access. ``lunaNMR.utils.parameter_manager`` imports
``lunaNMR.core.ps2d_config`` while ``lunaNMR.utils`` is still initialising; an eager
import of ``core_integrator`` here would reach back into ``lunaNMR.integrators`` and
``lunaNMR.utils`` half-loaded.
"""

import importlib

_LAZY_EXPORTS = {
    'EnhancedVoigtIntegrator': '.core_integrator',
    'EnhancedVoigtFitter': '.enhanced_voigt_fitter',
    'EnhancedPeakPicker': '.enhanced_peak_picker',
}


def __getattr__(name):
    module_path = _LAZY_EXPORTS.get(name)
    if module_path is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(importlib.import_module(module_path, __name__), name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(list(globals().keys()) + list(_LAZY_EXPORTS.keys()))


__all__ = [
    'EnhancedVoigtIntegrator',
    'EnhancedVoigtFitter',
    'EnhancedPeakPicker'
]
