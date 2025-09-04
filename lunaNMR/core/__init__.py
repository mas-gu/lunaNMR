"""
Core Processing Engines for LunaNMR

This module contains the core processing engines including
integrators, fitters, and peak pickers.

Author: Guillaume Mas
Date: 2025
"""

try:
    from .core_integrator import CoreIntegrator
except ImportError:
    pass

try:
    from .enhanced_voigt_fitter import EnhancedVoigtFitter
except ImportError:
    pass

try:
    from .enhanced_peak_picker import EnhancedPeakPicker
except ImportError:
    pass

try:
    from .integrated_detection_fitter import IntegratedDetectionFitter
except ImportError:
    pass

__all__ = [
    'CoreIntegrator',
    'EnhancedVoigtFitter',
    'EnhancedPeakPicker',
    'IntegratedDetectionFitter'
]
