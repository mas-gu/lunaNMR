"""
Specialized Integrators for LunaNMR

This module contains specialized integration classes for
different types of NMR data processing workflows.

Author: Guillaume Mas
Date: 2025
"""

try:
    from .inplace_advanced_nmr_integrator import InplaceAdvancedNMRIntegrator
except ImportError:
    pass

try:
    from .inplace_series_nmr_integrator import InplaceSeriesNMRIntegrator
except ImportError:
    pass

try:
    from .simple_pattern_matcher import SimplePatternMatcher
except ImportError:
    pass

__all__ = [
    'InplaceAdvancedNMRIntegrator',
    'InplaceSeriesNMRIntegrator',
    'SimplePatternMatcher'
]
