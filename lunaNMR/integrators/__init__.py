"""
Specialized Integrators for LunaNMR

This module contains specialized integration classes for
different types of NMR data processing workflows.

Author: Guillaume Mas
Date: 2025
"""

try:
    from .inplace_advanced_nmr_integrator import InPlaceAdvancedNMRIntegrator
except ImportError:
    pass

__all__ = [
    'InPlaceAdvancedNMRIntegrator'
]
