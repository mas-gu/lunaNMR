"""
Spectrum Processors for LunaNMR

This module contains various spectrum processing classes for
different types of NMR data analysis workflows.
"""

try:
    from .multi_spectrum_processor import MultiSpectrumProcessor
except ImportError:
    pass

try:
    from .single_spectrum_processor import SingleSpectrumProcessor
except ImportError:
    pass

__all__ = [
    'MultiSpectrumProcessor',
    'SingleSpectrumProcessor'
]
