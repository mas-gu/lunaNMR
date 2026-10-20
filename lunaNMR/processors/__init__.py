"""
Spectrum Processors for LunaNMR

This module contains various spectrum processing classes for
different types of NMR data analysis workflows.
"""

try:
    from .series_processor import SeriesProcessor
except ImportError:
    pass

try:
    from .multi_spectrum_processor import MultiSpectrumProcessor
except ImportError:
    pass

try:
    from .single_spectrum_processor import SingleSpectrumProcessor
except ImportError:
    pass

try:
    from .parallel_fitting import ParallelFitting
except ImportError:
    pass

__all__ = [
    'SeriesProcessor',
    'MultiSpectrumProcessor',
    'SingleSpectrumProcessor',
    'ParallelFitting'
]
