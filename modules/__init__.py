"""
Optional Modules for LunaNMR

This package contains optional modules that extend LunaNMR functionality.
These modules may have additional dependencies and are loaded conditionally.
"""

try:
    from .dynamiXs import *
except ImportError:
    pass

__all__ = []