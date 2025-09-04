"""
Validation and Testing Tools for LunaNMR

This module contains validation tools for checking installation
and system requirements.

Author: Guillaume Mas
Date: 2025
"""

try:
    from .verify_installation import verify_dependencies, check_installation
except ImportError:
    pass

__all__ = [
    'verify_dependencies',
    'check_installation'
]
