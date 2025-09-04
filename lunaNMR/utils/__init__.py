"""
Utility Classes and Managers for LunaNMR

This module contains utility classes for configuration management,
file handling, parameter management, and optimization.

Author: Guillaume Mas
Date: 2025
"""

try:
    from .config_manager import ConfigManager
except ImportError:
    pass

try:
    from .file_manager import FileManager
except ImportError:
    pass

try:
    from .parameter_manager import ParameterManager
except ImportError:
    pass

try:
    from .global_optimization_manager import GlobalOptimizationManager
except ImportError:
    pass

__all__ = [
    'ConfigManager',
    'FileManager',
    'ParameterManager',
    'GlobalOptimizationManager'
]
