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

try:
    from .font_config import FontManager, get_font_manager, configure_emoji_support, get_display_text
except ImportError:
    pass

try:
    from .delay_extractor import DelayExtractor
except ImportError:
    pass

try:
    from .relaxation_data_bridge import RelaxationDataBridge
except ImportError:
    pass

try:
    from .output_manager import (
        OutputManager,
        OutputLevel,
        log_progress,
        log_info,
        log_warning,
        log_error
    )
except ImportError:
    pass

try:
    from .project_manager import ProjectManager
except ImportError:
    pass

try:
    from .undo_manager import (
        PeakUndoManager,
        DeletePeakCommand,
        DeleteMultiplePeaksCommand,
        MovePeakCommand,
        AddPeakCommand
    )
except ImportError:
    pass

__all__ = [
    'ConfigManager',
    'FileManager',
    'ParameterManager',
    'GlobalOptimizationManager',
    'FontManager',
    'get_font_manager',
    'configure_emoji_support',
    'get_display_text',
    'DelayExtractor',
    'RelaxationDataBridge',
    'OutputManager',
    'OutputLevel',
    'log_progress',
    'log_info',
    'log_warning',
    'log_error',
    'ProjectManager',
    'PeakUndoManager',
    'DeletePeakCommand',
    'DeleteMultiplePeaksCommand',
    'MovePeakCommand',
    'AddPeakCommand',
]
