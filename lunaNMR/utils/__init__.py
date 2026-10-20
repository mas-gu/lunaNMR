# ABOUTME: Utility classes and managers for LunaNMR: output routing, delays, projects, undo.
# ABOUTME: Exported lazily (PEP 562) so importing one utility never loads the others, Qt among them.

"""
Utility Classes and Managers for LunaNMR

This module contains utility classes for configuration management,
file handling and parameter management.

A name is imported on first access. Importing ``lunaNMR.utils.output_manager`` loads that
module alone; an eager import here would load the project manager and the Qt undo stack
with it, in every process that logs a line.
"""

import importlib

_LAZY_EXPORTS = {
    'DelayExtractor': '.delay_extractor',
    'RelaxationDataBridge': '.relaxation_data_bridge',
    'OutputManager': '.output_manager',
    'OutputLevel': '.output_manager',
    'log_progress': '.output_manager',
    'log_info': '.output_manager',
    'log_warning': '.output_manager',
    'log_error': '.output_manager',
    'ProjectManager': '.project_manager',
    'PeakUndoManager': '.undo_manager',
    'DeletePeakCommand': '.undo_manager',
    'DeleteMultiplePeaksCommand': '.undo_manager',
    'MovePeakCommand': '.undo_manager',
    'AddPeakCommand': '.undo_manager',
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


__all__ = list(_LAZY_EXPORTS)
